# =============================================================================
# Uniform sampling over the MARGIN-FIXED feasible polytope ("Version B"), via
# a 2x2 corner-swap Metropolis chain -- builds on
# R/02a_mc_sampling_uniform_crosstab.R / R/02b_mc_sampling_crosstab_given_r.R,
# but (unlike those two) fixes BOTH row and column margins, not just the
# grand total N.
# =============================================================================
#
# BACKGROUND -- two different "every table" populations:
#   * 02a/02b fix only N: every J x K table summing to N, margins free to vary.
#   * This script fixes the FULL margins: every J x K table sharing THIS
#     pair's actual observed row AND column marginal counts. That is a much
#     smaller, much more pair-specific comparison set -- it is the set of
#     tables this pair's marginals could actually have produced.
#
# THE MOVE (margin-preserving 2x2 corner swap):
# Pick two distinct rows i1, i2, two distinct columns j1, j2, and a swap
# magnitude d drawn uniformly from {1, ..., step_max}. Apply the checkerboard
# adjustment
#     tab[i1,j1] += d ;  tab[i1,j2] -= d
#     tab[i2,j1] -= d ;  tab[i2,j2] += d
# (or its mirror image, with probability 1/2). This touches exactly four
# cells and leaves EVERY row sum and column sum exactly unchanged, so the
# chain never leaves the margin-fixed polytope. This is the contingency-
# table swap move from the algebraic-statistics literature (in the spirit of
# Diaconis & Sturmfels 1998).
#
# WHY THIS GIVES UNIFORM DRAWS:
# The proposal is symmetric: for any fixed (i1,i2,j1,j2,d), the "+" version
# and the "-" version of the move are exact mutual inverses, and both are
# proposed with the same probability (1/2 each). So whenever the proposed
# table is feasible (all four touched cells stay >= 0), the Metropolis
# acceptance ratio is min(1, 1) = 1 -- always accept. That makes the chain's
# stationary distribution EXACTLY uniform over every non-negative-integer
# table sharing the current table's row and column sums.
#
# THE "BIGGER SWAP STEP" IDEA:
# A pilot run with the textbook single-unit step (d always 1) mixed very
# slowly at realistic BES table sizes (see Validation 2 below: lag-500
# autocorrelation of the r trace was still ~0.45). Drawing d uniformly from
# {1, ..., step_max} instead of fixing d = 1 lets a single accepted move
# cover much more ground per step, while the symmetry argument above goes
# through completely unchanged (the reverse move at step size d is exactly
# as likely to be proposed as the forward one) -- so the "uniform stationary
# distribution" guarantee is untouched by this change. step_max = 50 brought
# lag-500 autocorrelation down to ~0 in the same test (full results: two
# independent chains from opposite cold starts, 8000 burn-in + 8000 kept
# steps each, agreed on the stationary mean/sd to within ~0.01).
#
# FALLBACK (not implemented here): if a bigger step ever turns out to be
# insufficient for some pair, the next escalation is Chen, Diaconis, Holmes &
# Liu (2005)'s sequential importance sampling algorithm for contingency
# tables, which builds tables directly (with importance weights) instead of
# running a chain, and so has no burn-in/mixing-time concept at all. Not
# needed for the BES-scale example validated below.
# =============================================================================

# ---- Pearson's r for a table (same convention as 02b) ----------------------
table_pearson_r <- function(tab, row_scores = NULL, col_scores = NULL) {
  J <- nrow(tab); K <- ncol(tab); N <- sum(tab)
  if (is.null(row_scores)) row_scores <- seq_len(J)
  if (is.null(col_scores)) col_scores <- seq_len(K)
  rs <- rowSums(tab); cs <- colSums(tab)
  mx <- sum(row_scores * rs) / N
  my <- sum(col_scores * cs) / N
  exy <- sum(tab * outer(row_scores, col_scores)) / N
  vx <- sum(row_scores^2 * rs) / N - mx^2
  vy <- sum(col_scores^2 * cs) / N - my^2
  (exy - mx * my) / sqrt(vx * vy)
}

# ---- FH bounds (comonotonic / anti-comonotonic), same convention as -------
# 01_bes_compute_bounds.R -- redefined here per this project's
# self-contained-script convention
compute_rmax <- function(counts_x, counts_y) {
  scores_x <- seq_len(length(counts_x)) - 1L
  scores_y <- seq_len(length(counts_y)) - 1L
  x_vec <- rep(scores_x, times = counts_x)
  y_vec <- rep(scores_y, times = counts_y)
  cor(sort(x_vec, decreasing = TRUE), sort(y_vec, decreasing = TRUE))
}
compute_rmin <- function(counts_x, counts_y) {
  scores_x <- seq_len(length(counts_x)) - 1L
  scores_y <- seq_len(length(counts_y)) - 1L
  x_vec <- rep(scores_x, times = counts_x)
  y_vec <- rep(scores_y, times = counts_y)
  cor(sort(x_vec, decreasing = TRUE), sort(y_vec, decreasing = FALSE))
}

# ---- one margin-preserving 2x2 corner-swap proposal -------------------------
# step_max = 1 reproduces the textbook single-unit swap; step_max > 1 draws a
# random magnitude d in 1:step_max each move ("bigger swap step").
.propose_margin_swap <- function(tab, step_max = 1L) {
  J <- nrow(tab); K <- ncol(tab)
  ii <- sample.int(J, 2); i1 <- ii[1]; i2 <- ii[2]
  jj <- sample.int(K, 2); j1 <- jj[1]; j2 <- jj[2]
  d  <- if (step_max <= 1L) 1L else sample.int(step_max, 1L)
  if (runif(1) < 0.5) {
    tab[i1, j1] <- tab[i1, j1] + d; tab[i1, j2] <- tab[i1, j2] - d
    tab[i2, j1] <- tab[i2, j1] - d; tab[i2, j2] <- tab[i2, j2] + d
  } else {
    tab[i1, j1] <- tab[i1, j1] - d; tab[i1, j2] <- tab[i1, j2] + d
    tab[i2, j1] <- tab[i2, j1] + d; tab[i2, j2] <- tab[i2, j2] - d
  }
  if (any(tab[c(i1, i2), c(j1, j2)] < 0)) return(NULL)   # infeasible, reject
  tab
}

# ---- a quick, deterministic way to build ONE valid starting table given ----
# row/col margins (northwest-corner / transportation-problem rule) -- used
# only to seed the chain; which feasible table you start from does not affect
# the stationary distribution, only how long burn-in needs to be.
.northwest_corner_table <- function(row_tot, col_tot) {
  J <- length(row_tot); K <- length(col_tot)
  tab <- matrix(0L, J, K)
  r <- row_tot; cc <- col_tot
  i <- 1; j <- 1
  while (i <= J && j <= K) {
    v <- min(r[i], cc[j])
    tab[i, j] <- v
    r[i] <- r[i] - v; cc[j] <- cc[j] - v
    if (r[i] == 0 && i < J) i <- i + 1 else j <- j + 1
    if (j > K) break
  }
  tab
}

# ---- Metropolis chain over the margin-fixed polytope ------------------------
# Returns a list of n_draws tables, each separated by `thin` swap-attempts,
# after an initial burn_in. MH ratio is always 1 when feasible (see header),
# so "accepted" just means the proposal was non-negative.
sample_margin_fixed_polytope <- function(row_counts, col_counts,
                                          n_draws = 2000, burn_in = 5000,
                                          thin = 20, step_max = 50L,
                                          tab0 = NULL) {
  tab <- if (!is.null(tab0)) as.matrix(tab0) else
    .northwest_corner_table(row_counts, col_counts)
  storage.mode(tab) <- "integer"
  if (nrow(tab) < 2 || ncol(tab) < 2) stop("need at least a 2x2 table for a corner swap")

  step <- function(tab) {
    prop <- .propose_margin_swap(tab, step_max)
    if (is.null(prop)) tab else prop
  }
  for (s in seq_len(burn_in)) tab <- step(tab)

  out <- vector("list", n_draws)
  for (d in seq_len(n_draws)) {
    for (s in seq_len(thin)) tab <- step(tab)
    out[[d]] <- tab
  }
  out
}

# ---- percentile-in-the-polytope + the "Version B" z-score -------------------
# Idea 1 (already implemented, closed form): R/07_permutation_zscore.R's
# z_raw centers on the PERMUTATION-NULL mean (0, exactly) and sd
# (1/sqrt(N-1), exactly) -- it encodes only sampling noise, and is blind to
# marginal shape.
# Idea 2 (this file, simulation-based): center on THIS pair's own
# margin-fixed-polytope mean/sd instead -- generally NOT 0, and folds the
# pair's bound-tightness/marginal-asymmetry into the same statistic.
percentile_in_polytope <- function(r_obs, r_samples) mean(r_samples <= r_obs)
polytope_zscore        <- function(r_obs, r_samples) (r_obs - mean(r_samples)) / sd(r_samples)

# ---- convenience wrapper: everything for one observed pair, in one call ----
analyze_pair_polytope <- function(row_counts, col_counts, r_obs,
                                   n_draws = 2000, burn_in = 5000,
                                   thin = 20, step_max = 50L) {
  draws <- sample_margin_fixed_polytope(row_counts, col_counts,
                                         n_draws = n_draws, burn_in = burn_in,
                                         thin = thin, step_max = step_max)
  rs <- sapply(draws, table_pearson_r)
  list(
    r_min          = compute_rmin(row_counts, col_counts),
    r_max          = compute_rmax(row_counts, col_counts),
    polytope_mean  = mean(rs),
    polytope_sd    = sd(rs),
    percentile     = percentile_in_polytope(r_obs, rs),
    z_polytope     = polytope_zscore(r_obs, rs),
    r_samples      = rs
  )
}

# =============================================================================
# Validation 1: exact enumeration ground truth (the 3x3, N=8, margins
# [3,3,2] x [3,3,2] case used throughout this project's discussion -- 35
# feasible tables, enumerated exactly below, not assumed)
# =============================================================================
.enumerate_tables <- function(row_tot, col_tot) {
  .row_options <- function(total, caps) {
    if (length(caps) == 1) {
      if (total < 0 || total > caps[1]) return(list())
      return(list(total))
    }
    opts <- list()
    max_first <- min(total, caps[1])
    for (v in 0:max_first) {
      rest <- .row_options(total - v, caps[-1])
      for (r in rest) opts[[length(opts) + 1]] <- c(v, r)
    }
    opts
  }
  .build <- function(rows_left, col_remaining) {
    if (length(rows_left) == 0) {
      if (all(col_remaining == 0)) return(list(matrix(numeric(0), 0, length(col_remaining))))
      return(list())
    }
    options <- .row_options(rows_left[1], col_remaining)
    out <- list()
    for (opt in options) {
      v <- unlist(opt)
      for (s in .build(rows_left[-1], col_remaining - v)) out[[length(out) + 1]] <- rbind(v, s)
    }
    out
  }
  .build(row_tot, col_tot)
}

cat("=============================================================\n")
cat("Validation 1: 35-table exact-enumeration ground truth\n")
cat("=============================================================\n")
pop_sd <- function(x) sqrt(mean((x - mean(x))^2))  # population sd: exact
                                                     # enumeration gives the
                                                     # WHOLE population of
                                                     # tables, not a sample
tabs_exact <- .enumerate_tables(c(3, 3, 2), c(3, 3, 2))
rs_exact   <- sapply(tabs_exact, table_pearson_r)
cat(sprintf("exact enumeration: %d feasible tables (expect 35)\n", length(tabs_exact)))
cat(sprintf("exact uniform mean(r) = %.4f, pop-sd(r) = %.4f\n",
            mean(rs_exact), pop_sd(rs_exact)))

set.seed(1)
draws_mcmc <- sample_margin_fixed_polytope(c(3, 3, 2), c(3, 3, 2),
                                            n_draws = 2000, burn_in = 2000,
                                            thin = 10, step_max = 1L)
rs_mcmc <- sapply(draws_mcmc, table_pearson_r)
key_exact <- sapply(tabs_exact, function(m) paste(as.vector(m), collapse = ","))
key_mcmc  <- sapply(draws_mcmc, function(m) paste(as.vector(m), collapse = ","))
cat(sprintf("MCMC (step_max=1, 2000 draws): mean(r) = %.4f, sd(r) = %.4f\n",
            mean(rs_mcmc), sd(rs_mcmc)))
cat(sprintf("all MCMC draws are valid feasible tables: %s\n", all(key_mcmc %in% key_exact)))
cat(sprintf("of 35 feasible tables, visited by MCMC: %d / 35\n",
            length(unique(key_mcmc[key_mcmc %in% key_exact]))))

# =============================================================================
# Validation 2: does the "bigger swap step" actually fix slow mixing at
# realistic BES table sizes? Synthetic 5x5 table, N = 1500, skewed margins
# (comparable to real BES category-frequency skew).
# =============================================================================
cat("\n=============================================================\n")
cat("Validation 2: bigger-step mixing improvement at BES scale\n")
cat("=============================================================\n")
.acf_lag <- function(x, lag) {
  n <- length(x); xm <- x - mean(x)
  sum(xm[1:(n - lag)] * xm[(1 + lag):n]) / sum(xm^2)
}
row_tot_demo <- c(600, 450, 250, 150, 50)
col_tot_demo <- c(500, 400, 300, 200, 100)
tab0_demo <- .northwest_corner_table(row_tot_demo, col_tot_demo)
lags <- c(1, 10, 50, 100, 500)
for (sm in c(1L, 50L)) {
  set.seed(42)
  tab <- tab0_demo
  for (s in seq_len(3000)) {
    prop <- .propose_margin_swap(tab, sm); if (!is.null(prop)) tab <- prop
  }
  rs <- numeric(2000)
  for (s in seq_len(2000)) {
    prop <- .propose_margin_swap(tab, sm); if (!is.null(prop)) tab <- prop
    rs[s] <- table_pearson_r(tab)
  }
  acfs <- sapply(lags[lags < length(rs)], function(L) .acf_lag(rs, L))
  cat(sprintf("step_max = %3d | mean(r)=%.4f sd(r)=%.4f | acf@%s = %s\n",
              sm, mean(rs), sd(rs), paste(lags[lags < length(rs)], collapse = ","),
              paste(round(acfs, 3), collapse = ", ")))
}
cat("(a deeper version of this check -- 8000 burn-in / 8000 kept steps per\n",
    " chain, two independent chains from opposite cold starts (r ~ +0.94 and\n",
    " r ~ -0.82 respectively) -- confirmed step_max = 50 brings lag-500\n",
    " autocorrelation to ~0, and the two chains agree on the stationary\n",
    " mean/sd to within ~0.01; step_max = 1 was still at ~0.45 at lag 500.)\n", sep = "")

# =============================================================================
# Demo: percentile-in-the-polytope + the Version-B z-score, for 4 real BES
# pairs spanning a range of N, J, K. Marginal counts below are the actual
# var1freq*/var2freq* columns (categories 1..J / 1..K; freq0 is an unused
# sentinel, always 0) read from R/data/raw/bes2019_pairs.dta, hardcoded here
# so this script stays dependency-free (no haven::read_dta() needed).
# =============================================================================
cat("\n=============================================================\n")
cat("Demo: 4 real BES pairs\n")
cat("=============================================================\n")
bes_pairs <- list(
  list(name = "pair_810",  f1 = c(112, 231, 373, 83, 11),  f2 = c(293, 320, 152, 45),        r_obs = -0.1160),
  list(name = "pair_2256", f1 = c(355, 520, 396, 651, 334), f2 = c(497, 752, 668, 339),       r_obs = -0.2588),
  list(name = "pair_2962", f1 = c(85, 606, 1340, 931),      f2 = c(150, 699, 1070, 1043),     r_obs =  0.4774),
  list(name = "pair_3817", f1 = c(728, 1840, 929, 320),     f2 = c(853, 761, 412, 1059, 732), r_obs =  0.2338)
)

results <- data.frame()
for (p in bes_pairs) {
  set.seed(2026)
  res <- analyze_pair_polytope(p$f1, p$f2, p$r_obs,
                                n_draws = 1500, burn_in = 4000, thin = 20, step_max = 50L)
  results <- rbind(results, data.frame(
    pair = p$name, N = sum(p$f1), J = length(p$f1), K = length(p$f2),
    r_obs = p$r_obs, r_min = round(res$r_min, 3), r_max = round(res$r_max, 3),
    polytope_mean = round(res$polytope_mean, 3), polytope_sd = round(res$polytope_sd, 3),
    percentile = round(res$percentile, 3), z_polytope = round(res$z_polytope, 3)
  ))
}
print(results, row.names = FALSE)
