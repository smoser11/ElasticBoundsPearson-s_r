# =============================================================================
# Uniform sampling over the MARGIN-FIXED polytope ("Version B")
# Every J x K table sharing a SPECIFIC pair's actual row AND column marginals,
# weighted EQUALLY -- builds the "percentile-in-the-feasible-polytope" idea
# described in paper/PlainLanguageExplanations.md.
# =============================================================================
#
# HOW THIS DIFFERS FROM THE REST OF R/02*:
#   02a fixes only the grand total N (margins float freely across draws).
#   02b builds on 02a to hit a TARGET r, with margins still free.
#   This file (02c) is the one that actually fixes BOTH row and column
#   margins -- the polytope of tables this pair's marginals could have
#   produced. Margin-free (02a-style) sampling has no Markov chain at all,
#   since stars-and-bars draws each table directly and exactly; margin-FIXED
#   sampling does not have an equivalent closed-form draw, so this file uses
#   MCMC instead.
#
# HOW THIS DIFFERS FROM 07_permutation_zscore.R's z_raw:
#   z_raw is built on the PERMUTATION-NULL (multivariate hypergeometric)
#   distribution over the same margin-fixed tables -- closed form, mean = 0
#   exactly, sd = 1/sqrt(N-1) exactly, concentrated near "looks independent."
#   This file's distribution is UNIFORM over the same set of tables -- every
#   table equally likely regardless of how independent-looking it is. They
#   are genuinely different distributions over the same tables, not two
#   versions of the same thing. See the validation block at the bottom for a
#   worked numeric comparison (35-table example: permutation-null gives
#   mean=0.0000, sd=0.3780 = 1/sqrt(7) exactly; uniform-over-polytope gives
#   mean=0.0330, sd=0.4562).
#
# THE SAMPLER -- variable-magnitude 2x2 corner-swap MCMC:
#   Pick two rows and two columns at random, and a random step size
#   s in {1, ..., max_step} with a random sign. Apply the checkerboard
#   +s/-s/-s/+s pattern at the four intersection cells. This preserves both
#   touched row sums and both touched column sums exactly (every other row/
#   column is untouched), so every visited table stays in the margin-fixed
#   polytope. The proposal is symmetric (the reverse move has identical
#   proposal probability), so the Metropolis acceptance ratio is exactly 1
#   whenever the move keeps all four touched cells non-negative -- this is
#   what makes the stationary distribution exactly uniform, not just
#   approximately so. max_step = 1 reproduces the classic single-unit swap
#   (Diaconis & Sturmfels 1998-style); max_step > 1 is the "bigger swap step"
#   idea for improving mixing at realistic BES table sizes (N in the
#   hundreds to thousands) -- see the mixing comparison below, which is the
#   actual evidence behind choosing to try this before reaching for Chen,
#   Diaconis, Holmes & Liu (2005)'s sequential importance sampling. SIS is
#   the documented fallback if a pair's table still won't mix even with a
#   large max_step; it is NOT implemented in this file.
# =============================================================================

# ---- Pearson's r for a table (identical convention to 02b) -----------------
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

# ---- Feasible starting table (northwest-corner / transportation rule) ------
# Not random -- just ONE valid table satisfying row_marg/col_marg exactly, to
# seed the MCMC chain. Standard construction: fill cells left-to-right,
# top-to-bottom, each time taking as much as the smaller of the remaining row
# capacity and remaining column capacity.
nw_corner_table <- function(row_marg, col_marg) {
  J <- length(row_marg); K <- length(col_marg)
  stopifnot(sum(row_marg) == sum(col_marg))
  r <- as.integer(row_marg); cc <- as.integer(col_marg)
  tab <- matrix(0L, nrow = J, ncol = K)
  i <- 1L; j <- 1L
  while (i <= J && j <= K) {
    val <- min(r[i], cc[j])
    tab[i, j] <- val
    r[i] <- r[i] - val
    cc[j] <- cc[j] - val
    if (r[i] == 0L && i < J) {
      i <- i + 1L
    } else if (cc[j] == 0L && j < K) {
      j <- j + 1L
    } else {
      break
    }
  }
  tab
}

# ---- Variable-magnitude 2x2 corner-swap proposal (symmetric) ---------------
# Returns NULL if the move would push any of the 4 touched cells negative
# (caller treats NULL as a rejection -- stay at the current table).
.propose_swap_var <- function(tab, max_step = 1L) {
  J <- nrow(tab); K <- ncol(tab)
  rows <- sample.int(J, 2)
  cols <- sample.int(K, 2)
  s <- sample.int(max_step, 1)
  if (runif(1) < 0.5) s <- -s
  i1 <- rows[1]; i2 <- rows[2]; j1 <- cols[1]; j2 <- cols[2]
  d11 <- tab[i1, j1] + s
  d12 <- tab[i1, j2] - s
  d21 <- tab[i2, j1] - s
  d22 <- tab[i2, j2] + s
  if (d11 < 0 || d12 < 0 || d21 < 0 || d22 < 0) return(NULL)
  tab[i1, j1] <- d11; tab[i1, j2] <- d12
  tab[i2, j1] <- d21; tab[i2, j2] <- d22
  tab
}

# =============================================================================
# Main driver: uniform draws over the margin-fixed polytope
# =============================================================================
# n_iter       total MCMC steps to run
# burn_in      steps discarded before any draw is kept
# thin         keep one draw every `thin` post-burn-in steps
# max_step     the "bigger swap step" knob -- see mixing comparison below
# track_full_chain  if TRUE, also return the un-thinned r trace (for
#                   diagnose_mixing()); costs O(n_iter) memory, off by default
sample_polytope_uniform <- function(row_marg, col_marg,
                                     n_iter = 20000, burn_in = 2000, thin = 10,
                                     max_step = 1L, init_tab = NULL,
                                     row_scores = NULL, col_scores = NULL,
                                     track_full_chain = FALSE) {
  stopifnot(sum(row_marg) == sum(col_marg))
  tab <- if (!is.null(init_tab)) init_tab else nw_corner_table(row_marg, col_marg)

  n_accept <- 0L
  n_draws_max <- max(0L, (n_iter - burn_in) %/% thin)
  draws_r <- numeric(n_draws_max)
  d_idx <- 0L
  full_chain_r <- if (track_full_chain) numeric(n_iter) else NULL

  for (it in seq_len(n_iter)) {
    prop <- .propose_swap_var(tab, max_step)
    if (!is.null(prop)) {
      tab <- prop
      n_accept <- n_accept + 1L
    }
    if (track_full_chain) {
      full_chain_r[it] <- table_pearson_r(tab, row_scores, col_scores)
    }
    if (it > burn_in && (it - burn_in) %% thin == 0) {
      d_idx <- d_idx + 1L
      draws_r[d_idx] <- table_pearson_r(tab, row_scores, col_scores)
    }
  }

  list(
    draws_r = draws_r[seq_len(d_idx)],
    final_table = tab,
    acceptance_rate = n_accept / n_iter,
    full_chain_r = full_chain_r
  )
}

# =============================================================================
# The two rescaling statistics this file exists to produce
# =============================================================================

# Where r_obs falls within the empirical Version-B (uniform-over-polytope)
# distribution.
polytope_percentile <- function(r_obs, draws_r) {
  mean(draws_r <= r_obs)
}

# Idea-2 z-score: rescale r_obs against the polytope's OWN mean/sd, instead of
# z_raw's permutation-null mean=0, sd=1/sqrt(N-1) (07_permutation_zscore.R).
# This is NOT centered at 0 in general -- the polytope mean tracks this
# pair's bound asymmetry (abs_asym / C1 / C2 in 01_bes_compute_bounds.R), so
# this statistic blends bound-tightness information into the same number
# that z_raw deliberately keeps separate. See the validation block for why.
polytope_zscore <- function(r_obs, draws_r) {
  (r_obs - mean(draws_r)) / sd(draws_r)
}

# =============================================================================
# Mixing diagnostics for a tracked full chain (track_full_chain = TRUE)
# =============================================================================
# Reports lag-k autocorrelation of the r trace at the requested lags, plus a
# simple effective-sample-size estimate (Geyer's initial-positive-sequence
# idea): ESS = n / (1 + 2 * sum of the leading run of positive autocorrelations).
diagnose_mixing <- function(full_chain_r, lags = c(1, 10, 50, 100, 500)) {
  n <- length(full_chain_r)
  ac <- acf(full_chain_r, lag.max = max(lags), plot = FALSE)$acf[, 1, 1]
  rho_at_lags <- ac[lags + 1]  # acf()'s first entry is lag 0
  names(rho_at_lags) <- paste0("lag", lags)

  rho_all <- ac[-1]
  keep <- cumsum(rho_all <= 0) == 0  # leading run of positive autocorrelations
  ess <- n / (1 + 2 * sum(rho_all[keep]))

  list(rho_at_lags = rho_at_lags, ess = ess, n = n)
}

# =============================================================================
# EXACT enumeration -- ground truth for small cases ONLY (combinatorial
# blowup for real BES-sized tables; this is a validation tool, not meant to
# run on production-scale pairs).
# =============================================================================
.compositions_bounded <- function(total, upper) {
  K <- length(upper)
  if (K == 1) {
    if (total >= 0 && total <= upper[1]) return(list(total)) else return(list())
  }
  out <- list()
  max_first <- min(total, upper[1])
  for (v1 in 0:max_first) {
    for (rest in .compositions_bounded(total - v1, upper[-1])) {
      out[[length(out) + 1]] <- c(v1, rest)
    }
  }
  out
}

enumerate_tables <- function(row_marg, col_marg) {
  J <- length(row_marg); K <- length(col_marg)
  stopifnot(sum(row_marg) == sum(col_marg))

  rec <- function(i, col_rem) {
    if (i == J) return(list(matrix(col_rem, nrow = 1)))  # last row is forced
    out <- list()
    for (v in .compositions_bounded(row_marg[i], col_rem)) {
      for (sub in rec(i + 1, col_rem - v)) {
        out[[length(out) + 1]] <- rbind(matrix(v, nrow = 1), sub)
      }
    }
    out
  }

  lapply(rec(1, col_marg), function(m) matrix(as.integer(m), nrow = J, ncol = K))
}

# =============================================================================
# VALIDATION 1 -- the familiar 35-table example (margins [3,3,2] x [3,3,2],
# N = 8), checked against EXACT enumeration computed right here in R (not
# borrowed from any earlier scratch work). Population-convention sd is used
# throughout this validation (divide by the count of tables / draws, not
# count-1) since the exact side is a fully known finite population, not a
# sample of some larger population -- using R's default sd() (which divides
# by n-1) against the exact side would silently introduce a ~1.5% mismatch
# here that has nothing to do with whether the sampler is correct.
# =============================================================================
cat("=== VALIDATION 1: 35-table example, margins [3,3,2] x [3,3,2], N=8 ===\n")

pop_sd <- function(x) sqrt(mean((x - mean(x))^2))

row_marg1 <- c(3, 3, 2); col_marg1 <- c(3, 3, 2)
exact_tables <- enumerate_tables(row_marg1, col_marg1)
cat(sprintf("Exact enumeration found %d tables (expect 35)\n", length(exact_tables)))

exact_rs <- sapply(exact_tables, table_pearson_r)
cat(sprintf("Exact, uniform-weighted:   mean(r) = %.4f   sd(r) = %.4f\n",
            mean(exact_rs), pop_sd(exact_rs)))
cat("  (expect mean=0.0330, sd=0.4562 -- matches the permutation-null comparison\n")
cat("   in PlainLanguageExplanations.md, which used mean=0.0000, sd=0.3780=1/sqrt(7))\n")

set.seed(20260626)
mc1 <- sample_polytope_uniform(row_marg1, col_marg1, n_iter = 200000, burn_in = 5000,
                                thin = 20, max_step = 2)
cat(sprintf("\nMCMC (max_step=2, n_draws=%d, accept_rate=%.3f):\n",
            length(mc1$draws_r), mc1$acceptance_rate))
cat(sprintf("  mean(r) = %.4f   sd(r) = %.4f\n", mean(mc1$draws_r), pop_sd(mc1$draws_r)))
cat(sprintf("  deviation from exact: dmean = %+.4f   dsd = %+.4f\n",
            mean(mc1$draws_r) - mean(exact_rs), pop_sd(mc1$draws_r) - pop_sd(exact_rs)))

# Coverage check: does the chain actually visit every one of the 35 exact
# tables, not just hover near a subset?
key_of <- function(tab) paste(as.vector(tab), collapse = ",")
exact_keys <- vapply(exact_tables, key_of, character(1))

set.seed(99)
tab <- nw_corner_table(row_marg1, col_marg1)
visited <- character(50000)
for (it in seq_len(50000)) {
  prop <- .propose_swap_var(tab, max_step = 2)
  if (!is.null(prop)) tab <- prop
  visited[it] <- key_of(tab)
}
cat(sprintf("\nCoverage: visited %d / %d distinct exact tables over 50000 iterations\n",
            length(intersect(unique(visited), exact_keys)), length(exact_keys)))

# =============================================================================
# VALIDATION 2 -- does the "bigger swap step" plan actually fix mixing at
# realistic BES table sizes? Same synthetic 5x5, N=1500, skewed-marginal
# case used earlier to FIND the mixing problem under max_step=1. Compares
# autocorrelation decay across several step sizes so the choice of max_step
# is an empirical one, not a guess.
# =============================================================================
cat("\n\n=== VALIDATION 2: mixing comparison, 5x5 / N=1500, skewed marginals ===\n")

row_marg2 <- c(600, 450, 250, 150, 50)
col_marg2 <- c(500, 400, 300, 200, 100)

for (max_step in c(1, 5, 15, 30, 60)) {
  set.seed(7)
  mc2 <- sample_polytope_uniform(row_marg2, col_marg2, n_iter = 100000, burn_in = 0,
                                  thin = 1, max_step = max_step, track_full_chain = TRUE)
  diag <- diagnose_mixing(mc2$full_chain_r, lags = c(1, 10, 50, 100, 500))
  cat(sprintf("max_step=%3d  accept=%.3f  ", max_step, mc2$acceptance_rate))
  cat(paste(sprintf("rho(lag%d)=%.3f", c(1, 10, 50, 100, 500), diag$rho_at_lags),
            collapse = "  "))
  cat(sprintf("  ESS=%.0f\n", diag$ess))
}
cat("(expect: lag-500 autocorrelation falling from ~0.95 at max_step=1 to\n")
cat(" well under 0.15 by max_step=30-60 -- if so, the bigger-step plan is\n")
cat(" working and Chen/Diaconis/Holmes/Liu (2005) SIS is not yet needed.)\n")

# =============================================================================
# VALIDATION 3 -- a real BES pair (not synthetic). Smallest table size (4x4)
# and smallest N (1891) among all 7,503 pairs, so the demo runs fast; the
# marginals and r_obs below are this pair's ACTUAL values from
# R/data/raw/bes2019_pairs.dta, not invented numbers.
# =============================================================================
cat("\n\n=== VALIDATION 3: a real BES pair (4x4, N=1891) ===\n")

library(haven)
bes <- read_dta("R/data/raw/bes2019_pairs.dta")
bes$cells <- bes$var1cats * bes$var2cats
row_pick <- bes[order(bes$cells, bes$nobs), ][1, ]

get_counts <- function(row, var_prefix, n_cats) {
  col0 <- paste0(var_prefix, "freq0")
  start <- if (col0 %in% names(row) && !is.na(row[[col0]]) && row[[col0]] > 0) 0 else 1
  vapply(start:(start + n_cats - 1), function(i) {
    col <- paste0(var_prefix, "freq", i)
    if (col %in% names(row) && !is.na(row[[col]])) as.numeric(row[[col]]) else 0
  }, numeric(1))
}

row_marg3 <- get_counts(row_pick, "var1", as.integer(row_pick$var1cats))
col_marg3 <- get_counts(row_pick, "var2", as.integer(row_pick$var2cats))
r_obs3 <- as.numeric(row_pick$corr)

cat(sprintf("row marginal: %s\n", paste(row_marg3, collapse = ", ")))
cat(sprintf("col marginal: %s\n", paste(col_marg3, collapse = ", ")))
cat(sprintf("N = %d   r_obs = %.4f\n", sum(row_marg3), r_obs3))

set.seed(2026)
mc3 <- sample_polytope_uniform(row_marg3, col_marg3, n_iter = 300000, burn_in = 10000,
                                thin = 30, max_step = 40)
cat(sprintf("\nMCMC (max_step=40, n_draws=%d, accept_rate=%.3f):\n",
            length(mc3$draws_r), mc3$acceptance_rate))
cat(sprintf("  polytope mean(r) = %.4f   sd(r) = %.4f\n",
            mean(mc3$draws_r), sd(mc3$draws_r)))
cat(sprintf("  percentile-in-polytope of r_obs: %.1f%%\n",
            100 * polytope_percentile(r_obs3, mc3$draws_r)))
cat(sprintf("  polytope z-score: %.3f\n",
            polytope_zscore(r_obs3, mc3$draws_r)))
cat("(expect roughly: polytope mean(r) near 0, sd(r) around 0.19-0.20,\n")
cat(" percentile in the mid/high-90s, z-score around 1.7-1.8 -- this pair's\n")
cat(" r_obs=0.338 sits unusually high even among ALL tables sharing its own\n")
cat(" marginals, not just relative to the two FH extremes.)\n")
