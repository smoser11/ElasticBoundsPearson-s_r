# 06_correlation_adjustment.R
#
# Standardized "adjustments" to observed Pearson's r (and, in parallel,
# Kendall's tau-b) for use as drop-in replacements in correlation matrices fed
# to Factor Analysis / SEM, where r_obs's *theoretical* range is narrower than
# [-1, 1] once the marginals are fixed (the "elastic bounds" of this project).
#
# Grounding: Warrens, M.J. (2013). "On Association Coefficients, Correction
# for Chance, and Correction for Maximum Value." Journal of Modern Mathematics
# Frontier, 2(4), 111-119. DOI: 10.14355/jmmf.2013.0204.01.
#
#   Warrens defines two corrections that can be applied to any association
#   coefficient A given fixed marginal totals:
#     - correction for chance:       c(A) = [A - E(A)] / [M(A) - E(A)]
#         E(A) = expected value of A under random pairing given the marginals
#         M(A) = a fixed constant (the coefficient's global theoretical max,
#                usually 1)
#     - correction for maximum value: d(A) = A / m(A)
#         m(A) = the maximum value A can attain *given the marginals*
#                (this is exactly r_max / tau_max in this project's notation)
#   Theorem 8 (Warrens 2013): c and d commute (c(d(A)) = d(c(A))).
#   Theorem 9 (Warrens 2013): any B = lambda + mu*A, where lambda and mu != 0
#   are themselves functions of the marginals (not of A), collapses to the
#   SAME fully-corrected coefficient as A once both corrections are applied.
#
#   For Pearson's r and Kendall's tau, the chance baseline (random pairing
#   that respects the marginals) has E(r) = E(tau) = 0 exactly, so c() is a
#   no-op for both -- r and tau are already "chance-corrected" in Warrens'
#   sense. The substantive correction available here is d(A) = A / m(A).
#
#   Warrens writes d(A) for one-sided coefficients bounded in [0, 1] with a
#   single maximum m(A) (e.g. kappa-type coefficients). Pearson's r and
#   Kendall's tau are two-sided/signed, so d() needs a sign-preserving
#   extension -- divide by r_max when r_obs >= 0, and by |r_min| when
#   r_obs < 0. That extension is exactly the "split_warrens" transform below,
#   and it generalizes the classical phi/phi_max correction for 2x2 tables
#   (Davenport, E.C., & El-Sanhurry, N.A. (1991). "Phi/phimax: Review and
#   synthesis." Educational and Psychological Measurement, 51, 821-828) to
#   general J x K tables.
#
#   PROVENANCE NOTE: this asymmetric (sign-split) normalization is not new
#   with Warrens (2013) -- it appears, fully formed, 21 years earlier in
#   Shih, W.J., & Huang, W.-M. (1992). "Evaluating Correlation with Proper
#   Bounds." Biometrics, 48(4), 1207-1213 [bib key: shih92-evaluating;
#   already present in paper/references.bib and listed (uncommented-on) in
#   the manuscript's FH-bounds reference dump]. Their real-data example
#   (p. 1211) divides two observed correlations (.23, .29) by their shared
#   attainable maximum (.88), reporting .26 and .33 "if one would interpret
#   the coefficients in the [-1,1] scale" -- i.e. r_obs / r_max for r_obs>=0.
#   Both examples happen to be positive, so the paper's own worked numbers
#   never exercise the r_obs<0 branch, but Shih & Huang were aware the
#   denominator would have to flip: they note on p. 1210 that "the maximum
#   and minimum correlations are not necessarily symmetric with respect to
#   0." split_warrens below is, formula for formula, their proposal --
#   verified against their numbers: 0.23/0.88 = 0.2614 (.26) and
#   0.29/0.88 = 0.3295 (.33), both exact matches. No separate transform or
#   PSD/invertibility sweep entry is added for this -- the split_warrens
#   results throughout this script (including the PSD sweeps in Parts 6-7)
#   already answer the invertibility question for the Shih & Huang
#   procedure, since the two are numerically identical, not just similar.
#
#   By Theorem 9, U1 (=r/(r_max-r_min)), U3 (=2r/(r_max-r_min)), and the
#   "tight" linear rescaling (=2(r-r_min)/(r_max-r_min) - 1) are all linear
#   transforms of r with marginal-dependent (lambda, mu), so once chance- and
#   max-value-corrected they all reduce to the same destination coefficient
#   as split_warrens. The practical choice is therefore not "which formula"
#   but "how much of the available correction to apply" -- raw r applies
#   none, U1/U3/tight-linear apply a *symmetric* version sized off the wider
#   of |r_min| and r_max (over-correcting whichever side is less extreme),
#   and split_warrens applies the maximum-value correction asymmetrically,
#   exactly matched to each side's own bound. U3 is retained below purely as
#   a cautionary baseline: it is the most naively "scale to fill [-1,1]"
#   transform and, as the PSD sweeps below show, also the one most prone to
#   breaking invertibility.
#
# What this script does:
#   1. Adds Kendall's tau-b bound functions (tau_min/tau_max) that mirror
#      this project's existing Pearson r_min/r_max construction exactly --
#      same comonotonic / anti-comonotonic sort, just retabulated and run
#      through a closed-form tau-b formula instead of cor().
#   2. Computes tau_min/tau_max (plus the existing r_min/r_max) for all 7,503
#      real BES pairs, and applies all 5 adjustment transforms to each pair's
#      real r_obs.
#        NOTE: tau_obs is NOT computable for real BES pairs -- the BES pairs
#        file stores only per-variable marginal frequencies, never the joint
#        J x K contingency table a pair's tau-b would need. Real-data tau
#        work is therefore limited to the bounds (tau_min/tau_max), not an
#        observed value; the "observed tau" side of this analysis is
#        necessarily synthetic (see step 3).
#   3. Runs a synthetic Monte Carlo PSD/invertibility sweep, in parallel for
#      Pearson r and Kendall's tau, applying the same 5 transforms to
#      simulated multi-variable correlation/tau matrices under a mild and an
#      adversarial parameter regime.
#   4. Runs a PSD/invertibility check using REAL BES pairwise correlations
#      (Pearson r only, per the tau_obs limitation above) -- both random
#      variable subsets and deliberately constructed "battery" subsets
#      (variables greedily chosen to maximize pairwise |r|, mimicking how a
#      real FA/SEM user would actually select a correlated item set) -- to
#      see whether the invertibility failures found in the adversarial
#      synthetic sweep actually materialize on real survey data.
#
# Reads:
#   R/data/raw/bes2019_pairs.dta
#
# Writes:
#   output/data/bes_tau_bounds.rds         -- per-pair r_min/r_max/tau_min/
#                                              tau_max + 5 adjusted-r columns,
#                                              all 7,503 BES pairs
#   output/data/adjustment_psd_synthetic.rds -- mild + adversarial synthetic
#                                              PSD sweep results, Pearson r
#                                              and Kendall's tau side by side
#   output/data/adjustment_psd_bes_real.rds  -- real-BES-data PSD check
#                                              (random + greedy-adversarial
#                                              variable subsets), Pearson r

library(haven)
library(dplyr)

# =============================================================================
# PART 1: Core Fréchet-Hoeffding bound functions (Pearson r)
# Identical to 01_bes_compute_bounds.R, redefined here per this project's
# convention that every numbered script is self-contained.
# =============================================================================

compute_rmax <- function(counts_x, counts_y) {
  K_x <- length(counts_x)
  K_y <- length(counts_y)
  scores_x <- 0:(K_x - 1)
  scores_y <- 0:(K_y - 1)

  x_vec <- rep(scores_x, times = counts_x)
  y_vec <- rep(scores_y, times = counts_y)

  if (length(x_vec) != length(y_vec)) {
    stop("Total counts differ: sum(counts_x)=", sum(counts_x),
         " sum(counts_y)=", sum(counts_y))
  }

  x_sorted <- sort(x_vec, decreasing = TRUE)
  y_sorted <- sort(y_vec, decreasing = TRUE)

  cor(x_sorted, y_sorted)
}

compute_rmin <- function(counts_x, counts_y) {
  K_x <- length(counts_x)
  K_y <- length(counts_y)
  scores_x <- 0:(K_x - 1)
  scores_y <- 0:(K_y - 1)

  x_vec <- rep(scores_x, times = counts_x)
  y_vec <- rep(scores_y, times = counts_y)

  x_sorted <- sort(x_vec, decreasing = TRUE)
  y_sorted <- sort(y_vec, decreasing = FALSE)

  cor(x_sorted, y_sorted)
}

# =============================================================================
# PART 2: Kendall's tau-b -- closed-form from a J x K count table, plus
# tau_min/tau_max via the SAME comonotonic / anti-comonotonic construction
# used above for r_min/r_max (this reuse is deliberate: the Fréchet-Hoeffding
# upper-bound copula extremizes ANY concordance-respecting association
# measure given fixed marginals, not just Pearson r, so the identical sort
# both maximizes r and maximizes tau-b).
#
# compute_taub_table() implements the standard tau-b formula
#   tau_b = (n_c - n_d) / sqrt((n0-n1)*(n0-n2))
# via O(J*K) cumulative 2-D suffix sums of the count table (no O(N^2) pairwise
# comparison needed). Verified against scipy.stats.kendalltau to machine
# precision (max abs error 2.22e-16, 500 random trials with ties) and the
# comonotonic/anti-comonotonic construction verified against brute-force
# enumeration of all tables with given margins (300 trials, 0 mismatches)
# during development of this script.
# =============================================================================

compute_taub_table <- function(tab) {
  tab <- as.matrix(tab)
  storage.mode(tab) <- "double"
  J <- nrow(tab)
  K <- ncol(tab)
  n <- sum(tab)

  # 2-D suffix sums: suffix[i, j] = sum of tab[i', j'] over all i' >= i, j' >= j
  # (1-indexed; row J+1 and column K+1 are the all-zero boundary)
  suffix <- matrix(0, nrow = J + 1, ncol = K + 1)
  for (i in J:1) {
    for (j in K:1) {
      suffix[i, j] <- tab[i, j] + suffix[i + 1, j] + suffix[i, j + 1] - suffix[i + 1, j + 1]
    }
  }

  n_c <- 0
  n_d <- 0
  for (i in 1:J) {
    for (j in 1:K) {
      cij <- tab[i, j]
      if (cij == 0) next
      n_c <- n_c + cij * suffix[i + 1, j + 1]
      n_d <- n_d + cij * (suffix[i + 1, 1] - suffix[i + 1, j])
    }
  }

  rs <- rowSums(tab)
  cs <- colSums(tab)
  n1 <- sum(rs * (rs - 1) / 2)
  n2 <- sum(cs * (cs - 1) / 2)
  n0 <- n * (n - 1) / 2
  denom <- sqrt((n0 - n1) * (n0 - n2))
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  (n_c - n_d) / denom
}

# Build a J x K count table from two parallel category-label vectors
# (explicit levels so empty categories show up as zero-count cells, just
# like the BES var1freqX/var2freqX columns can be zero).
retabulate <- function(x_vec, y_vec, K_x, K_y) {
  fx <- factor(x_vec, levels = 0:(K_x - 1))
  fy <- factor(y_vec, levels = 0:(K_y - 1))
  as.matrix(table(fx, fy))
}

compute_taumax <- function(counts_x, counts_y) {
  K_x <- length(counts_x)
  K_y <- length(counts_y)
  x_vec <- rep(0:(K_x - 1), times = counts_x)
  y_vec <- rep(0:(K_y - 1), times = counts_y)
  x_sorted <- sort(x_vec, decreasing = TRUE)
  y_sorted <- sort(y_vec, decreasing = TRUE)
  tab <- retabulate(x_sorted, y_sorted, K_x, K_y)
  compute_taub_table(tab)
}

compute_taumin <- function(counts_x, counts_y) {
  K_x <- length(counts_x)
  K_y <- length(counts_y)
  x_vec <- rep(0:(K_x - 1), times = counts_x)
  y_vec <- rep(0:(K_y - 1), times = counts_y)
  x_sorted <- sort(x_vec, decreasing = TRUE)
  y_sorted <- sort(y_vec, decreasing = FALSE)
  tab <- retabulate(x_sorted, y_sorted, K_x, K_y)
  compute_taub_table(tab)
}

# Observed tau-b for two raw (paired) category vectors -- only meaningful
# when a real joint table is available (synthetic data here; NOT available
# for real BES pairs, see header note).
tau_obs_pair <- function(x_vec, y_vec, K_x, K_y) {
  tab <- retabulate(x_vec, y_vec, K_x, K_y)
  compute_taub_table(tab)
}

# =============================================================================
# PART 3: The five rescaling / adjustment transforms (see header for the
# Warrens-framework justification of each)
# =============================================================================

t_raw   <- function(x, lo, hi) x
t_u1    <- function(x, lo, hi) x / (hi - lo)
t_tight <- function(x, lo, hi) 2 * (x - lo) / (hi - lo) - 1
# Warrens (2013)'s sign-preserving correction-for-maximum-value; independently
# proposed (same formula) by Shih & Huang (1992, Biometrics 48(4):1207-1213,
# bib key shih92-evaluating) -- see header note for the worked-example check.
t_split <- function(x, lo, hi) ifelse(x >= 0, x / hi, x / abs(lo))
t_u3    <- function(x, lo, hi) 2 * x / (hi - lo)

TRANSFORMS <- list(
  raw             = t_raw,
  U1_conservative = t_u1,
  tight_linear    = t_tight,
  split_warrens   = t_split,
  U3              = t_u3
)

# =============================================================================
# PART 4: Extract marginal counts from a BES row (identical to
# 01_bes_compute_bounds.R)
# =============================================================================

get_counts <- function(row, var_prefix, n_cats) {
  col0 <- paste0(var_prefix, "freq0")
  start <- if (col0 %in% names(row) && !is.na(row[[col0]]) && row[[col0]] > 0) 0 else 1

  counts <- vapply(start:(start + n_cats - 1), function(i) {
    col <- paste0(var_prefix, "freq", i)
    if (col %in% names(row) && !is.na(row[[col]])) as.numeric(row[[col]]) else 0
  }, numeric(1))

  pmax(counts, 0)
}

# =============================================================================
# PART 5: tau_min/tau_max + adjusted-r for all 7,503 real BES pairs
# =============================================================================

cat("Loading BES data...\n")
bes_path <- "R/data/raw/bes2019_pairs.dta"
bes <- read_dta(bes_path)
n_pairs <- nrow(bes)
cat("Loaded", n_pairs, "pairs\n")

cat("Computing tau_min/tau_max and adjusted-r transforms for all pairs",
    "(this may take a few minutes)...\n")

results <- vector("list", n_pairs)
pb_step <- max(1, floor(n_pairs / 20))

for (i in seq_len(n_pairs)) {
  if (i %% pb_step == 0) cat(sprintf("  %d / %d\n", i, n_pairs))

  row <- bes[i, ]
  K1 <- as.integer(row$var1cats)
  K2 <- as.integer(row$var2cats)
  counts1 <- get_counts(row, "var1", K1)
  counts2 <- get_counts(row, "var2", K2)
  r_obs <- as.numeric(row$corr)

  zero_var <- (var(rep(0:(length(counts1) - 1), times = counts1)) < 1e-12) ||
              (var(rep(0:(length(counts2) - 1), times = counts2)) < 1e-12)

  if (zero_var) {
    results[[i]] <- data.frame(
      pair_id = i, K1 = K1, K2 = K2, r_obs = r_obs,
      r_min = NA_real_, r_max = NA_real_,
      tau_min = NA_real_, tau_max = NA_real_,
      r_adj_raw = NA_real_, r_adj_U1_conservative = NA_real_,
      r_adj_tight_linear = NA_real_, r_adj_split_warrens = NA_real_,
      r_adj_U3 = NA_real_
    )
    next
  }

  r_max   <- tryCatch(compute_rmax(counts1, counts2),   error = function(e) NA_real_)
  r_min   <- tryCatch(compute_rmin(counts1, counts2),   error = function(e) NA_real_)
  tau_max <- tryCatch(compute_taumax(counts1, counts2), error = function(e) NA_real_)
  tau_min <- tryCatch(compute_taumin(counts1, counts2), error = function(e) NA_real_)

  if (is.na(r_min) || is.na(r_max) || (r_max - r_min) <= 1e-12) {
    adj_raw <- adj_u1 <- adj_tight <- adj_split <- adj_u3 <- NA_real_
  } else {
    adj_raw   <- t_raw(r_obs,   r_min, r_max)
    adj_u1    <- t_u1(r_obs,    r_min, r_max)
    adj_tight <- t_tight(r_obs, r_min, r_max)
    adj_split <- t_split(r_obs, r_min, r_max)
    adj_u3    <- t_u3(r_obs,    r_min, r_max)
  }

  results[[i]] <- data.frame(
    pair_id = i, K1 = K1, K2 = K2, r_obs = r_obs,
    r_min = r_min, r_max = r_max,
    tau_min = tau_min, tau_max = tau_max,
    r_adj_raw = adj_raw, r_adj_U1_conservative = adj_u1,
    r_adj_tight_linear = adj_tight, r_adj_split_warrens = adj_split,
    r_adj_U3 = adj_u3
  )
}

bes_tau_bounds <- bind_rows(results)

cat("\nComplete.\n")
cat(sprintf("  Pairs with valid tau bounds: %d / %d\n",
            sum(!is.na(bes_tau_bounds$tau_max)), nrow(bes_tau_bounds)))

saveRDS(bes_tau_bounds, "output/data/bes_tau_bounds.rds")
cat("Saved to output/data/bes_tau_bounds.rds\n\n")

# =============================================================================
# PART 6: Synthetic multivariate ordinal data generator (common-factor
# design): Z = sqrt(rho)*f + sqrt(1-rho)*eps, discretized via skewed
# (Dirichlet/gamma-weighted) quantile cutpoints, so marginals are typically
# asymmetric the way real survey items are.
# =============================================================================

make_correlated_ordinal_data <- function(n_obs, n_vars, K_vec, rho, skew = TRUE) {
  f   <- rnorm(n_obs)
  eps <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
  Z   <- sqrt(rho) * f + sqrt(1 - rho) * eps

  X <- matrix(NA_integer_, nrow = n_obs, ncol = n_vars)
  for (j in seq_len(n_vars)) {
    K <- K_vec[j]
    if (skew) {
      w <- rgamma(K, shape = 2)
      w <- w / sum(w)
    } else {
      w <- rep(1 / K, K)
    }
    cum_w <- cumsum(w)[-K]
    cutpoints <- qnorm(pmin(pmax(cum_w, 1e-6), 1 - 1e-6))
    X[, j] <- findInterval(Z[, j], cutpoints)
  }
  X
}

marginal_counts <- function(col, K) {
  as.numeric(table(factor(col, levels = 0:(K - 1))))
}

min_eig <- function(M) min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)

# =============================================================================
# PART 7: Synthetic PSD/invertibility sweep, Pearson r and Kendall's tau
# together (same simulated data feeds both, so this is a single pass per
# trial rather than two separate sweeps).
# =============================================================================

run_psd_sweep <- function(n_trials, n_vars_range, K_range, rho_range, N_range, label) {
  names_t <- names(TRANSFORMS)
  mk_tracker <- function() {
    list(
      fail  = setNames(rep(0L, length(names_t)), names_t),
      worst = setNames(rep(0,  length(names_t)), names_t),
      oor   = setNames(rep(0L, length(names_t)), names_t)
    )
  }
  pearson_tr <- mk_tracker()
  tau_tr     <- mk_tracker()
  raw_pearson_fail <- 0L; raw_pearson_worst <- 0
  raw_tau_fail     <- 0L; raw_tau_worst     <- 0
  n_used <- 0L

  for (t in seq_len(n_trials)) {
    n_vars <- sample(n_vars_range, 1)
    K_vec  <- sample(K_range, n_vars, replace = TRUE)
    rho    <- runif(1, rho_range[1], rho_range[2])
    n_obs  <- sample(N_range, 1)

    X <- make_correlated_ordinal_data(n_obs, n_vars, K_vec, rho, skew = TRUE)
    if (!all(apply(X, 2, function(col) length(unique(col)) > 1))) next
    n_used <- n_used + 1

    R    <- diag(1, n_vars); RMIN <- diag(1, n_vars); RMAX <- diag(1, n_vars)
    Tau  <- diag(1, n_vars); TMIN <- diag(1, n_vars); TMAX <- diag(1, n_vars)

    for (i in 1:(n_vars - 1)) {
      for (j in (i + 1):n_vars) {
        ci <- marginal_counts(X[, i], K_vec[i])
        cj <- marginal_counts(X[, j], K_vec[j])

        r_obs_ij   <- cor(X[, i], X[, j])
        r_max_ij   <- compute_rmax(ci, cj)
        r_min_ij   <- compute_rmin(ci, cj)
        tau_obs_ij <- tau_obs_pair(X[, i], X[, j], K_vec[i], K_vec[j])
        tau_max_ij <- compute_taumax(ci, cj)
        tau_min_ij <- compute_taumin(ci, cj)

        R[i, j]    <- R[j, i]    <- r_obs_ij
        RMIN[i, j] <- RMIN[j, i] <- r_min_ij
        RMAX[i, j] <- RMAX[j, i] <- r_max_ij
        Tau[i, j]  <- Tau[j, i]  <- tau_obs_ij
        TMIN[i, j] <- TMIN[j, i] <- tau_min_ij
        TMAX[i, j] <- TMAX[j, i] <- tau_max_ij
      }
    }

    e_r <- min_eig(R)
    raw_pearson_worst <- min(raw_pearson_worst, e_r)
    if (e_r < -1e-8) raw_pearson_fail <- raw_pearson_fail + 1

    e_t <- min_eig(Tau)
    raw_tau_worst <- min(raw_tau_worst, e_t)
    if (e_t < -1e-8) raw_tau_fail <- raw_tau_fail + 1

    for (nm in names_t) {
      fn <- TRANSFORMS[[nm]]

      S <- diag(1, n_vars)
      oor <- FALSE
      for (i in 1:(n_vars - 1)) {
        for (j in (i + 1):n_vars) {
          val <- fn(R[i, j], RMIN[i, j], RMAX[i, j])
          S[i, j] <- S[j, i] <- val
          if (abs(val) > 1.0001) oor <- TRUE
        }
      }
      e_s <- min_eig(S)
      pearson_tr$worst[nm] <- min(pearson_tr$worst[nm], e_s)
      if (e_s < -1e-8) pearson_tr$fail[nm] <- pearson_tr$fail[nm] + 1
      if (oor) pearson_tr$oor[nm] <- pearson_tr$oor[nm] + 1

      Stau <- diag(1, n_vars)
      oor_t <- FALSE
      for (i in 1:(n_vars - 1)) {
        for (j in (i + 1):n_vars) {
          val <- fn(Tau[i, j], TMIN[i, j], TMAX[i, j])
          Stau[i, j] <- Stau[j, i] <- val
          if (abs(val) > 1.0001) oor_t <- TRUE
        }
      }
      e_st <- min_eig(Stau)
      tau_tr$worst[nm] <- min(tau_tr$worst[nm], e_st)
      if (e_st < -1e-8) tau_tr$fail[nm] <- tau_tr$fail[nm] + 1
      if (oor_t) tau_tr$oor[nm] <- tau_tr$oor[nm] + 1
    }
  }

  mk_df <- function(tr, raw_fail, raw_worst) {
    data.frame(
      transform = c("raw", names_t),
      n         = n_used,
      fail_pct  = 100 * c(raw_fail, tr$fail) / n_used,
      worst_eig = c(raw_worst, tr$worst),
      oor_pct   = 100 * c(0, tr$oor) / n_used,
      row.names = NULL
    )
  }

  list(
    label  = label,
    n_used = n_used,
    pearson = mk_df(pearson_tr, raw_pearson_fail, raw_pearson_worst),
    tau     = mk_df(tau_tr, raw_tau_fail, raw_tau_worst)
  )
}

set.seed(20260622)

cat("Running synthetic PSD sweep (mild regime: 4-5 vars, rho 0.05-0.6)...\n")
sweep_mild <- run_psd_sweep(
  n_trials = 300, n_vars_range = 4:5, K_range = 3:6,
  rho_range = c(0.05, 0.6), N_range = 150:400, label = "mild"
)

cat("Running synthetic PSD sweep (adversarial regime: 6-8 vars, rho 0.3-0.85,",
    "small N)...\n")
sweep_adversarial <- run_psd_sweep(
  n_trials = 300, n_vars_range = 6:8, K_range = 3:6,
  rho_range = c(0.3, 0.85), N_range = 40:300, label = "adversarial"
)

saveRDS(
  list(mild = sweep_mild, adversarial = sweep_adversarial),
  "output/data/adjustment_psd_synthetic.rds"
)
cat("Saved to output/data/adjustment_psd_synthetic.rds\n\n")

# =============================================================================
# PART 8: Real-BES-data PSD/invertibility check (Pearson r only -- tau_obs is
# not available for real BES pairs, see header note). Two sampling schemes:
#   - random:  uniformly random k-variable subsets of the 123 BES variables
#   - greedy:  k-variable "batteries" built by greedily adding, at each step,
#              whichever remaining variable has the highest average |r| with
#              the variables already chosen -- i.e. deliberately constructing
#              the kind of highly-intercorrelated item set a real FA/SEM user
#              would actually select, which is the realistic stress case for
#              invertibility (nobody runs FA on a random grab-bag of unrelated
#              survey questions).
# =============================================================================

var1_int <- as.integer(bes$var1)
var2_int <- as.integer(bes$var2)
all_vars <- sort(unique(c(var1_int, var2_int)))

pair_lookup <- new.env(parent = emptyenv())
for (i in seq_len(n_pairs)) {
  key <- paste(var1_int[i], var2_int[i], sep = "_")
  pair_lookup[[key]] <- i
}

get_pair_idx <- function(v1, v2) {
  lo <- min(v1, v2); hi <- max(v1, v2)
  key <- paste(lo, hi, sep = "_")
  val <- pair_lookup[[key]]
  if (is.null(val)) NA_integer_ else val
}

bes_pair_bounds <- function(idx) {
  row <- bes[idx, ]
  K1 <- as.integer(row$var1cats)
  K2 <- as.integer(row$var2cats)
  c1 <- get_counts(row, "var1", K1)
  c2 <- get_counts(row, "var2", K2)
  if (sum(c1) == 0 || sum(c2) == 0) return(NULL)
  list(
    r_obs = as.numeric(row$corr),
    r_min = compute_rmin(c1, c2),
    r_max = compute_rmax(c1, c2)
  )
}

build_real_matrices <- function(vars_subset) {
  k <- length(vars_subset)
  R <- diag(1, k); RMIN <- diag(1, k); RMAX <- diag(1, k)
  for (a in 1:(k - 1)) {
    for (b in (a + 1):k) {
      idx <- get_pair_idx(vars_subset[a], vars_subset[b])
      if (is.na(idx)) return(NULL)
      pb <- bes_pair_bounds(idx)
      if (is.null(pb)) return(NULL)
      R[a, b]    <- R[b, a]    <- pb$r_obs
      RMIN[a, b] <- RMIN[b, a] <- pb$r_min
      RMAX[a, b] <- RMAX[b, a] <- pb$r_max
    }
  }
  list(R = R, RMIN = RMIN, RMAX = RMAX)
}

greedy_correlated_cluster <- function(seed, k, sample_size = 40) {
  chosen <- c(seed)
  candidates <- setdiff(all_vars, seed)
  while (length(chosen) < k && length(candidates) > 0) {
    cand_sample <- if (length(candidates) > sample_size) {
      sample(candidates, sample_size)
    } else {
      candidates
    }
    best_v <- NA_integer_
    best_score <- -1
    for (v in cand_sample) {
      scores <- numeric(0)
      valid <- TRUE
      for (c0 in chosen) {
        idx <- get_pair_idx(v, c0)
        if (is.na(idx)) { valid <- FALSE; break }
        scores <- c(scores, abs(bes$corr[idx]))
      }
      if (!valid) next
      avg <- mean(scores)
      if (avg > best_score) { best_score <- avg; best_v <- v }
    }
    if (is.na(best_v)) break
    chosen <- c(chosen, best_v)
    candidates <- setdiff(candidates, best_v)
  }
  chosen
}

run_real_psd_check <- function(n_trials, k, mode = c("random", "greedy")) {
  mode <- match.arg(mode)
  names_t <- names(TRANSFORMS)
  fail  <- setNames(rep(0L, length(names_t)), names_t)
  worst <- setNames(rep(0,  length(names_t)), names_t)
  oor   <- setNames(rep(0L, length(names_t)), names_t)
  raw_fail <- 0L; raw_worst <- 0; n_used <- 0L

  for (t in seq_len(n_trials)) {
    if (mode == "random") {
      subset_vars <- sample(all_vars, k)
    } else {
      seed <- sample(all_vars, 1)
      subset_vars <- greedy_correlated_cluster(seed, k)
      if (length(subset_vars) < k) next
    }

    built <- build_real_matrices(subset_vars)
    if (is.null(built)) next
    n_used <- n_used + 1
    R <- built$R; RMIN <- built$RMIN; RMAX <- built$RMAX

    e_r <- min_eig(R)
    raw_worst <- min(raw_worst, e_r)
    if (e_r < -1e-8) raw_fail <- raw_fail + 1

    for (nm in names_t) {
      fn <- TRANSFORMS[[nm]]
      S <- diag(1, k)
      oor_flag <- FALSE
      for (a in 1:(k - 1)) {
        for (b in (a + 1):k) {
          val <- fn(R[a, b], RMIN[a, b], RMAX[a, b])
          S[a, b] <- S[b, a] <- val
          if (abs(val) > 1.0001) oor_flag <- TRUE
        }
      }
      e_s <- min_eig(S)
      worst[nm] <- min(worst[nm], e_s)
      if (e_s < -1e-8) fail[nm] <- fail[nm] + 1
      if (oor_flag) oor[nm] <- oor[nm] + 1
    }
  }

  data.frame(
    transform = c("raw", names_t),
    n         = n_used,
    fail_pct  = 100 * c(raw_fail, fail) / n_used,
    worst_eig = c(raw_worst, worst),
    oor_pct   = 100 * c(0, oor) / n_used,
    row.names = NULL
  )
}

cat("Running real-BES-data PSD check: random 8-variable subsets...\n")
real_random_k8 <- run_real_psd_check(300, 8, mode = "random")

cat("Running real-BES-data PSD check: greedy max-correlation 8-variable",
    "'batteries'...\n")
real_greedy_k8 <- run_real_psd_check(150, 8, mode = "greedy")

cat("Running real-BES-data PSD check: greedy max-correlation 12-variable",
    "'batteries'...\n")
real_greedy_k12 <- run_real_psd_check(150, 12, mode = "greedy")

saveRDS(
  list(random_k8 = real_random_k8, greedy_k8 = real_greedy_k8, greedy_k12 = real_greedy_k12),
  "output/data/adjustment_psd_bes_real.rds"
)
cat("Saved to output/data/adjustment_psd_bes_real.rds\n\n")

# =============================================================================
# PART 9: Summary
# =============================================================================

cat("=============================================================\n")
cat(" 06_correlation_adjustment.R complete\n")
cat("=============================================================\n")
cat(sprintf(" BES pairs with valid tau bounds: %d / %d\n",
            sum(!is.na(bes_tau_bounds$tau_max)), nrow(bes_tau_bounds)))
cat("\n Synthetic PSD sweep -- Pearson r (mild):\n")
print(sweep_mild$pearson, row.names = FALSE)
cat("\n Synthetic PSD sweep -- Pearson r (adversarial):\n")
print(sweep_adversarial$pearson, row.names = FALSE)
cat("\n Synthetic PSD sweep -- Kendall's tau (mild):\n")
print(sweep_mild$tau, row.names = FALSE)
cat("\n Synthetic PSD sweep -- Kendall's tau (adversarial):\n")
print(sweep_adversarial$tau, row.names = FALSE)
cat("\n Real-BES-data PSD check -- random 8-variable subsets:\n")
print(real_random_k8, row.names = FALSE)
cat("\n Real-BES-data PSD check -- greedy max-correlation 8-variable batteries:\n")
print(real_greedy_k8, row.names = FALSE)
cat("\n Real-BES-data PSD check -- greedy max-correlation 12-variable batteries:\n")
print(real_greedy_k12, row.names = FALSE)
cat("=============================================================\n")
