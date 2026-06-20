# =============================================================================
# Uniform sampling of J x K cross-tab tables, CONDITIONED on Pearson's r
# (total N fixed, row/column margins free) -- builds on sample_uniform_crosstab.R
# =============================================================================
#
# WHY THIS NEEDS MORE THAN REJECTION SAMPLING:
# r is continuous-valued, so the probability that a uniformly drawn table hits
# an exact target r0 is ~0. We instead target the conditional UNIFORM
# distribution over the level set { tables T : |r(T) - r0| <= eps }, for a
# tolerance eps you choose.
#
# METHOD:
#   Phase 1 (search): simulated annealing using single-unit "swap" moves
#     (move 1 count from one cell to another) to steer an arbitrary starting
#     table toward r0. This is a search heuristic only -- it does not need to
#     be exactly correct, just to get the chain into the target band.
#   Phase 2 (exact uniform sampling): once inside the band, run a Metropolis
#     chain using the SAME swap moves, accepting a proposal iff it stays in
#     the band. Because the proposal (pick two distinct cells uniformly,
#     move one unit between them) has selection probability that does NOT
#     depend on the current table, it is exactly symmetric -- so this MH
#     ratio is always 1 when feasible, giving genuinely uniform draws over
#     the band (not just "tables with roughly the right r").
# =============================================================================

# ---- Pearson's r for a table, treating rows/cols as ordinal -----------------
# Default scores are 1:J and 1:K. Supply custom scores if your categories
# aren't naturally 1..J / 1..K spaced.
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

# ---- single swap proposal (symmetric; tab stored as a length-m vector) -----
.propose_swap <- function(tab, m) {
  idx <- sample.int(m, 2)
  from <- idx[1]; to <- idx[2]
  if (tab[from] == 0L) return(NULL)   # infeasible (can't remove from an empty cell)
  tab[from] <- tab[from] - 1L
  tab[to]   <- tab[to] + 1L
  tab
}

# ---- Phase 1: simulated annealing to steer toward r0 -----------------------
.anneal_to_target <- function(tab, r0, J, K, row_scores, col_scores,
                               n_iter = 20000, temp_start = 1, temp_end = 0.01) {
  m <- length(tab)
  r_cur <- table_pearson_r(matrix(tab, J, K), row_scores, col_scores)
  for (it in seq_len(n_iter)) {
    temp <- temp_start * (temp_end / temp_start)^(it / n_iter)
    prop <- .propose_swap(tab, m)
    if (is.null(prop)) next
    r_prop <- table_pearson_r(matrix(prop, J, K), row_scores, col_scores)
    d_cur  <- (r_cur  - r0)^2
    d_prop <- (r_prop - r0)^2
    if (d_prop <= d_cur || runif(1) < exp(-(d_prop - d_cur) / (2 * temp^2))) {
      tab <- prop; r_cur <- r_prop
    }
  }
  list(tab = tab, r = r_cur)
}

# ---- Main function: n_draws uniform tables with r in [r0-eps, r0+eps] ------
sample_uniform_crosstab_given_r <- function(N, J, K, r0, eps,
                                             row_scores = NULL, col_scores = NULL,
                                             n_draws = 1000, thin = 50,
                                             anneal_iter = 20000, max_anneal_tries = 5,
                                             init_tab = NULL) {
  m <- J * K
  if (is.null(row_scores)) row_scores <- seq_len(J)
  if (is.null(col_scores)) col_scores <- seq_len(K)
  in_band <- function(r) abs(r - r0) <= eps

  tab <- if (!is.null(init_tab)) as.integer(init_tab) else
    as.integer(sample_uniform_crosstab(N, J, K))  # from sample_uniform_crosstab.R
  r_cur <- table_pearson_r(matrix(tab, J, K), row_scores, col_scores)

  tries <- 0
  while (!in_band(r_cur) && tries < max_anneal_tries) {
    res <- .anneal_to_target(tab, r0, J, K, row_scores, col_scores, n_iter = anneal_iter)
    tab <- res$tab; r_cur <- res$r
    tries <- tries + 1
  }
  if (!in_band(r_cur)) {
    stop(sprintf(
      "Could not reach r in [%.4f, %.4f] after %d annealing attempt(s); reached r = %.4f.\n  Try: wider eps, larger anneal_iter, or check r0 is attainable for this N, J, K.",
      r0 - eps, r0 + eps, max_anneal_tries, r_cur))
  }

  out <- vector("list", n_draws)
  for (d in seq_len(n_draws)) {
    for (s in seq_len(thin)) {
      prop <- .propose_swap(tab, m)
      if (is.null(prop)) next
      r_prop <- table_pearson_r(matrix(prop, J, K), row_scores, col_scores)
      if (in_band(r_prop)) tab <- prop   # MH ratio = 1 whenever feasible & in-band
    }
    out[[d]] <- matrix(tab, nrow = J, ncol = K)
  }
  out
}

# =============================================================================
# Demo / sanity check
# =============================================================================
# (Requires sample_uniform_crosstab() from sample_uniform_crosstab.R to be
#  loaded first, e.g. source("sample_uniform_crosstab.R"))
# source("sample_uniform_crosstab.R"))

source("sample_uniform_crosstab.R")
set.seed(1)
draws <- sample_uniform_crosstab_given_r(
  N = 60, J = 4, K = 4, r0 = 0.6, eps = 0.01,
  n_draws = 300, thin = 60
)
rs <- sapply(draws, table_pearson_r)
cat("range of r across draws:", range(rs), "\n")
cat("all within band:", all(abs(rs - 0.6) <= 0.01), "\n")
cat("unique tables among draws:", length(unique(draws)), "/", length(draws), "\n")

