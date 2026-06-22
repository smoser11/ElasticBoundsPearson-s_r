# 02_mc_simulation.R
#
# Monte Carlo simulation of correlation bounds across orcat pairs.
# For each (K1, K2) category-count combination AND each sample size N_obs,
# draws cross-tab tables UNIFORMLY AT RANDOM from the space of all K1 x K2
# non-negative-integer tables whose cells sum to N_obs (stars-and-bars
# method -- see R/02a_mc_sampling_uniform_crosstab.R for the full derivation
# and a sanity check against the closed-form marginal distribution).
#
# Note this means marginals are no longer drawn independently: a table's
# row/column sums are simply whatever falls out of that uniform draw over
# the joint table space. For each draw we record:
#   - r_obs: the empirical Pearson's r of the actual simulated joint table
#   - r_min, r_max, and the usual constraint/marginal-shape metrics,
#     computed from that table's REALIZED row/column sums
#
# Output: output/data/mc_bounds.rds

library(dplyr)

# ---------------------------------------------------------------------------
# Reuse bound functions (inline to keep script self-contained)
# ---------------------------------------------------------------------------

compute_rmax <- function(counts_x, counts_y) {
  K_x <- length(counts_x); K_y <- length(counts_y)
  scores_x <- 0:(K_x - 1); scores_y <- 0:(K_y - 1)
  x_vec <- rep(scores_x, times = counts_x)
  y_vec <- rep(scores_y, times = counts_y)
  cor(sort(x_vec, decreasing = TRUE), sort(y_vec, decreasing = TRUE))
}

compute_rmin <- function(counts_x, counts_y) {
  K_x <- length(counts_x); K_y <- length(counts_y)
  scores_x <- 0:(K_x - 1); scores_y <- 0:(K_y - 1)
  x_vec <- rep(scores_x, times = counts_x)
  y_vec <- rep(scores_y, times = counts_y)
  cor(sort(x_vec, decreasing = TRUE), sort(y_vec, decreasing = FALSE))
}

shannon_entropy <- function(counts) {
  p <- counts / sum(counts); p <- p[p > 0]
  -sum(p * log(p))
}

tv_distance <- function(counts) {
  p <- counts / sum(counts)
  0.5 * sum(abs(p - rev(p)))
}

bhattacharyya_coef <- function(counts) {
  p <- counts / sum(counts)
  sum(sqrt(p * rev(p)))
}

overlap_coef <- function(counts) {
  p <- counts / sum(counts)
  sum(pmin(p, rev(p)))
}

# ---------------------------------------------------------------------------
# Uniform sampling of a J x K cross-tab with fixed grand total N (row/column
# margins are NOT fixed -- only the grand total is). Stars-and-bars method:
# lay out N "stars" and m-1 "bars" (m = J*K) in N+m-1 slots; choosing which
# m-1 slots are bars, uniformly at random, induces a uniform distribution
# over all compositions of N into m parts (cell counts = gaps between
# consecutive bars). See R/02a_mc_sampling_uniform_crosstab.R for the
# derivation and a closed-form sanity check against empirical frequencies.
# ---------------------------------------------------------------------------
sample_uniform_crosstab <- function(N, J, K) {
  m <- J * K
  total_slots <- N + m - 1
  bars <- sort(sample.int(total_slots, m - 1))
  parts <- diff(c(0L, bars, total_slots + 1L)) - 1L
  matrix(parts, nrow = J, ncol = K)
}

# ---------------------------------------------------------------------------
# Empirical Pearson's r of a sampled J x K table, scoring rows/columns
# 0:(J-1) / 0:(K-1) -- the same convention compute_rmax/compute_rmin use.
# Computed from the table's cell-count moments directly (closed form),
# rather than expanding all N_obs observations into score vectors, since
# this runs inside the innermost loop below.
# ---------------------------------------------------------------------------
compute_robs <- function(tab) {
  J <- nrow(tab); K <- ncol(tab); N <- sum(tab)
  row_scores <- 0:(J - 1); col_scores <- 0:(K - 1)
  rs <- rowSums(tab); cs <- colSums(tab)
  mx <- sum(row_scores * rs) / N
  my <- sum(col_scores * cs) / N
  exy <- sum(tab * outer(row_scores, col_scores)) / N
  vx <- sum(row_scores^2 * rs) / N - mx^2
  vy <- sum(col_scores^2 * cs) / N - my^2
  (exy - mx * my) / sqrt(vx * vy)
}

# ---------------------------------------------------------------------------
# Main simulation
# ---------------------------------------------------------------------------
set.seed(2025)

K_values     <- c(3, 4, 5, 6, 7, 9, 10, 11)
N_obs_values <- c(20, 50, 100, 200, 500, 1000, 2000, 5000)  # sample-size sweep
N_per_K      <- 500   # number of tables per (K1, K2, N_obs) combination.
                       # NOTE: was 2000 before the N_obs sweep was added below;
                       # scaled down here because adding 8 N_obs values
                       # multiplies the total table count by 8x. With
                       # N_per_K=500: 8 x 8 x 8 x 500 = 256,000 tables total
                       # (vs. 128,000 before, at the old N_per_K=2000 with no
                       # N_obs dimension). Raise this back toward 2000 for
                       # denser per-cell coverage if runtime allows.

cat("Running MC simulation...\n")
cat("K values:", paste(K_values, collapse=", "), "\n")
cat("N_obs values:", paste(N_obs_values, collapse=", "), "\n")
cat("Tables per (K1, K2, N_obs) combination:", N_per_K, "\n\n")

all_results <- list()
chunk_id <- 1

for (K1 in K_values) {
  for (K2 in K_values) {
    for (N_obs in N_obs_values) {
      cat(sprintf("  K1=%d, K2=%d, N_obs=%d ...", K1, K2, N_obs))
      chunk <- vector("list", N_per_K)

      for (sim in seq_len(N_per_K)) {
        tab <- sample_uniform_crosstab(N_obs, K1, K2)
        counts1 <- rowSums(tab)   # realized row marginal (length K1)
        counts2 <- colSums(tab)   # realized column marginal (length K2)

        # Skip zero-variance cases (all mass landed in a single row/column)
        if (var(rep(0:(K1-1), times=counts1)) < 1e-12 ||
            var(rep(0:(K2-1), times=counts2)) < 1e-12) {
          chunk[[sim]] <- NULL
          next
        }

        r_obs <- tryCatch(compute_robs(tab), error = function(e) NA)
        r_max <- tryCatch(compute_rmax(counts1, counts2), error = function(e) NA)
        r_min <- tryCatch(compute_rmin(counts1, counts2), error = function(e) NA)

        chunk[[sim]] <- list(
          K1 = K1, K2 = K2, N_obs = N_obs,
          r_obs = r_obs,
          r_min = r_min, r_max = r_max,
          C1 = if (!is.na(r_max) && !is.na(r_min)) (abs(r_min) + r_max) / 2 else NA,
          C2 = if (!is.na(r_max) && !is.na(r_min)) (r_min^2 + r_max^2) / 2 else NA,
          H1  = shannon_entropy(counts1),
          H2  = shannon_entropy(counts2),
          TV1 = tv_distance(counts1),
          TV2 = tv_distance(counts2),
          BC1 = bhattacharyya_coef(counts1),
          BC2 = bhattacharyya_coef(counts2),
          OVL1 = overlap_coef(counts1),
          OVL2 = overlap_coef(counts2)
        )
      }

      chunk_df <- bind_rows(lapply(Filter(Negate(is.null), chunk), as.data.frame))
      all_results[[chunk_id]] <- chunk_df
      chunk_id <- chunk_id + 1
      cat(sprintf(" %d rows\n", nrow(chunk_df)))
    }
  }
}

mc <- bind_rows(all_results)
mc$H_mean  <- (mc$H1 + mc$H2) / 2
mc$TV_mean <- (mc$TV1 + mc$TV2) / 2
mc$asymmetry <- abs(mc$r_min) - mc$r_max
mc$K_label <- paste0(mc$K1, "×", mc$K2)
mc$K_max <- pmax(mc$K1, mc$K2)
mc$K_min <- pmin(mc$K1, mc$K2)

cat(sprintf("\nTotal MC pairs: %d\n", nrow(mc)))
cat(sprintf("C1 range: [%.3f, %.3f]\n", min(mc$C1, na.rm=TRUE), max(mc$C1, na.rm=TRUE)))
cat(sprintf("Asymmetric (|r_min| > r_max): %d (%.1f%%)\n",
            sum(mc$asymmetry > 0, na.rm=TRUE),
            100 * mean(mc$asymmetry > 0, na.rm=TRUE)))

saveRDS(mc, "output/data/mc_bounds.rds")
cat("Saved to output/data/mc_bounds.rds\n")
