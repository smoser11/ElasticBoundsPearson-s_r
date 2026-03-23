# 02_mc_simulation.R
#
# Monte Carlo simulation of correlation bounds across orcat pairs.
# Generates random marginal distributions for K = 4, 5, 6, 7, 10, 11
# and computes rmin, rmax, constraint metrics for each pair.
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
# Draw random marginal distribution as integer counts
# N = effective sample size to expand to
# ---------------------------------------------------------------------------
random_marginal <- function(K, N = 1000) {
  # Draw from Dirichlet(1,...,1) = uniform on simplex
  # via Gamma(1,1) trick
  raw <- rgamma(K, shape = 1, rate = 1)
  p <- raw / sum(raw)
  counts <- as.integer(round(p * N))
  # Ensure sum matches N
  diff <- N - sum(counts)
  if (diff != 0) counts[which.max(counts)] <- counts[which.max(counts)] + diff
  pmax(counts, 0L)
}

# ---------------------------------------------------------------------------
# Main simulation
# ---------------------------------------------------------------------------
set.seed(2025)

K_values <- c(4, 5, 6, 7, 10, 11)
N_per_K <- 2000    # number of pairs per K combination (K1, K2)
N_obs    <- 1000   # sample size for expanding marginals

cat("Running MC simulation...\n")
cat("K values:", paste(K_values, collapse=", "), "\n")
cat("Pairs per K combination:", N_per_K, "\n\n")

all_results <- list()
chunk_id <- 1

for (K1 in K_values) {
  for (K2 in K_values) {
    cat(sprintf("  K1=%d, K2=%d ...", K1, K2))
    chunk <- vector("list", N_per_K)

    for (sim in seq_len(N_per_K)) {
      counts1 <- random_marginal(K1, N_obs)
      counts2 <- random_marginal(K2, N_obs)

      # Skip zero-variance cases
      if (var(rep(0:(K1-1), times=counts1)) < 1e-12 ||
          var(rep(0:(K2-1), times=counts2)) < 1e-12) {
        chunk[[sim]] <- NULL
        next
      }

      r_max <- tryCatch(compute_rmax(counts1, counts2), error = function(e) NA)
      r_min <- tryCatch(compute_rmin(counts1, counts2), error = function(e) NA)

      chunk[[sim]] <- list(
        K1 = K1, K2 = K2,
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
