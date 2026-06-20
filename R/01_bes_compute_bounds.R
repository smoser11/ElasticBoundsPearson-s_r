# 01_bes_compute_bounds.R
#
# Compute rmin, rmax, constraint metrics, and asymmetry measures
# for all 7,503 BES 2019 variable pairs.
#
# Outputs:
#   output/data/bes_bounds.rds  — data frame with all computed quantities

library(haven)
library(dplyr)

# ---------------------------------------------------------------------------
# Core FH bounds functions (exact, works with integer counts)
# ---------------------------------------------------------------------------

# Compute rmax via comonotonic coupling (sort both ascending, compute r)
compute_rmax <- function(counts_x, counts_y) {
  K_x <- length(counts_x)
  K_y <- length(counts_y)
  scores_x <- 0:(K_x - 1)
  scores_y <- 0:(K_y - 1)

  # Expand to observation vectors
  x_vec <- rep(scores_x, times = counts_x)
  y_vec <- rep(scores_y, times = counts_y)

  if (length(x_vec) != length(y_vec)) {
    stop("Total counts differ: sum(counts_x)=", sum(counts_x),
         " sum(counts_y)=", sum(counts_y))
  }

  # Comonotonic: both sorted descending
  x_sorted <- sort(x_vec, decreasing = TRUE)
  y_sorted <- sort(y_vec, decreasing = TRUE)

  cor(x_sorted, y_sorted)
}

# Compute rmin via anti-comonotonic coupling (one asc, one desc)
compute_rmin <- function(counts_x, counts_y) {
  K_x <- length(counts_x)
  K_y <- length(counts_y)
  scores_x <- 0:(K_x - 1)
  scores_y <- 0:(K_y - 1)

  x_vec <- rep(scores_x, times = counts_x)
  y_vec <- rep(scores_y, times = counts_y)

  # Anti-comonotonic: x descending, y ascending
  x_sorted <- sort(x_vec, decreasing = TRUE)
  y_sorted <- sort(y_vec, decreasing = FALSE)

  cor(x_sorted, y_sorted)
}

# Shannon entropy: -sum(p * log(p)), 0*log(0) = 0
shannon_entropy <- function(counts) {
  p <- counts / sum(counts)
  p <- p[p > 0]
  -sum(p * log(p))
}

# Total variation distance between p and its reversal
tv_distance <- function(counts) {
  p <- counts / sum(counts)
  p_rev <- rev(p)
  0.5 * sum(abs(p - p_rev))
}

# Bhattacharyya coefficient between p and its reversal
bhattacharyya_coef <- function(counts) {
  p <- counts / sum(counts)
  p_rev <- rev(p)
  sum(sqrt(p * p_rev))
}

# Overlap coefficient between p and its reversal
overlap_coef <- function(counts) {
  p <- counts / sum(counts)
  p_rev <- rev(p)
  sum(pmin(p, p_rev))
}

# ---------------------------------------------------------------------------
# Extract marginal counts from a BES row
# ---------------------------------------------------------------------------
get_counts <- function(row, var_prefix, n_cats) {
  # Determine start index (0 or 1)
  col0 <- paste0(var_prefix, "freq0")
  start <- if (col0 %in% names(row) && !is.na(row[[col0]]) && row[[col0]] > 0) 0 else 1

  counts <- vapply(start:(start + n_cats - 1), function(i) {
    col <- paste0(var_prefix, "freq", i)
    if (col %in% names(row) && !is.na(row[[col]])) as.numeric(row[[col]]) else 0
  }, numeric(1))

  pmax(counts, 0)
}

# ---------------------------------------------------------------------------
# Load BES data
# ---------------------------------------------------------------------------
cat("Loading BES data...\n")
bes_path <- "R/data/raw/bes2019_pairs.dta"
bes <- read_dta(bes_path)
cat("Loaded", nrow(bes), "pairs\n")

# ---------------------------------------------------------------------------
# Compute bounds for all pairs
# ---------------------------------------------------------------------------
cat("Computing bounds (this may take a minute)...\n")

n_pairs <- nrow(bes)
results <- vector("list", n_pairs)

pb_step <- floor(n_pairs / 20)

for (i in seq_len(n_pairs)) {
  if (i %% pb_step == 0) cat(sprintf("  %d / %d\n", i, n_pairs))

  row <- bes[i, ]

  counts1 <- get_counts(row, "var1", as.integer(row$var1cats))
  counts2 <- get_counts(row, "var2", as.integer(row$var2cats))

  # Skip if zero variance in either variable
  if (var(rep(0:(length(counts1)-1), times=counts1)) < 1e-12 ||
      var(rep(0:(length(counts2)-1), times=counts2)) < 1e-12) {
    results[[i]] <- list(
      r_min = NA, r_max = NA, r_obs = row$corr,
      C1 = NA, C2 = NA,
      H1 = NA, H2 = NA, H_mean = NA,
      TV1 = NA, TV2 = NA, TV_mean = NA,
      BC1 = NA, BC2 = NA,
      OVL1 = NA, OVL2 = NA,
      asymmetry = NA,
      K1 = row$var1cats, K2 = row$var2cats
    )
    next
  }

  r_max <- tryCatch(compute_rmax(counts1, counts2), error = function(e) NA)
  r_min <- tryCatch(compute_rmin(counts1, counts2), error = function(e) NA)

  C1 <- if (!is.na(r_max) && !is.na(r_min)) (abs(r_min) + r_max) / 2 else NA
  C2 <- if (!is.na(r_max) && !is.na(r_min)) (r_min^2 + r_max^2) / 2 else NA

  H1  <- shannon_entropy(counts1)
  H2  <- shannon_entropy(counts2)
  TV1 <- tv_distance(counts1)
  TV2 <- tv_distance(counts2)
  BC1 <- bhattacharyya_coef(counts1)
  BC2 <- bhattacharyya_coef(counts2)
  OVL1 <- overlap_coef(counts1)
  OVL2 <- overlap_coef(counts2)

  # Bound asymmetry: |r_min| - r_max (positive means |r_min| > r_max)
  asym <- if (!is.na(r_min) && !is.na(r_max)) abs(r_min) - r_max else NA

  results[[i]] <- list(
    r_min = r_min, r_max = r_max, r_obs = as.numeric(row$corr),
    C1 = C1, C2 = C2,
    H1 = H1, H2 = H2, H_mean = (H1 + H2) / 2,
    TV1 = TV1, TV2 = TV2, TV_mean = (TV1 + TV2) / 2,
    BC1 = BC1, BC2 = BC2,
    OVL1 = OVL1, OVL2 = OVL2,
    asymmetry = asym,
    K1 = as.integer(row$var1cats), K2 = as.integer(row$var2cats),
    skew1 = as.numeric(row$var1skewness),
    skew2 = as.numeric(row$var2skewness)
  )
}

# ---------------------------------------------------------------------------
# Combine into data frame
# ---------------------------------------------------------------------------
bounds_df <- bind_rows(lapply(results, as.data.frame))
bounds_df$pair_id <- seq_len(nrow(bounds_df))
bounds_df$K_max <- pmax(bounds_df$K1, bounds_df$K2)
bounds_df$K_min <- pmin(bounds_df$K1, bounds_df$K2)
bounds_df$K_label <- paste0(bounds_df$K1, "×", bounds_df$K2)

cat("\nComplete. Summary:\n")
cat(sprintf("  Pairs with valid bounds: %d / %d\n",
            sum(!is.na(bounds_df$r_max)), nrow(bounds_df)))
cat(sprintf("  C1 range: [%.3f, %.3f]\n",
            min(bounds_df$C1, na.rm=TRUE), max(bounds_df$C1, na.rm=TRUE)))
cat(sprintf("  Asymmetric (|r_min| > r_max): %d pairs\n",
            sum(bounds_df$asymmetry > 0, na.rm=TRUE)))

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
saveRDS(bounds_df, "output/data/bes_bounds.rds")
cat("Saved to output/data/bes_bounds.rds\n")
