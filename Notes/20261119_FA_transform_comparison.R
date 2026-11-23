# 20261119_FA_transform_comparison.R
#
# Compares every rescaling of Pearson's r in R/06_correlation_adjustment.R, plus
# polychoric correlations as a benchmark, in the factor-analysis simulation design of
# R/10_fa_ceiling_diagnostic_pilot.R.
#
# Methods (lo = r_min, hi = r_max for the pair, computed from the two marginals):
#   raw         r
#   split       r / hi if r >= 0, r / |lo| if r < 0   (Cureton's phi/phi_max; Shih & Huang;
#                                                      Shojima's rho*; "split_warrens" in R/06)
#   U3          2r / (hi - lo) = r / C1               (stretch by the share of [-1, 1] covered)
#   U1          r / (hi - lo)                         (values lie in an interval of width 1)
#   tight       2(r - lo)/(hi - lo) - 1               (maps [lo, hi] linearly onto [-1, 1])
#   polychoric  psych::polychoric (tetrachoric for binary items), no smoothing
#
# For each simulated data set we record:
#   - number of factors retained by Horn's parallel analysis (95th percentile of
#     random-data eigenvalues, n = 1000, p = 10), as in R/10;
#   - whether the matrix is a valid correlation matrix (positive semi-definite, PSD);
#   - RMSE of the off-diagonal entries against the TRUE latent correlations
#     (lambda^2 within a factor, lambda^2 * phi across factors).
#
# Run from the repository root:
#   Rscript Notes/20261119_FA_transform_comparison.R > Notes/20261119_FA_transform_comparison_output.txt 2>&1

suppressMessages({library(data.table); library(psych)})
set.seed(20261119)

rmax_from_props <- function(sx, px, sy, py) {
  cx <- cumsum(px); cy <- cumsum(py); cx[length(cx)] <- 1; cy[length(cy)] <- 1
  br <- sort(unique(c(0, cx, cy))); mid <- (head(br, -1) + tail(br, -1)) / 2; w <- diff(br)
  qx <- sx[pmin(findInterval(mid, cx) + 1, length(sx))]; qy <- sy[pmin(findInterval(mid, cy) + 1, length(sy))]
  mx <- sum(px * sx); my <- sum(py * sy)
  (sum(w * qx * qy) - mx * my) / (sqrt(sum(px * sx^2) - mx^2) * sqrt(sum(py * sy^2) - my^2))
}
bound_matrices <- function(X) {
  p <- ncol(X); Mx <- diag(p); Mn <- -diag(p)
  pr <- lapply(seq_len(p), function(j) { t <- table(X[, j]); list(s = as.numeric(names(t)), p = as.numeric(t) / sum(t)) })
  for (i in 1:(p - 1)) for (j in (i + 1):p) {
    Mx[i, j] <- Mx[j, i] <- rmax_from_props(pr[[i]]$s, pr[[i]]$p, pr[[j]]$s, pr[[j]]$p)
    Mn[i, j] <- Mn[j, i] <- -rmax_from_props(pr[[i]]$s, pr[[i]]$p, -rev(pr[[j]]$s), rev(pr[[j]]$p))
  }
  list(hi = Mx, lo = Mn)
}
gen_data <- function(n, p, lambda, K, shifts, nfac = 1, phi = 0.3) {
  if (nfac == 1) {
    f <- rnorm(n); Z <- sapply(seq_len(p), function(j) lambda * f + sqrt(1 - lambda^2) * rnorm(n))
  } else {
    f <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, phi, phi, 1), 2)); g <- rep(1:2, each = p / 2)
    Z <- sapply(seq_len(p), function(j) lambda * f[, g[j]] + sqrt(1 - lambda^2) * rnorm(n))
  }
  base <- if (K == 2) 0 else seq(-1.5, 1.5, length.out = K - 1)
  sapply(seq_len(p), function(j) findInterval(Z[, j], base + shifts[j]) + 1)
}
true_cor <- function(p, lambda, nfac, phi = 0.3) {
  if (nfac == 1) { S <- matrix(lambda^2, p, p) } else {
    g <- rep(1:2, each = p / 2); S <- outer(g, g, function(a, b) ifelse(a == b, lambda^2, lambda^2 * phi))
  }
  diag(S) <- 1; S
}
pa_threshold <- function(n, p, reps = 200)
  apply(replicate(reps, eigen(cor(matrix(rnorm(n * p), n, p)), only.values = TRUE)$values), 1, quantile, 0.95)
pa_count <- function(e, thr) { k <- 0; while (k < length(e) && e[k + 1] > thr[k + 1]) k <- k + 1; k }

transforms <- list(
  raw   = function(r, lo, hi) r,
  split = function(r, lo, hi) ifelse(r >= 0, r / hi, r / abs(lo)),
  U3    = function(r, lo, hi) 2 * r / (hi - lo),
  U1    = function(r, lo, hi) r / (hi - lo),
  tight = function(r, lo, hi) 2 * (r - lo) / (hi - lo) - 1
)

n <- 1000; p <- 10; thr <- pa_threshold(n, p)
patterns <- list(symmetric = rep(0, p), same_strong = rep(1.0, p),
                 opposite_mild = rep(c(0.5, -0.5), p / 2), opposite_strong = rep(c(1.0, -1.0), p / 2),
                 opposite_vstrong = rep(c(1.5, -1.5), p / 2))

run_one <- function(K, pat, lambda, nfac) {
  X <- gen_data(n, p, lambda, K, patterns[[pat]], nfac)
  if (any(apply(X, 2, function(x) length(unique(x))) < 2)) return(NULL)
  R <- cor(X); B <- bound_matrices(X); Strue <- true_cor(p, lambda, nfac)
  off <- upper.tri(R)
  mats <- lapply(transforms, function(f) { M <- f(R, B$lo, B$hi); diag(M) <- 1; M })
  pc <- tryCatch(suppressWarnings(suppressMessages(polychoric(X, smooth = FALSE, correct = 0)$rho)), error = function(e) NULL)
  if (!is.null(pc)) mats$polychoric <- pc
  rbindlist(lapply(names(mats), function(m) {
    e <- eigen(mats[[m]], symmetric = TRUE, only.values = TRUE)$values
    data.table(K = K, pattern = pat, lambda = lambda, nfac = nfac, method = m,
               k = pa_count(e, thr), not_psd = min(e) < -1e-8,
               rmse = sqrt(mean((mats[[m]][off] - Strue[off])^2)))
  }))
}

g1 <- CJ(K = c(2, 3, 5), pattern = names(patterns), lambda = c(0.5, 0.7, 0.9), rep = 1:20)
r1 <- rbindlist(lapply(seq_len(nrow(g1)), function(i) run_one(g1$K[i], g1$pattern[i], g1$lambda[i], 1)))
g2 <- CJ(K = c(2, 3, 5), pattern = names(patterns), lambda = c(0.7, 0.9), rep = 1:20)
r2 <- rbindlist(lapply(seq_len(nrow(g2)), function(i) run_one(g2$K[i], g2$pattern[i], g2$lambda[i], 2)))

lev <- c("raw", "split", "U3", "U1", "tight", "polychoric")
r1[, method := factor(method, lev)]; r2[, method := factor(method, lev)]
r1[, skew := fifelse(grepl("opposite", pattern), "opposite skew", "symmetric or same skew")]
r2[, skew := fifelse(grepl("opposite", pattern), "opposite skew", "symmetric or same skew")]

cat("A. ONE TRUE FACTOR (", nrow(r1) / length(lev), "data sets): overall\n")
print(r1[, .(spurious_2plus = round(mean(k >= 2), 3), missed_factor = round(mean(k == 0), 3),
             not_psd = round(mean(not_psd), 3), rmse_vs_latent = round(mean(rmse), 3)), by = method][order(method)])

cat("\nB. ONE TRUE FACTOR, by skew pattern and number of categories (spurious rate / not-PSD rate)\n")
tab <- r1[, .(sp = mean(k >= 2), np = mean(not_psd)), by = .(method, skew, K)]
tab[, cell := sprintf("%.2f / %.2f", sp, np)]
print(dcast(tab, skew + K ~ method, value.var = "cell"))

cat("\nC. ONE TRUE FACTOR, conditions where raw r gives spurious factors in at least half the data sets\n")
bad <- r1[method == "raw", .(raw_sp = mean(k >= 2)), by = .(K, pattern, lambda)][raw_sp >= 0.5]
tabc <- merge(r1, bad[, .(K, pattern, lambda)], by = c("K", "pattern", "lambda"))
print(tabc[, .(spurious = round(mean(k >= 2), 2), not_psd = round(mean(not_psd), 2)), by = method][order(method)])
cat("conditions:", nrow(bad), "\n")

cat("\nD. TWO TRUE FACTORS (", nrow(r2) / length(lev), "data sets): share finding at least 2 factors / exactly 2\n")
print(r2[, .(found_2plus = round(mean(k >= 2), 3), exactly_2 = round(mean(k == 2), 3),
             not_psd = round(mean(not_psd), 3), rmse_vs_latent = round(mean(rmse), 3)), by = method][order(method)])
print(dcast(r2[, .(cell = sprintf("%.2f", mean(k == 2))), by = .(method, skew, K)], skew + K ~ method, value.var = "cell"))
