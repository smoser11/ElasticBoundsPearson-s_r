# 20261118_FA_rescaling_check.R
#
# Tests the central claim of Cees's Option 1 (memo of 11 Nov 2026): that correcting each
# correlation for its bound (r / r_max if r >= 0, r / |r_min| otherwise) "would eliminate
# or diminish" the biases in factor analysis. Same simulation design as
# R/10_fa_ceiling_diagnostic_pilot.R (latent-normal factor model, thresholds shifted to
# create skew, Horn's parallel analysis with the 95th percentile of random-data eigenvalues),
# with fewer replications.
#
# Run from the repository root:
#   Rscript Notes/20261118_FA_rescaling_check.R > Notes/20261118_FA_rescaling_check_output.txt 2>&1

suppressMessages(library(data.table))
set.seed(20261118)

rmax_from_props <- function(sx, px, sy, py) {
  cx <- cumsum(px); cy <- cumsum(py); cx[length(cx)] <- 1; cy[length(cy)] <- 1
  br <- sort(unique(c(0, cx, cy))); mid <- (head(br, -1) + tail(br, -1)) / 2; w <- diff(br)
  qx <- sx[pmin(findInterval(mid, cx) + 1, length(sx))]; qy <- sy[pmin(findInterval(mid, cy) + 1, length(sy))]
  mx <- sum(px * sx); my <- sum(py * sy)
  (sum(w * qx * qy) - mx * my) / (sqrt(sum(px * sx^2) - mx^2) * sqrt(sum(py * sy^2) - my^2))
}
bound_matrices <- function(X) {
  p <- ncol(X); Mx <- diag(p); Mn <- diag(p)
  pr <- lapply(seq_len(p), function(j) { t <- table(X[, j]); list(s = as.numeric(names(t)), p = as.numeric(t) / sum(t)) })
  for (i in 1:(p - 1)) for (j in (i + 1):p) {
    Mx[i, j] <- Mx[j, i] <- rmax_from_props(pr[[i]]$s, pr[[i]]$p, pr[[j]]$s, pr[[j]]$p)
    Mn[i, j] <- Mn[j, i] <- -rmax_from_props(pr[[i]]$s, pr[[i]]$p, -rev(pr[[j]]$s), rev(pr[[j]]$p))
  }
  list(max = Mx, min = Mn)
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
pa_threshold <- function(n, p, reps = 200)
  apply(replicate(reps, eigen(cor(matrix(rnorm(n * p), n, p)), only.values = TRUE)$values), 1, quantile, 0.95)
pa_count <- function(R, thr) { e <- eigen(R, symmetric = TRUE, only.values = TRUE)$values; k <- 0; while (k < length(e) && e[k + 1] > thr[k + 1]) k <- k + 1; k }

n <- 1000; p <- 10; thr <- pa_threshold(n, p)
patterns <- list(symmetric = rep(0, p), same_strong = rep(1.0, p),
                 opposite_mild = rep(c(0.5, -0.5), p / 2), opposite_strong = rep(c(1.0, -1.0), p / 2),
                 opposite_vstrong = rep(c(1.5, -1.5), p / 2))
run_one <- function(K, pat, lambda, nfac) {
  X <- gen_data(n, p, lambda, K, patterns[[pat]], nfac)
  if (any(apply(X, 2, function(x) length(unique(x))) < 2)) return(NULL)
  R <- cor(X); B <- bound_matrices(X)
  S <- ifelse(R >= 0, R / B$max, R / abs(B$min)); diag(S) <- 1      # split rescaling
  data.table(K = K, pattern = pat, lambda = lambda, nfac = nfac,
             k_raw = pa_count(R, thr), k_rescaled = pa_count(S, thr),
             min_eig_rescaled = min(eigen(S, symmetric = TRUE, only.values = TRUE)$values))
}
g1 <- CJ(K = c(2, 3, 5), pattern = names(patterns), lambda = c(0.5, 0.7, 0.9), rep = 1:20)
r1 <- rbindlist(lapply(seq_len(nrow(g1)), function(i) run_one(g1$K[i], g1$pattern[i], g1$lambda[i], 1)))
g2 <- CJ(K = c(2, 5), pattern = c("opposite_mild", "opposite_strong", "opposite_vstrong"), lambda = c(0.7, 0.9), rep = 1:20)
r2 <- rbindlist(lapply(seq_len(nrow(g2)), function(i) run_one(g2$K[i], g2$pattern[i], g2$lambda[i], 2)))

cat("ONE TRUE FACTOR: share of data sets where parallel analysis retains 2+ factors (spurious)\n")
print(r1[, .(n = .N, raw = round(mean(k_raw >= 2), 2), rescaled = round(mean(k_rescaled >= 2), 2),
             rescaled_not_psd = round(mean(min_eig_rescaled < -1e-8), 2)), by = .(pattern, K, lambda)][order(pattern, K, lambda)])
cat(sprintf("\noverall: spurious with raw r %.3f; with rescaled r %.3f; rescaled matrix not PSD %.3f (of %d data sets)\n",
            mean(r1$k_raw >= 2), mean(r1$k_rescaled >= 2), mean(r1$min_eig_rescaled < -1e-8), nrow(r1)))

cat("\nTWO TRUE FACTORS: share of data sets where parallel analysis finds both (2+ factors)\n")
print(r2[, .(n = .N, raw = round(mean(k_raw >= 2), 2), rescaled = round(mean(k_rescaled >= 2), 2),
             rescaled_not_psd = round(mean(min_eig_rescaled < -1e-8), 2)), by = .(pattern, K, lambda)][order(pattern, K, lambda)])
