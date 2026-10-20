# 10_fa_ceiling_diagnostic_pilot.R
#
# PILOT. Does the "ceiling matrix" predict when Pearson-based factor retention
# finds spurious factors in ordinal data?
#
# The ceiling matrix R_max holds r_max for every item pair. It is the Pearson
# correlation matrix of perfectly one-dimensional data with the observed
# response distributions (every item a monotone function of one latent
# variable), so it needs nothing but the item frequencies. Two uses:
#   (1) Screening: how far is R_max from one-dimensional? (D1 = second
#       eigenvalue of R_max; D2 = off-diagonal residual after the first
#       component.) Compared with the usual rules of thumb (number of
#       categories, item skewness).
#   (2) Attribution: when Pearson parallel analysis retains a second factor, is
#       it aligned with R_max's second component (|cosine| of the eigenvectors)?
#       Compared between spurious second factors (one-factor truth) and real
#       ones (two-factor truth).
#
# Data: latent normal factor model, discretised by thresholds that are shifted
# to create skew (as in Rhemtulla et al. 2012; Garrido et al. 2013).
# Writes: output/data/fa_ceiling_pilot.rds

suppressMessages(library(data.table))
set.seed(20261020)

# r_max from proportions (as in 08_anes_pairs_bounds.R)
rmax_from_props <- function(sx, px, sy, py) {
  cx <- cumsum(px); cy <- cumsum(py); cx[length(cx)] <- 1; cy[length(cy)] <- 1
  br <- sort(unique(c(0, cx, cy))); mid <- (head(br, -1) + tail(br, -1)) / 2; w <- diff(br)
  qx <- sx[pmin(findInterval(mid, cx) + 1, length(sx))]; qy <- sy[pmin(findInterval(mid, cy) + 1, length(sy))]
  mx <- sum(px * sx); my <- sum(py * sy)
  (sum(w * qx * qy) - mx * my) / (sqrt(sum(px * sx^2) - mx^2) * sqrt(sum(py * sy^2) - my^2))
}
ceiling_matrix <- function(X) {
  p <- ncol(X); M <- diag(p)
  props <- lapply(seq_len(p), function(j) { t <- table(X[, j]); list(s = as.numeric(names(t)), p = as.numeric(t) / sum(t)) })
  for (i in 1:(p - 1)) for (j in (i + 1):p)
    M[i, j] <- M[j, i] <- rmax_from_props(props[[i]]$s, props[[i]]$p, props[[j]]$s, props[[j]]$p)
  M
}
skewness <- function(x) { m <- mean(x); mean((x - m)^3) / mean((x - m)^2)^1.5 }

# Data generation
gen_data <- function(n, p, lambda, K, shifts, nfac = 1, phi = 0.3) {
  if (nfac == 1) {
    f <- rnorm(n); Z <- sapply(seq_len(p), function(j) lambda * f + sqrt(1 - lambda^2) * rnorm(n))
  } else {
    f <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, phi, phi, 1), 2))
    g <- rep(1:2, each = p / 2)       # items 1..p/2 on factor 1, rest on factor 2
    Z <- sapply(seq_len(p), function(j) lambda * f[, g[j]] + sqrt(1 - lambda^2) * rnorm(n))
  }
  base <- if (K == 2) 0 else seq(-1.5, 1.5, length.out = K - 1)
  sapply(seq_len(p), function(j) findInterval(Z[, j], base + shifts[j]) + 1)
}

# Horn's parallel analysis on PCA eigenvalues (95th percentile of random-normal eigenvalues)
pa_threshold <- function(n, p, reps = 200) {
  apply(replicate(reps, eigen(cor(matrix(rnorm(n * p), n, p)), only.values = TRUE)$values), 1, quantile, 0.95)
}
pa_count <- function(R, thr) { e <- eigen(R, only.values = TRUE)$values; k <- 0; while (k < length(e) && e[k + 1] > thr[k + 1]) k <- k + 1; k }

n <- 1000; p <- 10
thr <- pa_threshold(n, p)

# Skew patterns: threshold shifts (positive shift = more mass in low categories)
patterns <- list(
  symmetric     = rep(0, p),
  same_mild     = rep(0.5, p),
  same_strong   = rep(1.0, p),
  opposite_mild   = rep(c(0.5, -0.5), p / 2),     # alternating so skew split is not the factor split
  opposite_strong = rep(c(1.0, -1.0), p / 2),
  opposite_vstrong = rep(c(1.5, -1.5), p / 2)
)

run_one <- function(K, pat, lambda, nfac) {
  X <- gen_data(n, p, lambda, K, patterns[[pat]], nfac)
  if (any(apply(X, 2, function(x) length(unique(x))) < 2)) return(NULL)
  R <- cor(X); Rm <- ceiling_matrix(X)
  eR <- eigen(R, symmetric = TRUE); eM <- eigen(Rm, symmetric = TRUE)
  resid <- Rm - eM$values[1] * tcrossprod(eM$vectors[, 1]); diag(resid) <- 0
  sk <- apply(X, 2, skewness)
  opp <- outer(sk, sk); opp <- mean(pmax(0, -opp[upper.tri(opp)]))
  data.table(K = K, pattern = pat, lambda = lambda, nfac = nfac,
             n_factors_pearson = pa_count(R, thr),
             D1_eig2_ceiling = eM$values[2],
             D2_resid_ceiling = sqrt(mean(resid[upper.tri(resid)]^2)),
             max_abs_skew = max(abs(sk)), opp_skew_index = opp,
             align_v2 = abs(sum(eR$vectors[, 2] * eM$vectors[, 2])),
             align_skew = abs(cor(eR$vectors[, 2], sk)),
             mean_rmax = mean(Rm[upper.tri(Rm)]), mean_r = mean(R[upper.tri(R)]))
}

grid1 <- CJ(K = c(2, 3, 5, 7), pattern = names(patterns), lambda = c(0.5, 0.7, 0.9), rep = 1:40)
cat("One-factor truth:", nrow(grid1), "datasets\n")
res1 <- rbindlist(lapply(seq_len(nrow(grid1)), function(i) run_one(grid1$K[i], grid1$pattern[i], grid1$lambda[i], 1)))
grid2 <- CJ(K = c(2, 5), pattern = c("opposite_mild", "opposite_strong", "opposite_vstrong"), lambda = c(0.7, 0.9), rep = 1:40)
cat("Two-factor truth:", nrow(grid2), "datasets\n")
res2 <- rbindlist(lapply(seq_len(nrow(grid2)), function(i) run_one(grid2$K[i], grid2$pattern[i], grid2$lambda[i], 2)))

saveRDS(list(one = res1, two = res2), "output/data/fa_ceiling_pilot.rds")

# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------
res1[, over := n_factors_pearson >= 2]
cat("\n1. Over-extraction rate (Pearson parallel analysis retains >= 2 factors; truth = 1):\n")
print(dcast(res1[, .(rate = round(mean(over), 2)), by = .(pattern, K, lambda)], pattern + lambda ~ K, value.var = "rate"))

auc <- function(score, y) { r <- rank(score); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
cat("\n2. How well does each marginals-only screen predict over-extraction? (AUC; 0.5 = useless, 1 = perfect)\n")
preds <- c("D1_eig2_ceiling", "D2_resid_ceiling", "max_abs_skew", "opp_skew_index")
print(round(sapply(preds, function(v) auc(res1[[v]], res1$over)), 3))
cat("   ...and within loading levels (lambda):\n")
print(round(t(sapply(c(0.5, 0.7, 0.9), function(l) { s <- res1[lambda == l]; sapply(preds, function(v) if (length(unique(s$over)) > 1) auc(s[[v]], s$over) else NA) })), 3))
cat("   Number of categories alone:", round(auc(-res1$K, res1$over), 3), "\n")

cat("\n3. Attribution: alignment of the retained second factor with R_max's second component\n")
spur <- res1[over == TRUE]
real <- res2[n_factors_pearson >= 2]
cat(sprintf("   spurious 2nd factors (1-factor truth, n=%d): median |cos| = %.2f; median |cor with skew| = %.2f\n",
            nrow(spur), median(spur$align_v2), median(spur$align_skew)))
cat(sprintf("   real 2nd factors     (2-factor truth, n=%d): median |cos| = %.2f; median |cor with skew| = %.2f\n",
            nrow(real), median(real$align_v2), median(real$align_skew)))
if (nrow(spur) > 0 && nrow(real) > 0)
  cat(sprintf("   AUC separating spurious from real by |cos| alignment: %.3f\n",
              auc(c(spur$align_v2, real$align_v2), c(rep(TRUE, nrow(spur)), rep(FALSE, nrow(real))))))
