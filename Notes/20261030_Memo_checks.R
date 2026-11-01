# 20261030_Memo_checks.R
#
# Recomputes the numbers from our own analyses that are quoted in the memo to Cees
# (Notes/20261030_Memo_Scott_Cees_Meeting20261031.qmd), from the data and saved outputs
# in this repository. Numbers quoted from published work are NOT recomputed: we have no
# access to other authors' data. The one exception is Shojima's (2025) bounds (section 3),
# which depend only on the marginals printed in his tables. Section 4 re-estimates
# Baldassarri & Gelman's model from the public ANES file, not from their own data.
# Run from the repository root:
#   Rscript Notes/20261030_Memo_checks.R > Notes/20261030_Memo_checks_output.txt 2>&1
#
# Reads: R/data/raw/bes2019_pairs.dta, output/data/bes_bounds.rds,
#        output/data/anes_pairs_bounds.rds, output/data/anes_item_marginals.rds,
#        output/data/anes_alignment_decomposition.rds, output/data/fa_ceiling_pilot.rds,
#        output/data/adjustment_psd_synthetic.rds, output/data/adjustment_psd_bes_real.rds
# (the ANES files are produced by R/08 and R/09 and need the ANES Cumulative Data File).

suppressMessages({library(haven); library(data.table); library(lme4)})
options(width = 120)
section <- function(s) cat("\n==============================\n", s, "\n==============================\n")

# Exact r_max from two marginals (scores ascending, proportions): comonotone coupling
rmax_from_props <- function(sx, px, sy, py) {
  cx <- cumsum(px); cy <- cumsum(py); cx[length(cx)] <- 1; cy[length(cy)] <- 1
  br <- sort(unique(c(0, cx, cy))); mid <- (head(br, -1) + tail(br, -1)) / 2; w <- diff(br)
  qx <- sx[pmin(findInterval(mid, cx) + 1, length(sx))]; qy <- sy[pmin(findInterval(mid, cy) + 1, length(sy))]
  mx <- sum(px * sx); my <- sum(py * sy)
  (sum(w * qx * qy) - mx * my) / (sqrt(sum(px * sx^2) - mx^2) * sqrt(sum(py * sy^2) - my^2))
}
rmin_from_props <- function(sx, px, sy, py) -rmax_from_props(sx, px, -rev(sy), rev(py))
r_of_table <- function(P, wa, wb) {
  pa <- rowSums(P); pb <- colSums(P); ma <- sum(pa * wa); mb <- sum(pb * wb)
  (sum(outer(wa, wb) * P) - ma * mb) / sqrt((sum(pa * wa^2) - ma^2) * (sum(pb * wb^2) - mb^2))
}

# ---------------------------------------------------------------------------
section("1. BES: how far correlations are from their ceilings (memo 4.1)")
b <- readRDS("output/data/bes_bounds.rds")
d <- read_dta("R/data/raw/bes2019_pairs.dta")
stopifnot(nrow(b) == nrow(d), all(b$pair_id == seq_len(nrow(d))))
rstar <- ifelse(b$r_obs >= 0, b$r_obs / b$r_max, b$r_obs / abs(b$r_min))
cat(sprintf("pairs: %d\nmedian |r| = %.3f\nmedian r_max = %.3f\nmedian |r*| (r over the bound on its own side) = %.3f\nshare of pairs with |r*| > .5 = %.3f\n",
            nrow(b), median(abs(b$r_obs)), median(b$r_max), median(abs(rstar)), mean(abs(rstar) > .5)))

# ---------------------------------------------------------------------------
section("2. Wasserstein formulas for r_max and r_min, all BES pairs (memo 7.1)")
# Marginals as R/01 extracts them: codes start..start+K-1 (start = 0 if freq0 > 0)
get_counts <- function(row, v, K) {
  start <- if (!is.na(row[[paste0(v, "freq0")]]) && row[[paste0(v, "freq0")]] > 0) 0 else 1
  sapply(start:(start + K - 1), function(i) as.numeric(row[[paste0(v, "freq", i)]]))
}
w2sq <- function(zx, px, zy, py) {   # squared 2-Wasserstein distance via quantile functions
  cx <- cumsum(px); cy <- cumsum(py); cx[length(cx)] <- 1; cy[length(cy)] <- 1
  br <- sort(unique(c(0, cx, cy))); mid <- (head(br, -1) + tail(br, -1)) / 2; w <- diff(br)
  qx <- zx[pmin(findInterval(mid, cx) + 1, length(zx))]; qy <- zy[pmin(findInterval(mid, cy) + 1, length(zy))]
  sum(w * (qx - qy)^2)
}
stdz <- function(s, p) { m <- sum(p * s); (s - m) / sqrt(sum(p * (s - m)^2)) }
skew <- function(s, p) { m <- sum(p * s); sum(p * (s - m)^3) / sum(p * (s - m)^2)^1.5 }
res <- t(sapply(seq_len(nrow(d)), function(i) {
  row <- d[i, ]
  c1 <- get_counts(row, "var1", row$var1cats); c2 <- get_counts(row, "var2", row$var2cats)
  s1 <- seq_along(c1) - 1; s2 <- seq_along(c2) - 1; p1 <- c1 / sum(c1); p2 <- c2 / sum(c2)
  k1 <- p1 > 0; k2 <- p2 > 0; s1 <- s1[k1]; p1 <- p1[k1]; s2 <- s2[k2]; p2 <- p2[k2]
  z1 <- stdz(s1, p1); z2 <- stdz(s2, p2)
  c(rmax_W = 1 - w2sq(z1, p1, z2, p2) / 2,
    rmin_W = -1 + w2sq(z1, p1, -rev(z2), rev(p2)) / 2,
    sk1 = skew(s1, p1), sk2 = skew(s2, p2))
}))
cat(sprintf("max |r_max(formula) - r_max(R/01 sorting)| = %.2e\nmax |r_min(formula) - r_min(R/01 sorting)| = %.2e\n",
            max(abs(res[, "rmax_W"] - b$r_max)), max(abs(res[, "rmin_W"] - b$r_min))))
asym <- abs(b$r_min) - b$r_max; prod_sk <- -res[, "sk1"] * res[, "sk2"]
cat(sprintf("cor(|r_min| - r_max, -s1*s2) = %.3f\nsign agreement: %.3f overall; %.3f where |s1*s2| > .1\n",
            cor(asym, prod_sk), mean(sign(asym) == sign(prod_sk)),
            mean((sign(asym) == sign(prod_sk))[abs(prod_sk) > .1])))
p <- c(.6, .3, .1); s <- 0:2
cat(sprintf("two identical items p = (.6,.3,.1): r_max = %.3f, r_min = %.3f\n",
            rmax_from_props(s, p, s, p), rmin_from_props(s, p, s, p)))

# ---------------------------------------------------------------------------
section("3. Shojima (2025) examples and the ordering point (memo 3, row 3)")
piA <- c(.12, .36, .52); piB <- c(.16, .16, .28, .26, .14)          # Shojima Table 3, p. 512
cat(sprintf("Table 3 (3x5): r_max = %.3f (Shojima p. 513: .858); r_min = %.3f (p. 518: -.845)\n",
            rmax_from_props(1:3, piA, 1:5, piB), rmin_from_props(1:3, piA, 1:5, piB)))
cat(sprintf("Table 1 (2x3): r_max = %.3f, r_min = %.3f (Shojima p. 509: .795, -.767)\n",
            rmax_from_props(1:2, c(.30, .70), 1:3, c(.25, .40, .35)),
            rmin_from_props(1:2, c(.30, .70), 1:3, c(.25, .40, .35))))
nwcr <- function(a, b) {             # northwest-corner rule in the listed category order
  P <- matrix(0, length(a), length(b)); i <- 1; j <- 1
  while (i <= length(a) && j <= length(b)) {
    x <- min(a[i], b[j]); P[i, j] <- x; a[i] <- a[i] - x; b[j] <- b[j] - x
    if (a[i] < 1e-12) i <- i + 1 else j <- j + 1
  }
  P
}
P <- nwcr(piA, piB)
cat(sprintf("NWCR table, scores 1,2,3: r = %.3f (equals the maximum)\n", r_of_table(P, 1:3, 1:5)))
wA <- c(1, 3, 2); o <- order(wA)
cat(sprintf("Same NWCR table, item-A scores 1,3,2: r = %.3f; true maximum with these scores = %.3f\n",
            r_of_table(P, wA, 1:5), rmax_from_props(wA[o], piA[o], 1:5, piB)))
cat("(Appendix 5's argument never uses the order of the scores, so if it were complete it would\n prove the NWCR optimal here too.)\n")

# ---------------------------------------------------------------------------
section("4. ANES: replication checks (memo 4.2)")
pr <- readRDS("output/data/anes_pairs_bounds.rds")
cat(sprintf("pair-years: %d; surveys: %d (%d-%d); bound violations: %d\n", nrow(pr), uniqueN(pr$year),
            min(pr$year), max(pr$year), sum(pr$r_w < pr$r_min - 1e-9 | pr$r_w > pr$r_max + 1e-9)))
m <- readRDS("output/data/anes_item_marginals.rds")
cat("number of categories per item (items x K):\n"); print(table(m[, .(K = uniqueN(score)), by = item]$K))
bg <- pr[type == "issue_party" & year <= 2004][, pair := paste(var1, var2)]
for (v in c("r_wgcna", "r_w")) for (c0 in c(1980, 1988)) {
  bg[, t := (year - c0) / 10]
  f <- fixef(suppressMessages(lmer(as.formula(paste(v, "~ t + (1 + t | pair)")), data = bg)))
  cat(sprintf("B&G Table 1A model, %s, t zeroed at %d: intercept %.3f, slope %.3f per decade (%d pair-years, %d issues)\n",
              v, c0, f[1], f[2], nrow(bg), uniqueN(bg$pair)))
}
x <- readRDS("output/data/anes_alignment_decomposition.rds")
cat("2016 levels, K&M weighting (compare K&M Fig. 3, p. 7):\n")
print(x$levels[measure == "r_wgcna" & year == 2016, .(group, r = round(value, 2))])

# ---------------------------------------------------------------------------
section("5. ANES: trend decomposition (memo 4.3 table)")
print(x$trends[, .(group, period, trend_r = round(slope_r, 3), bounds = round(bounds_part, 3),
                   within = round(position_part, 3), trend_rmax = round(slope_rmax, 3))])
cat("\nstrong party identifiers (codes 1 and 7):\n")
print(m[item == "partyid" & year %in% c(1972, 2004, 2016), .(strong = round(sum(prop[score %in% c(1, 7)]), 2)), by = year])
for (it in c("gays_adoptkids", "gayer_military")) {
  cat("\n", it, "(share in category 1, i.e. in favour):\n")
  print(m[item == it & score == 1, .(year, prop = round(prop, 2))])
}
p2 <- pr[year <= 2016][, rs := r_w / r_max][, pair := paste(var1, var2)]
setorder(p2, pair, year)
p2[, `:=`(dr = r_w - shift(r_w), dbound = (rs + shift(rs)) / 2 * (r_max - shift(r_max))), by = pair]
q <- p2[!is.na(dr)]
cat(sprintf("\nwave-to-wave changes in individual pairs: %d\nmechanical part: median |.| = %.3f, 90th pct = %.3f, share of var(dr) = %.3f, share > .05 = %.4f\n",
            nrow(q), median(abs(q$dbound)), quantile(abs(q$dbound), .9), var(q$dbound) / var(q$dr), mean(abs(q$dbound) > .05)))
cat(sprintf("between-wave total-variation distance of item marginals: median %.3f\n", median(x$tv$tv)))
cat(sprintf("share of pair-years (1972-2016) with r < 0: %.3f\n", mean(p2$r_w < 0)))
p2[type == "issue_issue", grp := fifelse(dom1 == dom2, dom1, "cross")]
print(p2[type == "issue_issue", .(share_negative = round(mean(r_w < 0), 3)), by = grp])

# ---------------------------------------------------------------------------
section("6. Factor-analysis pilot, R/10 (memo 4.4)")
fa <- readRDS("output/data/fa_ceiling_pilot.rds"); r1 <- fa$one; r2 <- fa$two
r1[, over := n_factors_pearson >= 2]
auc <- function(score, y) { r <- rank(score); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
cat(sprintf("one-factor data sets: %d; two-factor: %d; spurious 2nd factor in %.3f\n", nrow(r1), nrow(r2), mean(r1$over)))
print(dcast(r1[, .(rate = round(mean(over), 2)), by = .(pattern, K, lambda)], pattern + lambda ~ K, value.var = "rate"))
for (v in c("D1_eig2_ceiling", "opp_skew_index", "max_abs_skew")) cat(sprintf("AUC %s = %.3f\n", v, auc(r1[[v]], r1$over)))
cat(sprintf("AUC number of categories = %.3f\n", auc(-r1$K, r1$over)))
cat("AUC of R_max 2nd eigenvalue within loading levels:", round(sapply(c(.5, .7, .9), function(l) auc(r1[lambda == l]$D1_eig2_ceiling, r1[lambda == l]$over)), 3), "\n")
sp <- r1[over == TRUE]; re <- r2[n_factors_pearson >= 2]
cat(sprintf("alignment with R_max 2nd component: spurious median |cos| %.2f (n=%d), real %.2f (n=%d), AUC %.3f\n",
            median(sp$align_v2), nrow(sp), median(re$align_v2), nrow(re),
            auc(c(sp$align_v2, re$align_v2), c(rep(TRUE, nrow(sp)), rep(FALSE, nrow(re))))))
v <- 25:36; s <- d$var1 %in% v & d$var2 %in% v
M <- diag(12); R <- diag(12); i <- match(d$var1[s], v); j <- match(d$var2[s], v)
M[cbind(i, j)] <- M[cbind(j, i)] <- b$r_max[s]; R[cbind(i, j)] <- R[cbind(j, i)] <- b$r_obs[s]
cat(sprintf("BES items v25-v36 (%d pairs): eigenvalues of R: %s; of R_max: %s\n", sum(s),
            paste(round(eigen(R)$values[1:3], 2), collapse = ", "), paste(round(eigen(M)$values[1:3], 2), collapse = ", ")))

# ---------------------------------------------------------------------------
section("7. Invertibility (PSD) checks, R/06 (memo 7.4)")
syn <- readRDS("output/data/adjustment_psd_synthetic.rds"); real <- readRDS("output/data/adjustment_psd_bes_real.rds")
cat("adversarial synthetic, Pearson (failure %):\n"); print(syn$adversarial$pearson[, c("transform", "n", "fail_pct")])
cat("mild synthetic, Pearson (failure %):\n"); print(syn$mild$pearson[, c("transform", "n", "fail_pct")])
for (k in names(real)) { cat("BES", k, "(failure %):\n"); print(real[[k]][, c("transform", "n", "fail_pct")]) }

# ---------------------------------------------------------------------------
section("8. Margin standardization keeps odds ratios but changes r (memo 7.7)")
P <- matrix(c(30, 15, 5, 2,  15, 25, 12, 5,  6, 14, 24, 14,  2, 6, 15, 30), 4, byrow = TRUE); P <- P / sum(P)
ipf <- function(P, a, b, it = 500) { for (k in 1:it) { P <- P * (a / rowSums(P)); P <- t(t(P) * (b / colSums(P))) }; P }
lor <- function(P) log(P[-1, -1] * P[-nrow(P), -ncol(P)] / (P[-1, -ncol(P)] * P[-nrow(P), -1]))
Q <- ipf(P, c(.70, .15, .10, .05), c(.05, .10, .15, .70))
cat(sprintf("original table: r = %.3f; after giving it lopsided, opposite-skewed margins: r = %.3f\nmax change in any local log odds ratio = %.1e\n",
            r_of_table(P, 1:4, 1:4), r_of_table(Q, 1:4, 1:4), max(abs(lor(P) - lor(Q)))))
