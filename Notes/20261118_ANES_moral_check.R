# 20261118_ANES_moral_check.R
#
# A closer look at the ANES moral-issue alignment result (memo of 31 Oct, Section 4.3).
# Uses our own estimates from the public ANES Cumulative Data File (output of
# R/08_anes_pairs_bounds.R), not Kozlowski & Murphy's data or models.
#
# Moral domain = pairs of moral items (Kozlowski & Murphy's moral issues). Before 1986
# only ONE moral pair is observed (women_equals x abortion), so the long-run period
# here is 1986-2016; 2004-2016 is K&M's headline period. All moral correlations are
# positive, so r* = r / r_max throughout (no r_min needed in this domain).
#
# Run from the repository root:
#   Rscript Notes/20261118_ANES_moral_check.R > Notes/20261118_ANES_moral_check_output.txt 2>&1

suppressMessages({library(data.table); library(lme4)})
options(width = 130)
p <- readRDS("output/data/anes_pairs_bounds.rds")[year <= 2016]
m <- p[type == "issue_issue" & dom1 == "moral" & dom2 == "moral"]
m[, `:=`(pair = paste(var1, var2, sep = " x "), rstar = r_w / r_max)]
gay <- c("protect_homo", "gayer_military", "gays_adoptkids")
m[, gaypair := var1 %in% gay | var2 %in% gay]
cat(sprintf("moral pair-years: %d; pairs: %d; negative r: %d\n", nrow(m), uniqueN(m$pair), sum(m$r_w < 0)))

section <- function(s) cat("\n==============================\n", s, "\n==============================\n")

section("1. Trends per decade with standard errors (mixed models, random intercept and slope per pair)")
trend <- function(dt, v) {
  fit <- suppressMessages(lmer(as.formula(paste(v, "~ t + (1 + t | pair)")), data = dt))
  co <- summary(fit)$coefficients
  c(est = co["t", "Estimate"], se = co["t", "Std. Error"])
}
for (per in list(c(1986, 2016), c(2004, 2016))) {
  for (sub in c("all moral pairs", "pairs with a gay-rights item", "other moral pairs")) {
    dt <- m[year >= per[1] & year <= per[2]]
    if (sub == "pairs with a gay-rights item") dt <- dt[gaypair == TRUE]
    if (sub == "other moral pairs") dt <- dt[gaypair == FALSE]
    dt[, t := (year - 2000) / 10]
    out <- sapply(c("r_w", "r_max", "rstar"), function(v) trend(dt, v))
    cat(sprintf("%d-%d, %-29s (%3d pair-years, %2d pairs): trend r %+.3f (se %.3f) | r_max %+.3f (se %.3f) | r* %+.3f (se %.3f)\n",
                per[1], per[2], sub, nrow(dt), uniqueN(dt$pair),
                out["est", "r_w"], out["se", "r_w"], out["est", "r_max"], out["se", "r_max"],
                out["est", "rstar"], out["se", "rstar"]))
  }
}

section("2. Year-by-year means over the same set of pairs (pairs observed in all of 2004, 2008, 2012, 2016)")
keep <- m[year %in% c(2004, 2008, 2012, 2016), .N, by = pair][N == 4, pair]
cat("pairs observed in all four years:", length(keep), "\n")
print(m[pair %in% keep & year >= 2004, .(n_pairs = .N, mean_r = round(mean(r_w), 3),
        mean_rmax = round(mean(r_max), 3), mean_rstar = round(mean(rstar), 3)), by = year][order(year)])

section("3. 2004 vs 2016, pairs observed in both years")
a <- m[year == 2004, .(pair, gaypair, r04 = r_w, mx04 = r_max, rs04 = rstar)]
b <- m[year == 2016, .(pair, r16 = r_w, mx16 = r_max, rs16 = rstar)]
ab <- merge(a, b, by = "pair")
ab[, cf16 := rs16 * mx04]     # 2016 position within the bounds, but 2004 ceiling
print(ab[, .(pair, gaypair, r04 = round(r04, 3), r16 = round(r16, 3), rmax04 = round(mx04, 3),
             rmax16 = round(mx16, 3), rstar04 = round(rs04, 3), rstar16 = round(rs16, 3))][order(-gaypair, pair)])
cat(sprintf("\nmeans over %d pairs: r 2004 %.3f -> 2016 %.3f (change %+.3f); r_max %.3f -> %.3f; r* %.3f -> %.3f\n",
            nrow(ab), mean(ab$r04), mean(ab$r16), mean(ab$r16 - ab$r04), mean(ab$mx04), mean(ab$mx16),
            mean(ab$rs04), mean(ab$rs16)))
cat(sprintf("2016 r if the ceiling had stayed at its 2004 level (r*_2016 x r_max_2004): %.3f (change from 2004 %+.3f)\n",
            mean(ab$cf16), mean(ab$cf16 - ab$r04)))
for (g in c(TRUE, FALSE)) {
  s <- ab[gaypair == g]
  cat(sprintf("  %-28s (%d pairs): r %+.3f, r_max %+.3f, r* %+.3f\n", ifelse(g, "pairs with a gay-rights item", "other moral pairs"),
              nrow(s), mean(s$r16 - s$r04), mean(s$mx16 - s$mx04), mean(s$rs16 - s$rs04)))
}

section("4. The marginals behind the falling ceiling (weighted share in category 1)")
mi <- readRDS("output/data/anes_item_marginals.rds")[item %in% gay & year >= 1992 & year <= 2016]
print(dcast(mi[score == 1, .(item, year, p1 = round(prop, 2))], item ~ year, value.var = "p1"))
cat("(gayer_military and gays_adoptkids have 2 categories; protect_homo has 5.)\n")
