# 09_anes_alignment_decomposition.R
#
# How much of the measured change in ANES issue alignment (Kozlowski & Murphy
# 2021, Fig. 3) and issue partisanship (Baldassarri & Gelman 2008) comes from
# the attainable bounds moving, and how much from tighter association given
# the marginals?
#
# Every correlation factors exactly as r = r_star * r_max, with r_star = r / r_max
# (the share of the attainable ceiling that is used). For one item pair observed
# in waves a and b, the change splits exactly (two-factor Shapley split):
#     dr = mean(r_star) * d(r_max)   [bounds moving]
#        + mean(r_max) * d(r_star)   [position within bounds]
# Changes are chained over consecutive waves using the pairs observed in both
# waves, so the changing item battery does not create spurious change.
#
# Reads:  output/data/anes_pairs_bounds.rds, output/data/anes_item_marginals.rds
#         (both from 08_anes_pairs_bounds.R)
# Writes: output/data/anes_alignment_decomposition.rds
#         output/figures/16_anes_alignment_levels.pdf
#         output/figures/17_anes_alignment_decomposition.pdf

suppressMessages({library(data.table); library(lme4); library(ggplot2)})

pairs <- readRDS("output/data/anes_pairs_bounds.rds")[year <= 2016]
marg  <- readRDS("output/data/anes_item_marginals.rds")[year <= 2016]
pairs[, r_star := r_w / r_max]
pairs[, pair := paste(var1, var2, sep = "__")]

dom_label <- c(econ = "Economic", civr = "Civil rights", moral = "Moral", frgn = "Security")
pairs[type == "issue_issue", group := fifelse(dom1 == dom2, dom_label[dom1], "Cross-domain")]
pairs[type == "issue_party", group := "Issue x party ID"]
pairs[type == "issue_ideo",  group := "Issue x ideology"]
grp_levels <- c("Economic", "Civil rights", "Moral", "Security", "Cross-domain",
                "Issue x party ID", "Issue x ideology")
pa <- pairs[!is.na(group)]
pa[, group := factor(group, levels = grp_levels)]

# ---------------------------------------------------------------------------
# 1. Levels by year: mixed models with year effects and pair random intercepts
#    (adjusts for which pairs were asked in which year, as K&M do)
# ---------------------------------------------------------------------------
year_means <- function(dt, v) {
  m <- lmer(as.formula(paste(v, "~ 0 + factor(year) + (1 | pair)")), data = dt)
  f <- fixef(m)
  data.table(year = as.integer(sub("factor\\(year\\)", "", names(f))), value = unname(f))
}
levels_dt <- rbindlist(lapply(grp_levels, function(g) {
  dt <- pa[group == g]
  rbindlist(lapply(c(r_wgcna = "r_wgcna", r_w = "r_w", r_max = "r_max", r_star = "r_star"),
                   function(v) year_means(dt, v)[, `:=`(measure = v, group = g)]))
}))

# ---------------------------------------------------------------------------
# 2. Chained decomposition over consecutive waves
# ---------------------------------------------------------------------------
waves <- sort(unique(pa$year))
steps <- rbindlist(lapply(seq_along(waves)[-1], function(i) {
  a <- pa[year == waves[i - 1], .(pair, group, r_a = r_w, m_a = r_max, s_a = r_star)]
  b <- pa[year == waves[i],     .(pair, r_b = r_w, m_b = r_max, s_b = r_star)]
  ab <- merge(a, b, by = "pair")
  ab[, .(n_pairs = .N,
         d_r      = mean(r_b - r_a),
         d_bounds = mean((s_a + s_b) / 2 * (m_b - m_a)),
         d_pos    = mean((m_a + m_b) / 2 * (s_b - s_a))),
     by = group][, `:=`(from = waves[i - 1], to = waves[i])]
}))
steps <- steps[order(group, to)]
steps[, `:=`(cum_r = cumsum(d_r), cum_bounds = cumsum(d_bounds), cum_pos = cumsum(d_pos)), by = group]

# ---------------------------------------------------------------------------
# 3. K&M's headline period: 2004 -> 2016, and B&G's period: 1972 -> 2004
# ---------------------------------------------------------------------------
period_split <- function(y0, y1) {
  s <- steps[from >= y0 & to <= y1]
  s[, .(d_r = sum(d_r), d_bounds = sum(d_bounds), d_pos = sum(d_pos),
        min_pairs = min(n_pairs)), by = group][
    , share_bounds := d_bounds / d_r][, period := paste0(y0, "-", y1)]
}
periods <- rbind(period_split(1972, 2004), period_split(2004, 2016))
# NB: chaining is only as good as the overlap between consecutive waves. Before
# 2004 some within-domain links rest on a single item pair (see min_pairs), so
# the trend models below are the main estimates; the chain is kept as a check.

# ---------------------------------------------------------------------------
# 3b. Trend decomposition: mixed models with pair random intercepts and slopes
#     (time in decades). slope(r) ~= mean(r_star) * slope(r_max) + mean(r_max) * slope(r_star)
# ---------------------------------------------------------------------------
trend_slope <- function(dt, v) {
  m <- suppressMessages(lmer(as.formula(paste(v, "~ t + (1 + t | pair)")), data = dt))
  fixef(m)[["t"]]
}
trends <- rbindlist(lapply(grp_levels, function(g) rbindlist(lapply(list(c(1972, 2016), c(2004, 2016)), function(p) {
  dt <- pa[group == g & year >= p[1] & year <= p[2]][, t := (year - 1990) / 10]
  if (uniqueN(dt$year) < 3) return(NULL)
  b_m <- trend_slope(dt, "r_max"); b_s <- trend_slope(dt, "r_star")
  data.table(group = g, period = paste(p, collapse = "-"), pair_years = nrow(dt),
             slope_r = trend_slope(dt, "r_w"), slope_rmax = b_m, slope_rstar = b_s,
             bounds_part = mean(dt$r_star) * b_m, position_part = mean(dt$r_max) * b_s)
}))))

# ---------------------------------------------------------------------------
# 4. How big are the between-wave marginal shifts? (total variation distance)
# ---------------------------------------------------------------------------
tv <- rbindlist(lapply(seq_along(waves)[-1], function(i) {
  a <- marg[year == waves[i - 1], .(item, score, pa = prop)]
  b <- marg[year == waves[i],     .(item, score, pb = prop)]
  ab <- merge(a, b, by = c("item", "score"), all = TRUE)
  ab[is.na(pa), pa := 0][is.na(pb), pb := 0]
  ab <- ab[item %in% intersect(a$item, b$item)]
  ab[, .(tv = 0.5 * sum(abs(pb - pa))), by = item][, `:=`(from = waves[i - 1], to = waves[i])]
}))

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
cat("\nTrend decomposition (per decade; design-weighted r):\n")
print(trends[, .(group, period, pair_years, slope_r = round(slope_r, 3), bounds_part = round(bounds_part, 3),
                 position_part = round(position_part, 3), slope_rmax = round(slope_rmax, 3))])
cat("\nChained decomposition (check only; fragile where min_pairs is small):\n")
print(periods[, .(period, group, d_r = round(d_r, 3), d_bounds = round(d_bounds, 3),
                  d_pos = round(d_pos, 3), min_pairs)])
cat("\nBetween-wave marginal shifts (total variation distance, items asked in consecutive waves):\n")
print(round(quantile(tv$tv, c(.1, .25, .5, .75, .9, 1)), 3))

saveRDS(list(levels = levels_dt, trends = trends, steps = steps, periods = periods, tv = tv),
        "output/data/anes_alignment_decomposition.rds")

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
lv <- levels_dt[measure %in% c("r_wgcna", "r_max", "r_star")]
lv[, measure := factor(measure, levels = c("r_wgcna", "r_max", "r_star"),
                       labels = c("Correlation r (K&M weighting)", "Ceiling r_max",
                                  "Share of ceiling used, r* = r / r_max"))]
p1 <- ggplot(lv, aes(year, value, colour = group)) +
  geom_line() + geom_point(size = 0.8) +
  facet_wrap(~ measure, scales = "free_y", ncol = 1) +
  annotate("rect", xmin = 2004, xmax = 2016, ymin = -Inf, ymax = Inf, alpha = 0.06) +
  labs(x = NULL, y = NULL, colour = NULL,
       title = "ANES issue alignment and partisanship, face-to-face respondents, 1972-2016",
       subtitle = "Year effects from mixed models with pair random intercepts; shaded = post-2004 (K&M)") +
  theme_bw(base_size = 10)
ggsave("output/figures/16_anes_alignment_levels.pdf", p1, width = 8, height = 9)

tr <- melt(trends[, .(group, period, `Total trend in r` = slope_r, `Bounds moving` = bounds_part,
                      `Position within bounds` = position_part)],
            id.vars = c("group", "period"), variable.name = "part")
tr[, group := factor(group, levels = rev(grp_levels))]
p2 <- ggplot(tr, aes(value, group, fill = part)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  facet_wrap(~ period) +
  scale_fill_manual(values = c("grey30", "#D55E00", "#0072B2")) +
  labs(x = "Change in correlation per decade", y = NULL, fill = NULL,
       title = "Decomposing trends in ANES correlations: bounds moving vs tighter association",
       subtitle = "Mixed models with pair random intercepts and slopes; r = r* x r_max, so trend(r) ~ mean(r*) x trend(r_max) + mean(r_max) x trend(r*)") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave("output/figures/17_anes_alignment_decomposition.pdf", p2, width = 10, height = 6)
cat("\nSaved figures 16 and 17 and output/data/anes_alignment_decomposition.rds\n")
