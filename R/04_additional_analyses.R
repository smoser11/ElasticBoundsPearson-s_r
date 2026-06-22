# 04_additional_analyses.R
#
# Two additional analyses:
#   (A) The "spur" structure in the BES cloud: what drives pairs where r_max
#       is high but r_min is much less negative than expected?
#   (B) Asymmetry magnitude vs K: does |r_min - r_max| vary with category count
#       even when direction is unpredictable?

library(dplyr)
library(ggplot2)
library(tidyr)

theme_set(theme_bw(base_size = 12))

bes <- readRDS("output/data/bes_bounds.rds")
mc  <- readRDS("output/data/mc_bounds.rds")

save_fig <- function(p, name, w=7, h=5) {
  path <- paste0("output/figures/", name, ".pdf")
  ggsave(path, p, width=w, height=h)
  cat("Saved:", path, "\n")
  invisible(p)
}

# ===========================================================================
# (A) SPUR STRUCTURE: identify and characterise "spur" pairs
# A spur pair = r_max high BUT r_min less negative than expected
# Define: spur = r_max > 0.75  AND  r_min > -0.70  (right of -0.70)
# Compare to "main cloud" = r_max > 0.75 AND r_min <= -0.70
# ===========================================================================

bes_classified <- bes %>%
  filter(!is.na(r_max), !is.na(r_min)) %>%
  mutate(
    region = case_when(
      r_max > 0.75 & r_min > -0.70  ~ "Spur (high r_max, weak r_min)",
      r_max > 0.75 & r_min <= -0.70 ~ "Main cloud",
      TRUE                           ~ "Low r_max"
    ),
    region = factor(region, levels = c("Main cloud",
                                       "Spur (high r_max, weak r_min)",
                                       "Low r_max"))
  )

cat("=== SPUR ANALYSIS ===\n")
cat("Region counts:\n")
print(table(bes_classified$region))

# Compare characteristics across regions
spur_summary <- bes_classified %>%
  group_by(region) %>%
  summarise(
    n          = n(),
    med_K1     = median(K1),
    med_K2     = median(K2),
    med_C1     = round(median(C1, na.rm=TRUE), 3),
    med_H_mean = round(median(H_mean, na.rm=TRUE), 3),
    med_TV1    = round(median(TV1, na.rm=TRUE), 3),
    med_TV2    = round(median(TV2, na.rm=TRUE), 3),
    med_skew1  = round(median(skew1, na.rm=TRUE), 3),
    med_skew2  = round(median(skew2, na.rm=TRUE), 3),
    pct_asym   = round(100 * mean(asymmetry > 0, na.rm=TRUE), 1),
    .groups = "drop"
  )
cat("\nSpur summary:\n")
print(as.data.frame(spur_summary))

# Which K combinations produce spur pairs?
cat("\nSpur pairs K1 x K2 distribution:\n")
spur_K <- bes_classified %>%
  filter(region == "Spur (high r_max, weak r_min)") %>%
  count(K1, K2) %>% arrange(desc(n))
print(head(spur_K, 15))

# Spur vs main: TV distance comparison
# Key hypothesis: spur pairs have one variable with high TV (asymmetric marginal)
cat("\nTV distances (one variable much more skewed in spur?):\n")
spur_TV <- bes_classified %>%
  filter(region %in% c("Main cloud", "Spur (high r_max, weak r_min)")) %>%
  mutate(max_TV = pmax(TV1, TV2), min_TV = pmin(TV1, TV2))
print(
  spur_TV %>%
    group_by(region) %>%
    summarise(med_maxTV = round(median(max_TV), 3),
              med_minTV = round(median(min_TV), 3),
              med_diffTV = round(median(max_TV - min_TV), 3),
              .groups="drop")
)

# ----- Figure A1: cloud with regions highlighted -----
# TRANSFORM NOTE: this is the same (r_min, r_max) cloud as Figures 1-2, so it
# has the same problem -- on a linear scale the ~7,000 "main cloud" pairs that
# hug both FH extremes are crushed into one overplotted corner, which is
# exactly the region the spur (red) needs to be visually separated from. Apply
# the identical log10(distance-to-bound) transform used in Figs 1-2: x = 1 +
# r_min, y = 1 - r_max. The region classification itself is computed upstream
# on the raw r_min/r_max thresholds (r_max > 0.75, r_min > -0.70) and is
# unchanged here -- only the plotting coordinates differ -- and because the
# transform is monotonic per axis, those thresholds still draw as a straight
# vertical/horizontal split (just relocated to dist_min = 0.30 and dist_max =
# 0.25). The symmetry diagonal |r_min| = r_max again becomes the plain line
# y = x in distance space (see Fig 1 note), so spur pairs -- close to their
# r_max bound but comparatively far from their r_min bound -- now show up as
# a distinct band in the lower-right rather than smearing into the main cloud.
p_spur1 <- bes_classified %>%
  mutate(
    dist_min = 1 + r_min,   # distance of r_min from the lower extreme (-1)
    dist_max = 1 - r_max    # distance of r_max from the upper extreme (+1)
  ) %>%
  ggplot(aes(x = dist_min, y = dist_max, colour = region)) +
  geom_point(data = . %>% filter(region == "Main cloud"),
             alpha = 0.15, size = 0.5) +
  geom_point(data = . %>% filter(region == "Low r_max"),
             alpha = 0.25, size = 0.6) +
  geom_point(data = . %>% filter(region == "Spur (high r_max, weak r_min)"),
             alpha = 0.5, size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40",
              linewidth = 0.5) +
  scale_colour_manual(
    values = c("Main cloud" = "steelblue",
               "Spur (high r_max, weak r_min)" = "tomato",
               "Low r_max" = "grey60"),
    name = NULL
  ) +
  scale_x_log10(breaks = c(0.005, 0.01, 0.03, 0.1, 0.3, 0.7)) +
  scale_y_log10(breaks = c(0.005, 0.01, 0.03, 0.1, 0.3, 0.7)) +
  labs(
    title = "BES 2019: Identifying the 'spur' region",
    subtitle = paste(
      "Log scale = distance from each FH extreme (see Fig 1);",
      "red = spur (close to r_max bound, far from r_min bound);",
      "blue = main cloud",
      sep = "\n"
    ),
    x = expression(1 + r[min] ~ "  (log scale)"),
    y = expression(1 - r[max] ~ "  (log scale)")
  ) +
  theme(legend.position = "bottom")
save_fig(p_spur1, "10_bes_spur_highlighted")

# ----- Figure A2: TV distance comparison (spur vs main cloud) -----
spur_TV_long <- spur_TV %>%
  filter(region %in% c("Main cloud", "Spur (high r_max, weak r_min)")) %>%
  select(region, TV1, TV2) %>%
  pivot_longer(c(TV1, TV2), names_to = "variable", values_to = "TV") %>%
  mutate(variable = recode(variable, TV1 = "Variable 1", TV2 = "Variable 2"))

p_spur2 <- spur_TV_long %>%
  ggplot(aes(x = TV, fill = region, colour = region)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  facet_wrap(~variable) +
  scale_fill_manual(
    values = c("Main cloud" = "steelblue",
               "Spur (high r_max, weak r_min)" = "tomato"),
    name = NULL
  ) +
  scale_colour_manual(
    values = c("Main cloud" = "steelblue",
               "Spur (high r_max, weak r_min)" = "tomato"),
    name = NULL
  ) +
  labs(
    title = "BES 2019: TV distance distribution — spur vs main cloud",
    subtitle = "Spur pairs have one variable with markedly higher TV (more asymmetric marginal)",
    x = "TV distance (p vs reversed p)",
    y = "Density"
  ) +
  theme(legend.position = "bottom")
save_fig(p_spur2, "11_spur_tv_comparison", w=8, h=5)

# ----- Figure A3: asymmetry |r_min| - r_max within spur region -----
# Spur pairs by definition have r_min > -0.7 and r_max > 0.75
# so asymmetry = |r_min| - r_max = |r_min| - r_max
# For spur: r_max high, |r_min| low => asymmetry strongly NEGATIVE (r_max > |r_min|)
cat("\nAsymmetry in spur vs main cloud:\n")
bes_classified %>%
  filter(region != "Low r_max") %>%
  group_by(region) %>%
  summarise(med_asym = round(median(asymmetry, na.rm=TRUE), 3),
            mean_asym = round(mean(asymmetry, na.rm=TRUE), 3),
            .groups="drop") %>%
  print()

# ===========================================================================
# (B) ASYMMETRY MAGNITUDE vs K
# Even though direction is ~50/50, does |asymmetry| vary with K?
# Prediction: more categories → smaller |asymmetry| (bounds closer to symmetric)
# ===========================================================================

cat("\n=== ASYMMETRY MAGNITUDE vs K ===\n")

# BES: |asymmetry| by K
bes_abs <- bes %>%
  filter(!is.na(asymmetry)) %>%
  mutate(abs_asym = abs(asymmetry),
         K_min_f = factor(K_min),
         K_max_f = factor(K_max))

cat("BES |asymmetry| by K_min:\n")
bes_abs %>%
  group_by(K_min) %>%
  summarise(n=n(), med=round(median(abs_asym),3),
            q25=round(quantile(abs_asym,0.25),3),
            q75=round(quantile(abs_asym,0.75),3),
            .groups="drop") %>%
  print()

# MC: |asymmetry| by K (same K pairs)
mc_abs <- mc %>%
  filter(K1 == K2, !is.na(asymmetry)) %>%
  mutate(abs_asym = abs(asymmetry), K_f = factor(K1))

cat("\nMC |asymmetry| by K (K×K pairs):\n")
mc_abs %>%
  group_by(K1) %>%
  summarise(n=n(), med=round(median(abs_asym),3),
            q25=round(quantile(abs_asym,0.25),3),
            q75=round(quantile(abs_asym,0.75),3),
            .groups="drop") %>%
  print()

# MC: full K1 x K2 x N_obs grid -- not just the K1==K2 "square pairs" subset
# used above/in Figure 13. Now that 02_mc_simulation.R sweeps N_obs and
# records each draw's empirical Pearson's r (r_obs) alongside r_min/r_max,
# condition the |asymmetry| breakdown on all three of K1, K2, and N_obs, and
# report r_obs descriptive stats (mean, median, sd) for the same groups.
# This is a separate summary table -- it doesn't touch mc_abs, bes_abs, or
# any figure above/below.
cat("\nMC |asymmetry| and r_obs by K1, K2, N_obs (full grid):\n")
mc_grid <- mc %>%
  filter(!is.na(asymmetry), !is.na(r_obs)) %>%
  mutate(abs_asym = abs(asymmetry))

mc_grid_summary <- mc_grid %>%
  group_by(K1, K2, N_obs) %>%
  summarise(
    n          = n(),
    med_asym   = round(median(abs_asym), 3),
    q25_asym   = round(quantile(abs_asym, 0.25), 3),
    q75_asym   = round(quantile(abs_asym, 0.75), 3),
    mean_r_obs = round(mean(r_obs), 3),
    med_r_obs  = round(median(r_obs), 3),
    sd_r_obs   = round(sd(r_obs), 3),
    .groups = "drop"
  ) %>%
  arrange(K1, K2, N_obs)

print(as.data.frame(mc_grid_summary))

# ----- Figure B1: |asymmetry| by K — BES boxplot -----
# TRANSFORM NOTE: same abs_asym variable, same right-skew-against-zero shape
# as Figure 13 below (BES min ~1e-6, median ~0.028, max ~0.62, 24% of pairs
# < 0.01) -- on a linear y-axis the violins are squashed teardrops with
# almost all the visible shape sitting in the bottom few pixels. Apply the
# identical scales::pseudo_log_trans() fix used for the MC version in Fig 13
# so the low-asymmetry detail (most of the data) isn't buried.
K_cols <- c("4"="#E41A1C", "5"="#FF7F00", "6"="#D4C700",
            "7"="#4DAF4A", "10"="#377EB8", "11"="#984EA3")

p_asym_K_bes <- bes_abs %>%
  mutate(K_min_f = factor(K_min)) %>%
  ggplot(aes(x = K_min_f, y = abs_asym, fill = K_min_f)) +
  geom_violin(alpha = 0.6, linewidth = 0.4, draw_quantiles = 0.5) +
  scale_fill_manual(values = K_cols, guide = "none") +
  scale_y_continuous(
    trans  = scales::pseudo_log_trans(sigma = 0.001, base = 10),
    breaks = c(0, 0.001, 0.01, 0.1, 0.5)
  ) +
  labs(
    title = "BES 2019: Bound asymmetry magnitude by K",
    subtitle = paste(
      "Distribution of ||r_min| - r_max| by minimum number of",
      "categories; pseudo-log y-axis (linear near 0, log beyond)",
      sep = "\n"
    ),
    x = "Minimum K in pair",
    y = expression("|"*"|"*r[min]*"|"*" - "*r[max]*"|" ~ "  (pseudo-log scale)")
  )
save_fig(p_asym_K_bes, "12_bes_abs_asymmetry_by_K")

# ----- Figure B2: |asymmetry| by K — MC density -----
# TRANSFORM NOTE: abs_asym is heavily right-skewed (median ~0.045, max ~0.77)
# and -- unlike r_min/r_max in Figures 1-2 -- it can legitimately equal 0 (one
# exact 0 and ~100 near-zero floating-point cases out of 72,000 MC draws,
# where the sampled marginals happen to be almost perfectly symmetric). A
# literal log10(abs_asym) is undefined at 0 and would send those near-zero
# draws to roughly -15 on the log axis, wrecking the KDE bandwidth for
# everyone else. scales::pseudo_log_trans() is the standard fix: linear near
# zero, logarithmic once abs_asym exceeds sigma, so it still compresses the
# long right tail without blowing up at the boundary.
K_cols_v <- c("4"="#E41A1C", "5"="#FF7F00", "6"="#D4C700",
              "7"="#4DAF4A", "10"="#377EB8", "11"="#984EA3")

p_asym_K_mc <- mc_abs %>%
  ggplot(aes(x = abs_asym, fill = K_f, colour = K_f)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  scale_fill_manual(values = K_cols_v, name = "K") +
  scale_colour_manual(values = K_cols_v, name = "K") +
  scale_x_continuous(
    trans  = scales::pseudo_log_trans(sigma = 0.001, base = 10),
    breaks = c(0, 0.001, 0.01, 0.1, 0.5)
  ) +
  labs(
    title = "MC simulation: Bound asymmetry magnitude by K",
    subtitle = paste(
      "Distribution of ||r_min| - r_max| for K×K pairs (same K",
      "both variables); pseudo-log x-axis (linear near 0, log beyond)",
      sep = "\n"
    ),
    x = expression("|"*"|"*r[min]*"|"*" - "*r[max]*"|" ~ "  (pseudo-log scale)"),
    y = "Density"
  )
save_fig(p_asym_K_mc, "13_mc_abs_asymmetry_by_K_density")

# ----- Figure B3: |asymmetry| vs entropy, faceted by K (BES) -----
# TRANSFORM NOTE: same abs_asym boundary-at-zero issue as Figs 12-13 -- on a
# linear y-axis nearly all 7,503 points are crushed into a thin band near 0,
# which is exactly the region the loess trend needs to resolve. Same
# scales::pseudo_log_trans() fix.
p_asym_entropy <- bes_abs %>%
  mutate(K_max_f = factor(K_max)) %>%
  ggplot(aes(x = H_mean, y = abs_asym)) +
  geom_point(aes(colour = K_max_f), alpha = 0.2, size = 0.6) +
  geom_smooth(method = "loess", se = FALSE, colour = "black", linewidth = 0.9) +
  scale_colour_manual(values = K_cols, name = "K (max)") +
  scale_y_continuous(
    trans  = scales::pseudo_log_trans(sigma = 0.001, base = 10),
    breaks = c(0, 0.001, 0.01, 0.1, 0.5)
  ) +
  labs(
    title = "BES 2019: Entropy vs asymmetry magnitude",
    subtitle = paste(
      "Higher entropy (less concentrated) corresponds to",
      "smaller asymmetry magnitude; pseudo-log y-axis",
      "(linear near 0, log beyond)",
      sep = "\n"
    ),
    x = expression(bar(H)(X,Y)),
    y = expression("|"*"|"*r[min]*"|"*" - "*r[max]*"|" ~ "  (pseudo-log scale)")
  )
save_fig(p_asym_entropy, "14_bes_entropy_vs_abs_asymmetry")

cat("\nAll additional figures saved.\n")
