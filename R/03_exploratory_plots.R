# 03_exploratory_plots.R
#
# Exploratory visualizations of BES and MC bounds data.
# Produces a set of figures saved to output/figures/.

library(dplyr)
library(ggplot2)
library(tidyr)

theme_set(theme_bw(base_size = 12))

bes <- readRDS("output/data/bes_bounds.rds")
mc  <- readRDS("output/data/mc_bounds.rds")

# ---- convenience ----
K_cols <- c("4"="#E41A1C", "5"="#FF7F00", "6"="#FFFF33",
            "7"="#4DAF4A", "10"="#377EB8", "11"="#984EA3")

save_fig <- function(p, name, w=7, h=5) {
  path <- paste0("output/figures/", name, ".pdf")
  ggsave(path, p, width=w, height=h)
  cat("Saved:", path, "\n")
  invisible(p)
}

# ===========================================================================
# 1. BES: distribution of (r_min, r_max) — scatter cloud
#
# TRANSFORM NOTE: r_min in [-1,0) and r_max in (0,1] both pile up hard against
# their Frechet-Hoeffding extreme (r_min near -1, r_max near +1 for ~7,000 of
# the 7,503 pairs) -- a linear scale buries almost the whole dataset in one
# corner. Empirically (1-r_max) and (1+r_min) are bounded strictly away from
# zero (min ~0.005, confirmed via output/data/bes_bounds.rds), so log10 of
# each "distance to its own extreme" is well-defined for every pair and
# spreads that corner out instead of overplotting it. The symmetry condition
# |r_min| = r_max becomes (1+r_min) = (1-r_max) algebraically, so the dashed
# diagonal is still the line y=x in this distance space -- it just flips which
# side counts as "more constrained" (see subtitle).
# ===========================================================================
p1 <- bes %>%
  mutate(
    K_pair   = factor(paste0(K1, "×", K2)),
    dist_min = 1 + r_min,   # distance of r_min from the lower extreme (-1)
    dist_max = 1 - r_max    # distance of r_max from the upper extreme (+1)
  ) %>%
  ggplot(aes(x = dist_min, y = dist_max)) +
  geom_point(alpha = 0.2, size = 0.6, colour = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40",
              linewidth = 0.5) +
  scale_x_log10(breaks = c(0.005, 0.01, 0.03, 0.1, 0.3, 0.7)) +
  scale_y_log10(breaks = c(0.005, 0.01, 0.03, 0.1, 0.3, 0.7)) +
  labs(
    title = "BES 2019: Attainable bounds for all 7,503 variable pairs",
    subtitle = paste(
      "Log scale = distance from each FH extreme (most pairs hug",
      "the boundary, so linear axes overplot that corner).",
      "Dashed = symmetry (|r_min|=r_max); points above =",
      "|r_min| > r_max (r_min sits closer to its bound)",
      sep = "\n"
    ),
    x = expression(1 + r[min] ~ "  (log scale)"),
    y = expression(1 - r[max] ~ "  (log scale)")
  )
save_fig(p1, "01_bes_rmin_rmax_cloud")

# ===========================================================================
# 2. BES: cloud coloured by K_max
# Same log-distance-to-bound transform as Figure 1 (see note above).
# ===========================================================================
p2 <- bes %>%
  mutate(
    K_max_f  = factor(K_max),
    dist_min = 1 + r_min,
    dist_max = 1 - r_max
  ) %>%
  ggplot(aes(x = dist_min, y = dist_max, colour = K_max_f)) +
  geom_point(alpha = 0.3, size = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey30",
              linewidth = 0.5) +
  scale_colour_manual(values = K_cols, name = "K (max)") +
  scale_x_log10(breaks = c(0.005, 0.01, 0.03, 0.1, 0.3, 0.7)) +
  scale_y_log10(breaks = c(0.005, 0.01, 0.03, 0.1, 0.3, 0.7)) +
  labs(
    title = "BES 2019: Bounds coloured by maximum number of categories",
    subtitle = "Log scale = distance from each FH extreme bound (see Fig 1)",
    x = expression(1 + r[min] ~ "  (log scale)"),
    y = expression(1 - r[max] ~ "  (log scale)")
  )
save_fig(p2, "02_bes_rmin_rmax_by_K")

# ===========================================================================
# 3. BES: distribution of C1 overall and by K
# ===========================================================================
p3a <- bes %>%
  ggplot(aes(x = C1)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
  labs(
    title = "BES 2019: Distribution of constraint severity (C1)",
    subtitle = "C1 = (|r_min| + r_max)/2; C1=1 means unconstrained",
    x = expression(C[1]),
    y = "Count"
  )
save_fig(p3a, "03a_bes_C1_histogram")

# C1 by K combination
bes_K_summary <- bes %>%
  mutate(K_pair = factor(paste0(pmin(K1,K2), "×", pmax(K1,K2)))) %>%
  group_by(K1, K2) %>%
  summarise(med_C1 = median(C1, na.rm=TRUE),
            q25 = quantile(C1, 0.25, na.rm=TRUE),
            q75 = quantile(C1, 0.75, na.rm=TRUE),
            n = n(), .groups="drop") %>%
  mutate(K_label = paste0(K1, "×", K2))

p3b <- bes %>%
  mutate(K_label = paste0(K1, "×", K2)) %>%
  ggplot(aes(x = reorder(K_label, K1 + K2), y = C1)) +
  geom_boxplot(fill = "steelblue", alpha = 0.6, outlier.size = 0.3,
               outlier.alpha = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  coord_flip() +
  labs(
    title = "BES 2019: Constraint severity by category combination",
    x = "K1 × K2",
    y = expression(C[1])
  )
save_fig(p3b, "03b_bes_C1_by_K", w=6, h=8)

# ===========================================================================
# 4. MC: distribution of C1 by K (same K pairs for comparability)
# ===========================================================================
mc_sym <- mc %>% filter(K1 == K2) %>%
  mutate(K_f = factor(K1, levels = c(4,5,6,7,10,11)))

p4 <- mc_sym %>%
  ggplot(aes(x = C1, fill = K_f, colour = K_f)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  scale_fill_manual(values = K_cols, name = "K") +
  scale_colour_manual(values = K_cols, name = "K") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  labs(
    title = "MC simulation: Distribution of C1 by number of categories (K×K pairs)",
    subtitle = "More categories → distribution shifts toward C1=1 (less constraint)",
    x = expression(C[1]),
    y = "Density"
  )
save_fig(p4, "04_mc_C1_by_K_density")

# ===========================================================================
# 5. Entropy vs C1 (MC data)
# ===========================================================================
p5a <- mc %>%
  sample_n(min(nrow(mc), 5000)) %>%
  mutate(K_max_f = factor(K_max)) %>%
  ggplot(aes(x = H_mean, y = C1, colour = K_max_f)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, colour = "black") +
  scale_colour_manual(values = K_cols, name = "K (max)") +
  labs(
    title = "MC: Mean entropy predicts constraint severity (C1)",
    subtitle = "Higher entropy (less concentrated) → larger C1 (less constrained)",
    x = expression(bar(H)(X,Y)),
    y = expression(C[1])
  )
save_fig(p5a, "05a_mc_entropy_vs_C1")

# 5b. Entropy vs asymmetry
p5b <- mc %>%
  sample_n(min(nrow(mc), 5000)) %>%
  mutate(K_max_f = factor(K_max)) %>%
  ggplot(aes(x = H_mean, y = asymmetry, colour = K_max_f)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, colour = "black") +
  scale_colour_manual(values = K_cols, name = "K (max)") +
  labs(
    title = "MC: Entropy vs bound asymmetry",
    subtitle = expression(paste("Asymmetry = |", r[min], "| - ", r[max],
                                "; positive = |r_min| > r_max")),
    x = expression(bar(H)(X,Y)),
    y = expression("|"*r[min]*"| - "*r[max])
  )
save_fig(p5b, "05b_mc_entropy_vs_asymmetry")

# 5c. TV distance vs asymmetry (mean TV across both marginals)
p5c <- mc %>%
  sample_n(min(nrow(mc), 5000)) %>%
  mutate(K_max_f = factor(K_max)) %>%
  ggplot(aes(x = TV_mean, y = asymmetry, colour = K_max_f)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, colour = "black") +
  scale_colour_manual(values = K_cols, name = "K (max)") +
  labs(
    title = "MC: Mean TV distance vs bound asymmetry",
    subtitle = "TV measures skewness of marginal relative to its reversal",
    x = expression(bar(TV)(X,Y)),
    y = expression("|"*r[min]*"| - "*r[max])
  )
save_fig(p5c, "05c_mc_TV_vs_asymmetry")

# ===========================================================================
# 6. BES vs MC: compare C1 distributions
# ===========================================================================
mc_C1   <- data.frame(C1 = mc$C1[!is.na(mc$C1)],   source = "MC simulation")
bes_C1  <- data.frame(C1 = bes$C1[!is.na(bes$C1)],  source = "BES 2019")
combined <- bind_rows(mc_C1, bes_C1)

p6 <- combined %>%
  ggplot(aes(x = C1, fill = source, colour = source)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("MC simulation"="steelblue", "BES 2019"="tomato"),
                    name = NULL) +
  scale_colour_manual(values = c("MC simulation"="steelblue", "BES 2019"="tomato"),
                      name = NULL) +
  labs(
    title = "C1 distribution: MC simulation vs BES 2019",
    subtitle = "BES distribution skews higher (less constrained) than random marginals",
    x = expression(C[1]),
    y = "Density"
  )
save_fig(p6, "06_C1_bes_vs_mc")

# ===========================================================================
# 7. BES: asymmetry distribution
# ===========================================================================
p7 <- bes %>%
  ggplot(aes(x = asymmetry)) +
  geom_histogram(bins = 60, fill = "steelblue", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
  annotate("text", x = 0.04, y = Inf, vjust = 1.5, hjust = 0,
           label = sprintf("|r_min|>r_max: %d (%.0f%%)",
                           sum(bes$asymmetry > 0, na.rm=TRUE),
                           100*mean(bes$asymmetry > 0, na.rm=TRUE)),
           colour = "red", size = 3.5) +
  annotate("text", x = -0.35, y = Inf, vjust = 1.5, hjust = 0,
           label = sprintf("r_max>|r_min|: %d (%.0f%%)",
                           sum(bes$asymmetry < 0, na.rm=TRUE),
                           100*mean(bes$asymmetry < 0, na.rm=TRUE)),
           colour = "steelblue4", size = 3.5) +
  labs(
    title = "BES 2019: Distribution of bound asymmetry",
    subtitle = expression(paste("Asymmetry = |", r[min], "| - ", r[max])),
    x = expression("|"*r[min]*"| - "*r[max]),
    y = "Count"
  )
save_fig(p7, "07_bes_asymmetry_histogram")

# ===========================================================================
# 8. BES: r_obs relative to bounds (where does observed r sit?)
# ===========================================================================
bes_pos <- bes %>%
  filter(!is.na(r_min), !is.na(r_max), !is.na(r_obs)) %>%
  mutate(
    # Normalise r_obs within [r_min, r_max] → [0, 1]
    r_norm = (r_obs - r_min) / (r_max - r_min),
    # Rescaled to [-1, 1]
    r_rescaled = 2 * r_norm - 1
  )

p8a <- bes_pos %>%
  ggplot(aes(x = r_obs, y = r_rescaled)) +
  geom_point(alpha = 0.2, size = 0.6, colour = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  labs(
    title = "BES 2019: Observed r vs rescaled r*",
    subtitle = paste(
      "r* = (r_obs - r_min)/(r_max - r_min) × 2 - 1;",
      "dashed = no rescaling effect",
      sep = "\n"
    ),
    x = expression(r[obs]),
    y = expression(r^"*"~"(rescaled to [-1,1])")
  )
save_fig(p8a, "08a_bes_r_vs_rescaled")

p8b <- bes_pos %>%
  ggplot(aes(x = r_rescaled - r_obs)) +
  geom_histogram(bins = 60, fill = "steelblue", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
  labs(
    title = "BES 2019: Rescaling adjustment (r* - r_obs)",
    subtitle = "Positive = rescaling inflates; negative = rescaling deflates",
    x = expression(r^"*" - r[obs]),
    y = "Count"
  )
save_fig(p8b, "08b_bes_rescaling_adjustment")

# ===========================================================================
# 9. BES: r_obs as fraction of r_max (for positive correlations)
# ===========================================================================
bes_pos_only <- bes_pos %>% filter(r_obs > 0)

p9 <- bes_pos_only %>%
  mutate(frac_of_max = r_obs / r_max) %>%
  ggplot(aes(x = frac_of_max)) +
  geom_histogram(bins = 50, fill = "tomato", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30") +
  labs(
    title = "BES 2019: Observed r as fraction of r_max (positive pairs only)",
    subtitle = "Many pairs have r_obs near their theoretical maximum",
    x = expression(r[obs] / r[max]),
    y = "Count"
  )
save_fig(p9, "09_bes_r_fraction_of_rmax")

# ===========================================================================
# Summary statistics
# ===========================================================================
cat("\n--- BES Summary ---\n")
cat(sprintf("N pairs: %d\n", nrow(bes)))
cat(sprintf("C1: mean=%.3f, median=%.3f, sd=%.3f\n",
            mean(bes$C1, na.rm=TRUE), median(bes$C1, na.rm=TRUE), sd(bes$C1, na.rm=TRUE)))
cat(sprintf("r_max: mean=%.3f, min=%.3f, max=%.3f\n",
            mean(bes$r_max, na.rm=TRUE), min(bes$r_max, na.rm=TRUE), max(bes$r_max, na.rm=TRUE)))
cat(sprintf("r_min: mean=%.3f, min=%.3f, max=%.3f\n",
            mean(bes$r_min, na.rm=TRUE), min(bes$r_min, na.rm=TRUE), max(bes$r_min, na.rm=TRUE)))
cat(sprintf("|r_min| > r_max: %d (%.1f%%)\n",
            sum(bes$asymmetry > 0, na.rm=TRUE),
            100*mean(bes$asymmetry > 0, na.rm=TRUE)))
cat(sprintf("Pairs where r_obs > 0.5 * r_max: %d (%.1f%%)\n",
            sum(bes$r_obs > 0 & bes$r_obs > 0.5 * bes$r_max, na.rm=TRUE),
            100 * mean(bes$r_obs > 0 & bes$r_obs > 0.5 * bes$r_max, na.rm=TRUE)))
cat(sprintf("Pairs where r_obs/r_max > 0.9: %d (%.1f%%)\n",
            sum(bes$r_obs > 0 & (bes$r_obs / bes$r_max) > 0.9, na.rm=TRUE),
            100 * mean(bes$r_obs > 0 & (bes$r_obs / bes$r_max) > 0.9, na.rm=TRUE)))

cat("\n--- MC Summary ---\n")
cat(sprintf("N pairs: %d\n", nrow(mc)))
cat(sprintf("C1: mean=%.3f, median=%.3f, sd=%.3f\n",
            mean(mc$C1, na.rm=TRUE), median(mc$C1, na.rm=TRUE), sd(mc$C1, na.rm=TRUE)))
cat(sprintf("|r_min| > r_max: %d (%.1f%%)\n",
            sum(mc$asymmetry > 0, na.rm=TRUE),
            100*mean(mc$asymmetry > 0, na.rm=TRUE)))
