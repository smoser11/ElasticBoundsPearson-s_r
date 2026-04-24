# 05_hypothesis_testing.R
#
# Figure for §7.3: Hypothesis Testing
#
# Two-panel figure contrasting permutation null distributions for:
#   Pair A — constrained (small C1, concentrated marginals)
#   Pair B — unconstrained (large C1, spread marginals)
#
# Both pairs have the same observed r, and the same t-test p-value.
# The figure shows how a "significant" r means very different things
# depending on where it falls within the pair-specific attainable range.

library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(2024)
n <- 1000

theme_set(theme_bw(base_size = 11))

save_fig <- function(p, name, w = 8, h = 4.5) {
  path <- paste0("output/figures/", name, ".pdf")
  ggsave(path, p, width = w, height = h)
  cat("Saved:", path, "\n")
  invisible(p)
}

# ---------------------------------------------------------------------------
# Helper: compute rmax and rmin from count vectors
# ---------------------------------------------------------------------------
compute_rmax <- function(cx, cy) {
  x <- rep(seq_along(cx) - 1L, times = cx)
  y <- rep(seq_along(cy) - 1L, times = cy)
  cor(sort(x, decreasing = TRUE), sort(y, decreasing = TRUE))
}

compute_rmin <- function(cx, cy) {
  x <- rep(seq_along(cx) - 1L, times = cx)
  y <- rep(seq_along(cy) - 1L, times = cy)
  cor(sort(x, decreasing = TRUE), sort(y, decreasing = FALSE))
}

# ---------------------------------------------------------------------------
# Helper: generate individual-level data from count vectors with target r
# Uses a mixing approach: fraction alpha from comonotonic ordering, rest random
# ---------------------------------------------------------------------------
make_data <- function(cx, cy, alpha = 0.3) {
  x_full <- rep(seq_along(cx) - 1L, times = cx)
  y_full <- rep(seq_along(cy) - 1L, times = cy)
  # Comonotonic arrangement
  x_co <- sort(x_full, decreasing = TRUE)
  y_co <- sort(y_full, decreasing = TRUE)
  # Independently shuffled (null) arrangement
  x_null <- sample(x_full)
  y_null <- sample(y_full)
  # Mix: alpha fraction from comonotonic, rest from null
  use_co <- runif(n) < alpha
  x_obs <- ifelse(use_co, x_co, x_null)
  y_obs <- ifelse(use_co, y_co, sample(y_full))
  list(x = x_obs, y = y_obs)
}

# ---------------------------------------------------------------------------
# Helper: permutation test
# ---------------------------------------------------------------------------
run_permtest <- function(x, y, n_perm = 5000) {
  r_obs <- cor(x, y)
  r_perm <- replicate(n_perm, cor(x, sample(y)))
  p_perm <- mean(abs(r_perm) >= abs(r_obs))
  list(r_obs = r_obs, r_perm = r_perm, p_perm = p_perm)
}

# ---------------------------------------------------------------------------
# Pair A: CONSTRAINED
# K=4 x K=4, both with very skewed (concentrated) marginals
# One left-skewed, one right-skewed → rmax well below 1
# ---------------------------------------------------------------------------
cx_A <- c(30L, 70L, 200L, 700L)   # K=4, heavy mass at category 3 (right-skewed)
cy_A <- c(700L, 200L, 70L, 30L)   # K=4, heavy mass at category 0 (left-skewed)

rmax_A <- compute_rmax(cx_A, cy_A)
rmin_A <- compute_rmin(cx_A, cy_A)
C1_A   <- (abs(rmin_A) + rmax_A) / 2

cat(sprintf("Pair A: rmin=%.3f  rmax=%.3f  C1=%.3f\n", rmin_A, rmax_A, C1_A))

dat_A  <- make_data(cx_A, cy_A, alpha = 0.55)
res_A  <- run_permtest(dat_A$x, dat_A$y)

t_A <- res_A$r_obs * sqrt(n - 2) / sqrt(1 - res_A$r_obs^2)
p_t_A <- 2 * pt(-abs(t_A), df = n - 2)

cat(sprintf("Pair A: r_obs=%.3f  t=%.2f  p(t)=%.4f  p(perm)=%.4f\n",
            res_A$r_obs, t_A, p_t_A, res_A$p_perm))
cat(sprintf("Pair A: r_obs/rmax = %.1f%%\n", 100 * res_A$r_obs / rmax_A))

# ---------------------------------------------------------------------------
# Pair B: UNCONSTRAINED
# K=7 x K=7, roughly uniform marginals → rmax close to 1
# ---------------------------------------------------------------------------
cx_B <- c(120L, 150L, 160L, 160L, 160L, 130L, 120L)   # K=7, roughly uniform
cy_B <- c(115L, 145L, 165L, 165L, 165L, 125L, 120L)
# Adjust to exactly n=1000
cx_B[4] <- cx_B[4] + (n - sum(cx_B))
cy_B[4] <- cy_B[4] + (n - sum(cy_B))

rmax_B <- compute_rmax(cx_B, cy_B)
rmin_B <- compute_rmin(cx_B, cy_B)
C1_B   <- (abs(rmin_B) + rmax_B) / 2

cat(sprintf("Pair B: rmin=%.3f  rmax=%.3f  C1=%.3f\n", rmin_B, rmax_B, C1_B))

# Tune alpha so r_obs ≈ r_obs_A
dat_B  <- make_data(cx_B, cy_B, alpha = 0.14)
res_B  <- run_permtest(dat_B$x, dat_B$y)

t_B <- res_B$r_obs * sqrt(n - 2) / sqrt(1 - res_B$r_obs^2)
p_t_B <- 2 * pt(-abs(t_B), df = n - 2)

cat(sprintf("Pair B: r_obs=%.3f  t=%.2f  p(t)=%.4f  p(perm)=%.4f\n",
            res_B$r_obs, t_B, p_t_B, res_B$p_perm))
cat(sprintf("Pair B: r_obs/rmax = %.1f%%\n", 100 * res_B$r_obs / rmax_B))

# ---------------------------------------------------------------------------
# Build figure
# ---------------------------------------------------------------------------
make_panel <- function(res, rmin, rmax, r_obs, label, C1) {
  df_perm <- data.frame(r = res$r_perm)
  p_val   <- res$p_perm

  ggplot(df_perm, aes(x = r)) +
    geom_histogram(aes(y = after_stat(density)), bins = 60,
                   fill = "grey70", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = rmax, linetype = "dashed", colour = "#D55E00", linewidth = 0.8) +
    geom_vline(xintercept = rmin, linetype = "dashed", colour = "#0072B2", linewidth = 0.8) +
    geom_vline(xintercept = r_obs, colour = "black",   linewidth = 1.1) +
    annotate("text", x = rmax + 0.01, y = Inf, label = expression(r[max]),
             hjust = 0, vjust = 1.5, colour = "#D55E00", size = 3.2) +
    annotate("text", x = rmin - 0.01, y = Inf, label = expression(r[min]),
             hjust = 1, vjust = 1.5, colour = "#0072B2", size = 3.2) +
    annotate("text", x = r_obs + 0.01, y = Inf,
             label = sprintf("r[obs] == %.2f", r_obs),
             hjust = 0, vjust = 3.2, size = 3.2, parse = TRUE) +
    annotate("text", x = (rmin + rmax) / 2, y = -Inf,
             label = sprintf("italic(p)[perm] == '%.3f'", p_val),
             hjust = 0.5, vjust = -0.5, size = 3.2, parse = TRUE) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
    labs(
      title = label,
      subtitle = bquote(r[min] == .(round(rmin,2)) ~~ r[max] == .(round(rmax,2)) ~~ C[1] == .(round(C1,2))),
      x = expression(r ~ "(permutation null)"),
      y = "Density"
    ) +
    theme(plot.subtitle = element_text(size = 9))
}

p_A <- make_panel(res_A, rmin_A, rmax_A, res_A$r_obs,
                  "Pair A: Constrained marginals", C1_A)
p_B <- make_panel(res_B, rmin_B, rmax_B, res_B$r_obs,
                  "Pair B: Unconstrained marginals", C1_B)

fig <- p_A + p_B +
  plot_annotation(
    title    = "Permutation null distributions for two contrasting variable pairs",
    subtitle = paste0(
      "Vertical lines: r_min (blue dashed), r_max (orange dashed), r_obs (black solid).\n",
      "Both pairs have similar r_obs and similar t-test p-values, ",
      "yet r_obs occupies very different positions within the attainable range."
    ),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9,  colour = "grey30")
    )
  )

save_fig(fig, "15_permutation_null_comparison", w = 10, h = 5)
cat("\nDone.\n")
