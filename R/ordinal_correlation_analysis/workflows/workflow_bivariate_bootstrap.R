# ================================================================
# BOOTSTRAP UNCERTAINTY ANALYSIS WORKFLOW  
# Modular sections - run independently or sequentially
# 
# Purpose: Quantify uncertainty in correlation bounds estimates using
#          bootstrap methods and generate confidence regions
# ================================================================

library(here)
library(dplyr)
library(ggplot2)

# Load cache management utilities
source(here("R", "ordinal_correlation_analysis", "utilities", "cache_management.R"))

cat("🔄 BOOTSTRAP UNCERTAINTY ANALYSIS WORKFLOW\n")
cat("===========================================\n\n")

# ================================================================
# CONFIGURATION PARAMETERS
# ================================================================

params <- list(
  B = 1000,                    # Number of bootstrap samples
  confidence_levels = c(0.90, 0.95, 0.99),  # Confidence levels to compute
  store_full_samples = TRUE,   # Whether to store all bootstrap samples  
  force_regenerate = FALSE,    # Set TRUE to ignore existing cache
  
  # Analysis parameters
  subsample_size = 500,        # Subsample size for computational efficiency
  parallel = FALSE             # Use parallel processing (if available)
)

cat("📋 Configuration:\n")
cat("   Bootstrap samples:", params$B, "\n")
cat("   Confidence levels:", paste(params$confidence_levels, collapse = ", "), "\n")
cat("   Store full samples:", params$store_full_samples, "\n")
cat("   Subsample size:", params$subsample_size, "\n")
cat("   Force regenerate:", params$force_regenerate, "\n\n")

# ================================================================
# SECTION 1: LOAD REQUIRED DATA
# ================================================================

cat("=== SECTION 1: LOAD SIMULATION DATA ===\n")

# Check for Monte Carlo simulation data
mc_files <- list.files(here("R", "ordinal_correlation_analysis", "data", "raw"), 
                      pattern = "MCsim.*\\.rds", full.names = TRUE)

if (length(mc_files) == 0) {
  stop("No Monte Carlo simulation data found. Please run workflow_bivariate_mc.R first.")
}

# Use the most recent MC simulation file
mc_file <- mc_files[which.max(file.mtime(mc_files))]
cat("📂 Loading MC data:", basename(mc_file), "\n")
mc_data <- readRDS(mc_file)

cat("✅ Simulation data loaded\n")
cat("   Configurations:", length(mc_data), "\n")
if (length(mc_data) > 0 && !is.null(mc_data[[1]]$contingency_tables)) {
  total_tables <- sum(sapply(mc_data, function(x) length(x$contingency_tables)))
  cat("   Total tables available:", total_tables, "\n")
}
cat("\n")

# ================================================================
# SECTION 2: BOOTSTRAP BOUNDS ESTIMATION
# ================================================================

cat("=== SECTION 2: BOOTSTRAP BOUNDS ESTIMATION ===\n")

# Cache file for bootstrap results
bootstrap_cache_file <- here("R", "ordinal_correlation_analysis", "data", "processed",
                            generate_cache_filename("bootstrap_results", 
                                                   params[c("B", "subsample_size")]))

bootstrap_results <- cache_or_compute(
  cache_file = bootstrap_cache_file,
  force_regenerate = params$force_regenerate,
  compute_func = function() {
    cat("Running bootstrap analysis on simulation data...\n")
    
    # Source required bootstrap functions
    source(here("R", "ordinal_correlation_analysis", "1_bivariate_ordcats_correlation", 
               "1_rmin_rmax_rhat", "bootstrapJoint_r_minmax.R"))
    source(here("R", "correlation_bounds_core.R"))
    
    # Select subset of simulation data for computational efficiency
    if (params$subsample_size < length(mc_data)) {
      set.seed(123)  # For reproducibility
      selected_configs <- sample(seq_along(mc_data), params$subsample_size)
      mc_subset <- mc_data[selected_configs]
      cat("Using", length(mc_subset), "configurations (subsampled for efficiency)\n")
    } else {
      mc_subset <- mc_data
      cat("Using all", length(mc_subset), "configurations\n")
    }
    
    # Run bootstrap analysis
    full_bootstrap_results <- bootstrap_all_configurations(
      mc_subset, 
      B = params$B, 
      store_full_samples = params$store_full_samples
    )
    
    cat("Bootstrap analysis completed\n")
    return(full_bootstrap_results)
  }
)

cat("✅ Bootstrap estimation ready\n")
cat("   Bootstrap samples per config:", params$B, "\n")
cat("   Configurations analyzed:", length(bootstrap_results), "\n\n")

# ================================================================
# SECTION 3: CONFIDENCE INTERVAL COMPUTATION
# ================================================================

cat("=== SECTION 3: CONFIDENCE INTERVALS ===\n")

# Compute confidence intervals for each configuration
ci_results <- list()

for (i in seq_along(bootstrap_results)) {
  config_result <- bootstrap_results[[i]]
  
  # Extract bootstrap samples
  if (!is.null(config_result$bootstrap_samples)) {
    r_min_samples <- config_result$bootstrap_samples$r_min
    r_max_samples <- config_result$bootstrap_samples$r_max
    
    # Compute confidence intervals for each level
    ci_data <- data.frame()
    
    for (level in params$confidence_levels) {
      alpha <- 1 - level
      
      # Compute quantiles
      r_min_ci <- quantile(r_min_samples, c(alpha/2, 1-alpha/2), na.rm = TRUE)
      r_max_ci <- quantile(r_max_samples, c(alpha/2, 1-alpha/2), na.rm = TRUE)
      
      ci_data <- rbind(ci_data, data.frame(
        config_id = i,
        confidence_level = level,
        r_min_lower = r_min_ci[1],
        r_min_upper = r_min_ci[2], 
        r_max_lower = r_max_ci[1],
        r_max_upper = r_max_ci[2],
        r_min_mean = mean(r_min_samples, na.rm = TRUE),
        r_max_mean = mean(r_max_samples, na.rm = TRUE)
      ))
    }
    
    ci_results[[i]] <- ci_data
  }
}

# Combine all CI results
all_ci_results <- do.call(rbind, ci_results)

cat("📊 Confidence intervals computed:\n")
cat("   Configurations with CIs:", length(unique(all_ci_results$config_id)), "\n") 
cat("   Confidence levels:", paste(unique(all_ci_results$confidence_level), collapse = ", "), "\n")
cat("\n")

# Save CI results
ci_cache_file <- here("R", "ordinal_correlation_analysis", "output", "reports",
                     generate_cache_filename("confidence_intervals", params[c("B")], "rds"))
saveRDS(all_ci_results, ci_cache_file)
cat("💾 Confidence intervals saved:", basename(ci_cache_file), "\n\n")

# ================================================================
# SECTION 4: UNCERTAINTY VISUALIZATIONS
# ================================================================

cat("=== SECTION 4: UNCERTAINTY VISUALIZATIONS ===\n")

figures_dir <- here("R", "ordinal_correlation_analysis", "output", "figures") 
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Plot 1: Bootstrap confidence intervals (95% level)
ci_95 <- all_ci_results %>% filter(confidence_level == 0.95)

plot1 <- ggplot(ci_95, aes(x = config_id)) +
  geom_ribbon(aes(ymin = r_min_lower, ymax = r_min_upper), alpha = 0.3, fill = "red") +
  geom_ribbon(aes(ymin = r_max_lower, ymax = r_max_upper), alpha = 0.3, fill = "blue") +
  geom_line(aes(y = r_min_mean), color = "red", size = 0.8) +
  geom_line(aes(y = r_max_mean), color = "blue", size = 0.8) +
  labs(
    title = "Bootstrap Confidence Intervals for Correlation Bounds",
    subtitle = "95% confidence intervals with bootstrap means (Red: r_min, Blue: r_max)",
    x = "Configuration ID",
    y = "Correlation Value"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

plot1_file <- file.path(figures_dir, "bootstrap_confidence_intervals.png")
ggsave(plot1_file, plot1, width = 12, height = 8, dpi = 300)

# Plot 2: Uncertainty magnitude analysis
ci_95 <- ci_95 %>%
  mutate(
    r_min_uncertainty = r_min_upper - r_min_lower,
    r_max_uncertainty = r_max_upper - r_max_lower,
    total_uncertainty = r_min_uncertainty + r_max_uncertainty
  )

plot2 <- ggplot(ci_95) +
  geom_point(aes(x = r_min_uncertainty, y = r_max_uncertainty), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Bootstrap Uncertainty: r_min vs r_max",
    subtitle = "Width of 95% confidence intervals",
    x = "r_min Uncertainty (CI Width)",  
    y = "r_max Uncertainty (CI Width)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

plot2_file <- file.path(figures_dir, "bootstrap_uncertainty_comparison.png")
ggsave(plot2_file, plot2, width = 10, height = 8, dpi = 300)

# Plot 3: Distribution of uncertainty magnitudes
plot3 <- ggplot(ci_95, aes(x = total_uncertainty)) +
  geom_histogram(bins = 30, fill = "lightblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = mean(ci_95$total_uncertainty), linetype = "dashed", color = "red") +
  labs(
    title = "Distribution of Total Bootstrap Uncertainty",
    subtitle = paste("Mean total uncertainty:", round(mean(ci_95$total_uncertainty), 3)),
    x = "Total Uncertainty (Sum of r_min and r_max CI widths)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

plot3_file <- file.path(figures_dir, "bootstrap_uncertainty_distribution.png")
ggsave(plot3_file, plot3, width = 10, height = 6, dpi = 300)

cat("📈 Uncertainty visualizations generated:\n")
cat("   1. Confidence intervals:", basename(plot1_file), "\n")
cat("   2. Uncertainty comparison:", basename(plot2_file), "\n") 
cat("   3. Uncertainty distribution:", basename(plot3_file), "\n\n")

# ================================================================
# SECTION 5: BOOTSTRAP SUMMARY TABLES
# ================================================================

cat("=== SECTION 5: BOOTSTRAP SUMMARY TABLES ===\n")

tables_dir <- here("R", "ordinal_correlation_analysis", "output", "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# Create uncertainty summary by confidence level
uncertainty_summary <- all_ci_results %>%
  mutate(
    r_min_width = r_min_upper - r_min_lower,
    r_max_width = r_max_upper - r_max_lower
  ) %>%
  group_by(confidence_level) %>%
  summarise(
    n_configs = n(),
    mean_r_min_width = round(mean(r_min_width, na.rm = TRUE), 4),
    mean_r_max_width = round(mean(r_max_width, na.rm = TRUE), 4),
    median_r_min_width = round(median(r_min_width, na.rm = TRUE), 4),
    median_r_max_width = round(median(r_max_width, na.rm = TRUE), 4),
    .groups = "drop"
  )

# Save uncertainty summary
uncertainty_file <- file.path(tables_dir, "bootstrap_uncertainty_summary.csv")
write.csv(uncertainty_summary, uncertainty_file, row.names = FALSE)

# Create detailed bootstrap statistics
bootstrap_stats <- data.frame(
  Statistic = c(
    "Number of Bootstrap Samples",
    "Configurations Analyzed",
    "Mean r_min CI Width (95%)", 
    "Mean r_max CI Width (95%)",
    "Median r_min CI Width (95%)",
    "Median r_max CI Width (95%)",
    "Mean Total Uncertainty (95%)",
    "Max Total Uncertainty (95%)"
  ),
  Value = c(
    params$B,
    length(unique(all_ci_results$config_id)),
    round(mean(ci_95$r_min_uncertainty, na.rm = TRUE), 4),
    round(mean(ci_95$r_max_uncertainty, na.rm = TRUE), 4), 
    round(median(ci_95$r_min_uncertainty, na.rm = TRUE), 4),
    round(median(ci_95$r_max_uncertainty, na.rm = TRUE), 4),
    round(mean(ci_95$total_uncertainty, na.rm = TRUE), 4),
    round(max(ci_95$total_uncertainty, na.rm = TRUE), 4)
  )
)

bootstrap_stats_file <- file.path(tables_dir, "bootstrap_statistics.csv")
write.csv(bootstrap_stats, bootstrap_stats_file, row.names = FALSE)

cat("📊 Bootstrap summary tables generated:\n")
cat("   1. Uncertainty by confidence level:", basename(uncertainty_file), "\n")
cat("   2. Bootstrap statistics:", basename(bootstrap_stats_file), "\n\n")

# ================================================================
# WORKFLOW COMPLETION SUMMARY
# ================================================================

cat("🎉 BOOTSTRAP UNCERTAINTY WORKFLOW COMPLETE!\n")
cat("============================================\n")
cat("📁 Generated files:\n")
cat("   Bootstrap results:", basename(bootstrap_cache_file), "\n")
cat("   Confidence intervals:", basename(ci_cache_file), "\n")
cat("   Figures: 3 uncertainty visualization files\n")
cat("   Tables: 2 uncertainty summary files\n\n")

cat("📊 Uncertainty analysis summary:\n")
cat("   Bootstrap samples:", params$B, "\n")
cat("   Configurations analyzed:", length(unique(all_ci_results$config_id)), "\n")
cat("   Mean 95% CI width (r_min):", round(mean(ci_95$r_min_uncertainty, na.rm = TRUE), 4), "\n")
cat("   Mean 95% CI width (r_max):", round(mean(ci_95$r_max_uncertainty, na.rm = TRUE), 4), "\n\n")

cat("🔄 Integration with other workflows:\n")
cat("   - Results can be used to assess bounds reliability\n")
cat("   - Uncertainty estimates inform interpretation of BES analysis\n")
cat("   - Bootstrap samples available for advanced confidence region analysis\n\n")