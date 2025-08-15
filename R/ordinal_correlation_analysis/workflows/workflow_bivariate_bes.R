# ================================================================
# BES DATA ANALYSIS WORKFLOW
# Modular sections - run independently or sequentially
# 
# Purpose: Analyze British Election Study 2019 data to compute theoretical
#          correlation bounds for real-world ordinal variable pairs
# ================================================================

library(here)
library(dplyr)
library(ggplot2)
library(readstata13)

# Load cache management utilities
source(here("R", "ordinal_correlation_analysis", "utilities", "cache_management.R"))

cat("🇬🇧 BES DATA ANALYSIS WORKFLOW\n")
cat("==============================\n\n")

# ================================================================
# CONFIGURATION PARAMETERS
# ================================================================

params <- list(
  nsim = 2000,                 # Permutation simulations for bounds
  confidence_level = 0.95,     # Confidence level for intervals
  force_regenerate = FALSE,    # Set TRUE to ignore existing cache
  
  # Data filtering parameters
  min_nobs = 100,              # Minimum sample size for inclusion
  max_missing_prop = 0.1       # Maximum proportion missing data
)

cat("📋 Configuration:\n")
cat("   Permutation simulations:", params$nsim, "\n")
cat("   Confidence level:", params$confidence_level, "\n")
cat("   Minimum sample size:", params$min_nobs, "\n")
cat("   Force regenerate:", params$force_regenerate, "\n\n")

# ================================================================
# SECTION 1: LOAD AND VALIDATE BES DATA
# ================================================================

cat("=== SECTION 1: BES DATA LOADING ===\n")

# Load BES data
bes_data_file <- here("R", "ordinal_correlation_analysis", "data", "processed", 
                     "correlation and other data about pairs of BES2019 variables.dta")

if (!file.exists(bes_data_file)) {
  stop("BES data file not found: ", bes_data_file, 
       "\nPlease ensure the data file is in the correct location.")
}

cat("📂 Loading BES data:", basename(bes_data_file), "\n")
bes_data <- read.dta13(bes_data_file)

cat("✅ BES data loaded successfully\n")
cat("   Total variable pairs:", nrow(bes_data), "\n")
cat("   Columns:", ncol(bes_data), "\n")

# Validate required columns
required_cols <- c("var1", "var2", "corr", "nobs", "var1cats", "var2cats")
missing_cols <- setdiff(required_cols, names(bes_data))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Data quality checks
cat("📊 Data quality summary:\n")
cat("   Mean correlation:", round(mean(bes_data$corr, na.rm = TRUE), 3), "\n")
cat("   Mean sample size:", round(mean(bes_data$nobs, na.rm = TRUE), 0), "\n")
cat("   Variable pairs with nobs >=", params$min_nobs, ":", 
    sum(bes_data$nobs >= params$min_nobs, na.rm = TRUE), "\n")
cat("\n")

# ================================================================
# SECTION 2: BES BOUNDS COMPUTATION
# ================================================================

cat("=== SECTION 2: BES BOUNDS COMPUTATION ===\n")

# Cache file for bounds computation
bounds_cache_file <- here("R", "ordinal_correlation_analysis", "output", "reports",
                         generate_cache_filename("bes_bounds", params[c("nsim", "min_nobs")]))

bes_bounds_data <- cache_or_compute(
  cache_file = bounds_cache_file,
  force_regenerate = params$force_regenerate,
  compute_func = function() {
    cat("Computing correlation bounds for BES data...\n")
    
    # Source required functions
    source(here("R", "correlation_bounds_core.R"))
    source(here("R", "ordinal_correlation_analysis", "1_bivariate_ordcats_correlation", 
               "1_rmin_rmax_rhat", "2_bes_illustrative_example", "bes_data_analysis.R"))
    
    # Filter data based on quality criteria
    filtered_data <- bes_data %>%
      filter(nobs >= params$min_nobs,
             !is.na(corr),
             !is.na(var1cats),
             !is.na(var2cats))
    
    cat("Analyzing", nrow(filtered_data), "variable pairs (after filtering)\n")
    
    # Compute bounds for all pairs
    bounds_results <- analyze_all_bes_bounds(
      filtered_data, 
      nsim = params$nsim,
      progress = TRUE
    )
    
    cat("Bounds computation completed\n")
    return(bounds_results)
  }
)

cat("✅ BES bounds computation ready\n")
cat("   Variable pairs analyzed:", nrow(bes_bounds_data), "\n")
if ("r_min" %in% names(bes_bounds_data)) {
  cat("   Mean r_min:", round(mean(bes_bounds_data$r_min, na.rm = TRUE), 3), "\n")
  cat("   Mean r_max:", round(mean(bes_bounds_data$r_max, na.rm = TRUE), 3), "\n")
  cat("   Mean bounds range:", round(mean(bes_bounds_data$r_max - bes_bounds_data$r_min, na.rm = TRUE), 3), "\n")
}
cat("\n")

# ================================================================
# SECTION 3: BES BOUNDS ANALYSIS AND SUMMARY
# ================================================================

cat("=== SECTION 3: BES BOUNDS ANALYSIS ===\n")

# Add derived variables for analysis
bes_analysis <- bes_bounds_data %>%
  mutate(
    bounds_range = r_max - r_min,
    bounds_asymmetry = r_max + r_min,
    rescaled_r = (observed_r - r_min) / (r_max - r_min) * 2 - 1,
    bounds_utilization = abs(observed_r) / pmax(abs(r_min), abs(r_max)),
    is_symmetric = abs(bounds_asymmetry) < 0.01
  )

# Generate summary statistics
cat("📊 BES Bounds Analysis Summary:\n")
cat("   Variable pairs:", nrow(bes_analysis), "\n")
cat("   Mean observed correlation:", round(mean(bes_analysis$observed_r, na.rm = TRUE), 3), "\n")
cat("   Mean theoretical range:", round(mean(bes_analysis$bounds_range, na.rm = TRUE), 3), "\n")
cat("   Proportion with symmetric bounds:", round(mean(bes_analysis$is_symmetric, na.rm = TRUE), 3), "\n")
cat("   Mean bounds utilization:", round(mean(bes_analysis$bounds_utilization, na.rm = TRUE), 3), "\n")
cat("\n")

# Save analysis results
analysis_cache_file <- here("R", "ordinal_correlation_analysis", "output", "reports",
                           generate_cache_filename("bes_analysis", params[c("nsim")], "rds"))
saveRDS(bes_analysis, analysis_cache_file)
cat("💾 BES analysis saved:", basename(analysis_cache_file), "\n\n")

# ================================================================
# SECTION 4: BES VISUALIZATIONS
# ================================================================

cat("=== SECTION 4: BES VISUALIZATIONS ===\n")

figures_dir <- here("R", "ordinal_correlation_analysis", "output", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Plot 1: BES Bounds landscape
plot1 <- ggplot(bes_analysis, aes(x = r_min, y = r_max)) +
  geom_point(alpha = 0.5, color = "darkblue", size = 0.8) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray", alpha = 0.7) +
  labs(
    title = "Correlation Bounds Landscape: British Election Study 2019",
    subtitle = paste("Analysis of", nrow(bes_analysis), "ordinal variable pairs"),
    x = "Minimum Possible Correlation (r_min)",
    y = "Maximum Possible Correlation (r_max)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

plot1_file <- file.path(figures_dir, "bes_bounds_landscape.png")
ggsave(plot1_file, plot1, width = 12, height = 8, dpi = 300)

# Plot 2: Observed vs theoretical bounds comparison
plot2 <- ggplot(bes_analysis, aes(x = bounds_range, y = abs(observed_r))) +
  geom_point(alpha = 0.5, color = "forestgreen") +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  labs(
    title = "Observed Correlation vs Theoretical Range",
    subtitle = "Relationship between bounds width and observed correlation magnitude",
    x = "Theoretical Bounds Range (r_max - r_min)",
    y = "Absolute Observed Correlation"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

plot2_file <- file.path(figures_dir, "bes_observed_vs_bounds.png") 
ggsave(plot2_file, plot2, width = 10, height = 8, dpi = 300)

# Plot 3: Asymmetry analysis
plot3 <- ggplot(bes_analysis, aes(x = bounds_asymmetry)) +
  geom_histogram(bins = 60, fill = "orange", alpha = 0.7, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = mean(bes_analysis$bounds_asymmetry, na.rm = TRUE), 
             linetype = "dashed", color = "blue", size = 1) +
  labs(
    title = "Distribution of Bounds Asymmetry: BES 2019",
    subtitle = paste("Asymmetry = r_max + r_min;", 
                     round(100 * mean(bes_analysis$is_symmetric, na.rm = TRUE), 1), 
                     "% have symmetric bounds"),
    x = "Bounds Asymmetry",
    y = "Number of Variable Pairs"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

plot3_file <- file.path(figures_dir, "bes_bounds_asymmetry.png")
ggsave(plot3_file, plot3, width = 10, height = 6, dpi = 300)

# Plot 4: Rescaled correlations comparison  
plot4 <- ggplot(bes_analysis, aes(x = observed_r, y = rescaled_r)) +
  geom_point(alpha = 0.4, color = "purple") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Original vs Rescaled Correlations",
    subtitle = "Linear rescaling to [-1, 1] range using theoretical bounds",
    x = "Original Observed Correlation",
    y = "Rescaled Correlation"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

plot4_file <- file.path(figures_dir, "bes_original_vs_rescaled.png")
ggsave(plot4_file, plot4, width = 10, height = 8, dpi = 300)

cat("📈 BES visualizations generated:\n")
cat("   1. Bounds landscape:", basename(plot1_file), "\n")
cat("   2. Observed vs bounds:", basename(plot2_file), "\n")
cat("   3. Asymmetry distribution:", basename(plot3_file), "\n")
cat("   4. Original vs rescaled:", basename(plot4_file), "\n\n")

# ================================================================
# SECTION 5: BES SUMMARY TABLES
# ================================================================

cat("=== SECTION 5: BES SUMMARY TABLES ===\n")

tables_dir <- here("R", "ordinal_correlation_analysis", "output", "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# Create comprehensive summary table
bes_summary_stats <- data.frame(
  Statistic = c(
    "Number of Variable Pairs",
    "Mean Observed Correlation",
    "SD Observed Correlation", 
    "Mean r_min",
    "SD r_min",
    "Mean r_max", 
    "SD r_max",
    "Mean Bounds Range",
    "SD Bounds Range",
    "Mean Bounds Asymmetry",
    "Proportion Symmetric Bounds",
    "Mean Sample Size",
    "Median Sample Size"
  ),
  Value = c(
    nrow(bes_analysis),
    round(mean(bes_analysis$observed_r, na.rm = TRUE), 4),
    round(sd(bes_analysis$observed_r, na.rm = TRUE), 4),
    round(mean(bes_analysis$r_min, na.rm = TRUE), 4),
    round(sd(bes_analysis$r_min, na.rm = TRUE), 4),
    round(mean(bes_analysis$r_max, na.rm = TRUE), 4),
    round(sd(bes_analysis$r_max, na.rm = TRUE), 4),
    round(mean(bes_analysis$bounds_range, na.rm = TRUE), 4),
    round(sd(bes_analysis$bounds_range, na.rm = TRUE), 4),
    round(mean(bes_analysis$bounds_asymmetry, na.rm = TRUE), 4),
    round(mean(bes_analysis$is_symmetric, na.rm = TRUE), 4),
    round(mean(bes_analysis$n_obs, na.rm = TRUE), 0),
    round(median(bes_analysis$n_obs, na.rm = TRUE), 0)
  )
)

# Save summary statistics
summary_file <- file.path(tables_dir, "bes_summary_statistics.csv")
write.csv(bes_summary_stats, summary_file, row.names = FALSE)

# Create detailed results table (subset for inspection)
detailed_results <- bes_analysis %>%
  select(var1, var2, observed_r, r_min, r_max, bounds_range, 
         bounds_asymmetry, rescaled_r, n_obs) %>%
  arrange(desc(bounds_range)) %>%
  head(20)

detailed_file <- file.path(tables_dir, "bes_detailed_results_top20.csv")
write.csv(detailed_results, detailed_file, row.names = FALSE)

cat("📊 BES summary tables generated:\n")
cat("   1. Summary statistics:", basename(summary_file), "\n")  
cat("   2. Top 20 detailed results:", basename(detailed_file), "\n\n")

# ================================================================
# WORKFLOW COMPLETION SUMMARY  
# ================================================================

cat("🎉 BES DATA ANALYSIS WORKFLOW COMPLETE!\n")
cat("========================================\n")
cat("📁 Generated files:\n")
cat("   Bounds data:", basename(bounds_cache_file), "\n")
cat("   Analysis data:", basename(analysis_cache_file), "\n") 
cat("   Figures: 4 visualization files\n")
cat("   Tables: 2 summary files\n\n")

cat("📊 Key findings:\n")
cat("   Variable pairs analyzed:", nrow(bes_analysis), "\n")
cat("   Average bounds range:", round(mean(bes_analysis$bounds_range, na.rm = TRUE), 3), "\n")
cat("   Symmetric bounds:", round(100 * mean(bes_analysis$is_symmetric, na.rm = TRUE), 1), "%\n")
cat("   Mean rescaling effect:", round(mean(abs(bes_analysis$rescaled_r) - abs(bes_analysis$observed_r), na.rm = TRUE), 3), "\n\n")

cat("🔄 Next steps:\n")
cat("   - Run workflow_bivariate_asymmetry.R for detailed asymmetry analysis\n")
cat("   - Run workflow_matrices.R for matrix property analysis\n") 
cat("   - Compare with Monte Carlo results using comparison workflows\n\n")