# ================================================================
# MONTE CARLO SIMULATION WORKFLOW
# Modular sections - run independently or sequentially
# 
# Purpose: Generate simulated contingency tables, compute theoretical bounds,
#          and analyze Monte Carlo simulation results
# ================================================================

library(here)
library(dplyr)
library(ggplot2)

# Load cache management utilities
source(here("R", "ordinal_correlation_analysis", "utilities", "cache_management.R"))

cat("🎲 MONTE CARLO SIMULATION WORKFLOW\n")
cat("===================================\n\n")

# ================================================================
# CONFIGURATION PARAMETERS
# Modify these as needed for your analysis
# ================================================================

params <- list(
  numsims = 10,              # Number of simulated contingency tables 1000
  seed = 42,                   # Random seed for reproducibility  
  force_regenerate = FALSE,    # Set TRUE to ignore existing cache
  
  # Analysis parameters
  nsim_bounds = 10,          # Simulations for bounds computation 1000
  confidence_level = 0.95,     # Confidence level for intervals
  
  # Progress reporting
  verbose = TRUE,              # Show detailed progress messages
  progress_every = 10          # Report progress every N configurations
)

cat("📋 Configuration:\n")
cat("   Number of simulations:", params$numsims, "\n")
cat("   Random seed:", params$seed, "\n")
cat("   Force regenerate:", params$force_regenerate, "\n")
cat("   Bounds simulations:", params$nsim_bounds, "\n\n")

# ================================================================
# SECTION 1: CHECK/GENERATE MC SIMULATION DATA
# 
# This section generates the fundamental Monte Carlo simulated
# contingency tables with known marginal distributions
# ================================================================

cat("=== SECTION 1: MC SIMULATION DATA ===\n")

# Define cache file path
mc_cache_file <- here("R", "ordinal_correlation_analysis", "data", "raw", 
                     generate_cache_filename("MCsim", params[c("numsims", "seed")]))

mc_data <- cache_or_compute(
  cache_file = mc_cache_file,
  force_regenerate = params$force_regenerate,
  compute_func = function() {
    cat("Generating new Monte Carlo simulation data...\n")
    
    # Source the MC simulation functions
    source(here("R", "ordinal_correlation_analysis", "1_bivariate_ordcats_correlation", 
               "1_rmin_rmax_rhat", "1_monte_carlo_simulation", "make_MCsimulated_data.R"))
    
    # Set seed for reproducibility
    set.seed(params$seed)
    
    # Run Monte Carlo simulation
    sim_result <- run_mc_simulation(numsims = params$numsims)
    
    cat("Generated", length(sim_result), "simulation configurations\n")
    return(sim_result)
  }
)

cat("✅ MC simulation data ready\n")
cat("   Configurations:", length(mc_data$tables), "\n")
if (length(mc_data$tables) > 0 && !is.null(mc_data$tables[[1]])) {
  cat("   Tables per config:", length(mc_data$tables[[1]]), "\n")
}
cat("\n")

# ================================================================
# SECTION 2: CHECK/GENERATE MC BOUNDS ANALYSIS
# 
# This section computes theoretical correlation bounds for each
# simulated contingency table and analyzes the results
# ================================================================

cat("=== SECTION 2: MC BOUNDS ANALYSIS ===\n")

# Define cache file for analysis results
analysis_cache_file <- here("R", "ordinal_correlation_analysis", "output", "reports",
                           generate_cache_filename("mc_bounds_analysis", 
                                                  params[c("numsims", "seed", "nsim_bounds")]))

mc_bounds_analysis <- cache_or_compute(
  cache_file = analysis_cache_file,
  force_regenerate = params$force_regenerate,
  compute_func = function() {
    cat("Analyzing Monte Carlo bounds...\n")
    
    # Source required functions
    source(here("R", "correlation_bounds_core.R"))
    
    # Initialize results list
    analysis_results <- list()
    
    # Process each configuration
    for (config_id in seq_along(mc_data$tables)) {
      config_tables <- mc_data$tables[[config_id]]
      
      cat("Processing configuration", config_id, "of", length(mc_data$tables), "\n")
      
      # Analyze each contingency table in this configuration
      table_results <- list()
      
      for (table_id in seq_along(config_tables)) {
        table <- config_tables[[table_id]]
        
        # Extract marginal distributions
        marginal_x <- rowSums(table)
        marginal_y <- colSums(table)
        
        # Compute theoretical bounds
        r_max <- max_corr_bound(marginal_x, marginal_y)
        r_min <- min_corr_bound(marginal_x, marginal_y)
        
        # Store results
        table_results[[table_id]] <- list(
          config_id = config_id,
          table_id = table_id,
          table = table,
          marginal_x = marginal_x,
          marginal_y = marginal_y,
          r_min = r_min,
          r_max = r_max,
          bounds_range = r_max - r_min,
          bounds_asymmetry = r_max + r_min
        )
      }
      
      analysis_results[[config_id]] <- list(
        config_id = config_id,
        config_info = mc_data$config_info[config_id, ],
        table_results = table_results
      )
    }
    
    cat("Bounds analysis completed for all configurations\n")
    return(analysis_results)
  }
)

cat("✅ MC bounds analysis ready\n")
cat("   Analyzed configurations:", length(mc_bounds_analysis), "\n")
if (length(mc_bounds_analysis) > 0 && !is.null(mc_bounds_analysis[[1]]$table_results)) {
  total_tables <- sum(sapply(mc_bounds_analysis, function(x) length(x$table_results)))
  cat("   Total tables analyzed:", total_tables, "\n")
}
cat("\n")

# ================================================================
# SECTION 3: GENERATE MC SUMMARY STATISTICS
# 
# This section creates summary statistics and data frames
# suitable for visualization and further analysis
# ================================================================

cat("=== SECTION 3: MC SUMMARY STATISTICS ===\n")

# Create summary data frame
mc_summary <- data.frame()
total_configs <- length(mc_bounds_analysis)

if (params$verbose) {
  cat("📊 Building summary data frame from", total_configs, "configurations...\n")
}

for (config_idx in seq_along(mc_bounds_analysis)) {
  config <- mc_bounds_analysis[[config_idx]]
  
  # Progress reporting
  if (params$verbose && (config_idx %% params$progress_every == 0 || config_idx == total_configs)) {
    tables_in_config <- length(config$table_results)
    cat(sprintf("   📈 Processing config %d/%d (%d tables)\n", 
                config_idx, total_configs, tables_in_config))
  }
  
  for (table_result in config$table_results) {
    mc_summary <- rbind(mc_summary, data.frame(
      config_id = table_result$config_id,
      table_id = table_result$table_id,
      r_min = table_result$r_min,
      r_max = table_result$r_max,
      bounds_range = table_result$bounds_range,
      bounds_asymmetry = table_result$bounds_asymmetry,
      n_rows = nrow(table_result$table),
      n_cols = ncol(table_result$table),
      total_n = sum(table_result$table)
    ))
  }
}

if (params$verbose) {
  cat("✅ Summary data frame complete with", nrow(mc_summary), "table results\n\n")
}

# Generate summary statistics
cat("📊 Monte Carlo Bounds Summary:\n")
cat("   Total simulated tables:", nrow(mc_summary), "\n")
cat("   Mean r_min:", round(mean(mc_summary$r_min, na.rm = TRUE), 4), "\n")
cat("   Mean r_max:", round(mean(mc_summary$r_max, na.rm = TRUE), 4), "\n")
cat("   Mean bounds range:", round(mean(mc_summary$bounds_range, na.rm = TRUE), 4), "\n")
cat("   Mean asymmetry:", round(mean(mc_summary$bounds_asymmetry, na.rm = TRUE), 4), "\n")
cat("\n")

# Save summary data
summary_cache_file <- here("R", "ordinal_correlation_analysis", "output", "reports",
                          generate_cache_filename("mc_summary_data", 
                                                 params[c("numsims", "seed")], "rds"))
saveRDS(mc_summary, summary_cache_file)
cat("💾 Summary data saved:", basename(summary_cache_file), "\n\n")

# ================================================================
# SECTION 4: GENERATE MC VISUALIZATIONS
# 
# This section creates publication-ready visualizations of the
# Monte Carlo simulation results
# ================================================================

cat("=== SECTION 4: MC VISUALIZATIONS ===\n")

# Ensure figures directory exists
figures_dir <- here("R", "ordinal_correlation_analysis", "output", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Plot 1: Bounds landscape (r_min vs r_max)
plot1 <- ggplot(mc_summary, aes(x = r_min, y = r_max)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Monte Carlo Simulation: Correlation Bounds Landscape",
    subtitle = paste("Based on", nrow(mc_summary), "simulated contingency tables"),
    x = "Minimum Possible Correlation (r_min)",
    y = "Maximum Possible Correlation (r_max)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

plot1

plot1_file <- file.path(figures_dir, "mc_bounds_landscape.png")
ggsave(plot1_file, plot1, width = 10, height = 8, dpi = 300)

# Plot 2: Bounds range distribution
plot2 <- ggplot(mc_summary, aes(x = bounds_range)) +
  geom_histogram(bins = 50, fill = "lightblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = mean(mc_summary$bounds_range), 
             linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "Distribution of Correlation Bounds Ranges",
    subtitle = paste("Mean range:", round(mean(mc_summary$bounds_range, na.rm = TRUE), 3)),
    x = "Bounds Range (r_max - r_min)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

plot2

plot2_file <- file.path(figures_dir, "mc_bounds_range_distribution.png")
ggsave(plot2_file, plot2, width = 10, height = 6, dpi = 300)

# Plot 3: Bounds asymmetry distribution
plot3 <- ggplot(mc_summary, aes(x = bounds_asymmetry)) +
  geom_histogram(bins = 50, fill = "lightgreen", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = mean(mc_summary$bounds_asymmetry), 
             linetype = "dashed", color = "blue", size = 1) +
  labs(
    title = "Distribution of Correlation Bounds Asymmetry",
    subtitle = "Asymmetry = r_max + r_min (perfect symmetry = 0)",
    x = "Bounds Asymmetry",
    y = "Count"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

plot3

plot3_file <- file.path(figures_dir, "mc_bounds_asymmetry_distribution.png")
ggsave(plot3_file, plot3, width = 10, height = 6, dpi = 300)

cat("📈 Generated visualizations:\n")
cat("   1. Bounds landscape:", basename(plot1_file), "\n")
cat("   2. Range distribution:", basename(plot2_file), "\n") 
cat("   3. Asymmetry distribution:", basename(plot3_file), "\n\n")

# ================================================================
# SECTION 5: GENERATE MC SUMMARY TABLES
# 
# This section creates summary tables in multiple formats
# ================================================================

cat("=== SECTION 5: MC SUMMARY TABLES ===\n")

# Ensure tables directory exists
tables_dir <- here("R", "ordinal_correlation_analysis", "output", "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# Create summary statistics table
summary_stats <- data.frame(
  Statistic = c(
    "Number of Simulated Tables",
    "Mean r_min", 
    "SD r_min",
    "Mean r_max",
    "SD r_max", 
    "Mean Bounds Range",
    "SD Bounds Range",
    "Mean Bounds Asymmetry",
    "SD Bounds Asymmetry"
  ),
  Value = c(
    nrow(mc_summary),
    round(mean(mc_summary$r_min, na.rm = TRUE), 4),
    round(sd(mc_summary$r_min, na.rm = TRUE), 4),
    round(mean(mc_summary$r_max, na.rm = TRUE), 4),
    round(sd(mc_summary$r_max, na.rm = TRUE), 4),
    round(mean(mc_summary$bounds_range, na.rm = TRUE), 4),
    round(sd(mc_summary$bounds_range, na.rm = TRUE), 4),
    round(mean(mc_summary$bounds_asymmetry, na.rm = TRUE), 4),
    round(sd(mc_summary$bounds_asymmetry, na.rm = TRUE), 4)
  )
)

# Save in multiple formats
csv_file <- file.path(tables_dir, "mc_summary_statistics.csv")
write.csv(summary_stats, csv_file, row.names = FALSE)

cat("📊 Summary statistics table saved:", basename(csv_file), "\n\n")

# ================================================================
# WORKFLOW COMPLETION SUMMARY
# ================================================================

cat("🎉 MONTE CARLO WORKFLOW COMPLETE!\n")
cat("==================================\n")
cat("📁 Generated files:\n")
cat("   Data: ", basename(mc_cache_file), "\n")
cat("   Analysis: ", basename(analysis_cache_file), "\n")
cat("   Summary: ", basename(summary_cache_file), "\n")
cat("   Figures: 3 visualization files\n") 
cat("   Tables: 1 summary statistics file\n\n")

cat("📊 Analysis summary:\n")
cat("   Simulated tables:", nrow(mc_summary), "\n")
cat("   Average bounds range:", round(mean(mc_summary$bounds_range), 3), "\n")
cat("   Proportion asymmetric:", round(mean(abs(mc_summary$bounds_asymmetry) > 0.01), 3), "\n\n")

cat("🔄 To rerun sections:\n")
cat("   - Modify params$force_regenerate = TRUE for fresh computation\n")
cat("   - Modify params$numsims for different simulation size\n") 
cat("   - Run individual sections by highlighting code blocks\n\n")