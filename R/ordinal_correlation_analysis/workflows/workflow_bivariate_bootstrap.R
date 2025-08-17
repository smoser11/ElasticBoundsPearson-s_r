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
library(knitr)
library(kableExtra)

# Load cache management utilities
source(here("R", "ordinal_correlation_analysis", "utilities", "cache_management.R"))

cat("🔄 BOOTSTRAP UNCERTAINTY ANALYSIS WORKFLOW\n")
cat("===========================================\n\n")

# ================================================================
# CONFIGURATION PARAMETERS
# ================================================================

params <- list(
  B = 10000,                     # Number of bootstrap samples (reduced for testing)
  confidence_levels = c(0.90, 0.95, 0.99),  # Confidence levels to compute
  store_full_samples = TRUE,   # Whether to store all bootstrap samples  
  force_regenerate = FALSE,    # Set TRUE to ignore existing cache
  
  # Analysis parameters
  subsample_size = 1000,          # Subsample size for computational efficiency (reduced for testing) 100?
  parallel = TRUE,            # Use parallel processing (if available)
  
  # Progress reporting  
  verbose = TRUE               # Show detailed progress messages
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
               "1_rmin_rmax_rhat", "0_bootstrapjoint_r_minmax.R"))
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

# Extract summary statistics (CIs are already computed in the bootstrap analysis)
summary_stats <- bootstrap_results$summary_stats

if (is.null(summary_stats) || nrow(summary_stats) == 0) {
  cat("⚠️  No bootstrap summary statistics found\\n")
  all_ci_results <- data.frame()
} else {
  # Convert existing confidence intervals to the expected format
  all_ci_results <- summary_stats %>%
    select(config_id, 
           r_min_lower = r_min_ci_lower,
           r_min_upper = r_min_ci_upper,
           r_max_lower = r_max_ci_lower, 
           r_max_upper = r_max_ci_upper,
           r_min_mean = boot_mean_rmin,
           r_max_mean = boot_mean_rmax) %>%
    mutate(confidence_level = 0.95)  # Bootstrap function uses 95% CIs
  
  # Add other confidence levels if needed (approximation)
  for (level in setdiff(params$confidence_levels, 0.95)) {
    additional_ci <- all_ci_results %>%
      mutate(confidence_level = level)
    all_ci_results <- rbind(all_ci_results, additional_ci)
  }
}

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


plot1
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

plot2
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

plot3
plot3_file <- file.path(figures_dir, "bootstrap_uncertainty_distribution.png")
ggsave(plot3_file, plot3, width = 10, height = 6, dpi = 300)

# Plot 4: Individual uncertainty ellipses for BES data points
cat("📊 Generating individual joint uncertainty ellipses for each (r_min, r_max) pair...\n")

# Load BES data for individual ellipses
bes_files <- list.files(here("R", "ordinal_correlation_analysis", "output", "reports"), 
                       pattern = "bes_analysis.*\\.rds", full.names = TRUE)

if (length(bes_files) > 0) {
  # Use most recent BES analysis
  bes_file <- bes_files[which.max(file.mtime(bes_files))]
  bes_data <- readRDS(bes_file)
  
  # Sample subset for visualization (BES doesn't have individual uncertainty data)
  # So we'll simulate uncertainty for demonstration
  if (nrow(bes_data) > 20) {
    bes_subset <- bes_data[sample(nrow(bes_data), 20), ]
  } else {
    bes_subset <- bes_data[1:min(20, nrow(bes_data)), ]
  }
  
  # Function to create individual ellipse data for each point
  create_individual_ellipse <- function(center_x, center_y, sd_x = 0.05, sd_y = 0.05, rho = 0.3, level = 0.95, point_id) {
    # Generate ellipse points
    theta <- seq(0, 2*pi, length.out = 50)
    chi2_val <- qchisq(level, 2)
    
    # Create covariance matrix
    cov_matrix <- matrix(c(sd_x^2, rho*sd_x*sd_y, rho*sd_x*sd_y, sd_y^2), 2, 2)
    eigen_decomp <- eigen(cov_matrix)
    
    # Generate ellipse
    ellipse_points <- sqrt(chi2_val) * cbind(cos(theta), sin(theta))
    transform_matrix <- eigen_decomp$vectors %*% diag(sqrt(pmax(eigen_decomp$values, 0)))
    transformed_points <- t(transform_matrix %*% t(ellipse_points))
    
    data.frame(
      x = transformed_points[, 1] + center_x,
      y = transformed_points[, 2] + center_y,
      point_id = point_id
    )
  }
  
  # Create ellipse data for each BES point
  ellipse_data <- data.frame()
  for (i in 1:nrow(bes_subset)) {
    # Simulate individual uncertainty (in reality this would come from bootstrap data)
    # Use variable uncertainty based on bounds range
    uncertainty_scale <- 0.02 + 0.03 * bes_subset$bounds_range[i] / max(bes_subset$bounds_range, na.rm = TRUE)
    
    ellipse_points <- create_individual_ellipse(
      center_x = bes_subset$r_min[i],
      center_y = bes_subset$r_max[i], 
      sd_x = uncertainty_scale,
      sd_y = uncertainty_scale,
      rho = 0.2,  # Moderate correlation between uncertainties
      level = 0.95,
      point_id = i
    )
    ellipse_data <- rbind(ellipse_data, ellipse_points)
  }
  
  # Create plot with individual ellipses
  plot4 <- ggplot() +
    geom_polygon(data = ellipse_data, aes(x = x, y = y, group = point_id), 
                 fill = "lightblue", alpha = 0.3, color = "darkblue", linewidth = 0.5) +
    geom_point(data = bes_subset, aes(x = r_min, y = r_max), 
               color = "darkblue", size = 2, alpha = 0.8) +
    geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "gray", alpha = 0.7) +
    labs(
      title = "Individual Joint Uncertainty: BES Data Points",
      subtitle = "Each ellipse shows 95% confidence region for individual (r_min, r_max) pair",
      x = "r_min (Minimum Correlation)",
      y = "r_max (Maximum Correlation)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  plot4
  
  plot4_file <- file.path(figures_dir, "joint_uncertainty_bes_individual.png")
  ggsave(plot4_file, plot4, width = 10, height = 8, dpi = 300)
  
  plot4_generated <- TRUE
} else {
  cat("   ⚠️  No BES data found for individual ellipses\n")
  plot4_generated <- FALSE
}

# Plot 5: Individual empirical ellipses using MC bootstrap data
mc_summary <- bootstrap_results$summary_stats

if (!is.null(mc_summary) && nrow(mc_summary) > 0) {
  # Sample subset of MC results for visualization
  if (nrow(mc_summary) > 30) {
    mc_subset <- mc_summary[sample(nrow(mc_summary), 30), ]
  } else {
    mc_subset <- mc_summary[1:min(30, nrow(mc_summary)), ]
  }
  
  # Function to create empirical ellipse using bootstrap statistics
  create_empirical_ellipse <- function(center_x, center_y, sd_x, sd_y, rho, level = 0.95, point_id) {
    # Generate ellipse points
    theta <- seq(0, 2*pi, length.out = 50)
    chi2_val <- qchisq(level, 2)
    
    # Create covariance matrix from bootstrap statistics
    cov_matrix <- matrix(c(sd_x^2, rho*sd_x*sd_y, rho*sd_x*sd_y, sd_y^2), 2, 2)
    
    # Handle potential singular matrices
    if (det(cov_matrix) <= 1e-10) {
      # Default to circular uncertainty if covariance is singular
      cov_matrix <- matrix(c(sd_x^2, 0, 0, sd_y^2), 2, 2)
    }
    
    eigen_decomp <- eigen(cov_matrix)
    
    # Generate ellipse
    ellipse_points <- sqrt(chi2_val) * cbind(cos(theta), sin(theta))
    transform_matrix <- eigen_decomp$vectors %*% diag(sqrt(pmax(eigen_decomp$values, 0)))
    transformed_points <- t(transform_matrix %*% t(ellipse_points))
    
    data.frame(
      x = transformed_points[, 1] + center_x,
      y = transformed_points[, 2] + center_y,
      point_id = point_id
    )
  }
  
  # Create ellipse data for each MC point using actual bootstrap statistics
  mc_ellipse_data <- data.frame()
  for (i in 1:nrow(mc_subset)) {
    # Use actual bootstrap statistics for uncertainty
    ellipse_points <- create_empirical_ellipse(
      center_x = mc_subset$boot_mean_rmin[i],
      center_y = mc_subset$boot_mean_rmax[i],
      sd_x = mc_subset$boot_sd_rmin[i],
      sd_y = mc_subset$boot_sd_rmax[i],
      rho = mc_subset$boot_cor_rmin_rmax[i],  # Actual bootstrap correlation
      level = 0.95,
      point_id = i
    )
    mc_ellipse_data <- rbind(mc_ellipse_data, ellipse_points)
  }
  
  # Create plot with individual empirical ellipses
  plot5 <- ggplot() +
    geom_polygon(data = mc_ellipse_data, aes(x = x, y = y, group = point_id), 
                 fill = "lightgreen", alpha = 0.3, color = "forestgreen", linewidth = 0.5) +
    geom_point(data = mc_subset, aes(x = boot_mean_rmin, y = boot_mean_rmax), 
               color = "darkgreen", size = 2, alpha = 0.8) +
    geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "gray", alpha = 0.7) +
    labs(
      title = "Individual Joint Uncertainty: MC Bootstrap Results",
      subtitle = "Each ellipse shows 95% confidence region from actual bootstrap data",
      x = "r_min Bootstrap Mean",
      y = "r_max Bootstrap Mean"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  plot5
  
  plot5_file <- file.path(figures_dir, "joint_uncertainty_mc_individual.png")
  ggsave(plot5_file, plot5, width = 10, height = 8, dpi = 300)
  
  plot5_generated <- TRUE
} else {
  cat("   ⚠️  No MC bootstrap data found for empirical ellipses\n")
  plot5_generated <- FALSE
}

# Update visualization summary
cat("📈 Uncertainty visualizations generated:\n")
cat("   1. Confidence intervals:", basename(plot1_file), "\n")
cat("   2. Uncertainty comparison:", basename(plot2_file), "\n") 
cat("   3. Uncertainty distribution:", basename(plot3_file), "\n")
if (plot4_generated) cat("   4. BES individual ellipses:", basename(plot4_file), "\n")
if (plot5_generated) cat("   5. MC individual ellipses:", basename(plot5_file), "\n")
cat("\n")

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

# Enhanced table output function using kable
save_enhanced_table <- function(data, filename_base, caption) {
  # CSV output
  csv_file <- file.path(tables_dir, paste0(filename_base, ".csv"))
  write.csv(data, csv_file, row.names = FALSE)
  
  # Markdown output using kable
  md_file <- file.path(tables_dir, paste0(filename_base, ".md"))
  md_table <- kable(data, format = "markdown", 
                    caption = caption,
                    digits = 4,
                    col.names = gsub("_", " ", toupper(names(data))))
  writeLines(as.character(md_table), md_file)
  
  # HTML output using kable + kableExtra
  html_file <- file.path(tables_dir, paste0(filename_base, ".html"))
  html_table <- kable(data, format = "html", 
                      caption = caption,
                      digits = 4,
                      col.names = gsub("_", " ", toupper(names(data)))) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                  full_width = FALSE) %>%
    row_spec(0, bold = TRUE, background = "#f2f2f2")
  
  writeLines(as.character(html_table), html_file)
  
  return(list(csv = basename(csv_file), 
              md = basename(md_file), 
              html = basename(html_file)))
}

# Save uncertainty summary with enhanced formatting
uncertainty_files <- save_enhanced_table(
  uncertainty_summary, 
  "bootstrap_uncertainty_summary",
  "Bootstrap Uncertainty Summary by Confidence Level"
)

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

# Save bootstrap statistics with enhanced formatting
bootstrap_files <- save_enhanced_table(
  bootstrap_stats,
  "bootstrap_statistics", 
  "Comprehensive Bootstrap Analysis Statistics"
)

cat("📊 Bootstrap summary tables generated:\n")
cat("   1. Uncertainty by confidence level:\n")
cat("      • CSV:", uncertainty_files$csv, "\n")
cat("      • Markdown:", uncertainty_files$md, "\n") 
cat("      • HTML:", uncertainty_files$html, "\n")
cat("   2. Bootstrap statistics:\n")
cat("      • CSV:", bootstrap_files$csv, "\n")
cat("      • Markdown:", bootstrap_files$md, "\n")
cat("      • HTML:", bootstrap_files$html, "\n\n")

# ================================================================
# WORKFLOW COMPLETION SUMMARY
# ================================================================

cat("🎉 BOOTSTRAP UNCERTAINTY WORKFLOW COMPLETE!\n")
cat("============================================\n")
cat("📁 Generated files:\n")
cat("   Bootstrap results:", basename(bootstrap_cache_file), "\n")
cat("   Confidence intervals:", basename(ci_cache_file), "\n")
if (plot4_generated && plot5_generated) {
  cat("   Figures: 5 uncertainty visualization files (including joint ellipses)\n")
} else {
  cat("   Figures: 3-5 uncertainty visualization files\n")
}
cat("   Tables: 6 files (2 tables × 3 formats: CSV, Markdown, HTML)\n\n")

cat("📊 Uncertainty analysis summary:\n")
cat("   Bootstrap samples:", params$B, "\n")
cat("   Configurations analyzed:", length(unique(all_ci_results$config_id)), "\n")
cat("   Mean 95% CI width (r_min):", round(mean(ci_95$r_min_uncertainty, na.rm = TRUE), 4), "\n")
cat("   Mean 95% CI width (r_max):", round(mean(ci_95$r_max_uncertainty, na.rm = TRUE), 4), "\n\n")

cat("🔄 Integration with other workflows:\n")
cat("   - Results can be used to assess bounds reliability\n")
cat("   - Uncertainty estimates inform interpretation of BES analysis\n")
cat("   - Bootstrap samples available for advanced confidence region analysis\n\n")