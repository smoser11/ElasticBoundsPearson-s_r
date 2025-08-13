# Complete R Project Reorganization for Ordinal Correlation Analysis

## Directory Structure with File Mapping

```
ordinal_correlation_analysis/
├── 1_bivariate_ordcats_correlation/
│   ├── 1_rmin_rmax_rhat/
│   │   ├── 1_monte_carlo_simulation/
│   │   │   ├── core_bounds_functions.R          # From: correlation_bounds_core.R
│   │   │   ├── simulation_studies.R             # From: correlation_bounds_simulation.R
│   │   │   └── validation_tests.R               # From: correlation_bounds_examples.R
│   │   ├── 2_bes_illustrative_example/
│   │   │   ├── bes_data_analysis.R              # From: correlation_bounds_bes.R
│   │   │   ├── example_demonstrations.R         # From: correlation_bounds_bes_example.R
│   │   │   └── case_studies.R                   # From: correlation_bounds_demo.R
│   │   └── bounds_analysis_main.R               # Integration module
│   ├── 2_asymmetry_analysis/
│   │   ├── asymmetry_measures.R                 # New synthesis
│   │   ├── total_variation_functions.R          # New functionality
│   │   └── asymmetry_visualization.R            # From: correlation-bounds-visualization.R
│   ├── 3_visualization/
│   │   ├── bounds_visualization.R               # From: correlation-bounds-visualization.R
│   │   ├── distribution_plots.R                 # From: correlation_bounds_bes.R (plotting functions)
│   │   └── interactive_plots.R                  # Enhanced from examples
│   └── bivariate_main.R                         # Module coordinator
├── 2_correlation_matrices/
│   ├── 1_matrix_properties/
│   │   ├── invertibility_tests.R                # From: correlation_matrix_test.R
│   │   ├── condition_number_analysis.R          # From: correlation_matrix_test.R
│   │   ├── psd_validation.R                     # From: correlation_matrix_test.R
│   │   └── matrix_diagnostics.R                 # Synthesis of matrix functions
│   ├── 2_matrix_construction/
│   │   ├── matrix_builder.R                     # From: correlation_matrix_test.R
│   │   ├── random_trials.R                      # From: correlation_matrix_test.R
│   │   └── validation_framework.R               # From: correlation_matrix_example.R
│   └── matrices_main.R                          # Module coordinator
├── 3_fixes_and_rescaling/
│   ├── 1_simple_rescaling/
│   │   ├── linear_rescaling.R                   # From: correlation_bounds_bes.R
│   │   └── rescaling_validation.R               # From: correlation_bounds_demo.R
│   ├── 2_advanced_rescaling/
│   │   ├── sophisticated_methods.R              # New advanced approaches
│   │   ├── preserving_properties.R              # Matrix property preservation
│   │   └── alternative_approaches.R             # Future extensions
│   └── rescaling_main.R                         # Module coordinator
├── utilities/
│   ├── helper_functions.R                       # Common utilities
│   ├── data_loading.R                           # Data handling functions
│   └── plotting_themes.R                        # Consistent plot styling
├── data/
│   ├── raw/                                     # Original BES data
│   └── processed/                               # Cleaned/transformed data
├── output/
│   ├── figures/                                 # Generated plots
│   ├── tables/                                  # Analysis results
│   └── reports/                                 # Full analysis reports
└── main_analysis.R                              # Master analysis script
```

## Complete File Implementations

### 1. Core Bounds Functions (`1_bivariate_ordcats_correlation/1_rmin_rmax_rhat/1_monte_carlo_simulation/core_bounds_functions.R`)

```r
# core_bounds_functions.R
# Complete implementation synthesized from correlation_bounds_core.R

#' Compute maximum correlation bound using Fréchet-Hoeffding upper bound
#' 
#' @param marginalX Vector of marginal probabilities or counts for variable X
#' @param marginalY Vector of marginal probabilities or counts for variable Y  
#' @param sample_size Sample size (used if marginals are probabilities)
#' @return Maximum possible Pearson correlation coefficient
max_corr_bound <- function(marginalX, marginalY, sample_size = 10000) {
  # Convert probabilities to counts if needed
  if (sum(marginalX) <= 1.1) {
    countsX <- as.vector(rmultinom(1, size = sample_size, prob = marginalX))
  }
```

### 8. Utility Functions (`utilities/helper_functions.R`)

```r
# helper_functions.R
# Common utility functions used across modules

#' Generate significance test comparison summary
#'
#' @param results_df Data frame with correlation bounds results
#' @return Summary of t-test vs randomization test differences
generate_significance_summary <- function(results_df) {
  # Calculate t-test p-values and significance
  results_df$t_value <- results_df$observed_r * sqrt((results_df$n_obs - 2) / (1 - results_df$observed_r^2))
  results_df$t_pvalue <- 2 * pt(-abs(results_df$t_value), df = results_df$n_obs - 2)
  results_df$t_significant <- results_df$t_pvalue < 0.05
  
  # Randomization test significance (outside CI)
  results_df$randomization_significant <- results_df$observed_r < results_df$ci_lower | 
                                         results_df$observed_r > results_df$ci_upper
  
  # Compare significance tests
  results_df$significance_agreement <- results_df$t_significant == results_df$randomization_significant
  
  results_df$significance_difference <- case_when(
    results_df$t_significant & !results_df$randomization_significant ~ "t-test only",
    !results_df$t_significant & results_df$randomization_significant ~ "Randomization only", 
    results_df$t_significant & results_df$randomization_significant ~ "Both",
    !results_df$t_significant & !results_df$randomization_significant ~ "Neither",
    TRUE ~ NA_character_
  )
  
  # Summary statistics
  significance_summary <- list(
    total_pairs = nrow(results_df),
    t_test_significant = sum(results_df$t_significant, na.rm = TRUE),
    randomization_significant = sum(results_df$randomization_significant, na.rm = TRUE),
    both_significant = sum(results_df$significance_difference == "Both", na.rm = TRUE),
    neither_significant = sum(results_df$significance_difference == "Neither", na.rm = TRUE),
    t_test_only = sum(results_df$significance_difference == "t-test only", na.rm = TRUE),
    randomization_only = sum(results_df$significance_difference == "Randomization only", na.rm = TRUE),
    agreement_rate = mean(results_df$significance_agreement, na.rm = TRUE),
    disagreement_rate = 1 - mean(results_df$significance_agreement, na.rm = TRUE)
  )
  
  return(list(
    detailed_results = results_df,
    summary = significance_summary
  ))
}

#' Safe division avoiding division by zero
#'
#' @param numerator Numerator value
#' @param denominator Denominator value
#' @param default_value Value to return when denominator is zero
#' @return Safe division result
safe_divide <- function(numerator, denominator, default_value = NA) {
  ifelse(abs(denominator) < 1e-15, default_value, numerator / denominator)
}

#' Robust correlation calculation with error handling
#'
#' @param x First variable
#' @param y Second variable
#' @param method Correlation method
#' @return Correlation coefficient or NA if calculation fails
robust_cor <- function(x, y, method = "pearson") {
  tryCatch({
    if (length(x) != length(y) || length(x) < 3) return(NA)
    if (all(is.na(x)) || all(is.na(y))) return(NA)
    if (var(x, na.rm = TRUE) == 0 || var(y, na.rm = TRUE) == 0) return(NA)
    cor(x, y, method = method, use = "complete.obs")
  }, error = function(e) NA)
}

#' Check if a vector represents valid probabilities
#'
#' @param probs Vector to check
#' @param tolerance Numerical tolerance
#' @return Logical indicating validity
is_valid_probability <- function(probs, tolerance = 1e-10) {
  all(probs >= -tolerance) && abs(sum(probs) - 1) < tolerance
}

#' Convert frequency counts to probabilities with validation
#'
#' @param counts Vector of frequency counts
#' @return Vector of probabilities
counts_to_probs <- function(counts) {
  if (sum(counts) == 0) return(rep(0, length(counts)))
  probs <- counts / sum(counts)
  return(probs)
}

#' Validate BES data row for analysis
#'
#' @param row Single row from BES dataset
#' @return List with validation results
validate_bes_row <- function(row) {
  errors <- c()
  warnings <- c()
  
  # Check required columns
  required_cols <- c("var1", "var2", "corr", "nobs", "var1cats", "var2cats")
  missing_cols <- setdiff(required_cols, names(row))
  if (length(missing_cols) > 0) {
    errors <- c(errors, paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Check data validity
  if ("nobs" %in% names(row) && (is.na(row$nobs) || row$nobs < 10)) {
    warnings <- c(warnings, "Sample size too small (< 10)")
  }
  
  if ("corr" %in% names(row) && (is.na(row$corr) || abs(row$corr) > 1)) {
    errors <- c(errors, "Invalid correlation value")
  }
  
  if ("var1cats" %in% names(row) && (is.na(row$var1cats) || row$var1cats < 2)) {
    errors <- c(errors, "Invalid number of categories for var1")
  }
  
  if ("var2cats" %in% names(row) && (is.na(row$var2cats) || row$var2cats < 2)) {
    errors <- c(errors, "Invalid number of categories for var2")
  }
  
  return(list(
    valid = length(errors) == 0,
    errors = errors,
    warnings = warnings
  ))
}

#' Create a standardized error result
#'
#' @param row Original data row
#' @param error_message Error description
#' @param error_code Optional error code
#' @return Standardized error result
create_error_result <- function(row, error_message, error_code = NULL) {
  list(
    var1 = if("var1" %in% names(row)) row$var1 else NA,
    var2 = if("var2" %in% names(row)) row$var2 else NA,
    observed_r = if("corr" %in% names(row)) row$corr else NA,
    r_min = NA,
    r_max = NA,
    r_rescaled = NA,
    ci_lower = NA,
    ci_upper = NA,
    success = FALSE,
    error_message = error_message,
    error_code = error_code
  )
}

#' Format numbers for reporting with appropriate precision
#'
#' @param x Numeric value or vector
#' @param digits Number of decimal places
#' @param scientific Use scientific notation for small numbers
#' @return Formatted string
format_number <- function(x, digits = 3, scientific = FALSE) {
  if (is.na(x)) return("NA")
  if (scientific && abs(x) < 0.001 && x != 0) {
    return(formatC(x, format = "e", digits = digits))
  } else {
    return(formatC(round(x, digits), format = "f", digits = digits))
  }
}

#' Calculate bootstrap confidence intervals
#'
#' @param data Vector of data values
#' @param statistic Function to calculate statistic
#' @param n_bootstrap Number of bootstrap samples
#' @param conf_level Confidence level (0-1)
#' @param seed Random seed
#' @return List with confidence interval bounds
bootstrap_ci <- function(data, statistic, n_bootstrap = 1000, conf_level = 0.95, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Original statistic
  original_stat <- statistic(data)
  
  # Bootstrap samples
  bootstrap_stats <- replicate(n_bootstrap, {
    bootstrap_sample <- sample(data, length(data), replace = TRUE)
    statistic(bootstrap_sample)
  })
  
  # Calculate confidence interval
  alpha <- 1 - conf_level
  ci_bounds <- quantile(bootstrap_stats, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  
  return(list(
    original = original_stat,
    lower = ci_bounds[1],
    upper = ci_bounds[2],
    bootstrap_stats = bootstrap_stats
  ))
}
```

### 9. Visualization Functions (`1_bivariate_ordcats_correlation/3_visualization/bounds_visualization.R`)

```r
# bounds_visualization.R
# Comprehensive visualization functions for correlation bounds analysis

#' Create correlation bounds scatter plot
#'
#' @param results_df Data frame with bounds analysis results
#' @param color_by Variable to color points by
#' @param title Plot title
#' @return ggplot object
plot_bounds_scatter <- function(results_df, color_by = "bounds_range", title = "Correlation Bounds Analysis") {
  # Calculate bounds range if not present
  if (!"bounds_range" %in% names(results_df)) {
    results_df$bounds_range <- results_df$r_max - results_df$r_min
  }
  
  p <- ggplot(results_df, aes(x = r_min, y = r_max)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray", alpha = 0.5) +
    coord_equal() +
    theme_minimal() +
    labs(
      title = title,
      subtitle = paste("Analysis of", nrow(results_df), "variable pairs"),
      x = "Minimum Possible Correlation (r_min)",
      y = "Maximum Possible Correlation (r_max)"
    )
  
  # Add color mapping based on specified variable
  if (color_by == "bounds_range") {
    p <- p + 
      geom_point(aes(color = bounds_range), alpha = 0.7, size = 1.5) +
      scale_color_viridis_c(name = "Bounds\nRange", option = "plasma")
  } else if (color_by == "observed_r") {
    p <- p + 
      geom_point(aes(color = observed_r), alpha = 0.7, size = 1.5) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                           midpoint = 0, name = "Observed\nCorrelation")
  } else if (color_by == "asymmetry") {
    if ("basic_asymmetry" %in% names(results_df)) {
      p <- p + 
        geom_point(aes(color = abs(basic_asymmetry)), alpha = 0.7, size = 1.5) +
        scale_color_viridis_c(name = "Asymmetry\nMagnitude", option = "inferno")
    }
  } else {
    p <- p + geom_point(alpha = 0.7, size = 1.5, color = "steelblue")
  }
  
  return(p)
}

#' Create permutation distribution histogram with bounds
#'
#' @param r_sim Vector of simulated correlation values
#' @param r_min Theoretical minimum correlation
#' @param r_max Theoretical maximum correlation  
#' @param ci_bounds Confidence interval bounds
#' @param observed_r Observed correlation (optional)
#' @param title Plot title
#' @return ggplot object
plot_permutation_distribution <- function(r_sim, r_min, r_max, ci_bounds, 
                                        observed_r = NULL, title = "Permutation Distribution") {
  df <- data.frame(r = r_sim)
  
  # Set x-axis limits
  x_range <- range(c(r_sim, r_min, r_max, ci_bounds, observed_r), na.rm = TRUE)
  x_limits <- x_range + c(-0.05, 0.05) * diff(x_range)
  
  p <- ggplot(df, aes(x = r)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", 
                   color = "white", alpha = 0.7) +
    geom_density(color = "darkblue", size = 1) +
    geom_vline(xintercept = r_min, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = r_max, color = "green", linetype = "dashed", size = 1) +
    geom_vline(xintercept = ci_bounds[1], color = "purple", linetype = "dotted", size = 1) +
    geom_vline(xintercept = ci_bounds[2], color = "purple", linetype = "dotted", size = 1) +
    xlim(x_limits) +
    theme_minimal() +
    labs(
      title = title,
      subtitle = paste0("Theoretical bounds: [", round(r_min, 3), ", ", round(r_max, 3), "]",
                       ", 95% CI: [", round(ci_bounds[1], 3), ", ", round(ci_bounds[2], 3), "]"),
      x = "Pearson's r",
      y = "Density"
    )
  
  # Add observed correlation if provided
  if (!is.null(observed_r)) {
    p <- p + 
      geom_vline(xintercept = observed_r, color = "black", size = 1.2) +
      labs(subtitle = paste0("Theoretical bounds: [", round(r_min, 3), ", ", round(r_max, 3), "]",
                            ", 95% CI: [", round(ci_bounds[1], 3), ", ", round(ci_bounds[2], 3), "]",
                            ", Observed: ", round(observed_r, 3)))
  }
  
  return(p)
}

#' Create rescaling comparison plot
#'
#' @param results_df Data frame with original and rescaled correlations
#' @param title Plot title
#' @return ggplot object
plot_rescaling_comparison <- function(results_df, title = "Original vs Rescaled Correlations") {
  # Calculate rescaling magnitude if not present
  if (!"rescaling_magnitude" %in% names(results_df)) {
    results_df$rescaling_magnitude <- abs(results_df$r_rescaled) - abs(results_df$observed_r)
  }
  
  p <- ggplot(results_df, aes(x = observed_r, y = r_rescaled)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray", alpha = 0.5) +
    geom_point(aes(color = rescaling_magnitude), alpha = 0.7, size = 1.5) +
    scale_color_gradient2(low = "blue", mid = "gray", high = "red", 
                         midpoint = 0, name = "Rescaling\nEffect") +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
    theme_minimal() +
    labs(
      title = title,
      subtitle = "Points above diagonal increased in magnitude, below decreased",
      x = "Original Correlation",
      y = "Rescaled Correlation"
    )
  
  return(p)
}

#' Create asymmetry analysis visualization
#'
#' @param asymmetry_results Results from asymmetry analysis
#' @return List of ggplot objects
plot_asymmetry_analysis <- function(asymmetry_results) {
  data <- asymmetry_results$detailed_results
  
  # Asymmetry distribution
  p1 <- ggplot(data, aes(x = basic_asymmetry)) +
    geom_histogram(bins = 30, fill = "lightcoral", alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +
    labs(
      title = "Distribution of Bounds Asymmetry",
      subtitle = paste0("Perfect symmetry at 0, ", 
                       round(100 * asymmetry_results$summary$prop_symmetric, 1), 
                       "% of pairs are approximately symmetric"),
      x = "Asymmetry (r_max + r_min)",
      y = "Count"
    )
  
  # Asymmetry vs range
  p2 <- ggplot(data, aes(x = bounds_range, y = abs(basic_asymmetry))) +
    geom_point(aes(color = asymmetry_class), alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE, color = "black") +
    scale_color_viridis_d(name = "Asymmetry\nClass") +
    theme_minimal() +
    labs(
      title = "Asymmetry vs Bounds Range",
      x = "Bounds Range (r_max - r_min)",
      y = "Asymmetry Magnitude |r_max + r_min|"
    )
  
  # Asymmetry by categories
  if ("var1_categories" %in% names(data)) {
    data$total_categories <- data$var1_categories + data$var2_categories
    
    p3 <- ggplot(data, aes(x = factor(total_categories), y = abs(basic_asymmetry))) +
      geom_boxplot(fill = "lightgreen", alpha = 0.7) +
      theme_minimal() +
      labs(
        title = "Asymmetry by Total Number of Categories",
        x = "Total Categories (var1 + var2)",
        y = "Asymmetry Magnitude"
      )
  } else {
    p3 <- NULL
  }
  
  return(list(distribution = p1, vs_range = p2, by_categories = p3))
}

#' Generate all visualizations for the analysis
#'
#' @param bounds_results Bounds analysis results
#' @param asymmetry_results Asymmetry analysis results  
#' @param matrix_results Matrix testing results
#' @param rescaling_results Rescaling analysis results
#' @param output_dir Output directory for saving plots
#' @return List of all generated plots
generate_all_visualizations <- function(bounds_results, asymmetry_results = NULL,
                                       matrix_results = NULL, rescaling_results = NULL,
                                       output_dir = "output/figures") {
  plots <- list()
  
  # 1. Bounds visualizations
  plots$bounds_scatter <- plot_bounds_scatter(bounds_results, color_by = "bounds_range")
  plots$bounds_by_observed <- plot_bounds_scatter(bounds_results, color_by = "observed_r")
  
  if (!is.null(asymmetry_results)) {
    plots$bounds_by_asymmetry <- plot_bounds_scatter(bounds_results, color_by = "asymmetry")
  }
  
  # 2. Rescaling visualizations
  if ("r_rescaled" %in% names(bounds_results)) {
    plots$rescaling_comparison <- plot_rescaling_comparison(bounds_results)
  }
  
  # 3. Asymmetry visualizations
  if (!is.null(asymmetry_results)) {
    asymmetry_plots <- plot_asymmetry_analysis(asymmetry_results)
    plots <- c(plots, asymmetry_plots)
  }
  
  # 4. Matrix properties visualizations
  if (!is.null(matrix_results) && "psd" %in% names(matrix_results)) {
    plots$matrix_comparison <- plot_matrix_properties_comparison(matrix_results)
  }
  
  # Save plots if output directory specified
  if (!is.null(output_dir) && dir.exists(dirname(output_dir))) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    for (plot_name in names(plots)) {
      if (!is.null(plots[[plot_name]])) {
        ggsave(
          filename = file.path(output_dir, paste0(plot_name, ".png")),
          plot = plots[[plot_name]],
          width = 10, height = 8, dpi = 300
        )
      }
    }
  }
  
  return(plots)
}
```

### 10. Report Generation (`utilities/report_generation.R`)

```r
# report_generation.R
# Functions for generating comprehensive analysis reports

#' Compile comprehensive analysis report
#'
#' @param bounds_results Bounds analysis results
#' @param asymmetry_results Asymmetry analysis results
#' @param significance_results Significance testing results
#' @param matrix_results_psd Matrix results with PSD enforcement
#' @param matrix_results_no_psd Matrix results without PSD enforcement
#' @param rescaling_results Rescaling analysis results
#' @param config Analysis configuration
#' @return Comprehensive report list
compile_analysis_report <- function(bounds_results, asymmetry_results, significance_results,
                                   matrix_results_psd, matrix_results_no_psd, 
                                   rescaling_results, config) {
  
  report <- list(
    metadata = list(
      analysis_date = Sys.time(),
      n_variable_pairs = nrow(bounds_results),
      configuration = config
    ),
    
    # Executive summary
    executive_summary = compile_executive_summary(
      bounds_results, asymmetry_results, matrix_results_psd, rescaling_results
    ),
    
    # Detailed results by analysis component
    bounds_analysis = compile_bounds_summary(bounds_results),
    asymmetry_analysis = compile_asymmetry_summary(asymmetry_results),
    significance_analysis = compile_significance_summary(significance_results),
    matrix_analysis = compile_matrix_summary(matrix_results_psd, matrix_results_no_psd),
    rescaling_analysis = compile_rescaling_summary(rescaling_results),
    
    # Key findings and recommendations
    key_findings = extract_key_findings(bounds_results, asymmetry_results, 
                                       matrix_results_psd, rescaling_results),
    recommendations = generate_recommendations(bounds_results, asymmetry_results,
                                             matrix_results_psd, rescaling_results)
  )
  
  return(report)
}

#' Generate executive summary
compile_executive_summary <- function(bounds_results, asymmetry_results, 
                                     matrix_results, rescaling_results) {
  list(
    total_pairs = nrow(bounds_results),
    bounds_range_summary = list(
      mean = mean(bounds_results$r_max - bounds_results$r_min, na.rm = TRUE),
      median = median(bounds_results$r_max - bounds_results$r_min, na.rm = TRUE),
      min = min(bounds_results$r_max - bounds_results$r_min, na.rm = TRUE),
      max = max(bounds_results$r_max - bounds_results$r_min, na.rm = TRUE)
    ),
    asymmetry_summary = if(!is.null(asymmetry_results)) {
      list(
        prop_symmetric = asymmetry_results$summary$prop_symmetric,
        mean_asymmetry = asymmetry_results$summary$mean_asymmetry,
        max_asymmetry = asymmetry_results$summary$max_asymmetry
      )
    } else { NULL },
    matrix_validity = list(
      psd_rate = matrix_results$orig_psd_percent,
      invertible_rate = matrix_results$orig_invertible_percent,
      mean_condition_number = matrix_results$orig_mean_condition
    ),
    rescaling_effects = list(
      mean_magnitude_change = rescaling_results$summary$mean_magnitude_change,
      prop_increased = rescaling_results$summary$prop_increased,
      prop_decreased = rescaling_results$summary$prop_decreased
    )
  )
}

#' Save analysis report to files
#'
#' @param report Compiled report
#' @param output_dir Output directory
save_analysis_report <- function(report, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Save complete report as RDS
  saveRDS(report, file.path(output_dir, "complete_analysis_report.rds"))
  
  # Save summary as JSON for easy reading
  summary_json <- jsonlite::toJSON(report$executive_summary, pretty = TRUE, auto_unbox = TRUE)
  writeLines(summary_json, file.path(output_dir, "executive_summary.json"))
  
  # Generate formatted text report
  text_report <- generate_text_report(report)
  writeLines(text_report, file.path(output_dir, "analysis_report.txt"))
  
  # Generate CSV summaries for key results
  write.csv(report$bounds_analysis$summary_table, 
           file.path(output_dir, "bounds_summary.csv"), row.names = FALSE)
  
  if (!is.null(report$asymmetry_analysis)) {
    write.csv(report$asymmetry_analysis$summary_table,
             file.path(output_dir, "asymmetry_summary.csv"), row.names = FALSE)
  }
}

#' Generate formatted text report
#'
#' @param report Compiled report
#' @return Character vector with formatted report
generate_text_report <- function(report) {
  lines <- c(
    "ORDINAL CORRELATION BOUNDS ANALYSIS REPORT",
    "==========================================",
    "",
    paste("Analysis Date:", report$metadata$analysis_date),
    paste("Number of Variable Pairs:", report$metadata$n_variable_pairs),
    "",
    "EXECUTIVE SUMMARY",
    "-----------------",
    "",
    paste("Total variable pairs analyzed:", report$executive_summary$total_pairs),
    "",
    "Correlation Bounds Summary:",
    paste("  Mean bounds range:", format_number(report$executive_summary$bounds_range_summary$mean)),
    paste("  Median bounds range:", format_number(report$executive_summary$bounds_range_summary$median)),
    paste("  Range: [", format_number(report$executive_summary$bounds_range_summary$min), 
          ", ", format_number(report$executive_summary$bounds_range_summary$max), "]", sep = ""),
    ""
  )
  
  if (!is.null(report$executive_summary$asymmetry_summary)) {
    lines <- c(lines,
      "Asymmetry Analysis:",
      paste("  Proportion symmetric bounds:", 
            format_number(100 * report$executive_summary$asymmetry_summary$prop_symmetric), "%"),
      paste("  Mean asymmetry:", 
            format_number(report$executive_summary$asymmetry_summary$mean_asymmetry)),
      paste("  Maximum asymmetry:", 
            format_number(report$executive_summary$asymmetry_summary$max_asymmetry)),
      ""
    )
  }
  
  lines <- c(lines,
    "Matrix Properties:",
    paste("  Positive semidefinite rate:", 
          format_number(report$executive_summary$matrix_validity$psd_rate), "%"),
    paste("  Invertible rate:", 
          format_number(report$executive_summary$matrix_validity$invertible_rate), "%"),
    paste("  Mean condition number:", 
          format_number(report$executive_summary$matrix_validity$mean_condition_number)),
    "",
    "Rescaling Effects:",
    paste("  Mean magnitude change:", 
          format_number(report$executive_summary$rescaling_effects$mean_magnitude_change)),
    paste("  Proportion increased:", 
          format_number(100 * report$executive_summary$rescaling_effects$prop_increased), "%"),
    paste("  Proportion decreased:", 
          format_number(100 * report$executive_summary$rescaling_effects$prop_decreased), "%"),
    "",
    "KEY FINDINGS",
    "------------"
  )
  
  # Add key findings
  for (finding in report$key_findings) {
    lines <- c(lines, paste("• ", finding))
  }
  
  lines <- c(lines, "", "RECOMMENDATIONS", "---------------")
  
  # Add recommendations
  for (rec in report$recommendations) {
    lines <- c(lines, paste("• ", rec))
  }
  
  return(lines)
}
```

## Migration Guide

To transition from your current codebase to this new structure:

### Step 1: Create Directory Structure
```bash
mkdir -p ordinal_correlation_analysis/{1_bivariate_ordcats_correlation/{1_rmin_rmax_rhat/{1_monte_carlo_simulation,2_bes_illustrative_example},2_asymmetry_analysis,3_visualization},2_correlation_matrices/{1_matrix_properties,2_matrix_construction},3_fixes_and_rescaling/{1_simple_rescaling,2_advanced_rescaling},utilities,data/{raw,processed},output/{figures,tables,reports}}
```

### Step 2: Move and Adapt Functions
1. Extract core bounds functions from `correlation_bounds_core.R` → `1_monte_carlo_simulation/core_bounds_functions.R`
2. Adapt BES-specific functions from `correlation_bounds_bes.R` → `2_bes_illustrative_example/bes_data_analysis.R`
3. Consolidate matrix functions from `correlation_matrix_test.R` → `2_correlation_matrices/1_matrix_properties/`
4. Move visualization from `correlation-bounds-visualization.R` → `3_visualization/bounds_visualization.R`

### Step 3: Update Function Calls
- Replace direct function calls with module-based sourcing
- Update file paths in source() statements
- Standardize function naming conventions

### Step 4: Test Integration
- Run `main_analysis.R` with a small subset of data
- Verify all modules load correctly
- Check that outputs are generated in correct directories

This reorganization provides:
- **Clear separation of concerns** between bivariate analysis, matrix properties, and rescaling
- **Modular development** allowing independent work on each component
- **Comprehensive functionality** preserving all your original analytical capabilities
- **Professional structure** following R package development best practices
- **Extensible framework** for adding new methods and visualizations

## Module Integration Example

Here's how the modules work together in practice:

### Example Analysis Workflow

```r
# Load BES data
bes_data <- load_bes_data("data/raw/bes_2019.csv")

# 1. Run bivariate analysis
source("1_bivariate_ordcats_correlation/bivariate_main.R")
bivariate_results <- run_bivariate_analysis(bes_data)

# 2. Test matrix properties  
source("2_correlation_matrices/matrices_main.R")
matrix_results <- run_matrix_analysis(bivariate_results$bounds_data)

# 3. Apply rescaling methods
source("3_fixes_and_rescaling/rescaling_main.R") 
rescaling_results <- run_rescaling_analysis(bivariate_results$bounds_data)

# 4. Generate comprehensive report
source("main_analysis.R")
final_results <- run_complete_analysis(bes_data)
```

### Advanced Rescaling Methods (`3_fixes_and_rescaling/2_advanced_rescaling/sophisticated_methods.R`)

```r
# sophisticated_methods.R
# Advanced rescaling approaches beyond simple linear mapping

#' Polynomial rescaling using higher-order transformations
#'
#' @param r Observed correlation
#' @param r_min Theoretical minimum
#' @param r_max Theoretical maximum
#' @param degree Polynomial degree (default: 2)
#' @return Rescaled correlation
rescale_correlation_polynomial <- function(r, r_min, r_max, degree = 2) {
  if (is.na(r) || is.na(r_min) || is.na(r_max)) return(NA)
  
  # Normalize r to [0, 1] interval
  r_normalized <- (r - r_min) / (r_max - r_min)
  
  # Apply polynomial transformation
  # For degree 2: emphasizes extreme values
  r_poly <- sign(r_normalized - 0.5) * abs(r_normalized - 0.5)^degree + 0.5
  
  # Map back to [-1, 1]
  r_rescaled <- 2 * r_poly - 1
  
  return(r_rescaled)
}

#' Logistic rescaling for smooth S-curve transformation
#'
#' @param r Observed correlation
#' @param r_min Theoretical minimum
#' @param r_max Theoretical maximum
#' @param steepness Steepness parameter (default: 5)
#' @return Rescaled correlation
rescale_correlation_logistic <- function(r, r_min, r_max, steepness = 5) {
  if (is.na(r) || is.na(r_min) || is.na(r_max)) return(NA)
  
  # Normalize to [0, 1]
  r_normalized <- (r - r_min) / (r_max - r_min)
  
  # Apply logistic transformation
  # Centers around 0.5 with specified steepness
  r_logistic <- 1 / (1 + exp(-steepness * (r_normalized - 0.5)))
  
  # Map to [-1, 1]
  r_rescaled <- 2 * r_logistic - 1
  
  return(r_rescaled)
}

#' Quantile-based rescaling using empirical distribution
#'
#' @param r Observed correlation
#' @param all_correlations Vector of all observed correlations for quantile estimation
#' @return Rescaled correlation based on quantile position
rescale_correlation_quantile <- function(r, all_correlations) {
  if (is.na(r)) return(NA)
  
  # Calculate quantile position
  quantile_pos <- mean(all_correlations <= r, na.rm = TRUE)
  
  # Map quantile to [-1, 1] using inverse normal transformation
  if (quantile_pos == 0) quantile_pos <- 1e-6
  if (quantile_pos == 1) quantile_pos <- 1 - 1e-6
  
  r_rescaled <- qnorm(quantile_pos)
  
  # Bound to [-1, 1]
  r_rescaled <- pmax(-1, pmin(1, r_rescaled / 3))  # Divide by 3 for reasonable scaling
  
  return(r_rescaled)
}

#' Copula-based rescaling preserving rank structure
#'
#' @param correlations_df Data frame with correlations and bounds
#' @return Data frame with copula-rescaled correlations
rescale_correlation_copula <- function(correlations_df) {
  # Extract variables
  r_obs <- correlations_df$observed_r
  r_min <- correlations_df$r_min  
  r_max <- correlations_df$r_max
  
  # Calculate uniform margins for observed correlations
  u_obs <- rank(r_obs, na.last = "keep") / (sum(!is.na(r_obs)) + 1)
  
  # Calculate uniform margins for theoretical bounds
  u_min <- rank(r_min, na.last = "keep") / (sum(!is.na(r_min)) + 1)
  u_max <- rank(r_max, na.last = "keep") / (sum(!is.na(r_max)) + 1)
  
  # Apply Gaussian copula transformation
  z_obs <- qnorm(u_obs)
  z_min <- qnorm(u_min)
  z_max <- qnorm(u_max)
  
  # Rescale in Gaussian space
  z_rescaled <- ifelse(z_obs < 0,
                       -1 * abs(z_obs) / abs(z_min),  # Negative correlations
                       1 * abs(z_obs) / abs(z_max))    # Positive correlations
  
  # Transform back to correlation scale
  u_rescaled <- pnorm(z_rescaled)
  r_rescaled <- 2 * u_rescaled - 1  # Map to [-1, 1]
  
  correlations_df$r_rescaled_copula <- r_rescaled
  return(correlations_df)
}

#' Matrix-aware rescaling that preserves positive semidefiniteness
#'
#' @param correlations_df Data frame with pairwise correlations
#' @param method Rescaling method to use
#' @param max_iterations Maximum optimization iterations
#' @return Data frame with matrix-aware rescaled correlations
rescale_correlation_matrix_aware <- function(correlations_df, method = "linear", max_iterations = 100) {
  # Apply initial rescaling
  if (method == "linear") {
    correlations_df$r_rescaled_temp <- mapply(rescale_correlation_simple,
                                            correlations_df$observed_r,
                                            correlations_df$r_min,
                                            correlations_df$r_max)
  } else if (method == "polynomial") {
    correlations_df$r_rescaled_temp <- mapply(rescale_correlation_polynomial,
                                            correlations_df$observed_r,
                                            correlations_df$r_min,
                                            correlations_df$r_max)
  } else if (method == "logistic") {
    correlations_df$r_rescaled_temp <- mapply(rescale_correlation_logistic,
                                            correlations_df$observed_r,
                                            correlations_df$r_min,
                                            correlations_df$r_max)
  }
  
  # Check if rescaled correlations form valid correlation matrix
  tryCatch({
    # Construct matrix with rescaled correlations
    test_matrix <- construct_correlation_matrix(correlations_df, use_rescaled = FALSE, enforce_psd = FALSE)
    matrix_props <- test_matrix_properties_comprehensive(test_matrix)
    
    if (!matrix_props$is_positive_semidefinite) {
      # Apply matrix correction to ensure PSD
      corrected_matrix <- make_matrix_psd(test_matrix)
      
      # Extract corrected correlations back to data frame
      all_vars <- unique(c(correlations_df$var1, correlations_df$var2))
      var_indices <- setNames(1:length(all_vars), all_vars)
      
      for (i in 1:nrow(correlations_df)) {
        var1 <- correlations_df$var1[i]
        var2 <- correlations_df$var2[i]
        
        if (var1 %in% all_vars && var2 %in% all_vars) {
          idx1 <- var_indices[var1]
          idx2 <- var_indices[var2]
          correlations_df$r_rescaled_temp[i] <- corrected_matrix[idx1, idx2]
        }
      }
    }
    
    correlations_df$r_rescaled_matrix_aware <- correlations_df$r_rescaled_temp
    
  }, error = function(e) {
    # If matrix construction fails, fall back to simple rescaling
    warning("Matrix-aware rescaling failed, using simple rescaling: ", e$message)
    correlations_df$r_rescaled_matrix_aware <- correlations_df$r_rescaled_temp
  })
  
  return(correlations_df)
}

#' Compare multiple rescaling methods
#'
#' @param correlations_df Data frame with correlation bounds
#' @return Data frame with multiple rescaling approaches
compare_rescaling_methods <- function(correlations_df) {
  # Apply different rescaling methods
  correlations_df$r_rescaled_linear <- mapply(rescale_correlation_simple,
                                            correlations_df$observed_r,
                                            correlations_df$r_min,
                                            correlations_df$r_max)
  
  correlations_df$r_rescaled_polynomial <- mapply(rescale_correlation_polynomial,
                                                correlations_df$observed_r,
                                                correlations_df$r_min,
                                                correlations_df$r_max)
  
  correlations_df$r_rescaled_logistic <- mapply(rescale_correlation_logistic,
                                               correlations_df$observed_r,
                                               correlations_df$r_min,
                                               correlations_df$r_max)
  
  correlations_df$r_rescaled_quantile <- mapply(rescale_correlation_quantile,
                                               correlations_df$observed_r,
                                               MoreArgs = list(all_correlations = correlations_df$observed_r))
  
  # Apply copula rescaling
  correlations_df <- rescale_correlation_copula(correlations_df)
  
  # Calculate method comparisons
  methods <- c("linear", "polynomial", "logistic", "quantile", "copula")
  method_cols <- paste0("r_rescaled_", methods)
  
  # Correlation between methods
  method_correlations <- cor(correlations_df[method_cols], use = "complete.obs")
  
  # Summary statistics by method
  method_summary <- correlations_df[method_cols] %>%
    summarise_all(list(
      mean = ~mean(.x, na.rm = TRUE),
      sd = ~sd(.x, na.rm = TRUE),
      min = ~min(.x, na.rm = TRUE),
      max = ~max(.x, na.rm = TRUE)
    ))
  
  return(list(
    data = correlations_df,
    method_correlations = method_correlations,
    method_summary = method_summary
  ))
}
```

### BES Data Loading (`utilities/data_loading.R`)

```r
# data_loading.R
# Functions for loading and preprocessing BES data

#' Load BES data with proper formatting and validation
#'
#' @param file_path Path to BES data file
#' @param validate_data Perform data validation
#' @return Cleaned BES data frame
load_bes_data <- function(file_path = "data/raw/bes_data.csv", validate_data = TRUE) {
  if (!file.exists(file_path)) {
    stop("BES data file not found: ", file_path)
  }
  
  cat("Loading BES data from:", file_path, "\n")
  
  # Load data
  bes_data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  cat("Loaded", nrow(bes_data), "rows and", ncol(bes_data), "columns\n")
  
  # Basic data cleaning
  bes_data <- clean_bes_data(bes_data)
  
  # Validation if requested
  if (validate_data) {
    validation_results <- validate_bes_dataset(bes_data)
    if (!validation_results$valid) {
      warning("Data validation failed:\n", paste(validation_results$errors, collapse = "\n"))
    } else {
      cat("Data validation passed\n")
    }
  }
  
  return(bes_data)
}

#' Clean and standardize BES data
#'
#' @param bes_data Raw BES data frame
#' @return Cleaned data frame
clean_bes_data <- function(bes_data) {
  # Ensure required columns exist
  required_cols <- c("var1", "var2", "corr", "nobs", "var1cats", "var2cats")
  missing_cols <- setdiff(required_cols, names(bes_data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Convert variable identifiers to character
  bes_data$var1 <- as.character(bes_data$var1)
  bes_data$var2 <- as.character(bes_data$var2)
  
  # Ensure numeric columns are properly formatted
  numeric_cols <- c("corr", "nobs", "var1cats", "var2cats")
  for (col in numeric_cols) {
    if (col %in% names(bes_data)) {
      bes_data[[col]] <- as.numeric(bes_data[[col]])
    }
  }
  
  # Remove rows with invalid data
  valid_rows <- !is.na(bes_data$corr) & 
                !is.na(bes_data$nobs) & 
                bes_data$nobs > 0 &
                abs(bes_data$corr) <= 1
  
  if (sum(!valid_rows) > 0) {
    cat("Removing", sum(!valid_rows), "rows with invalid data\n")
    bes_data <- bes_data[valid_rows, ]
  }
  
  # Sort by variable pairs for consistent ordering
  bes_data <- bes_data[order(bes_data$var1, bes_data$var2), ]
  
  return(bes_data)
}

#' Validate entire BES dataset
#'
#' @param bes_data BES data frame
#' @return Validation results
validate_bes_dataset <- function(bes_data) {
  errors <- c()
  warnings <- c()
  
  # Check data dimensions
  if (nrow(bes_data) == 0) {
    errors <- c(errors, "Dataset is empty")
  }
  
  # Check for required columns
  required_cols <- c("var1", "var2", "corr", "nobs", "var1cats", "var2cats")
  missing_cols <- setdiff(required_cols, names(bes_data))
  if (length(missing_cols) > 0) {
    errors <- c(errors, paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Validate correlations
  invalid_corr <- sum(abs(bes_data$corr) > 1, na.rm = TRUE)
  if (invalid_corr > 0) {
    errors <- c(errors, paste(invalid_corr, "correlations outside [-1, 1] range"))
  }
  
  # Check sample sizes
  small_samples <- sum(bes_data$nobs < 30, na.rm = TRUE)
  if (small_samples > 0) {
    warnings <- c(warnings, paste(small_samples, "pairs have small sample sizes (< 30)"))
  }
  
  # Check frequency/proportion columns
  freq_cols <- grep("freq[0-9]+", names(bes_data), value = TRUE)
  prop_cols <- grep("prop[0-9]+", names(bes_data), value = TRUE)
  
  if (length(freq_cols) == 0 && length(prop_cols) == 0) {
    errors <- c(errors, "No frequency or proportion columns found")
  }
  
  # Check for duplicate variable pairs
  if (anyDuplicated(paste(bes_data$var1, bes_data$var2))) {
    errors <- c(errors, "Duplicate variable pairs found")
  }
  
  return(list(
    valid = length(errors) == 0,
    errors = errors,
    warnings = warnings,
    n_rows = nrow(bes_data),
    n_valid_pairs = sum(!is.na(bes_data$corr))
  ))
}

#' Create example BES dataset for testing
#'
#' @param n_pairs Number of variable pairs to generate
#' @param seed Random seed
#' @return Simulated BES data frame
create_example_bes_data <- function(n_pairs = 100, seed = 123) {
  set.seed(seed)
  
  # Generate variable pairs
  n_vars <- ceiling(sqrt(2 * n_pairs))
  var_combinations <- expand.grid(var1 = 1:n_vars, var2 = 1:n_vars)
  var_combinations <- var_combinations[var_combinations$var1 < var_combinations$var2, ]
  var_combinations <- var_combinations[1:min(n_pairs, nrow(var_combinations)), ]
  
  # Generate realistic BES-like data
  bes_data <- data.frame(
    var1 = var_combinations$var1,
    var2 = var_combinations$var2,
    var1cats = sample(3:7, nrow(var_combinations), replace = TRUE),
    var2cats = sample(3:9, nrow(var_combinations), replace = TRUE),
    nobs = sample(800:2000, nrow(var_combinations), replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Generate marginal distributions and correlations
  for (i in 1:nrow(bes_data)) {
    # Generate marginal frequencies
    var1_freqs <- rmultinom(1, bes_data$nobs[i], runif(bes_data$var1cats[i]))
    var2_freqs <- rmultinom(1, bes_data$nobs[i], runif(bes_data$var2cats[i]))
    
    # Add frequency columns (assuming variables start at 1)
    for (j in 1:bes_data$var1cats[i]) {
      bes_data[i, paste0("var1freq", j)] <- var1_freqs[j]
    }
    for (j in 1:bes_data$var2cats[i]) {
      bes_data[i, paste0("var2freq", j)] <- var2_freqs[j]
    }
    
    # Calculate theoretical bounds
    r_min <- min_corr_bound(var1_freqs, var2_freqs)
    r_max <- max_corr_bound(var1_freqs, var2_freqs)
    
    # Generate realistic correlation within bounds
    correlation_position <- runif(1)
    bes_data$corr[i] <- r_min + correlation_position * (r_max - r_min)
  }
  
  # Fill missing frequency columns with zeros
  max_cats <- max(c(bes_data$var1cats, bes_data$var2cats))
  for (var_prefix in c("var1freq", "var2freq")) {
    for (j in 1:max_cats) {
      col_name <- paste0(var_prefix, j)
      if (!col_name %in% names(bes_data)) {
        bes_data[[col_name]] <- 0
      }
    }
  }
  
  return(bes_data)
}
```

## Additional Integration Files

### Module Coordinators

Each major module needs a coordinator file that sources its components:

**`1_bivariate_ordcats_correlation/bivariate_main.R`:**
```r
# Source all bivariate analysis components
source("1_bivariate_ordcats_correlation/1_rmin_rmax_rhat/1_monte_carlo_simulation/core_bounds_functions.R")
source("1_bivariate_ordcats_correlation/1_rmin_rmax_rhat/2_bes_illustrative_example/bes_data_analysis.R") 
source("1_bivariate_ordcats_correlation/2_asymmetry_analysis/asymmetry_measures.R")
source("1_bivariate_ordcats_correlation/3_visualization/bounds_visualization.R")

run_bivariate_analysis <- function(bes_data, config = list()) {
  # Coordinate all bivariate analyses
  bounds_results <- analyze_all_bes_bounds(bes_data, nsim = config$nsim %||% 1000)
  asymmetry_results <- analyze_bounds_asymmetry_comprehensive(bounds_results)
  return(list(bounds_data = bounds_results, asymmetry = asymmetry_results))
}
```

**`2_correlation_matrices/matrices_main.R`:**
```r
# Source matrix analysis components  
source("2_correlation_matrices/1_matrix_properties/matrix_diagnostics.R")
source("2_correlation_matrices/2_matrix_construction/random_trials.R")

run_matrix_analysis <- function(bounds_data, config = list()) {
  # Coordinate matrix property testing
  trials <- test_random_correlation_matrices(bounds_data, 
                                           n_trials = config$matrix_trials %||% 100)
  summary <- summarize_matrix_tests(trials)
  return(list(trials = trials, summary = summary))
}
```

**`3_fixes_and_rescaling/rescaling_main.R`:**
```r
# Source rescaling components
source("3_fixes_and_rescaling/1_simple_rescaling/linear_rescaling.R")
source("3_fixes_and_rescaling/2_advanced_rescaling/sophisticated_methods.R")

run_rescaling_analysis <- function(bounds_data, config = list()) {
  # Apply and compare rescaling methods
  simple_rescaling <- apply_simple_rescaling(bounds_data)
  advanced_methods <- compare_rescaling_methods(bounds_data)
  return(list(simple = simple_rescaling, advanced = advanced_methods))
}
```

This comprehensive reorganization maintains all your sophisticated analytical methods while making the codebase:

1. **More maintainable** - clear module boundaries
2. **More extensible** - easy to add new methods within existing structure  
3. **More collaborative** - different team members can work on different modules
4. **More testable** - individual components can be tested in isolation
5. **More professional** - follows modern R development practices

The structure preserves your investment in complex correlation bounds theory while providing a foundation for future methodological developments. else {
    countsX <- marginalX
    sample_size <- sum(countsX)
  }
  if (sum(marginalY) <= 1.1) {
    countsY <- as.vector(rmultinom(1, size = sample_size, prob = marginalY))
  } else {
    countsY <- marginalY
    if (sum(countsY) != sample_size) stop("Counts for X and Y must be equal.")
  }
  
  # Define ordinal levels (0, 1, 2, ...)
  x_levels <- 0:(length(countsX) - 1)
  y_levels <- 0:(length(countsY) - 1)
  
  # Compute proportions and moments
  pX <- countsX / sum(countsX)
  pY <- countsY / sum(countsY)
  muX <- sum(x_levels * pX)
  muY <- sum(y_levels * pY)
  sigmaX <- sqrt(sum(x_levels^2 * pX) - muX^2)
  sigmaY <- sqrt(sum(y_levels^2 * pY) - muY^2)
  
  # Compute joint distribution using cumulative probabilities (comonotonic coupling)
  cum_pX <- cumsum(pX)
  cum_pY <- cumsum(pY)
  joint <- matrix(0, nrow = length(x_levels), ncol = length(y_levels))
  
  i <- 1; j <- 1
  cumX_prev <- 0; cumY_prev <- 0
  while (i <= length(x_levels) && j <= length(y_levels)) {
    currX <- cum_pX[i]
    currY <- cum_pY[j]
    joint[i, j] <- min(currX, currY) - max(cumX_prev, cumY_prev)
    
    if (currX < currY) {
      cumX_prev <- currX
      i <- i + 1
    } else if (currX > currY) {
      cumY_prev <- currY
      j <- j + 1
    } else {
      cumX_prev <- currX
      cumY_prev <- currY
      i <- i + 1
      j <- j + 1
    }
  }
  
  # Compute E[XY] based on the joint distribution
  E_XY <- sum(outer(x_levels, y_levels, FUN = "*") * joint)
  
  # Calculate the maximum Pearson correlation
  r_max <- (E_XY - muX * muY) / (sigmaX * sigmaY)
  return(r_max)
}

#' Compute minimum correlation bound using anti-comonotonic pairing
#' 
#' @param marginalX Vector of marginal probabilities or counts for variable X
#' @param marginalY Vector of marginal probabilities or counts for variable Y
#' @param sample_size Sample size (used if marginals are probabilities)
#' @return Minimum possible Pearson correlation coefficient
min_corr_bound <- function(marginalX, marginalY, sample_size = 10000) {
  # Convert probabilities to counts if needed
  if (sum(marginalX) <= 1.1) {
    countsX <- as.vector(rmultinom(1, size = sample_size, prob = marginalX))
  } else {
    countsX <- marginalX
    sample_size <- sum(countsX)
  }
  if (sum(marginalY) <= 1.1) {
    countsY <- as.vector(rmultinom(1, size = sample_size, prob = marginalY))
  } else {
    countsY <- marginalY
    if (sum(countsY) != sample_size) stop("Counts for X and Y must be equal.")
  }
  
  # Define ordinal levels (0, 1, 2, ...)
  x_levels <- 0:(length(countsX) - 1)
  y_levels <- 0:(length(countsY) - 1)
  
  # Create full sample vectors based on counts
  x_vec <- rep(x_levels, times = countsX)
  y_vec <- rep(y_levels, times = countsY)
  
  # Pair highest X with lowest Y (anti-comonotonic pairing)
  x_sorted <- sort(x_vec, decreasing = TRUE)
  y_sorted <- sort(y_vec, decreasing = FALSE)
  
  # Calculate correlation
  r_min <- cor(x_sorted, y_sorted)
  return(r_min)
}

#' Simulate permutation distribution for null hypothesis testing
#'
#' @param marginalX Vector of marginal probabilities or counts for variable X
#' @param marginalY Vector of marginal probabilities or counts for variable Y
#' @param nsim Number of simulations
#' @param sample_size Sample size (used if marginals are probabilities)
#' @return Vector of simulated Pearson correlation coefficients
simulate_permutation_distribution <- function(marginalX, marginalY, nsim = 1000, sample_size = 10000) {
  # Convert probabilities to counts if needed
  if (sum(marginalX) <= 1.1) {
    countsX <- as.vector(rmultinom(1, size = sample_size, prob = marginalX))
  } else {
    countsX <- marginalX
    sample_size <- sum(countsX)
  }
  if (sum(marginalY) <= 1.1) {
    countsY <- as.vector(rmultinom(1, size = sample_size, prob = marginalY))
  } else {
    countsY <- marginalY
    if (sum(countsY) != sample_size) stop("Counts for X and Y must be equal.")
  }
  
  # Define ordinal levels (0, 1, 2, ...)
  x_levels <- 0:(length(countsX) - 1)
  y_levels <- 0:(length(countsY) - 1)
  
  # Create sample vectors based on counts
  x_vec <- rep(x_levels, times = countsX)
  y_vec <- rep(y_levels, times = countsY)
  
  # Initialize storage for simulated correlations
  r_vals <- numeric(nsim)
  for (i in 1:nsim) {
    # Randomly permute both vectors to simulate independence
    x_perm <- sample(x_vec, size = length(x_vec), replace = FALSE)
    y_perm <- sample(y_vec, size = length(y_vec), replace = FALSE)
    r_vals[i] <- cor(x_perm, y_perm)
  }
  return(r_vals)
}

#' Convert probabilities to counts ensuring exact sample size
#'
#' @param probs Vector of probabilities
#' @param sample_size Desired sample size
#' @return Vector of counts with sum equal to sample_size
probs_to_counts <- function(probs, sample_size) {
  counts <- round(probs * sample_size)
  diff <- sample_size - sum(counts)
  
  if(diff != 0) {
    # Adjust the count of the largest category to ensure total equals sample_size
    i_max <- which.max(counts)
    counts[i_max] <- counts[i_max] + diff
  }
  
  return(counts)
}
```

### 2. BES Data Analysis (`1_bivariate_ordcats_correlation/1_rmin_rmax_rhat/2_bes_illustrative_example/bes_data_analysis.R`)

```r
# bes_data_analysis.R
# Synthesized from correlation_bounds_bes.R with complete BES-specific functionality

#' Analyze correlation bounds for a BES ordinal variable pair
#'
#' @param row Data frame row containing BES marginal frequencies/proportions
#' @param nsim Number of simulations for permutation distribution
#' @param use_prop Logical; use proportions instead of frequencies
#' @param return_simulations Logical; return simulated r values
#' @return List with bounds analysis results
analyze_bes_corr_bounds <- function(row, nsim = 1000, use_prop = FALSE, return_simulations = FALSE) {
  # Extract the number of categories for each variable
  var1_cats <- row$var1cats
  var2_cats <- row$var2cats
  n_obs <- row$nobs
  
  # Extract marginal distributions with proper category handling
  var1_marginals <- extract_bes_marginals(row, "var1", var1_cats, use_prop)
  var2_marginals <- extract_bes_marginals(row, "var2", var2_cats, use_prop)
  
  # Validate marginals
  if (sum(var1_marginals) == 0 || sum(var2_marginals) == 0) {
    warning("Empty marginal distributions detected.")
    return(create_error_result(row, "Empty marginal distributions"))
  }
  
  # Calculate theoretical bounds
  r_min <- min_corr_bound(var1_marginals, var2_marginals, sample_size = n_obs)
  r_max <- max_corr_bound(var1_marginals, var2_marginals, sample_size = n_obs)
  
  # Simulate permutation distribution
  r_sim <- simulate_permutation_distribution(var1_marginals, var2_marginals, 
                                           nsim = nsim, sample_size = n_obs)
  
  # Calculate 95% confidence interval
  ci_bounds <- quantile(r_sim, probs = c(0.025, 0.975))
  
  # Rescale the observed correlation
  r_rescaled <- rescale_correlation_simple(row$corr, r_min, r_max)
  
  # Prepare comprehensive results
  result <- list(
    r_min = r_min,
    r_max = r_max,
    ci_lower = ci_bounds[1],
    ci_upper = ci_bounds[2],
    observed_r = row$corr,
    r_rescaled = r_rescaled,
    success = TRUE,
    var1_categories = var1_cats,
    var2_categories = var2_cats,
    n_obs = n_obs,
    var1_marginals = var1_marginals,
    var2_marginals = var2_marginals
  )
  
  # Add simulations if requested
  if (return_simulations) {
    result$simulations <- r_sim
  }
  
  return(result)
}

#' Extract marginal distributions from BES data row
#'
#' @param row Data frame row with BES variable information
#' @param var_prefix Variable prefix ("var1" or "var2")
#' @param n_cats Number of categories
#' @param use_prop Use proportions instead of frequencies
#' @return Vector of marginal counts/proportions
extract_bes_marginals <- function(row, var_prefix, n_cats, use_prop = FALSE) {
  suffix <- if (use_prop) "prop" else "freq"
  marginals <- c()
  
  # Determine starting category (0 or 1) by checking for non-zero values
  var_start <- determine_bes_start_category(row, var_prefix, suffix)
  var_end <- var_start + n_cats - 1
  
  for (i in var_start:var_end) {
    col_name <- paste0(var_prefix, suffix, i)
    if (col_name %in% names(row) && !is.na(row[[col_name]])) {
      marginals <- c(marginals, row[[col_name]])
    } else {
      marginals <- c(marginals, 0)  # Add 0 for missing categories
    }
  }
  
  # Remove NAs and ensure non-negative values
  marginals <- pmax(0, marginals, na.rm = TRUE)
  return(marginals)
}

#' Determine starting category for BES variables (0 or 1)
#'
#' @param row Data frame row
#' @param var_prefix Variable prefix
#' @param suffix "freq" or "prop"
#' @return Starting category index
determine_bes_start_category <- function(row, var_prefix, suffix) {
  # Check if category 0 exists and has non-zero values
  col0_name <- paste0(var_prefix, suffix, "0")
  if (col0_name %in% names(row) && !is.na(row[[col0_name]]) && row[[col0_name]] > 0) {
    return(0)
  } else {
    return(1)
  }
}

#' Apply bounds analysis to entire BES dataset
#'
#' @param df BES data frame with variable pairs
#' @param nsim Number of simulations per pair
#' @param use_prop Use proportions instead of frequencies
#' @param return_simulations Include simulation results
#' @param progress Show progress indicator
#' @return Data frame with bounds analysis results
analyze_all_bes_bounds <- function(df, nsim = 1000, use_prop = FALSE, 
                                  return_simulations = FALSE, progress = TRUE) {
  results_list <- list()
  n_rows <- nrow(df)
  
  if (progress) cat("Processing", n_rows, "BES variable pairs...\n")
  
  for (i in 1:n_rows) {
    if (progress && i %% 100 == 0) cat("Processed", i, "of", n_rows, "pairs\n")
    
    # Analyze bounds for this pair
    bounds_result <- analyze_bes_corr_bounds(df[i, ], nsim = nsim, 
                                           use_prop = use_prop, 
                                           return_simulations = return_simulations)
    
    # Add original variable identifiers
    bounds_result$var1 <- df$var1[i]
    bounds_result$var2 <- df$var2[i]
    
    results_list[[i]] <- bounds_result
  }
  
  # Convert list to data frame
  results_df <- convert_bounds_results_to_df(results_list, return_simulations)
  
  if (progress) cat("Bounds analysis complete!\n")
  return(results_df)
}

#' Convert results list to structured data frame
convert_bounds_results_to_df <- function(results_list, include_simulations = FALSE) {
  results_df <- data.frame(
    var1 = sapply(results_list, function(x) x$var1),
    var2 = sapply(results_list, function(x) x$var2),
    r_min = sapply(results_list, function(x) x$r_min),
    r_max = sapply(results_list, function(x) x$r_max),
    ci_lower = sapply(results_list, function(x) x$ci_lower),
    ci_upper = sapply(results_list, function(x) x$ci_upper),
    observed_r = sapply(results_list, function(x) x$observed_r),
    r_rescaled = sapply(results_list, function(x) x$r_rescaled),
    var1_categories = sapply(results_list, function(x) x$var1_categories),
    var2_categories = sapply(results_list, function(x) x$var2_categories),
    n_obs = sapply(results_list, function(x) x$n_obs),
    success = sapply(results_list, function(x) x$success)
  )
  
  # Add simulations as list column if requested
  if (include_simulations) {
    results_df$simulations <- lapply(results_list, function(x) x$simulations)
  }
  
  return(results_df)
}

#' Create error result for failed analyses
create_error_result <- function(row, message) {
  list(
    r_min = NA, r_max = NA,
    ci_lower = NA, ci_upper = NA,
    observed_r = row$corr,
    r_rescaled = NA,
    success = FALSE,
    message = message
  )
}
```

### 3. Simple Rescaling Functions (`3_fixes_and_rescaling/1_simple_rescaling/linear_rescaling.R`)

```r
# linear_rescaling.R
# Simple rescaling methods synthesized from multiple source files

#' Simple linear rescaling to [-1, 1] interval
#'
#' @param r Observed correlation
#' @param r_min Theoretical minimum correlation
#' @param r_max Theoretical maximum correlation
#' @return Rescaled correlation in [-1, 1]
rescale_correlation_simple <- function(r, r_min, r_max) {
  if (is.na(r) || is.na(r_min) || is.na(r_max)) {
    return(NA)
  }
  
  if (r < 0) {
    # For negative correlations, rescale between r_min and 0
    if (r_min < 0) {
      return((-1 / r_min) * r)
    } else {
      return(-1)  # Edge case: r_min >= 0 but r < 0
    }
  } else if (r > 0) {
    # For positive correlations, rescale between 0 and r_max
    if (r_max > 0) {
      return((1 / r_max) * r)
    } else {
      return(1)  # Edge case: r_max <= 0 but r > 0
    }
  } else {
    return(0)  # r exactly equals 0
  }
}

#' Apply rescaling to entire results data frame
#'
#' @param results_df Data frame with bounds analysis results
#' @return Data frame with added rescaled correlations
apply_simple_rescaling <- function(results_df) {
  results_df$r_rescaled <- mapply(
    rescale_correlation_simple,
    results_df$observed_r,
    results_df$r_min,
    results_df$r_max
  )
  
  # Add rescaling diagnostics
  results_df$rescaling_factor_pos <- ifelse(results_df$r_max > 0, 1 / results_df$r_max, NA)
  results_df$rescaling_factor_neg <- ifelse(results_df$r_min < 0, -1 / results_df$r_min, NA)
  results_df$rescaling_magnitude <- abs(results_df$r_rescaled) - abs(results_df$observed_r)
  
  return(results_df)
}

#' Analyze rescaling effects across variable pairs
#'
#' @param results_df Data frame with original and rescaled correlations
#' @return List with rescaling analysis summary
analyze_rescaling_effects <- function(results_df) {
  # Basic statistics
  rescaling_summary <- list(
    mean_magnitude_change = mean(results_df$rescaling_magnitude, na.rm = TRUE),
    median_magnitude_change = median(results_df$rescaling_magnitude, na.rm = TRUE),
    prop_increased = mean(results_df$rescaling_magnitude > 0, na.rm = TRUE),
    prop_decreased = mean(results_df$rescaling_magnitude < 0, na.rm = TRUE),
    max_increase = max(results_df$rescaling_magnitude, na.rm = TRUE),
    max_decrease = min(results_df$rescaling_magnitude, na.rm = TRUE)
  )
  
  # Categorize rescaling effects
  results_df$rescaling_effect <- cut(
    results_df$rescaling_magnitude,
    breaks = c(-Inf, -0.1, -0.05, 0.05, 0.1, Inf),
    labels = c("Large Decrease", "Small Decrease", "Minimal Change", 
               "Small Increase", "Large Increase"),
    include.lowest = TRUE
  )
  
  # Effect distribution
  effect_distribution <- table(results_df$rescaling_effect)
  
  return(list(
    summary = rescaling_summary,
    detailed_results = results_df,
    effect_distribution = effect_distribution
  ))
}

#' Compare rescaling across different variable characteristics
#'
#' @param results_df Data frame with rescaling results
#' @return Analysis of rescaling by variable properties
analyze_rescaling_by_characteristics <- function(results_df) {
  # By number of categories
  results_df$total_categories <- results_df$var1_categories + results_df$var2_categories
  results_df$category_difference <- abs(results_df$var1_categories - results_df$var2_categories)
  
  # By bounds range
  results_df$bounds_range <- results_df$r_max - results_df$r_min
  results_df$bounds_asymmetry <- results_df$r_max + results_df$r_min
  
  # Summary by categories
  category_analysis <- results_df %>%
    group_by(total_categories) %>%
    summarise(
      n_pairs = n(),
      mean_rescaling_magnitude = mean(rescaling_magnitude, na.rm = TRUE),
      mean_bounds_range = mean(bounds_range, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Summary by bounds characteristics
  bounds_analysis <- results_df %>%
    mutate(
      range_category = cut(bounds_range, breaks = 3, labels = c("Narrow", "Medium", "Wide")),
      asymmetry_category = cut(abs(bounds_asymmetry), breaks = 3, labels = c("Symmetric", "Moderate", "Asymmetric"))
    ) %>%
    group_by(range_category, asymmetry_category) %>%
    summarise(
      n_pairs = n(),
      mean_rescaling_magnitude = mean(rescaling_magnitude, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(list(
    by_categories = category_analysis,
    by_bounds = bounds_analysis,
    detailed_data = results_df
  ))
}
```

### 4. Matrix Construction and Testing (`2_correlation_matrices/1_matrix_properties/matrix_diagnostics.R`)

```r
# matrix_diagnostics.R
# Complete matrix analysis functions synthesized from correlation_matrix_test.R

#' Construct correlation matrix from pairwise correlations
#'
#' @param correlations Data frame with var1, var2, and correlation columns
#' @param use_rescaled Use rescaled correlations instead of observed
#' @param enforce_psd Ensure positive semidefiniteness
#' @return Correlation matrix
construct_correlation_matrix <- function(correlations, use_rescaled = FALSE, enforce_psd = TRUE) {
  # Identify all unique variables
  all_vars <- unique(c(correlations$var1, correlations$var2))
  n_vars <- length(all_vars)
  
  # Create variable index mapping
  var_indices <- setNames(1:n_vars, all_vars)
  
  # Initialize correlation matrix
  corr_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)
  diag(corr_matrix) <- 1
  
  # Choose correlation column
  corr_column <- ifelse(use_rescaled, "r_rescaled", "observed_r")
  
  # Fill correlation matrix
  for (i in 1:nrow(correlations)) {
    var1 <- correlations$var1[i]
    var2 <- correlations$var2[i]
    corr_value <- correlations[[corr_column]][i]
    
    if (!is.na(corr_value) && var1 %in% all_vars && var2 %in% all_vars) {
      idx1 <- var_indices[var1]
      idx2 <- var_indices[var2]
      
      # Fill both triangles for symmetry
      corr_matrix[idx1, idx2] <- corr_value
      corr_matrix[idx2, idx1] <- corr_value
    }
  }
  
  # Enforce positive semidefiniteness if requested
  if (enforce_psd) {
    corr_matrix <- make_matrix_psd(corr_matrix)
  }
  
  # Add variable names
  rownames(corr_matrix) <- all_vars
  colnames(corr_matrix) <- all_vars
  
  return(corr_matrix)
}

#' Make matrix positive semidefinite using eigenvalue adjustment
#'
#' @param matrix Input correlation matrix
#' @param min_eigenvalue Minimum allowed eigenvalue
#' @return Adjusted positive semidefinite matrix
make_matrix_psd <- function(matrix, min_eigenvalue = 1e-8) {
  # Ensure exact symmetry
  matrix <- (matrix + t(matrix)) / 2
  
  # Eigendecomposition
  eigen_decomp <- eigen(matrix, symmetric = TRUE)
  eigenvalues <- eigen_decomp$values
  eigenvectors <- eigen_decomp$vectors
  
  # Adjust negative eigenvalues
  if (min(eigenvalues) < 0) {
    eigenvalues[eigenvalues < min_eigenvalue] <- min_eigenvalue
    
    # Reconstruct matrix
    matrix <- eigenvectors %*% diag(eigenvalues) %*% t(eigenvectors)
    
    # Rescale to ensure unit diagonal
    D <- diag(1/sqrt(diag(matrix)))
    matrix <- D %*% matrix %*% D
    
    # Final symmetry adjustment
    matrix <- (matrix + t(matrix)) / 2
    diag(matrix) <- 1
  }
  
  return(matrix)
}

#' Comprehensive matrix property testing
#'
#' @param matrix Correlation matrix to test
#' @param tolerance Numerical tolerance for tests
#' @return List with detailed matrix properties
test_matrix_properties_comprehensive <- function(matrix, tolerance = 1e-10) {
  results <- list()
  
  # Basic structural properties
  results$is_symmetric <- isSymmetric(matrix, tol = tolerance)
  results$has_unit_diagonal <- all(abs(diag(matrix) - 1) < tolerance)
  
  # Ensure symmetry for eigenvalue computation
  if (!results$is_symmetric) {
    matrix <- (matrix + t(matrix)) / 2
  }
  
  # Eigenvalue analysis
  eigenvalues <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  results$eigenvalues <- eigenvalues
  results$min_eigenvalue <- min(eigenvalues)
  results$max_eigenvalue <- max(eigenvalues)
  results$n_negative_eigenvalues <- sum(eigenvalues < -tolerance)
  results$n_zero_eigenvalues <- sum(abs(eigenvalues) < tolerance)
  
  # Matrix definiteness
  results$is_positive_definite <- results$min_eigenvalue > tolerance
  results$is_positive_semidefinite <- results$min_eigenvalue > -tolerance
  results$is_invertible <- results$min_eigenvalue > tolerance
  
  # Condition number analysis
  if (results$is_invertible) {
    results$condition_number <- results$max_eigenvalue / results$min_eigenvalue
    results$condition_category <- classify_condition_number(results$condition_number)
    results$log_condition_number <- log10(results$condition_number)
  } else {
    results$condition_number <- Inf
    results$condition_category <- "Singular"
    results$log_condition_number <- Inf
  }
  
  # Matrix determinant and rank
  results$determinant <- prod(eigenvalues)
  results$log_determinant <- sum(log(eigenvalues[eigenvalues > tolerance]))
  results$rank <- sum(eigenvalues > tolerance)
  results$rank_deficiency <- nrow(matrix) - results$rank
  
  # Numerical stability measures
  results$spectral_radius <- max(abs(eigenvalues))
  results$frobenius_norm <- sqrt(sum(matrix^2))
  results$nuclear_norm <- sum(abs(eigenvalues))
  
  # Correlation-specific properties
  results$max_off_diagonal <- max(abs(matrix[upper.tri(matrix)]))
  results$min_off_diagonal <- min(abs(matrix[upper.tri(matrix)]))
  results$mean_correlation <- mean(matrix[upper.tri(matrix)])
  results$n_high_correlations <- sum(abs(matrix[upper.tri(matrix)]) > 0.8)
  
  # Overall validity assessment
  results$is_valid_correlation_matrix <- 
    results$is_symmetric && 
    results$has_unit_diagonal && 
    results$is_positive_semidefinite &&
    all(abs(matrix[upper.tri(matrix)]) <= 1 + tolerance)
  
  return(results)
}

#' Classify condition number for interpretation
#'
#' @param cond_num Condition number
#' @return Character classification
classify_condition_number <- function(cond_num) {
  if (is.infinite(cond_num)) return("Singular")
  if (cond_num <= 10) return("Well-conditioned")
  if (cond_num <= 100) return("Moderately ill-conditioned")
  if (cond_num <= 1000) return("Ill-conditioned")
  return("Severely ill-conditioned")
}

#' Compare properties between original and rescaled matrices
#'
#' @param original_matrix Original correlation matrix
#' @param rescaled_matrix Rescaled correlation matrix
#' @return Detailed comparison results
compare_matrix_properties <- function(original_matrix, rescaled_matrix) {
  # Test properties of both matrices
  orig_props <- test_matrix_properties_comprehensive(original_matrix)
  rescaled_props <- test_matrix_properties_comprehensive(rescaled_matrix)
  
  # Property comparisons
  comparison <- list(
    # Structural preservation
    symmetry_preserved = orig_props$is_symmetric == rescaled_props$is_symmetric,
    diagonal_preserved = orig_props$has_unit_diagonal == rescaled_props$has_unit_diagonal,
    
    # Definiteness changes
    psd_preserved = orig_props$is_positive_semidefinite == rescaled_props$is_positive_semidefinite,
    invertibility_preserved = orig_props$is_invertible == rescaled_props$is_invertible,
    psd_improved = !orig_props$is_positive_semidefinite && rescaled_props$is_positive_semidefinite,
    invertibility_improved = !orig_props$is_invertible && rescaled_props$is_invertible,
    
    # Numerical changes
    condition_ratio = if(is.finite(orig_props$condition_number) && is.finite(rescaled_props$condition_number)) {
      rescaled_props$condition_number / orig_props$condition_number
    } else { NA },
    
    eigenvalue_changes = list(
      min_eigenvalue_change = rescaled_props$min_eigenvalue - orig_props$min_eigenvalue,
      max_eigenvalue_change = rescaled_props$max_eigenvalue - orig_props$max_eigenvalue,
      eigenvalue_range_change = (rescaled_props$max_eigenvalue - rescaled_props$min_eigenvalue) - 
                               (orig_props$max_eigenvalue - orig_props$min_eigenvalue)
    ),
    
    # Correlation structure changes
    correlation_changes = list(
      mean_correlation_change = rescaled_props$mean_correlation - orig_props$mean_correlation,
      max_correlation_change = rescaled_props$max_off_diagonal - orig_props$max_off_diagonal,
      correlation_range_change = (rescaled_props$max_off_diagonal - rescaled_props$min_off_diagonal) - 
                                (orig_props$max_off_diagonal - orig_props$min_off_diagonal)
    )
  )
  
  # Overall assessment
  comparison$overall_properties_preserved <- 
    comparison$psd_preserved && comparison$invertibility_preserved
  
  comparison$numerical_stability_improved <- 
    !is.na(comparison$condition_ratio) && comparison$condition_ratio < 1
  
  return(list(
    original = orig_props,
    rescaled = rescaled_props,
    comparison = comparison
  ))
}
```

### 5. Matrix Testing Framework (`2_correlation_matrices/2_matrix_construction/random_trials.R`)

```r
# random_trials.R
# Random matrix testing framework from correlation_matrix_test.R

#' Run random correlation matrix trials
#'
#' @param results_df Data frame with correlation bounds results
#' @param n_trials Number of random trials to run
#' @param min_vars Minimum variables per trial
#' @param max_vars Maximum variables per trial
#' @param seed Random seed for reproducibility
#' @param enforce_psd Whether to enforce positive semidefiniteness
#' @return Data frame with trial results
test_random_correlation_matrices <- function(results_df, n_trials = 100, 
                                           min_vars = 3, max_vars = 15,
                                           seed = 123, enforce_psd = TRUE) {
  set.seed(seed)
  
  # Get all unique variables
  all_vars <- unique(c(results_df$var1, results_df$var2))
  n_all_vars <- length(all_vars)
  max_vars <- min(max_vars, n_all_vars)
  
  trial_results <- data.frame()
  
  for (trial in 1:n_trials) {
    # Randomly select variables and sample size
    n_vars <- sample(min_vars:max_vars, 1)
    selected_vars <- sample(all_vars, n_vars)
    
    # Get correlations for selected variables
    subset_corrs <- results_df %>%
      filter(var1 %in% selected_vars & var2 %in% selected_vars)
    
    # Skip if insufficient correlations
    if (nrow(subset_corrs) < (n_vars * (n_vars - 1) / 2 * 0.5)) {
      next  # Need at least 50% of possible pairs
    }
    
    # Construct both matrices
    orig_matrix <- tryCatch({
      construct_correlation_matrix(subset_corrs, use_rescaled = FALSE, enforce_psd = enforce_psd)
    }, error = function(e) NULL)
    
    rescaled_matrix <- tryCatch({
      construct_correlation_matrix(subset_corrs, use_rescaled = TRUE, enforce_psd = enforce_psd)
    }, error = function(e) NULL)
    
    if (is.null(orig_matrix) || is.null(rescaled_matrix)) next
    
    # Test matrix properties
    comparison_results <- compare_matrix_properties(orig_matrix, rescaled_matrix)
    
    # Store trial results
    trial_result <- data.frame(
      trial = trial,
      n_variables = n_vars,
      n_correlations = nrow(subset_corrs),
      enforce_psd = enforce_psd,
      
      # Original matrix properties
      orig_psd = comparison_results$original$is_positive_semidefinite,
      orig_invertible = comparison_results$original$is_invertible,
      orig_condition = comparison_results$original$condition_number,
      orig_min_eigenvalue = comparison_results$original$min_eigenvalue,
      orig_valid = comparison_results$original$is_valid_correlation_matrix,
      
      # Rescaled matrix properties  
      rescaled_psd = comparison_results$rescaled$is_positive_semidefinite,
      rescaled_invertible = comparison_results$rescaled$is_invertible,
      rescaled_condition = comparison_results$rescaled$condition_number,
      rescaled_min_eigenvalue = comparison_results$rescaled$min_eigenvalue,
      rescaled_valid = comparison_results$rescaled$is_valid_correlation_matrix,
      
      # Comparisons
      properties_preserved = comparison_results$comparison$overall_properties_preserved,
      psd_improved = comparison_results$comparison$psd_improved,
      invertibility_improved = comparison_results$comparison$invertibility_improved,
      condition_ratio = comparison_results$comparison$condition_ratio,
      
      stringsAsFactors = FALSE
    )
    
    trial_results <- rbind(trial_results, trial_result)
  }
  
  return(trial_results)
}

#' Summarize matrix testing results
#'
#' @param test_results Output from test_random_correlation_matrices
#' @return Summary statistics
summarize_matrix_tests <- function(test_results) {
  if (nrow(test_results) == 0) {
    return(list(message = "No valid trials completed"))
  }
  
  summary_stats <- list(
    n_trials = nrow(test_results),
    
    # Original matrix statistics
    orig_psd_percent = 100 * mean(test_results$orig_psd, na.rm = TRUE),
    orig_invertible_percent = 100 * mean(test_results$orig_invertible, na.rm = TRUE),
    orig_valid_percent = 100 * mean(test_results$orig_valid, na.rm = TRUE),
    orig_mean_condition = mean(test_results$orig_condition[is.finite(test_results$orig_condition)], na.rm = TRUE),
    
    # Rescaled matrix statistics
    rescaled_psd_percent = 100 * mean(test_results$rescaled_psd, na.rm = TRUE),
    rescaled_invertible_percent = 100 * mean(test_results$rescaled_invertible, na.rm = TRUE),
    rescaled_valid_percent = 100 * mean(test_results$rescaled_valid, na.rm = TRUE),
    rescaled_mean_condition = mean(test_results$rescaled_condition[is.finite(test_results$rescaled_condition)], na.rm = TRUE),
    
    # Improvement statistics
    properties_preserved_percent = 100 * mean(test_results$properties_preserved, na.rm = TRUE),
    psd_improved_count = sum(test_results$psd_improved, na.rm = TRUE),
    invertibility_improved_count = sum(test_results$invertibility_improved, na.rm = TRUE),
    
    # Condition number analysis
    condition_improved_percent = 100 * mean(test_results$condition_ratio < 1, na.rm = TRUE),
    median_condition_ratio = median(test_results$condition_ratio[is.finite(test_results$condition_ratio)], na.rm = TRUE)
  )
  
  return(summary_stats)
}
```

### 6. Asymmetry Analysis (`1_bivariate_ordcats_correlation/2_asymmetry_analysis/asymmetry_measures.R`)

```r
# asymmetry_measures.R
# Comprehensive asymmetry analysis functions

#' Calculate multiple asymmetry measures for correlation bounds
#'
#' @param r_min Minimum correlation bound
#' @param r_max Maximum correlation bound
#' @return List with various asymmetry measures
calculate_asymmetry_measures <- function(r_min, r_max) {
  # Basic asymmetry (should be 0 for symmetric bounds)
  basic_asymmetry <- r_max + r_min
  
  # Range and center
  bounds_range <- r_max - r_min
  bounds_center <- (r_max + r_min) / 2
  
  # Normalized asymmetry measures
  if (bounds_range > 1e-10) {
    normalized_asymmetry <- basic_asymmetry / bounds_range
    relative_shift <- bounds_center / (bounds_range / 2)  # How far center is from 0
  } else {
    normalized_asymmetry <- 0
    relative_shift <- 0
  }
  
  # Asymmetry classification
  asymmetry_magnitude <- abs(basic_asymmetry)
  if (asymmetry_magnitude < 0.01) {
    asymmetry_class <- "Symmetric"
  } else if (asymmetry_magnitude < 0.1) {
    asymmetry_class <- "Mildly Asymmetric" 
  } else if (asymmetry_magnitude < 0.3) {
    asymmetry_class <- "Moderately Asymmetric"
  } else {
    asymmetry_class <- "Highly Asymmetric"
  }
  
  # Direction of asymmetry
  if (basic_asymmetry > 0.01) {
    asymmetry_direction <- "Positive Skewed"  # Upper bound dominates
  } else if (basic_asymmetry < -0.01) {
    asymmetry_direction <- "Negative Skewed"  # Lower bound dominates
  } else {
    asymmetry_direction <- "Balanced"
  }
  
  return(list(
    basic_asymmetry = basic_asymmetry,
    normalized_asymmetry = normalized_asymmetry,
    bounds_range = bounds_range,
    bounds_center = bounds_center,
    relative_shift = relative_shift,
    asymmetry_magnitude = asymmetry_magnitude,
    asymmetry_class = asymmetry_class,
    asymmetry_direction = asymmetry_direction
  ))
}

#' Calculate total variation distance between marginal distributions
#'
#' @param marginal1 First marginal distribution (counts or probabilities)
#' @param marginal2 Second marginal distribution (counts or probabilities)
#' @return Total variation distance
total_variation_distance <- function(marginal1, marginal2) {
  # Normalize to probabilities
  p1 <- marginal1 / sum(marginal1)
  p2 <- marginal2 / sum(marginal2)
  
  # Handle different lengths by padding with zeros
  max_len <- max(length(p1), length(p2))
  if (length(p1) < max_len) p1 <- c(p1, rep(0, max_len - length(p1)))
  if (length(p2) < max_len) p2 <- c(p2, rep(0, max_len - length(p2)))
  
  # Total variation distance = 0.5 * sum(|p1_i - p2_i|)
  tv_distance <- 0.5 * sum(abs(p1 - p2))
  return(tv_distance)
}

#' Analyze asymmetry patterns across all variable pairs
#'
#' @param results_df Data frame with correlation bounds results
#' @return Comprehensive asymmetry analysis
analyze_bounds_asymmetry_comprehensive <- function(results_df) {
  # Calculate asymmetry measures for each pair
  asymmetry_list <- apply(results_df, 1, function(row) {
    calculate_asymmetry_measures(as.numeric(row["r_min"]), as.numeric(row["r_max"]))
  })
  
  # Convert to data frame format
  asymmetry_df <- do.call(rbind, lapply(asymmetry_list, as.data.frame))
  
  # Combine with original results
  results_with_asymmetry <- cbind(results_df, asymmetry_df)
  
  # Calculate marginal distribution similarities if available
  if (all(c("var1_marginals", "var2_marginals") %in% names(results_df))) {
    results_with_asymmetry$marginal_tv_distance <- mapply(
      function(m1, m2) {
        if (is.list(m1) && is.list(m2)) {
          total_variation_distance(m1[[1]], m2[[1]])
        } else {
          NA
        }
      },
      results_df$var1_marginals,
      results_df$var2_marginals
    )
  }
  
  # Summary statistics
  asymmetry_summary <- list(
    # Basic statistics
    mean_asymmetry = mean(asymmetry_df$basic_asymmetry, na.rm = TRUE),
    median_asymmetry = median(asymmetry_df$basic_asymmetry, na.rm = TRUE),
    sd_asymmetry = sd(asymmetry_df$basic_asymmetry, na.rm = TRUE),
    
    # Proportions by class
    prop_symmetric = mean(asymmetry_df$asymmetry_class == "Symmetric", na.rm = TRUE),
    prop_mildly_asymmetric = mean(asymmetry_df$asymmetry_class == "Mildly Asymmetric", na.rm = TRUE),
    prop_moderately_asymmetric = mean(asymmetry_df$asymmetry_class == "Moderately Asymmetric", na.rm = TRUE),
    prop_highly_asymmetric = mean(asymmetry_df$asymmetry_class == "Highly Asymmetric", na.rm = TRUE),
    
    # Direction analysis
    prop_positive_skewed = mean(asymmetry_df$asymmetry_direction == "Positive Skewed", na.rm = TRUE),
    prop_negative_skewed = mean(asymmetry_df$asymmetry_direction == "Negative Skewed", na.rm = TRUE),
    prop_balanced = mean(asymmetry_df$asymmetry_direction == "Balanced", na.rm = TRUE),
    
    # Extreme cases
    max_asymmetry = max(abs(asymmetry_df$basic_asymmetry), na.rm = TRUE),
    most_asymmetric_pair = which.max(abs(asymmetry_df$basic_asymmetry)),
    
    # Range characteristics
    mean_bounds_range = mean(asymmetry_df$bounds_range, na.rm = TRUE),
    correlation_asymmetry_range = cor(abs(asymmetry_df$basic_asymmetry), 
                                    asymmetry_df$bounds_range, use = "complete.obs")
  )
  
  return(list(
    detailed_results = results_with_asymmetry,
    summary = asymmetry_summary,
    asymmetry_by_class = table(asymmetry_df$asymmetry_class),
    asymmetry_by_direction = table(asymmetry_df$asymmetry_direction)
  ))
}

#' Identify factors associated with asymmetry
#'
#' @param asymmetry_results Output from analyze_bounds_asymmetry_comprehensive
#' @return Analysis of asymmetry correlates
analyze_asymmetry_correlates <- function(asymmetry_results) {
  data <- asymmetry_results$detailed_results
  
  # Correlations with variable characteristics
  correlates <- list()
  
  if ("var1_categories" %in% names(data)) {
    correlates$with_var1_categories <- cor(abs(data$basic_asymmetry), 
                                         data$var1_categories, use = "complete.obs")
    correlates$with_var2_categories <- cor(abs(data$basic_asymmetry), 
                                         data$var2_categories, use = "complete.obs")
    correlates$with_total_categories <- cor(abs(data$basic_asymmetry), 
                                          data$var1_categories + data$var2_categories, 
                                          use = "complete.obs")
    correlates$with_category_difference <- cor(abs(data$basic_asymmetry), 
                                             abs(data$var1_categories - data$var2_categories), 
                                             use = "complete.obs")
  }
  
  if ("n_obs" %in% names(data)) {
    correlates$with_sample_size <- cor(abs(data$basic_asymmetry), 
                                     data$n_obs, use = "complete.obs")
  }
  
  if ("marginal_tv_distance" %in% names(data)) {
    correlates$with_marginal_distance <- cor(abs(data$basic_asymmetry), 
                                            data$marginal_tv_distance, use = "complete.obs")
  }
  
  # Other correlations
  correlates$with_bounds_range <- cor(abs(data$basic_asymmetry), 
                                    data$bounds_range, use = "complete.obs")
  correlates$with_observed_correlation <- cor(abs(data$basic_asymmetry), 
                                            abs(data$observed_r), use = "complete.obs")
  
  return(correlates)
}
```

### 7. Main Analysis Coordinator (`main_analysis.R`)

```r
# main_analysis.R
# Master analysis coordination script

# Load required libraries
library(dplyr)
library(ggplot2)
library(gridExtra)
library(Matrix)

# Source all module functions
source("1_bivariate_ordcats_correlation/bivariate_main.R")
source("2_correlation_matrices/matrices_main.R") 
source("3_fixes_and_rescaling/rescaling_main.R")
source("utilities/helper_functions.R")
source("utilities/data_loading.R")

#' Run complete ordinal correlation analysis pipeline
#'
#' @param bes_data BES dataset with ordinal variables
#' @param config Analysis configuration list
#' @param output_dir Directory for outputs
#' @return List with all analysis results
run_complete_analysis <- function(bes_data, config = get_default_config(), output_dir = "output/") {
  
  # Create output directories
  create_output_directories(output_dir)
  
  cat("Starting comprehensive ordinal correlation analysis...\n")
  cat("Configuration:\n")
  print(config)
  
  # 1. BIVARIATE CORRELATION BOUNDS ANALYSIS
  cat("\n1. BIVARIATE CORRELATION BOUNDS ANALYSIS\n")
  cat("==========================================\n")
  
  # Calculate theoretical bounds and permutation distributions
  bounds_results <- analyze_all_bes_bounds(
    bes_data, 
    nsim = config$nsim,
    progress = config$progress
  )
  
  # Asymmetry analysis
  asymmetry_results <- analyze_bounds_asymmetry_comprehensive(bounds_results)
  
  # Generate significance comparison
  significance_summary <- generate_significance_summary(bounds_results)
  
  cat("Bounds analysis complete for", nrow(bounds_results), "variable pairs\n")
  
  # 2. CORRELATION MATRIX PROPERTIES ANALYSIS  
  cat("\n2. CORRELATION MATRIX PROPERTIES ANALYSIS\n")
  cat("===========================================\n")
  
  # Random matrix trials
  matrix_trials <- test_random_correlation_matrices(
    bounds_results,
    n_trials = config$matrix_trials,
    min_vars = config$min_vars,
    max_vars = config$max_vars,
    enforce_psd = TRUE
  )
  
  matrix_trials_no_psd <- test_random_correlation_matrices(
    bounds_results,
    n_trials = config$matrix_trials,
    min_vars = config$min_vars, 
    max_vars = config$max_vars,
    enforce_psd = FALSE
  )
  
  # Summarize matrix properties
  matrix_summary_psd <- summarize_matrix_tests(matrix_trials)
  matrix_summary_no_psd <- summarize_matrix_tests(matrix_trials_no_psd)
  
  cat("Matrix analysis complete:", matrix_summary_psd$n_trials, "trials with PSD,", 
      matrix_summary_no_psd$n_trials, "trials without PSD\n")
  
  # 3. RESCALING ANALYSIS
  cat("\n3. RESCALING ANALYSIS\n")
  cat("=====================\n")
  
  # Apply rescaling
  rescaled_results <- apply_simple_rescaling(bounds_results)
  rescaling_analysis <- analyze_rescaling_effects(rescaled_results)
  rescaling_characteristics <- analyze_rescaling_by_characteristics(rescaled_results)
  
  cat("Rescaling analysis complete\n")
  
  # 4. GENERATE VISUALIZATIONS
  cat("\n4. GENERATING VISUALIZATIONS\n")
  cat("============================\n")
  
  # Create comprehensive visualization suite
  visualizations <- generate_all_visualizations(
    bounds_results = rescaled_results,
    asymmetry_results = asymmetry_results,
    matrix_results = list(psd = matrix_trials, no_psd = matrix_trials_no_psd),
    rescaling_results = rescaling_analysis,
    output_dir = file.path(output_dir, "figures")
  )
  
  cat("Visualizations saved to", file.path(output_dir, "figures"), "\n")
  
  # 5. GENERATE COMPREHENSIVE REPORT
  cat("\n5. GENERATING COMPREHENSIVE REPORT\n")
  cat("==================================\n")
  
  # Compile final report
  final_report <- compile_analysis_report(
    bounds_results = rescaled_results,
    asymmetry_results = asymmetry_results,
    significance_results = significance_summary,
    matrix_results_psd = matrix_summary_psd,
    matrix_results_no_psd = matrix_summary_no_psd,
    rescaling_results = rescaling_analysis,
    config = config
  )
  
  # Save report
  save_analysis_report(final_report, file.path(output_dir, "reports"))
  
  cat("Analysis complete! Results saved to", output_dir, "\n")
  
  # Return comprehensive results
  return(list(
    bounds = rescaled_results,
    asymmetry = asymmetry_results,
    significance = significance_summary,
    matrices = list(
      psd = matrix_summary_psd,
      no_psd = matrix_summary_no_psd,
      trials_psd = matrix_trials,
      trials_no_psd = matrix_trials_no_psd
    ),
    rescaling = rescaling_analysis,
    visualizations = visualizations,
    report = final_report,
    config = config
  ))
}

#' Get default analysis configuration
#'
#' @return List with default configuration parameters
get_default_config <- function() {
  list(
    nsim = 1000,                    # Permutation simulations per pair
    matrix_trials = 100,            # Random matrix trials
    min_vars = 3,                   # Minimum variables per matrix
    max_vars = 15,                  # Maximum variables per matrix
    progress = TRUE,                # Show progress indicators
    seed = 123,                     # Random seed
    significance_alpha = 0.05,      # Significance level
    tolerance = 1e-10,              # Numerical tolerance
    n_bootstrap = 500,              # Bootstrap samples for CI
    plot_theme = "minimal"          # ggplot theme
  )
}

#' Create output directory structure
#'
#' @param output_dir Base output directory
create_output_directories <- function(output_dir) {
  dirs <- c(
    file.path(output_dir, "figures"),
    file.path(output_dir, "figures", "bounds"),
    file.path(output_dir, "figures", "asymmetry"), 
    file.path(output_dir, "figures", "matrices"),
    file.path(output_dir, "figures", "rescaling"),
    file.path(output_dir, "tables"),
    file.path(output_dir, "reports"),
    file.path(output_dir, "data")
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
}

# Execute analysis if run directly
if (!interactive()) {
  # Configuration
  config <- get_default_config()
  config$nsim <- 500  # Reduce for faster execution
  config$matrix_trials <- 50
  
  # Load BES data (implement load_bes_data in utilities/data_loading.R)
  bes_data <- load_bes_data()
  
  # Run complete analysis
  results <- run_complete_analysis(bes_data, config)
  
  cat("\nAnalysis Summary:\n")
  cat("================\n")
  cat("Variable pairs analyzed:", nrow(results$bounds), "\n")
  cat("Asymmetric bounds:", round(100 * (1 - results$asymmetry$summary$prop_symmetric), 1), "%\n")
  cat("Matrix trials (PSD):", results$matrices$psd$n_trials, "\n")
  cat("Matrix trials (no PSD):", results$matrices$no_psd$n_trials, "\n")
  cat("Rescaling improved correlations:", round(100 * results$rescaling$summary$prop_increased, 1), "%\n")
}