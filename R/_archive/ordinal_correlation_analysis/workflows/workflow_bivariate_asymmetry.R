# ================================================================
# BIVARIATE ASYMMETRY ANALYSIS WORKFLOW
# ================================================================
# Purpose: Analyze asymmetry patterns in correlation bounds and their
#          implications for statistical inference
# ================================================================

library(here)
library(dplyr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(gridExtra)

# Load cache management utilities
source(here("R", "ordinal_correlation_analysis", "utilities", "cache_management.R"))

# ================================================================
# CONFIGURATION PARAMETERS
# ================================================================

params <- list(
  # Data sources
  use_mc_data = TRUE,              # Use Monte Carlo simulation data
  use_bes_data = TRUE,             # Use BES survey data
  force_regenerate = FALSE,        # Set TRUE to ignore existing cache
  
  # Analysis parameters
  asymmetry_threshold = 0.1,       # Threshold for significant asymmetry
  min_sample_size = 50,            # Minimum sample size for analysis
  confidence_level = 0.95,         # Confidence level for statistical tests
  
  # Visualization parameters
  max_points_plot = 200,           # Maximum points per plot for performance
  color_palette = "viridis",       # Color palette for plots
  
  # Progress reporting
  verbose = TRUE                   # Show detailed progress messages
)

cat("🔄 BIVARIATE ASYMMETRY ANALYSIS WORKFLOW\n")
cat("=========================================\n\n")

cat("📋 Configuration:\n")
cat("   Use MC data:", params$use_mc_data, "\n")
cat("   Use BES data:", params$use_bes_data, "\n")
cat("   Asymmetry threshold:", params$asymmetry_threshold, "\n")
cat("   Min sample size:", params$min_sample_size, "\n")
cat("   Force regenerate:", params$force_regenerate, "\n\n")

# ================================================================
# SECTION 1: LOAD AND PREPARE DATA
# ================================================================

cat("=== SECTION 1: DATA LOADING AND PREPARATION ===\n")

# Initialize data containers
asymmetry_data <- list()

# Load Monte Carlo data if requested
if (params$use_mc_data) {
  mc_files <- list.files(here("R", "ordinal_correlation_analysis", "output", "reports"), 
                        pattern = "mc_bounds_analysis.*\\.rds", full.names = TRUE)
  
  if (length(mc_files) > 0) {
    mc_file <- mc_files[which.max(file.mtime(mc_files))]
    mc_data <- readRDS(mc_file)
    
    if (params$verbose) {
      cat("📂 Loaded MC data:", basename(mc_file), "\n")
      cat("   Configurations:", length(mc_data), "\n")
    }
    
    # Convert MC data to asymmetry analysis format
    mc_asymmetry <- data.frame()
    for (config_id in seq_along(mc_data)) {
      config_results <- mc_data[[config_id]]
      if (!is.null(config_results$table_results)) {
        for (table_result in config_results$table_results) {
          mc_asymmetry <- rbind(mc_asymmetry, data.frame(
            source = "Monte Carlo",
            config_id = config_id,
            table_id = table_result$table_id,
            r_min = table_result$r_min,
            r_max = table_result$r_max,
            bounds_range = table_result$bounds_range,
            bounds_asymmetry = table_result$bounds_asymmetry,
            total_n = sum(table_result$table)
          ))
        }
      }
    }
    
    asymmetry_data$mc <- mc_asymmetry
    cat("✅ MC asymmetry data prepared:", nrow(mc_asymmetry), "observations\n")
  } else {
    cat("⚠️  No MC data found for asymmetry analysis\n")
  }
}

# Load BES data if requested
if (params$use_bes_data) {
  bes_files <- list.files(here("R", "ordinal_correlation_analysis", "output", "reports"), 
                         pattern = "bes_analysis.*\\.rds", full.names = TRUE)
  
  if (length(bes_files) > 0) {
    bes_file <- bes_files[which.max(file.mtime(bes_files))]
    bes_data <- readRDS(bes_file)
    
    if (params$verbose) {
      cat("📂 Loaded BES data:", basename(bes_file), "\n")
      cat("   Variable pairs:", nrow(bes_data), "\n")
    }
    
    # Convert BES data to asymmetry analysis format
    bes_asymmetry <- bes_data %>%
      filter(n_obs >= params$min_sample_size) %>%
      mutate(
        source = "BES Survey",
        config_id = row_number(),
        table_id = 1,
        total_n = n_obs
      ) %>%
      select(source, config_id, table_id, r_min, r_max, bounds_range, 
             bounds_asymmetry, total_n)
    
    asymmetry_data$bes <- bes_asymmetry
    cat("✅ BES asymmetry data prepared:", nrow(bes_asymmetry), "observations\n")
  } else {
    cat("⚠️  No BES data found for asymmetry analysis\n")
  }
}

# Combine all data sources
if (length(asymmetry_data) > 0) {
  combined_data <- do.call(rbind, asymmetry_data)
  cat("✅ Combined asymmetry dataset ready:", nrow(combined_data), "total observations\n\n")
} else {
  stop("❌ No data available for asymmetry analysis")
}

# ================================================================
# SECTION 2: ASYMMETRY PATTERN ANALYSIS
# ================================================================

cat("=== SECTION 2: ASYMMETRY PATTERN ANALYSIS ===\n")

# Cache file for asymmetry analysis
asymmetry_cache_file <- here("R", "ordinal_correlation_analysis", "output", "reports",
                           generate_cache_filename("asymmetry_analysis", 
                                                 params[c("asymmetry_threshold", "min_sample_size")]))

asymmetry_analysis <- cache_or_compute(
  cache_file = asymmetry_cache_file,
  force_regenerate = params$force_regenerate,
  compute_func = function() {
    cat("Computing asymmetry pattern analysis...\n")
    
    # Asymmetry classification
    analysis_results <- combined_data %>%
      mutate(
        # Classify asymmetry strength
        asymmetry_category = case_when(
          abs(bounds_asymmetry) < params$asymmetry_threshold ~ "Symmetric",
          bounds_asymmetry >= params$asymmetry_threshold ~ "Positive Asymmetry",
          bounds_asymmetry <= -params$asymmetry_threshold ~ "Negative Asymmetry",
          TRUE ~ "Mild Asymmetry"
        ),
        
        # Relative position of zero within bounds
        zero_position = ifelse(r_max != r_min, -r_min / (r_max - r_min), 0.5),
        
        # Asymmetry strength (absolute)
        asymmetry_strength = abs(bounds_asymmetry),
        
        # Range-adjusted asymmetry
        relative_asymmetry = ifelse(bounds_range > 0, bounds_asymmetry / bounds_range, 0)
      )
    
    # Summary statistics by source
    summary_by_source <- analysis_results %>%
      group_by(source) %>%
      summarise(
        n_observations = n(),
        mean_asymmetry = mean(bounds_asymmetry, na.rm = TRUE),
        median_asymmetry = median(bounds_asymmetry, na.rm = TRUE),
        sd_asymmetry = sd(bounds_asymmetry, na.rm = TRUE),
        prop_symmetric = mean(asymmetry_category == "Symmetric", na.rm = TRUE),
        prop_positive_asym = mean(asymmetry_category == "Positive Asymmetry", na.rm = TRUE),
        prop_negative_asym = mean(asymmetry_category == "Negative Asymmetry", na.rm = TRUE),
        mean_range = mean(bounds_range, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Summary statistics by asymmetry category
    summary_by_category <- analysis_results %>%
      group_by(asymmetry_category) %>%
      summarise(
        n_observations = n(),
        proportion = n() / nrow(analysis_results),
        mean_range = mean(bounds_range, na.rm = TRUE),
        median_range = median(bounds_range, na.rm = TRUE),
        mean_r_min = mean(r_min, na.rm = TRUE),
        mean_r_max = mean(r_max, na.rm = TRUE),
        .groups = "drop"
      )
    
    list(
      detailed_results = analysis_results,
      summary_by_source = summary_by_source,
      summary_by_category = summary_by_category,
      analysis_info = list(
        asymmetry_threshold = params$asymmetry_threshold,
        total_observations = nrow(analysis_results),
        computation_time = Sys.time()
      )
    )
  }
)

cat("✅ Asymmetry pattern analysis complete\n")
cat("   Total observations:", nrow(asymmetry_analysis$detailed_results), "\n")
cat("   Asymmetry categories:", length(unique(asymmetry_analysis$detailed_results$asymmetry_category)), "\n\n")

# ================================================================
# SECTION 3: ASYMMETRY VISUALIZATIONS
# ================================================================

cat("=== SECTION 3: ASYMMETRY VISUALIZATIONS ===\n")

# Create figures directory
figures_dir <- here("R", "ordinal_correlation_analysis", "output", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Sample data for visualization if too large
plot_data <- asymmetry_analysis$detailed_results
if (nrow(plot_data) > params$max_points_plot) {
  plot_data <- plot_data %>% 
    sample_n(params$max_points_plot)
  cat("📊 Sampled", params$max_points_plot, "points for visualization\n")
}

# Plot 1: Asymmetry distribution by source
plot1 <- ggplot(asymmetry_analysis$detailed_results, aes(x = bounds_asymmetry, fill = source)) +
  geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
  geom_vline(xintercept = c(-params$asymmetry_threshold, params$asymmetry_threshold), 
             linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.8) +
  facet_wrap(~source, scales = "free_y") +
  scale_fill_viridis_d(name = "Data Source") +
  labs(
    title = "Distribution of Bounds Asymmetry by Data Source",
    subtitle = paste("Red lines indicate asymmetry threshold (±", params$asymmetry_threshold, ")"),
    x = "Bounds Asymmetry",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

plot1

plot1_file <- file.path(figures_dir, "asymmetry_distribution_by_source.png")
ggsave(plot1_file, plot1, width = 12, height = 6, dpi = 300)

# Plot 2: Asymmetry vs Range relationship
plot2 <- ggplot(plot_data, aes(x = bounds_range, y = bounds_asymmetry, color = source)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = c(-params$asymmetry_threshold, params$asymmetry_threshold), 
             linetype = "dashed", color = "red", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
  scale_color_viridis_d(name = "Data Source") +
  labs(
    title = "Relationship Between Bounds Range and Asymmetry",
    subtitle = "LOESS smoothing shows trend patterns",
    x = "Bounds Range (r_max - r_min)",
    y = "Bounds Asymmetry"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

plot2

plot2_file <- file.path(figures_dir, "asymmetry_vs_range_relationship.png")
ggsave(plot2_file, plot2, width = 10, height = 8, dpi = 300)

# Plot 3: Bounds landscape colored by asymmetry
plot3 <- ggplot(plot_data, aes(x = r_min, y = r_max, color = bounds_asymmetry)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "gray", alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Asymmetry") +
  facet_wrap(~source) +
  labs(
    title = "Correlation Bounds Landscape Colored by Asymmetry",
    subtitle = "Blue = negative asymmetry, Red = positive asymmetry",
    x = "r_min (Minimum Correlation)",
    y = "r_max (Maximum Correlation)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

plot3

plot3_file <- file.path(figures_dir, "bounds_landscape_asymmetry.png")
ggsave(plot3_file, plot3, width = 12, height = 6, dpi = 300)

# Plot 4: Asymmetry category proportions
category_summary <- asymmetry_analysis$summary_by_category %>%
  mutate(asymmetry_category = factor(asymmetry_category, 
                                   levels = c("Negative Asymmetry", "Symmetric", 
                                            "Mild Asymmetry", "Positive Asymmetry")))

plot4 <- ggplot(category_summary, aes(x = asymmetry_category, y = proportion, 
                                     fill = asymmetry_category)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = paste0(round(proportion * 100, 1), "%")), 
            vjust = -0.5, fontface = "bold") +
  scale_fill_viridis_d(guide = "none") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Distribution of Asymmetry Categories",
    subtitle = paste("Based on threshold of ±", params$asymmetry_threshold),
    x = "Asymmetry Category",
    y = "Proportion of Observations"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

plot4

plot4_file <- file.path(figures_dir, "asymmetry_category_proportions.png")
ggsave(plot4_file, plot4, width = 10, height = 6, dpi = 300)

cat("📈 Asymmetry visualizations generated:\n")
cat("   1. Distribution by source:", basename(plot1_file), "\n")
cat("   2. Range vs asymmetry:", basename(plot2_file), "\n")
cat("   3. Bounds landscape:", basename(plot3_file), "\n")
cat("   4. Category proportions:", basename(plot4_file), "\n\n")

# ================================================================
# SECTION 4: ASYMMETRY SUMMARY TABLES
# ================================================================

cat("=== SECTION 4: ASYMMETRY SUMMARY TABLES ===\n")

# Create tables directory
tables_dir <- here("R", "ordinal_correlation_analysis", "output", "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

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

# Save summary by source
source_files <- save_enhanced_table(
  asymmetry_analysis$summary_by_source,
  "asymmetry_summary_by_source",
  "Asymmetry Analysis Summary by Data Source"
)

# Save summary by category
category_files <- save_enhanced_table(
  asymmetry_analysis$summary_by_category,
  "asymmetry_summary_by_category", 
  "Asymmetry Analysis Summary by Category"
)

# Create detailed statistical summary
detailed_stats <- data.frame(
  Statistic = c(
    "Total Observations",
    "Mean Asymmetry",
    "Median Asymmetry", 
    "Standard Deviation Asymmetry",
    "Proportion Symmetric",
    "Proportion Positive Asymmetry",
    "Proportion Negative Asymmetry",
    "Mean Bounds Range",
    "Asymmetry Threshold Used"
  ),
  Value = c(
    nrow(asymmetry_analysis$detailed_results),
    round(mean(asymmetry_analysis$detailed_results$bounds_asymmetry, na.rm = TRUE), 4),
    round(median(asymmetry_analysis$detailed_results$bounds_asymmetry, na.rm = TRUE), 4),
    round(sd(asymmetry_analysis$detailed_results$bounds_asymmetry, na.rm = TRUE), 4),
    round(mean(asymmetry_analysis$detailed_results$asymmetry_category == "Symmetric", na.rm = TRUE), 4),
    round(mean(asymmetry_analysis$detailed_results$asymmetry_category == "Positive Asymmetry", na.rm = TRUE), 4),
    round(mean(asymmetry_analysis$detailed_results$asymmetry_category == "Negative Asymmetry", na.rm = TRUE), 4),
    round(mean(asymmetry_analysis$detailed_results$bounds_range, na.rm = TRUE), 4),
    params$asymmetry_threshold
  )
)

# Save detailed statistics
stats_files <- save_enhanced_table(
  detailed_stats,
  "asymmetry_detailed_statistics",
  "Comprehensive Asymmetry Analysis Statistics"
)

cat("📊 Asymmetry summary tables generated:\n")
cat("   1. Summary by source:\n")
cat("      • CSV:", source_files$csv, "\n")
cat("      • Markdown:", source_files$md, "\n")
cat("      • HTML:", source_files$html, "\n")
cat("   2. Summary by category:\n")
cat("      • CSV:", category_files$csv, "\n")
cat("      • Markdown:", category_files$md, "\n")
cat("      • HTML:", category_files$html, "\n")
cat("   3. Detailed statistics:\n")
cat("      • CSV:", stats_files$csv, "\n")
cat("      • Markdown:", stats_files$md, "\n")
cat("      • HTML:", stats_files$html, "\n\n")

# ================================================================
# WORKFLOW COMPLETION SUMMARY
# ================================================================

cat("🎉 BIVARIATE ASYMMETRY ANALYSIS WORKFLOW COMPLETE!\n")
cat("==================================================\n")
cat("📁 Generated files:\n")
cat("   Asymmetry analysis:", basename(asymmetry_cache_file), "\n")
cat("   Figures: 4 asymmetry visualization files\n")
cat("   Tables: 9 files (3 tables × 3 formats: CSV, Markdown, HTML)\n\n")

cat("📊 Asymmetry analysis summary:\n")
cat("   Total observations:", nrow(asymmetry_analysis$detailed_results), "\n")
cat("   Data sources:", length(unique(asymmetry_analysis$detailed_results$source)), "\n")
cat("   Mean asymmetry:", round(mean(asymmetry_analysis$detailed_results$bounds_asymmetry, na.rm = TRUE), 4), "\n")
cat("   Proportion symmetric:", round(mean(asymmetry_analysis$detailed_results$asymmetry_category == "Symmetric", na.rm = TRUE), 3), "\n\n")

cat("🔄 Key findings:\n")
symmetric_prop <- mean(asymmetry_analysis$detailed_results$asymmetry_category == "Symmetric", na.rm = TRUE)
if (symmetric_prop > 0.5) {
  cat("   - Majority of bounds show symmetric patterns\n")
} else {
  cat("   - Asymmetric bounds are more common than symmetric ones\n")
}

pos_asym_prop <- mean(asymmetry_analysis$detailed_results$asymmetry_category == "Positive Asymmetry", na.rm = TRUE)
neg_asym_prop <- mean(asymmetry_analysis$detailed_results$asymmetry_category == "Negative Asymmetry", na.rm = TRUE)

if (pos_asym_prop > neg_asym_prop) {
  cat("   - Positive asymmetry is more prevalent than negative\n")
} else if (neg_asym_prop > pos_asym_prop) {
  cat("   - Negative asymmetry is more prevalent than positive\n")
} else {
  cat("   - Positive and negative asymmetries are balanced\n")
}

cat("\n🔄 Integration with other workflows:\n")
cat("   - Results inform interpretation of bounds reliability\n")
cat("   - Asymmetry patterns guide rescaling method selection\n")
cat("   - Findings relevant for matrix construction strategies\n")

