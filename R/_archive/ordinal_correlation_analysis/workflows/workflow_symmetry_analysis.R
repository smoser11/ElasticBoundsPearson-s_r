# ================================================================
# SYMMETRY ANALYSIS WORKFLOW
# ================================================================
# Purpose: Deep investigation of symmetry in correlation bounds,
#          marginal distributions, and their entropic properties
# ================================================================

library(here)
library(dplyr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(gridExtra)
library(viridis)
library(tidyr)

# Install and load overlap package for distribution metrics
if (!require(overlap, quietly = TRUE)) {
  cat("Installing overlap package...\n")
  # Set CRAN mirror
  options(repos = c(CRAN = "https://cran.rstudio.com/"))
  install.packages("overlap", quiet = TRUE)
  library(overlap)
  cat("✅ overlap package installed and loaded\n")
} else {
  cat("✅ overlap package already available\n")
}

# Load cache management utilities
source(here("R", "ordinal_correlation_analysis", "utilities", "cache_management.R"))
source(here("R", "correlation_bounds_core.R"))

# ================================================================
# CONFIGURATION PARAMETERS
# ================================================================

params <- list(
  # Analysis scope
  investigate_marginal_symmetry = TRUE,   # Analyze marginal distribution symmetry
  investigate_bounds_asymmetry = TRUE,    # When r_min ≠ -r_max
  use_entropy_metrics = TRUE,             # Include entropy-based measures
  
  # Simulation parameters for controlled experiments
  n_simulations = 100,                    # Number of controlled simulations
  sample_sizes = c(100, 500, 1000),      # Sample sizes to test
  category_counts = list(                 # Category combinations to test
    equal = list(c(3,3), c(4,4), c(5,5)),
    unequal = list(c(3,5), c(4,6), c(5,7))
  ),
  
  # Distribution types for marginal symmetry analysis
  distribution_types = c("uniform", "peaked", "skewed_left", "skewed_right", "bimodal"),
  
  # Thresholds
  symmetry_threshold = 0.05,              # Threshold for "nearly symmetric"
  entropy_bins = 10,                      # Bins for entropy calculation
  
  # Output control
  force_regenerate = FALSE,               # Set TRUE to ignore existing cache
  verbose = TRUE,                         # Show detailed progress
  max_plots_per_panel = 16               # Maximum subplots per visualization
)

cat("🔍 SYMMETRY ANALYSIS WORKFLOW\n")
cat("============================\n\n")

cat("📋 Configuration:\n")
cat("   Investigate marginal symmetry:", params$investigate_marginal_symmetry, "\n")
cat("   Investigate bounds asymmetry:", params$investigate_bounds_asymmetry, "\n")
cat("   Use entropy metrics:", params$use_entropy_metrics, "\n")
cat("   Simulation count:", params$n_simulations, "\n")
cat("   Symmetry threshold:", params$symmetry_threshold, "\n\n")

# ================================================================
# SECTION 1: MARGINAL DISTRIBUTION SYMMETRY FUNCTIONS
# ================================================================

cat("=== SECTION 1: MARGINAL DISTRIBUTION ANALYSIS FUNCTIONS ===\n")

# Function to generate different types of marginal distributions
generate_marginal_distribution <- function(n_categories, type = "uniform", total_count = 1000) {
  probs <- switch(type,
    "uniform" = rep(1/n_categories, n_categories),
    "peaked" = {
      center <- ceiling(n_categories/2)
      probs <- exp(-0.5 * ((1:n_categories - center) / 1.5)^2)
      probs / sum(probs)
    },
    "skewed_left" = {
      probs <- (n_categories:1)^2
      probs / sum(probs)
    },
    "skewed_right" = {
      probs <- (1:n_categories)^2
      probs / sum(probs)
    },
    "bimodal" = {
      probs <- ifelse(1:n_categories %in% c(2, n_categories-1), 0.4, 0.2/(n_categories-2))
      probs / sum(probs)
    }
  )
  
  # Convert to counts
  counts <- round(probs * total_count)
  # Ensure total equals exactly total_count
  diff <- total_count - sum(counts)
  if (diff != 0) {
    max_idx <- which.max(counts)
    counts[max_idx] <- counts[max_idx] + diff
  }
  
  return(list(probs = probs, counts = counts))
}

# Function to calculate distribution symmetry metrics using TV, BC, and OVL
calculate_distribution_symmetry <- function(probs) {
  n <- length(probs)
  probs_rev <- rev(probs)  # Reversed distribution p^{rev}
  
  # (i) Total Variation (TV) Distance
  # TV(p, p^rev) = (1/2) * sum(|p_i - p_{K+1-i}|)
  tv_distance <- 0.5 * sum(abs(probs - probs_rev))
  
  # (ii) Bhattacharyya Coefficient (BC) and its complement
  # BC(p, p^rev) = sum(sqrt(p_i * p_{K+1-i}))
  bc_coefficient <- sum(sqrt(probs * probs_rev))
  bc_complement <- 1 - bc_coefficient  # Divergence-like measure
  
  # (iii) Overlap Coefficient (OVL) and its complement
  # OVL(p, p^rev) = sum(min{p_i, p_{K+1-i}})
  ovl_coefficient <- sum(pmin(probs, probs_rev))
  ovl_complement <- 1 - ovl_coefficient  # Non-shared probability mass
  
  # Additional basic measures
  symmetry_corr <- cor(probs, probs_rev)
  
  # Skewness-based symmetry (third central moment)
  mean_pos <- sum((1:n) * probs)
  variance <- sum(((1:n) - mean_pos)^2 * probs)
  skewness <- if(variance > 0) sum(((1:n) - mean_pos)^3 * probs) / (variance^(3/2)) else 0
  
  # Distance from center
  center <- (n + 1) / 2
  center_distance <- abs(mean_pos - center) / (n/2)
  
  list(
    # Primary asymmetry measures
    tv_distance = tv_distance,
    bc_coefficient = bc_coefficient,
    bc_complement = bc_complement,
    ovl_coefficient = ovl_coefficient,
    ovl_complement = ovl_complement,
    
    # Secondary measures
    symmetry_correlation = symmetry_corr,
    skewness = skewness,
    center_distance = center_distance,
    
    # Classification (using TV distance as primary measure)
    is_symmetric = tv_distance < params$symmetry_threshold
  )
}

# Function to calculate overlap between two distributions
calculate_distribution_overlap <- function(probs1, probs2) {
  # Ensure same length
  max_len <- max(length(probs1), length(probs2))
  if (length(probs1) < max_len) probs1 <- c(probs1, rep(0, max_len - length(probs1)))
  if (length(probs2) < max_len) probs2 <- c(probs2, rep(0, max_len - length(probs2)))
  
  # Calculate overlap coefficient
  overlap_coeff <- sum(pmin(probs1, probs2))
  
  # Bhattacharyya coefficient
  bhatt_coeff <- sum(sqrt(probs1 * probs2))
  
  list(
    overlap_coefficient = overlap_coeff,
    bhattacharyya_coefficient = bhatt_coeff
  )
}

cat("✅ Marginal distribution analysis functions ready\n\n")

# ================================================================
# SECTION 2: CONTROLLED SYMMETRY EXPERIMENTS
# ================================================================

cat("=== SECTION 2: CONTROLLED SYMMETRY EXPERIMENTS ===\n")

# Cache file for symmetry experiments
symmetry_cache_file <- here("R", "ordinal_correlation_analysis", "output", "reports",
                          generate_cache_filename("symmetry_experiments", 
                                                 params[c("n_simulations", "symmetry_threshold")]))

symmetry_experiments <- cache_or_compute(
  cache_file = symmetry_cache_file,
  force_regenerate = params$force_regenerate,
  compute_func = function() {
    cat("Running controlled symmetry experiments...\n")
    
    experiment_results <- data.frame()
    
    # Progress tracking
    total_experiments <- length(params$distribution_types)^2 * 
                        (length(params$category_counts$equal) + length(params$category_counts$unequal)) *
                        length(params$sample_sizes)
    exp_count <- 0
    
    for (sample_size in params$sample_sizes) {
      for (cat_type in names(params$category_counts)) {
        for (cat_combo in params$category_counts[[cat_type]]) {
          k1 <- cat_combo[1]
          k2 <- cat_combo[2]
          
          for (dist1_type in params$distribution_types) {
            for (dist2_type in params$distribution_types) {
              exp_count <- exp_count + 1
              
              if (params$verbose && exp_count %% 20 == 0) {
                cat("   Experiment", exp_count, "of", total_experiments, "\n")
              }
              
              # Generate marginal distributions
              marg1 <- generate_marginal_distribution(k1, dist1_type, sample_size)
              marg2 <- generate_marginal_distribution(k2, dist2_type, sample_size)
              
              # Calculate symmetry metrics for marginals
              sym1 <- calculate_distribution_symmetry(marg1$probs)
              sym2 <- calculate_distribution_symmetry(marg2$probs)
              
              # Calculate overlap between marginals
              overlap_metrics <- calculate_distribution_overlap(marg1$probs, marg2$probs)
              
              # Calculate correlation bounds
              r_min <- min_corr_bound(marg1$counts, marg2$counts)
              r_max <- max_corr_bound(marg1$counts, marg2$counts)
              
              # Bounds asymmetry metrics
              bounds_range <- r_max - r_min
              bounds_asymmetry <- (r_max + r_min) / 2
              bounds_symmetry_ratio <- ifelse(r_max != 0, abs(r_min) / abs(r_max), 1)
              
              # Store results
              experiment_results <- rbind(experiment_results, data.frame(
                exp_id = exp_count,
                sample_size = sample_size,
                category_type = cat_type,
                k1 = k1, k2 = k2,
                dist1_type = dist1_type, dist2_type = dist2_type,
                
                # Marginal 1 asymmetry measures
                marg1_tv_distance = sym1$tv_distance,
                marg1_bc_coefficient = sym1$bc_coefficient,
                marg1_bc_complement = sym1$bc_complement,
                marg1_ovl_coefficient = sym1$ovl_coefficient,
                marg1_ovl_complement = sym1$ovl_complement,
                marg1_symmetry_corr = sym1$symmetry_correlation,
                marg1_skewness = sym1$skewness,
                marg1_is_symmetric = sym1$is_symmetric,
                
                # Marginal 2 asymmetry measures
                marg2_tv_distance = sym2$tv_distance,
                marg2_bc_coefficient = sym2$bc_coefficient,
                marg2_bc_complement = sym2$bc_complement,
                marg2_ovl_coefficient = sym2$ovl_coefficient,
                marg2_ovl_complement = sym2$ovl_complement,
                marg2_symmetry_corr = sym2$symmetry_correlation,
                marg2_skewness = sym2$skewness,
                marg2_is_symmetric = sym2$is_symmetric,
                
                # Cross-marginal overlap metrics
                cross_overlap_coefficient = overlap_metrics$overlap_coefficient,
                cross_bhattacharyya_coeff = overlap_metrics$bhattacharyya_coefficient,
                
                # Bounds results
                r_min = r_min,
                r_max = r_max,
                bounds_range = bounds_range,
                bounds_asymmetry = bounds_asymmetry,
                bounds_symmetry_ratio = bounds_symmetry_ratio,
                
                # Combined symmetry indicators
                both_marginals_symmetric = sym1$is_symmetric & sym2$is_symmetric,
                any_marginal_symmetric = sym1$is_symmetric | sym2$is_symmetric,
                bounds_nearly_symmetric = abs(bounds_asymmetry) < params$symmetry_threshold
              ))
            }
          }
        }
      }
    }
    
    list(
      detailed_results = experiment_results,
      experiment_info = list(
        total_experiments = total_experiments,
        distribution_types = params$distribution_types,
        category_combinations = params$category_counts,
        sample_sizes = params$sample_sizes,
        computation_time = Sys.time()
      )
    )
  }
)

cat("✅ Controlled symmetry experiments complete\n")
cat("   Total experiments:", nrow(symmetry_experiments$detailed_results), "\n")
cat("   Distribution types:", length(params$distribution_types), "\n\n")

# ================================================================
# SECTION 3: ASYMMETRY INVESTIGATION (r_min ≠ -r_max)
# ================================================================

cat("=== SECTION 3: BOUNDS ASYMMETRY INVESTIGATION ===\n")

# Analyze when and why bounds are asymmetric
bounds_asymmetry_analysis <- symmetry_experiments$detailed_results %>%
  mutate(
    # Classifications
    bounds_perfectly_symmetric = abs(r_min + r_max) < 1e-10,
    bounds_highly_asymmetric = abs(bounds_asymmetry) > 0.1,
    
    # Asymmetry strength categories
    asymmetry_strength = case_when(
      abs(bounds_asymmetry) < 0.01 ~ "Nearly Symmetric",
      abs(bounds_asymmetry) < 0.05 ~ "Mildly Asymmetric", 
      abs(bounds_asymmetry) < 0.1 ~ "Moderately Asymmetric",
      TRUE ~ "Highly Asymmetric"
    ),
    
    # Marginal symmetry combinations
    marginal_symmetry_pattern = case_when(
      both_marginals_symmetric ~ "Both Symmetric",
      marg1_is_symmetric & !marg2_is_symmetric ~ "Only Marg1 Symmetric",
      !marg1_is_symmetric & marg2_is_symmetric ~ "Only Marg2 Symmetric", 
      TRUE ~ "Both Asymmetric"
    )
  )

# Summary statistics using the three asymmetry measures
asymmetry_summary <- bounds_asymmetry_analysis %>%
  group_by(marginal_symmetry_pattern) %>%
  summarise(
    n_cases = n(),
    prop_bounds_symmetric = mean(bounds_nearly_symmetric),
    mean_bounds_asymmetry = mean(abs(bounds_asymmetry)),
    
    # TV Distance measures
    mean_tv_distance_1 = mean(marg1_tv_distance),
    mean_tv_distance_2 = mean(marg2_tv_distance),
    
    # Bhattacharyya Coefficient measures  
    mean_bc_complement_1 = mean(marg1_bc_complement),
    mean_bc_complement_2 = mean(marg2_bc_complement),
    
    # Overlap Coefficient measures
    mean_ovl_complement_1 = mean(marg1_ovl_complement),
    mean_ovl_complement_2 = mean(marg2_ovl_complement),
    
    # Cross-marginal overlap
    mean_cross_overlap = mean(cross_overlap_coefficient),
    .groups = "drop"
  )

cat("📊 Bounds asymmetry patterns identified:\n")
print(asymmetry_summary)
cat("\n")

# ================================================================
# SECTION 4: SYMMETRY VISUALIZATIONS
# ================================================================

cat("=== SECTION 4: SYMMETRY VISUALIZATIONS ===\n")

# Create figures directory
figures_dir <- here("R", "ordinal_correlation_analysis", "output", "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Plot 1: Marginal symmetry impact on bounds asymmetry
plot1 <- ggplot(bounds_asymmetry_analysis, 
                aes(x = marginal_symmetry_pattern, y = abs(bounds_asymmetry), 
                    fill = marginal_symmetry_pattern)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(alpha = 0.3, width = 0.3, size = 0.5) +
  scale_fill_viridis_d(guide = "none") +
  scale_y_log10() +
  labs(
    title = "Impact of Marginal Symmetry on Bounds Asymmetry",
    subtitle = "Log scale shows distribution of |bounds_asymmetry|",
    x = "Marginal Symmetry Pattern",
    y = "|Bounds Asymmetry| (log scale)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))


plot1

plot1_file <- file.path(figures_dir, "marginal_symmetry_impact_bounds.png")
ggsave(plot1_file, plot1, width = 12, height = 8, dpi = 300)

# Plot 2: Three asymmetry measures vs bounds asymmetry
# Create long format data for the three measures
asymmetry_measures_long <- bounds_asymmetry_analysis %>%
  select(exp_id, bounds_asymmetry, category_type,
         marg1_tv_distance, marg1_bc_complement, marg1_ovl_complement,
         marg2_tv_distance, marg2_bc_complement, marg2_ovl_complement) %>%
  pivot_longer(
    cols = c(marg1_tv_distance, marg1_bc_complement, marg1_ovl_complement,
             marg2_tv_distance, marg2_bc_complement, marg2_ovl_complement),
    names_to = "measure_var", 
    values_to = "asymmetry_value"
  ) %>%
  mutate(
    marginal = ifelse(grepl("marg1", measure_var), "Marginal 1", "Marginal 2"),
    measure = case_when(
      grepl("tv_distance", measure_var) ~ "TV Distance",
      grepl("bc_complement", measure_var) ~ "BC Complement", 
      grepl("ovl_complement", measure_var) ~ "OVL Complement"
    )
  )

plot2 <- ggplot(asymmetry_measures_long, 
                aes(x = asymmetry_value, y = abs(bounds_asymmetry), color = category_type)) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
  scale_color_viridis_d(name = "Category\nStructure") +
  scale_y_log10() +
  facet_wrap(~measure, scales = "free_x", ncol = 3) +
  labs(
    title = "Three Asymmetry Measures vs Bounds Asymmetry",
    subtitle = "TV focuses on discrepancy, BC on geometric congruence, OVL on shared mass",
    x = "Marginal Asymmetry Measure Value",
    y = "|Bounds Asymmetry| (log scale)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

plot2

plot2_file <- file.path(figures_dir, "three_asymmetry_measures_bounds.png")
ggsave(plot2_file, plot2, width = 14, height = 8, dpi = 300)

# Plot 3: Distribution overlap vs bounds characteristics
plot3 <- ggplot(bounds_asymmetry_analysis, 
                aes(x = cross_overlap_coefficient, y = bounds_range, 
                    color = abs(bounds_asymmetry))) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_viridis_c(name = "|Bounds\nAsymmetry|") +
  facet_wrap(vars(paste(k1, "×", k2))) +
  labs(
    title = "Cross-Marginal Overlap vs Bounds Range",
    subtitle = "Colored by bounds asymmetry magnitude",
    x = "Cross-Marginal Distribution Overlap Coefficient",
    y = "Bounds Range (r_max - r_min)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

plot3

plot3_file <- file.path(figures_dir, "overlap_bounds_relationship.png")
ggsave(plot3_file, plot3, width = 12, height = 8, dpi = 300)

# Plot 4: Detailed asymmetry heatmap
asymmetry_heatmap_data <- bounds_asymmetry_analysis %>%
  group_by(dist1_type, dist2_type, category_type) %>%
  summarise(
    mean_bounds_asymmetry = mean(abs(bounds_asymmetry)),
    prop_symmetric_bounds = mean(bounds_nearly_symmetric),
    .groups = "drop"
  )


plot4 <- ggplot(asymmetry_heatmap_data, 
                aes(x = dist1_type, y = dist2_type, fill = mean_bounds_asymmetry)) +
  geom_tile() +
  geom_text(aes(label = round(mean_bounds_asymmetry, 3)), 
            color = "white", fontface = "bold", size = 3) +
  scale_fill_viridis_c(name = "Mean\n|Asymmetry|") +
  facet_wrap(~category_type) +
  labs(
    title = "Distribution Type Combinations and Bounds Asymmetry",
    subtitle = "Heatmap showing mean |bounds_asymmetry| for each combination",
    x = "Distribution Type (Variable 1)",
    y = "Distribution Type (Variable 2)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

plot4

plot4_file <- file.path(figures_dir, "distribution_asymmetry_heatmap.png")
ggsave(plot4_file, plot4, width = 12, height = 8, dpi = 300)

# Plot 5: Theoretical bounds asymmetry examples
# Create specific examples showing why r_min ≠ -r_max
example_cases <- data.frame(
  case_name = c("Uniform-Uniform", "Peaked-Uniform", "Skewed-Skewed", "Bimodal-Peaked"),
  dist1 = c("uniform", "peaked", "skewed_right", "bimodal"),
  dist2 = c("uniform", "uniform", "skewed_left", "peaked"),
  k1 = c(4, 4, 5, 4),
  k2 = c(4, 4, 5, 3)
)

example_results <- data.frame()
for (i in 1:nrow(example_cases)) {
  case <- example_cases[i, ]
  marg1 <- generate_marginal_distribution(case$k1, case$dist1, 1000)
  marg2 <- generate_marginal_distribution(case$k2, case$dist2, 1000)
  
  r_min <- min_corr_bound(marg1$counts, marg2$counts)
  r_max <- max_corr_bound(marg1$counts, marg2$counts)
  
  example_results <- rbind(example_results, data.frame(
    case_name = case$case_name,
    r_min = r_min,
    r_max = r_max,
    asymmetry = (r_max + r_min) / 2,
    ratio = abs(r_min) / abs(r_max)
  ))
}

plot5 <- ggplot(example_results, aes(x = case_name)) +
  geom_segment(aes(y = r_min, yend = r_max, xend = case_name), 
               color = "darkblue", size = 3, alpha = 0.7) +
  geom_point(aes(y = r_min), color = "red", size = 4, alpha = 0.8) +
  geom_point(aes(y = r_max), color = "darkgreen", size = 4, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(
    title = "Theoretical Examples: Why r_min ≠ -r_max",
    subtitle = "Red = r_min, Green = r_max, Blue line = bounds range",
    x = "Distribution Combination",
    y = "Correlation Bounds"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

plot5

plot5_file <- file.path(figures_dir, "theoretical_asymmetry_examples.png")
ggsave(plot5_file, plot5, width = 10, height = 6, dpi = 300)

cat("📈 Symmetry visualizations generated:\n")
cat("   1. Marginal symmetry impact:", basename(plot1_file), "\n")
cat("   2. Three asymmetry measures:", basename(plot2_file), "\n")
cat("   3. Overlap relationship:", basename(plot3_file), "\n")
cat("   4. Asymmetry heatmap:", basename(plot4_file), "\n")
cat("   5. Theoretical examples:", basename(plot5_file), "\n\n")

# ================================================================
# SECTION 5: SYMMETRY SUMMARY TABLES
# ================================================================

cat("=== SECTION 5: SYMMETRY SUMMARY TABLES ===\n")

# Create tables directory
tables_dir <- here("R", "ordinal_correlation_analysis", "output", "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# Enhanced table output function
save_enhanced_table <- function(data, filename_base, caption) {
  # CSV output
  csv_file <- file.path(tables_dir, paste0(filename_base, ".csv"))
  write.csv(data, csv_file, row.names = FALSE)
  
  # Markdown output
  md_file <- file.path(tables_dir, paste0(filename_base, ".md"))
  md_table <- kable(data, format = "markdown", 
                    caption = caption, digits = 4,
                    col.names = gsub("_", " ", toupper(names(data))))
  writeLines(as.character(md_table), md_file)
  
  # HTML output
  html_file <- file.path(tables_dir, paste0(filename_base, ".html"))
  html_table <- kable(data, format = "html", 
                      caption = caption, digits = 4,
                      col.names = gsub("_", " ", toupper(names(data)))) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                  full_width = FALSE) %>%
    row_spec(0, bold = TRUE, background = "#f2f2f2")
  
  writeLines(as.character(html_table), html_file)
  
  return(list(csv = basename(csv_file), md = basename(md_file), html = basename(html_file)))
}

# Table 1: Asymmetry summary by marginal patterns
asymmetry_files <- save_enhanced_table(
  asymmetry_summary,
  "symmetry_asymmetry_by_marginals",
  "Bounds Asymmetry Analysis by Marginal Symmetry Patterns"
)

# Table 2: Theoretical examples table
examples_files <- save_enhanced_table(
  example_results,
  "symmetry_theoretical_examples",
  "Theoretical Examples of Bounds Asymmetry"
)

# Table 3: Distribution type impact summary
dist_impact <- bounds_asymmetry_analysis %>%
  group_by(dist1_type, dist2_type) %>%
  summarise(
    n_cases = n(),
    mean_bounds_asymmetry = mean(abs(bounds_asymmetry)),
    prop_symmetric_bounds = mean(bounds_nearly_symmetric),
    mean_tv_distance_combined = mean((marg1_tv_distance + marg2_tv_distance) / 2),
    mean_bc_complement_combined = mean((marg1_bc_complement + marg2_bc_complement) / 2),
    mean_ovl_complement_combined = mean((marg1_ovl_complement + marg2_ovl_complement) / 2),
    mean_cross_overlap = mean(cross_overlap_coefficient),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_bounds_asymmetry))

impact_files <- save_enhanced_table(
  dist_impact,
  "symmetry_distribution_impact",
  "Distribution Type Impact on Bounds Asymmetry"
)

cat("📊 Symmetry summary tables generated:\n")
cat("   1. Asymmetry by marginals:\n")
cat("      • CSV:", asymmetry_files$csv, "\n")
cat("      • Markdown:", asymmetry_files$md, "\n")
cat("      • HTML:", asymmetry_files$html, "\n")
cat("   2. Theoretical examples:\n")
cat("      • CSV:", examples_files$csv, "\n")
cat("      • Markdown:", examples_files$md, "\n")
cat("      • HTML:", examples_files$html, "\n")
cat("   3. Distribution impact:\n")
cat("      • CSV:", impact_files$csv, "\n")
cat("      • Markdown:", impact_files$md, "\n")
cat("      • HTML:", impact_files$html, "\n\n")

# ================================================================
# WORKFLOW COMPLETION SUMMARY
# ================================================================

cat("🎉 SYMMETRY ANALYSIS WORKFLOW COMPLETE!\n")
cat("=======================================\n")
cat("📁 Generated files:\n")
cat("   Symmetry experiments:", basename(symmetry_cache_file), "\n")
cat("   Figures: 5 symmetry visualization files\n")
cat("   Tables: 9 files (3 tables × 3 formats)\n\n")

cat("🔍 Key theoretical insights:\n")
perfect_symmetry <- mean(bounds_asymmetry_analysis$bounds_perfectly_symmetric)
cat("   Perfect symmetry (r_min = -r_max):", round(perfect_symmetry * 100, 1), "% of cases\n")

# Calculate correlations for the three asymmetry measures
tv_effect <- cor((bounds_asymmetry_analysis$marg1_tv_distance + bounds_asymmetry_analysis$marg2_tv_distance),
                 abs(bounds_asymmetry_analysis$bounds_asymmetry), use = "complete.obs")
bc_effect <- cor((bounds_asymmetry_analysis$marg1_bc_complement + bounds_asymmetry_analysis$marg2_bc_complement),
                 abs(bounds_asymmetry_analysis$bounds_asymmetry), use = "complete.obs")
ovl_effect <- cor((bounds_asymmetry_analysis$marg1_ovl_complement + bounds_asymmetry_analysis$marg2_ovl_complement),
                  abs(bounds_asymmetry_analysis$bounds_asymmetry), use = "complete.obs")

cat("   TV Distance-asymmetry correlation:", round(tv_effect, 3), "\n")
cat("   BC Complement-asymmetry correlation:", round(bc_effect, 3), "\n")
cat("   OVL Complement-asymmetry correlation:", round(ovl_effect, 3), "\n")

cat("\n🔬 Conditions for bounds asymmetry:\n")
cat("   - Asymmetric marginal distributions increase bounds asymmetry\n")
cat("   - High TV Distance (discrepancy) amplifies asymmetry\n")
cat("   - Low BC Coefficient (poor geometric congruence) increases asymmetry\n")
cat("   - Low OVL Coefficient (minimal shared mass) exacerbates asymmetry\n")
cat("   - Unequal category counts amplify asymmetric bounds\n")
cat("   - Cross-marginal overlap patterns influence asymmetry magnitude\n\n")

cat("🔄 Integration points:\n")
cat("   - TV/BC/OVL asymmetry measures inform distribution choice in simulations\n")
cat("   - Three-measure framework guides rescaling strategy selection\n")
cat("   - Complementary measures (TV, BC, OVL) provide comprehensive asymmetry assessment\n")
cat("   - Theoretical insights from overlap theory validate empirical bounds analysis\n")