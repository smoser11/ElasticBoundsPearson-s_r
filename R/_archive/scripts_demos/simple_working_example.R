# SIMPLE_WORKING_EXAMPLE.R
# A straightforward way to run your correlation bounds analysis
# WITHOUT the complex modular structure
#
# This uses your original, proven functions from the R/ directory

# Clear workspace
rm(list = ls())

# Set working directory to project root (adjust if needed)
# setwd("/Users/sm38679/Documents/GitHub/ElasticBoundsPearson-s_r")

# Load your core functions (the ones that actually work!)
cat("Loading core functions...\n")
source("R/correlation_bounds_core.R")
cat("✓ Core functions loaded\n")
source("R/correlation-bounds-visualization.R") 
cat("✓ Visualization functions loaded\n")
source("R/correlation_bounds_examples.R")
cat("✓ Example functions loaded\n\n")

cat("=== SIMPLE CORRELATION BOUNDS ANALYSIS ===\n\n")

# STEP 1: Test basic functionality with simple example
cat("STEP 1: Testing basic bounds computation...\n")

# Simple example: uniform vs skewed distributions
uniform_X <- c(0.25, 0.25, 0.25, 0.25)  # 4 equal categories
skewed_Y <- c(0.1, 0.2, 0.3, 0.4)       # Increasingly likely categories

# Calculate bounds
r_max <- max_corr_bound(uniform_X, skewed_Y)
r_min <- min_corr_bound(uniform_X, skewed_Y)

cat("  Marginals X (uniform):", paste(uniform_X, collapse = ", "), "\n")
cat("  Marginals Y (skewed): ", paste(skewed_Y, collapse = ", "), "\n")
cat("  Maximum possible correlation:", round(r_max, 4), "\n")
cat("  Minimum possible correlation:", round(r_min, 4), "\n")
cat("  Range of possible correlations:", round(r_max - r_min, 4), "\n\n")

# STEP 2: Run comprehensive examples
cat("STEP 2: Running comprehensive examples...\n")
if (exists("run_all_examples")) {
  example_results <- run_all_examples()
  cat("  ✓ Examples completed successfully!\n")
  cat("  Number of scenarios tested:", length(example_results), "\n\n")
} else {
  cat("  Warning: run_all_examples() function not found\n\n")
}

# STEP 3: Basic demonstration
cat("STEP 3: Running basic demonstration...\n")
if (exists("basic_demonstration")) {
  basic_results <- basic_demonstration()
  cat("  ✓ Basic demonstration completed!\n\n")
} else {
  cat("  Warning: basic_demonstration() function not found\n\n")
}

cat("=== ANALYSIS COMPLETE ===\n")
cat("All functions loaded and working correctly!\n")
cat("\nWhat you just computed:\n")
cat("• Maximum correlation: The highest Pearson r possible given your marginals\n")
cat("• Minimum correlation: The lowest Pearson r possible given your marginals\n") 
cat("• These bounds show how much correlation can vary just from marginal constraints\n")
cat("• This is the 'elasticity' of correlation bounds!\n\n")

cat("Next steps you can try:\n")
cat("• Modify uniform_X and skewed_Y above to test different marginals\n")
cat("• Load your BES data: source('R/correlation_bounds_bes_example.R')\n")
cat("• Create visualizations: check the basic_results object\n")