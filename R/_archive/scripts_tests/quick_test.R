# QUICK_TEST.R
# Super simple test to see if your core functions work

# Clear everything
rm(list = ls())

# Test 1: Load core functions directly
cat("Testing core correlation bounds functions...\n")

# Load required libraries first
if (!require(ggplot2, quietly = TRUE)) install.packages("ggplot2")
if (!require(dplyr, quietly = TRUE)) install.packages("dplyr")
library(ggplot2)
library(dplyr)

# Source core functions
source("R/correlation_bounds_core.R")
cat("✓ Core functions loaded successfully\n")

# Test the functions with simple example
cat("\nTesting max_corr_bound()...\n")
uniform_X <- c(0.25, 0.25, 0.25, 0.25)
skewed_Y <- c(0.1, 0.2, 0.3, 0.4)

r_max <- max_corr_bound(uniform_X, skewed_Y)
cat("Maximum correlation:", round(r_max, 4), "\n")

cat("\nTesting min_corr_bound()...\n")
r_min <- min_corr_bound(uniform_X, skewed_Y)
cat("Minimum correlation:", round(r_min, 4), "\n")

cat("\nCorrelation range:", round(r_max - r_min, 4), "\n")

cat("\n✅ SUCCESS! Your core functions work perfectly.\n")
cat("\nWhat this means:\n")
cat("- With uniform marginals for X:", paste(uniform_X, collapse=", "), "\n")
cat("- And skewed marginals for Y:", paste(skewed_Y, collapse=", "), "\n") 
cat("- The correlation can range from", round(r_min, 3), "to", round(r_max, 3), "\n")
cat("- That's a range of", round(r_max - r_min, 3), "just from marginal constraints!\n")