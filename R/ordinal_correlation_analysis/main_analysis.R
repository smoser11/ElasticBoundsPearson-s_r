# main_analysis.R
# Master analysis coordination script

# Required libraries
library(dplyr)
library(ggplot2)
library(gridExtra)

# Source all module coordinators
source("1_bivariate_ordcats_correlation/bivariate_main.R")
source("2_correlation_matrices/matrices_main.R") 
source("3_fixes_and_rescaling/rescaling_main.R")

#' Run complete ordinal correlation analysis pipeline
#'
#' @param bes_data BES dataset
#' @param config Analysis configuration  
#' @return Complete analysis results
run_complete_analysis <- function(bes_data, config = get_default_config()) {
	
	cat("Starting comprehensive ordinal correlation analysis...\n")
	cat("Dataset:", nrow(bes_data), "variable pairs\n")
	
	# 1. Bivariate bounds analysis
	cat("\n=== BIVARIATE ANALYSIS ===\n")
	bivariate_results <- run_bivariate_analysis(bes_data, config)
	
	# 2. Matrix properties analysis
	cat("\n=== MATRIX ANALYSIS ===\n") 
	matrix_results <- run_matrix_analysis(bivariate_results$bounds_data, config)
	
	# 3. Rescaling analysis
	cat("\n=== RESCALING ANALYSIS ===\n")
	rescaling_results <- run_rescaling_analysis(bivariate_results$bounds_data, config)
	
	cat("\n=== ANALYSIS COMPLETE ===\n")
	
	return(list(
		bivariate = bivariate_results,
		matrices = matrix_results,
		rescaling = rescaling_results,
		config = config
	))
}

#' Get default analysis configuration
get_default_config <- function() {
	list(
		nsim = 1000,           # Permutation simulations
		matrix_trials = 50,    # Matrix trials  
		min_vars = 3,          # Min variables per matrix
		max_vars = 10,         # Max variables per matrix
		progress = TRUE,       # Show progress
		seed = 123            # Random seed
	)
}

# Test if run directly
if (!interactive()) {
	# Load your BES data here
	# bes_data <- read.csv("path/to/your/bes_data.csv")
	
	# Run analysis
	# results <- run_complete_analysis(bes_data)
	
	cat("Run: results <- run_complete_analysis(your_bes_data)\n")
}