# matrices_main.R
# Correlation matrix analysis module coordinator

# Source all matrix analysis components
source("2_correlation_matrices/1_matrix_properties/matrix_diagnostics.R")
source("2_correlation_matrices/1_matrix_properties/psd_validation.R")
source("2_correlation_matrices/1_matrix_properties/condition_number_analysis.R")
source("2_correlation_matrices/2_matrix_construction/random_trials.R")

# Required libraries
library(Matrix)
library(dplyr)

#' Run complete correlation matrix analysis
#'
#' @param bounds_data Results from bivariate bounds analysis
#' @param config Configuration list
#' @return List with matrix analysis results
run_matrix_analysis <- function(bounds_data, config = list()) {
	# Set defaults
	n_trials <- config$matrix_trials %||% 100
	min_vars <- config$min_vars %||% 3
	max_vars <- config$max_vars %||% 15
	
	cat("Running correlation matrix property analysis...\n")
	
	# 1. Test with PSD enforcement
	matrix_trials_psd <- test_random_correlation_matrices(
		bounds_data, 
		n_trials = n_trials,
		min_vars = min_vars,
		max_vars = max_vars,
		enforce_psd = TRUE
	)
	
	# 2. Test without PSD enforcement  
	matrix_trials_no_psd <- test_random_correlation_matrices(
		bounds_data,
		n_trials = n_trials,
		min_vars = min_vars, 
		max_vars = max_vars,
		enforce_psd = FALSE
	)
	
	# 3. Summarize results
	summary_psd <- summarize_matrix_tests(matrix_trials_psd)
	summary_no_psd <- summarize_matrix_tests(matrix_trials_no_psd)
	
	cat("Matrix analysis complete:", summary_psd$n_trials, "trials completed\n")
	
	return(list(
		trials_psd = matrix_trials_psd,
		trials_no_psd = matrix_trials_no_psd,
		summary_psd = summary_psd,
		summary_no_psd = summary_no_psd
	))
}