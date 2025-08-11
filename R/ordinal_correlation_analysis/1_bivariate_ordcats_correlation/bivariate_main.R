# bivariate_main.R
# Bivariate correlation analysis module coordinator

# Source all bivariate analysis components
source("1_bivariate_ordcats_correlation/1_rmin_rmax_rhat/1_monte_carlo_simulation/core_bounds_functions.R")
source("1_bivariate_ordcats_correlation/1_rmin_rmax_rhat/2_bes_illustrative_example/bes_data_analysis.R")
source("1_bivariate_ordcats_correlation/2_asymmetry_analysis/asymmetry_measures.R")
source("1_bivariate_ordcats_correlation/3_visualization/bounds_visualization.R")

# Required libraries for this module
library(dplyr)
library(ggplot2)

#' Run complete bivariate correlation analysis
#'
#' @param bes_data BES dataset
#' @param config Configuration list
#' @return List with bounds and asymmetry results
run_bivariate_analysis <- function(bes_data, config = list()) {
	# Set defaults
	nsim <- config$nsim %||% 1000
	progress <- config$progress %||% TRUE
	
	cat("Running bivariate correlation bounds analysis...\n")
	
	# 1. Calculate correlation bounds
	bounds_results <- analyze_all_bes_bounds(
		bes_data, 
		nsim = nsim,
		progress = progress
	)
	
	# 2. Analyze asymmetry
	asymmetry_results <- analyze_bounds_asymmetry_comprehensive(bounds_results)
	
	# 3. Generate summary
	bounds_summary <- generate_bes_bounds_summary(bounds_results)
	
	cat("Bivariate analysis complete for", nrow(bounds_results), "variable pairs\n")
	
	return(list(
		bounds_data = bounds_results,
		bounds_summary = bounds_summary,
		asymmetry = asymmetry_results
	))
}

#' Default null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x