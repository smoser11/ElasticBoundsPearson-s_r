
# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

#' Simple rescaling function (temporary - will be moved to rescaling module)
rescale_correlation_simple <- function(r, r_min, r_max) {
	if (is.na(r) || is.na(r_min) || is.na(r_max)) return(NA)
	
	if (r < 0) {
		if (r_min < 0) {
			return((-1 / r_min) * r)
		} else {
			return(-1)
		}
	} else if (r > 0) {
		if (r_max > 0) {
			return((1 / r_max) * r)
		} else {
			return(1)
		}
	} else {
		return(0)
	}
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