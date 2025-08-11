# 3_fixes_and_rescaling/1_simple_rescaling/linear_rescaling.R
# Simple linear rescaling methods

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
	
	return(list(
		summary = rescaling_summary,
		detailed_results = results_df
	))
}

#' Analyze rescaling by variable characteristics
#'
#' @param results_df Data frame with rescaling results
#' @return Analysis by variable properties
analyze_rescaling_by_characteristics <- function(results_df) {
	# By number of categories
	results_df$total_categories <- results_df$var1_categories + results_df$var2_categories
	results_df$bounds_range <- results_df$r_max - results_df$r_min
	
	# Summary by categories
	if ("var1_categories" %in% names(results_df)) {
		category_analysis <- results_df %>%
			group_by(total_categories) %>%
			summarise(
				n_pairs = n(),
				mean_rescaling_magnitude = mean(rescaling_magnitude, na.rm = TRUE),
				mean_bounds_range = mean(bounds_range, na.rm = TRUE),
				.groups = "drop"
			)
	} else {
		category_analysis <- NULL
	}
	
	return(list(
		by_categories = category_analysis,
		detailed_data = results_df
	))
}