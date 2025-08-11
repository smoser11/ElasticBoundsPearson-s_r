# 1_bivariate_ordcats_correlation/2_asymmetry_analysis/asymmetry_measures.R
# Asymmetry analysis functions for correlation bounds

#' Calculate asymmetry measures for correlation bounds
#'
#' @param r_min Minimum correlation bound
#' @param r_max Maximum correlation bound
#' @return List with asymmetry measures
calculate_asymmetry_measures <- function(r_min, r_max) {
	# Basic asymmetry (should be 0 for symmetric bounds)
	basic_asymmetry <- r_max + r_min
	
	# Range and center
	bounds_range <- r_max - r_min
	bounds_center <- (r_max + r_min) / 2
	
	# Normalized asymmetry measures
	if (bounds_range > 1e-10) {
		normalized_asymmetry <- basic_asymmetry / bounds_range
		relative_shift <- bounds_center / (bounds_range / 2)
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
		asymmetry_direction <- "Positive Skewed"
	} else if (basic_asymmetry < -0.01) {
		asymmetry_direction <- "Negative Skewed"
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

#' Analyze asymmetry across all variable pairs
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
	
	# Summary statistics
	asymmetry_summary <- list(
		mean_asymmetry = mean(asymmetry_df$basic_asymmetry, na.rm = TRUE),
		median_asymmetry = median(asymmetry_df$basic_asymmetry, na.rm = TRUE),
		sd_asymmetry = sd(asymmetry_df$basic_asymmetry, na.rm = TRUE),
		prop_symmetric = mean(asymmetry_df$asymmetry_class == "Symmetric", na.rm = TRUE),
		prop_mildly_asymmetric = mean(asymmetry_df$asymmetry_class == "Mildly Asymmetric", na.rm = TRUE),
		prop_moderately_asymmetric = mean(asymmetry_df$asymmetry_class == "Moderately Asymmetric", na.rm = TRUE),
		prop_highly_asymmetric = mean(asymmetry_df$asymmetry_class == "Highly Asymmetric", na.rm = TRUE),
		max_asymmetry = max(abs(asymmetry_df$basic_asymmetry), na.rm = TRUE)
	)
	
	return(list(
		detailed_results = results_with_asymmetry,
		summary = asymmetry_summary,
		asymmetry_by_class = table(asymmetry_df$asymmetry_class),
		asymmetry_by_direction = table(asymmetry_df$asymmetry_direction)
	))
}