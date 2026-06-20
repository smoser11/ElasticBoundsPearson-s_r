# 1_bivariate_ordcats_correlation/3_visualization/bounds_visualization.R
# Visualization functions for correlation bounds analysis

# Required libraries
library(ggplot2)
library(gridExtra)

#' Create correlation bounds scatter plot
#'
#' @param results_df Data frame with bounds analysis results
#' @param color_by Variable to color points by
#' @return ggplot object
plot_bounds_scatter <- function(results_df, color_by = "bounds_range") {
	# Calculate bounds range if not present
	if (!"bounds_range" %in% names(results_df)) {
		results_df$bounds_range <- results_df$r_max - results_df$r_min
	}
	
	p <- ggplot(results_df, aes(x = r_min, y = r_max)) +
		geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
		geom_hline(yintercept = 0, linetype = "dotted", color = "gray") +
		geom_vline(xintercept = 0, linetype = "dotted", color = "gray") +
		coord_equal() +
		theme_minimal() +
		labs(
			title = "Correlation Bounds Analysis",
			x = "Minimum Possible Correlation (r_min)",
			y = "Maximum Possible Correlation (r_max)"
		)
	
	# Add color mapping
	if (color_by == "bounds_range") {
		p <- p + 
			geom_point(aes(color = bounds_range), alpha = 0.7, size = 1.5) +
			scale_color_viridis_c(name = "Bounds\nRange")
	} else {
		p <- p + geom_point(alpha = 0.7, size = 1.5, color = "steelblue")
	}
	
	return(p)
}

#' Create rescaling comparison plot
#'
#' @param results_df Data frame with original and rescaled correlations
#' @return ggplot object
plot_rescaling_comparison <- function(results_df) {
	if (!"rescaling_magnitude" %in% names(results_df)) {
		results_df$rescaling_magnitude <- abs(results_df$r_rescaled) - abs(results_df$observed_r)
	}
	
	p <- ggplot(results_df, aes(x = observed_r, y = r_rescaled)) +
		geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
		geom_hline(yintercept = 0, linetype = "dotted", color = "gray") +
		geom_vline(xintercept = 0, linetype = "dotted", color = "gray") +
		geom_point(aes(color = rescaling_magnitude), alpha = 0.7, size = 1.5) +
		scale_color_gradient2(low = "blue", mid = "gray", high = "red", 
							  midpoint = 0, name = "Rescaling\nEffect") +
		coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
		theme_minimal() +
		labs(
			title = "Original vs Rescaled Correlations",
			x = "Original Correlation",
			y = "Rescaled Correlation"
		)
	
	return(p)
}