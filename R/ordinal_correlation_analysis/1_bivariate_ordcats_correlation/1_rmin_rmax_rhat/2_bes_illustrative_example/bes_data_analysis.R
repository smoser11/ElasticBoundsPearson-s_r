
# bes_data_analysis.R
# BES-specific correlation bounds analysis functions
# Adapted from correlation_bounds_bes.R

# Source required dependencies
if (!exists("max_corr_bound")) {
	source("1_bivariate_ordcats_correlation/1_rmin_rmax_rhat/1_monte_carlo_simulation/core_bounds_functions.R")
}

# Required libraries
library(dplyr)  # If used for data manipulation



	# NEW (in bes_data_analysis.R):
analyze_bes_corr_bounds <- function(row, nsim = 1000, use_prop = FALSE, return_simulations = FALSE) {
		# Add dependency source at top of function
		# (since this will be in a different file structure)
		
		# Extract the number of categories for each variable
		var1_cats <- row$var1cats
		var2_cats <- row$var2cats
		n_obs <- row$nobs
		
		
		# Extract the marginal distributions, handling variables that start at either 0 or 1
		if (use_prop) {
			# Use proportions if requested
			var1_marginals <- c()
			var2_marginals <- c()
			
			# Determine start category for var1 (0 or 1)
			var1_start <- ifelse("var1prop0" %in% names(row) && !is.na(row$var1prop0) && row$var1prop0 > 0, 0, 1)
			var1_end <- var1_start + var1_cats - 1
			
			for (i in var1_start:var1_end) {
				prop_col <- paste0("var1prop", i)
				if (prop_col %in% names(row) && !is.na(row[[prop_col]])) {
					var1_marginals <- c(var1_marginals, row[[prop_col]])
				} else {
					var1_marginals <- c(var1_marginals, 0)  # Add 0 for missing categories
				}
			}
			
			# Determine start category for var2 (0 or 1)
			var2_start <- ifelse("var2prop0" %in% names(row) && !is.na(row$var2prop0) && row$var2prop0 > 0, 0, 1)
			var2_end <- var2_start + var2_cats - 1
			
			for (i in var2_start:var2_end) {
				prop_col <- paste0("var2prop", i)
				if (prop_col %in% names(row) && !is.na(row[[prop_col]])) {
					var2_marginals <- c(var2_marginals, row[[prop_col]])
				} else {
					var2_marginals <- c(var2_marginals, 0)  # Add 0 for missing categories
				}
			}
		} else {
			# Use frequencies (default)
			var1_marginals <- c()
			var2_marginals <- c()
			
			# Determine start category for var1 (0 or 1)
			var1_start <- ifelse("var1freq0" %in% names(row) && !is.na(row$var1freq0) && row$var1freq0 > 0, 0, 1)
			var1_end <- var1_start + var1_cats - 1
			
			for (i in var1_start:var1_end) {
				freq_col <- paste0("var1freq", i)
				if (freq_col %in% names(row) && !is.na(row[[freq_col]])) {
					var1_marginals <- c(var1_marginals, row[[freq_col]])
				} else {
					var1_marginals <- c(var1_marginals, 0)  # Add 0 for missing categories
				}
			}
			
			# Determine start category for var2 (0 or 1)
			var2_start <- ifelse("var2freq0" %in% names(row) && !is.na(row$var2freq0) && row$var2freq0 > 0, 0, 1)
			var2_end <- var2_start + var2_cats - 1
			
			for (i in var2_start:var2_end) {
				freq_col <- paste0("var2freq", i)
				if (freq_col %in% names(row) && !is.na(row[[freq_col]])) {
					var2_marginals <- c(var2_marginals, row[[freq_col]])
				} else {
					var2_marginals <- c(var2_marginals, 0)  # Add 0 for missing categories
				}
			}
		}
		
		# Remove any NAs and ensure non-negative values
		var1_marginals <- pmax(0, var1_marginals, na.rm = TRUE)
		var2_marginals <- pmax(0, var2_marginals, na.rm = TRUE)
		
		# Skip if we have empty marginals
		if (sum(var1_marginals) == 0 || sum(var2_marginals) == 0) {
			warning("Empty marginal distributions detected.")
			return(list(
				r_min = NA,
				r_max = NA,
				ci_lower = NA,
				ci_upper = NA,
				observed_r = row$corr,
				success = FALSE,
				message = "Empty marginal distributions"
			))
		}
		
		# Calculate theoretical bounds
		r_min <- min_corr_bound(var1_marginals, var2_marginals, sample_size = n_obs)
		r_max <- max_corr_bound(var1_marginals, var2_marginals, sample_size = n_obs)
		
		# Simulate permutation distribution
		r_sim <- simulate_permutation_distribution(var1_marginals, var2_marginals, 
												   nsim = nsim, sample_size = n_obs)
		
		# Calculate 95% confidence interval
		ci_bounds <- quantile(r_sim, probs = c(0.025, 0.975))
		
		# Rescale the observed correlation based on theoretical bounds
		# Maps r_min to -1, 0 to 0, and r_max to +1
		if (row$corr < 0) {
			# For negative correlations, rescale between r_min and 0
			if (r_min < 0) {  # Ensure we don't divide by zero
				r_rescaled <- (-1 / r_min) * row$corr
			} else {
				# If r_min >= 0, rescaling negative correlation isn't possible
				# (this should be rare since r_min is usually negative)
				r_rescaled <- -1  # Default to -1 in this edge case
			}
		} else if (row$corr > 0) {
			# For positive correlations, rescale between 0 and r_max
			if (r_max > 0) {  # Ensure we don't divide by zero
				r_rescaled <- (1 / r_max) * row$corr
			} else {
				# If r_max <= 0, rescaling positive correlation isn't possible
				# (this should be rare since r_max is usually positive)
				r_rescaled <- 1  # Default to 1 in this edge case
			}
		} else {
			# If correlation is exactly 0, keep it as 0
			r_rescaled <- 0
		}
		
		# Prepare results
		result <- list(
			r_min = r_min,
			r_max = r_max,
			ci_lower = ci_bounds[1],
			ci_upper = ci_bounds[2],
			observed_r = row$corr,
			r_rescaled = r_rescaled,
			success = TRUE,
			var1_categories = var1_cats,
			var2_categories = var2_cats,
			n_obs = n_obs,
			var1_start = var1_start,
			var2_start = var2_start,
			var1_marginals = var1_marginals,
			var2_marginals = var2_marginals
		)
		
		# Add simulations if requested
		if (return_simulations) {
			result$simulations <- r_sim
		}
		
		return(result)
		cat(class(result))
	}
		# [Copy rest of function body exactly, but update internal function calls]
		# Change: min_corr_bound() → min_corr_bound() (no change needed if properly sourced)
		# Change: max_corr_bound() → max_corr_bound() (no change needed if properly sourced)
	
	


#' Apply correlation bounds analysis to an entire data frame
#'
#' This function applies the analyze_corr_bounds function to each row of a data frame,
#' allowing for batch processing of many variable pairs.
#'
#' @param df A data frame where each row represents a pair of ordinal variables
#' @param nsim Number of simulations for the permutation distribution (default: 1000)
#' @param use_prop Logical; if TRUE, use proportion columns instead of frequency columns (default: FALSE)
#' @param return_simulations Logical; if TRUE, return the simulated r values (default: FALSE)
#' @param progress Logical; if TRUE, display a progress indicator (default: TRUE)
#'
#' @return A data frame with correlation bounds analysis results for each row
#'
#' @examples
#' # Assuming 'bes_data' is your data frame
#' all_results <- analyze_all_corr_bounds(bes_data)
analyze_all_bes_bounds <- function(df, nsim = 1000, use_prop = FALSE, 
									return_simulations = FALSE, progress = TRUE) {
	df <- bes_data
	# Initialize results list
	results_list <- list()
	
	# Progress tracking
	n_rows <- nrow(df)
	if (progress) cat("Processing", n_rows, "variable pairs...\n")
	
	# Process each row
	for (i in 1:n_rows) {
		#if (progress && i %% 100 == 0) cat("Processed", i, "of", n_rows, "pairs\n")
		
		# Extract row
		row <- df[i, ]
		
		# Analyze bounds
		bounds_result <- analyze_corr_bounds(row, nsim = nsim, use_prop = use_prop, 
											 return_simulations = return_simulations)
		
		# Add original variables
		bounds_result[var1] <- row$var1
		bounds_result$var2 <- row$var2
		
		# Store results
		results_list[[i]] <- bounds_result
	}
	
	# Convert list to data frame, managing nested structures
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
		var1_start = sapply(results_list, function(x) x$var1_start),
		var2_start = sapply(results_list, function(x) x$var2_start),
		n_obs = sapply(results_list, function(x) x$n_obs),
		success = sapply(results_list, function(x) x$success)
	)
	
	# Add simulations if requested (as a list column)
	if (return_simulations) {
		results_df$simulations <- lapply(results_list, function(x) x$simulations)
	}
	
	return(results_df)
}


#' Extract marginal distributions from BES data row
#'
#' @param row Data frame row with BES variable information
#' @param var_prefix Variable prefix ("var1" or "var2")
#' @param n_cats Number of categories
#' @param use_prop Use proportions instead of frequencies
#' @return Vector of marginal counts/proportions
extract_bes_marginals <- function(row, var_prefix, n_cats, use_prop = FALSE) {
	suffix <- if (use_prop) "prop" else "freq"
	marginals <- c()
	
	# Determine starting category (0 or 1)
	start_col_0 <- paste0(var_prefix, suffix, "0")
	start_col_1 <- paste0(var_prefix, suffix, "1")
	
	# Check if category 0 exists and has non-zero values
	if (start_col_0 %in% names(row) && !is.na(row[[start_col_0]]) && row[[start_col_0]] > 0) {
		var_start <- 0
	} else {
		var_start <- 1
	}
	
	var_end <- var_start + n_cats - 1
	
	# Extract marginals
	for (i in var_start:var_end) {
		col_name <- paste0(var_prefix, suffix, i)
		if (col_name %in% names(row) && !is.na(row[[col_name]])) {
			marginals <- c(marginals, row[[col_name]])
		} else {
			marginals <- c(marginals, 0)  # Add 0 for missing categories
		}
	}
	
	# Remove NAs and ensure non-negative values
	marginals <- pmax(0, marginals, na.rm = TRUE)
	return(marginals)
}

#' Determine starting category for BES variables (0 or 1)
#'
#' @param row Data frame row
#' @param var_prefix Variable prefix ("var1" or "var2")
#' @param suffix "freq" or "prop"
#' @return Starting category index (0 or 1)
determine_bes_start_category <- function(row, var_prefix, suffix) {
	# Check if category 0 exists and has non-zero values
	col0_name <- paste0(var_prefix, suffix, "0")
	if (col0_name %in% names(row) && !is.na(row[[col0_name]]) && row[[col0_name]] > 0) {
		return(0)
	} else {
		return(1)
	}
}


#' Validate BES data row for analysis
#'
#' @param row Single row from BES dataset
#' @return List with validation results
validate_bes_row <- function(row) {
	errors <- c()
	warnings <- c()
	
	# Check required columns
	required_cols <- c("var1", "var2", "corr", "nobs", "var1cats", "var2cats")
	missing_cols <- setdiff(required_cols, names(row))
	if (length(missing_cols) > 0) {
		errors <- c(errors, paste("Missing columns:", paste(missing_cols, collapse = ", ")))
	}
	
	# Check data validity
	if ("nobs" %in% names(row) && (is.na(row$nobs) || row$nobs < 10)) {
		warnings <- c(warnings, "Sample size too small (< 10)")
	}
	
	if ("corr" %in% names(row) && (is.na(row$corr) || abs(row$corr) > 1)) {
		errors <- c(errors, "Invalid correlation value")
	}
	
	if ("var1cats" %in% names(row) && (is.na(row$var1cats) || row$var1cats < 2)) {
		errors <- c(errors, "Invalid number of categories for var1")
	}
	
	if ("var2cats" %in% names(row) && (is.na(row$var2cats) || row$var2cats < 2)) {
		errors <- c(errors, "Invalid number of categories for var2")
	}
	
	return(list(
		valid = length(errors) == 0,
		errors = errors,
		warnings = warnings
	))
}

#' Create error result for failed analyses
create_error_result <- function(row, message) {
	list(
		r_min = NA, r_max = NA,
		ci_lower = NA, ci_upper = NA,
		observed_r = if("corr" %in% names(row)) row$corr else NA,
		r_rescaled = NA,
		success = FALSE,
		message = message
	)
}

#' Generate summary report for BES bounds analysis
#'
#' @param results_df The results data frame from analyze_all_bes_bounds
#' @return A data frame with summary statistics
generate_bes_bounds_summary <- function(results_df) {
	results_df <- results_df %>%
		filter(success == TRUE) %>%
		mutate(
			# Check if observed is outside empirical CI
			outside_ci = observed_r < ci_lower | observed_r > ci_upper,
			
			# Theoretical range
			theo_range = r_max - r_min,
			
			# Position within theoretical range
			rel_position = (observed_r - r_min) / (r_max - r_min)
		)
	
	return(results_df)
}
