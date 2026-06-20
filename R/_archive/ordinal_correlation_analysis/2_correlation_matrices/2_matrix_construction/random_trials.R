# 2_correlation_matrices/2_matrix_construction/random_trials.R
# Random correlation matrix testing framework

#' Run random correlation matrix trials
#'
#' @param results_df Data frame with correlation bounds results
#' @param n_trials Number of random trials to run
#' @param min_vars Minimum variables per trial
#' @param max_vars Maximum variables per trial
#' @param seed Random seed for reproducibility
#' @param enforce_psd Whether to enforce positive semidefiniteness
#' @return Data frame with trial results
test_random_correlation_matrices <- function(results_df, n_trials = 100, 
											 min_vars = 3, max_vars = 15,
											 seed = 123, enforce_psd = TRUE) {
	set.seed(seed)
	
	# Get all unique variables
	all_vars <- unique(c(results_df$var1, results_df$var2))
	n_all_vars <- length(all_vars)
	max_vars <- min(max_vars, n_all_vars)
	
	trial_results <- data.frame()
	
	for (trial in 1:n_trials) {
		# Randomly select variables and sample size
		n_vars <- sample(min_vars:max_vars, 1)
		selected_vars <- sample(all_vars, n_vars)
		
		# Get correlations for selected variables
		subset_corrs <- results_df %>%
			filter(var1 %in% selected_vars & var2 %in% selected_vars)
		
		# Skip if insufficient correlations
		if (nrow(subset_corrs) < (n_vars * (n_vars - 1) / 2 * 0.5)) {
			next  # Need at least 50% of possible pairs
		}
		
		# Construct both matrices
		orig_matrix <- tryCatch({
			construct_correlation_matrix(subset_corrs, use_rescaled = FALSE, enforce_psd = enforce_psd)
		}, error = function(e) NULL)
		
		rescaled_matrix <- tryCatch({
			construct_correlation_matrix(subset_corrs, use_rescaled = TRUE, enforce_psd = enforce_psd)
		}, error = function(e) NULL)
		
		if (is.null(orig_matrix) || is.null(rescaled_matrix)) next
		
		# Test matrix properties
		orig_props <- test_matrix_properties_comprehensive(orig_matrix)
		rescaled_props <- test_matrix_properties_comprehensive(rescaled_matrix)
		
		# Store trial results
		trial_result <- data.frame(
			trial = trial,
			n_variables = n_vars,
			n_correlations = nrow(subset_corrs),
			enforce_psd = enforce_psd,
			
			# Original matrix properties
			orig_psd = orig_props$is_positive_semidefinite,
			orig_invertible = orig_props$is_invertible,
			orig_condition = orig_props$condition_number,
			orig_min_eigenvalue = orig_props$min_eigenvalue,
			orig_valid = orig_props$is_valid_correlation_matrix,
			
			# Rescaled matrix properties  
			rescaled_psd = rescaled_props$is_positive_semidefinite,
			rescaled_invertible = rescaled_props$is_invertible,
			rescaled_condition = rescaled_props$condition_number,
			rescaled_min_eigenvalue = rescaled_props$min_eigenvalue,
			rescaled_valid = rescaled_props$is_valid_correlation_matrix,
			
			stringsAsFactors = FALSE
		)
		
		trial_results <- rbind(trial_results, trial_result)
	}
	
	return(trial_results)
}

#' Summarize matrix testing results
#'
#' @param test_results Output from test_random_correlation_matrices
#' @return Summary statistics
summarize_matrix_tests <- function(test_results) {
	if (nrow(test_results) == 0) {
		return(list(message = "No valid trials completed"))
	}
	
	summary_stats <- list(
		n_trials = nrow(test_results),
		
		# Original matrix statistics
		orig_psd_percent = 100 * mean(test_results$orig_psd, na.rm = TRUE),
		orig_invertible_percent = 100 * mean(test_results$orig_invertible, na.rm = TRUE),
		orig_valid_percent = 100 * mean(test_results$orig_valid, na.rm = TRUE),
		orig_mean_condition = mean(test_results$orig_condition[is.finite(test_results$orig_condition)], na.rm = TRUE),
		
		# Rescaled matrix statistics
		rescaled_psd_percent = 100 * mean(test_results$rescaled_psd, na.rm = TRUE),
		rescaled_invertible_percent = 100 * mean(test_results$rescaled_invertible, na.rm = TRUE),
		rescaled_valid_percent = 100 * mean(test_results$rescaled_valid, na.rm = TRUE),
		rescaled_mean_condition = mean(test_results$rescaled_condition[is.finite(test_results$rescaled_condition)], na.rm = TRUE)
	)
	
	return(summary_stats)
}