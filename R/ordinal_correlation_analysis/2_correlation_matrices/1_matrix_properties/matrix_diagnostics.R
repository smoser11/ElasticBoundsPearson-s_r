# 2_correlation_matrices/1_matrix_properties/matrix_diagnostics.R
# Core matrix construction and property testing

# Source required dependencies using here() for robust paths
library(here)
source(here("R", "ordinal_correlation_analysis", "2_correlation_matrices", "1_matrix_properties", "psd_validation.R"))

# Required libraries
library(Matrix)

#' Construct correlation matrix from pairwise correlations
#'
#' @param correlations Data frame with var1, var2, and correlation columns
#' @param use_rescaled Use rescaled correlations instead of observed
#' @param enforce_psd Ensure positive semidefiniteness
#' @return Correlation matrix
construct_correlation_matrix <- function(correlations, use_rescaled = FALSE, enforce_psd = TRUE) {
	# Input validation
	if (!is.data.frame(correlations)) {
		stop("correlations must be a data frame")
	}
	
	required_cols <- c("var1", "var2")
	if (!all(required_cols %in% names(correlations))) {
		stop("correlations must have var1 and var2 columns")
	}
	
	# Choose correlation column
	corr_column <- ifelse(use_rescaled, "r_rescaled", "observed_r")
	if (!corr_column %in% names(correlations)) {
		stop("Missing correlation column: ", corr_column)
	}
	
	# Identify all unique variables
	all_vars <- unique(c(correlations$var1, correlations$var2))
	n_vars <- length(all_vars)
	
	if (n_vars < 2) {
		stop("Need at least 2 variables to construct a matrix")
	}
	
	# Create variable index mapping
	var_indices <- setNames(1:n_vars, all_vars)
	
	# Initialize correlation matrix
	corr_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)
	diag(corr_matrix) <- 1
	
	# Fill correlation matrix
	for (i in 1:nrow(correlations)) {
		var1 <- correlations$var1[i]
		var2 <- correlations$var2[i]
		corr_value <- correlations[[corr_column]][i]
		
		if (!is.na(corr_value) && var1 %in% all_vars && var2 %in% all_vars) {
			idx1 <- var_indices[as.character(var1)]
			idx2 <- var_indices[as.character(var2)]
			
			# Fill both triangles for symmetry
			corr_matrix[idx1, idx2] <- corr_value
			corr_matrix[idx2, idx1] <- corr_value
		}
	}
	
	# Enforce positive semidefiniteness if requested
	if (enforce_psd) {
		corr_matrix <- make_matrix_psd(corr_matrix)
	}
	
	# Add variable names (ensure matrix still has dimensions)
	if (is.matrix(corr_matrix) && nrow(corr_matrix) == n_vars) {
		rownames(corr_matrix) <- as.character(all_vars)
		colnames(corr_matrix) <- as.character(all_vars)
	}
	
	return(corr_matrix)
}

#' Comprehensive matrix property testing (FIXED VERSION)
#'
#' @param matrix Correlation matrix to test
#' @param tolerance Numerical tolerance for tests
#' @return List with detailed matrix properties
test_matrix_properties_comprehensive <- function(matrix, tolerance = 1e-10) {
	# Initialize results list
	results <- list()
	
	# Input validation
	if (!is.matrix(matrix)) {
		stop("Input must be a matrix")
	}
	
	if (nrow(matrix) != ncol(matrix)) {
		stop("Input must be a square matrix")
	}
	
	# Basic structural properties
	results$is_symmetric <- isSymmetric(matrix, tol = tolerance)
	results$has_unit_diagonal <- all(abs(diag(matrix) - 1) < tolerance)
	
	# Ensure symmetry for eigenvalue computation
	work_matrix <- matrix
	if (!results$is_symmetric) {
		work_matrix <- (matrix + t(matrix)) / 2
	}
	
	# Eigenvalue analysis with error handling
	tryCatch({
		eigenvalues <- eigen(work_matrix, symmetric = TRUE, only.values = TRUE)$values
		results$eigenvalues <- eigenvalues
		results$min_eigenvalue <- min(eigenvalues)
		results$max_eigenvalue <- max(eigenvalues)
		results$n_negative_eigenvalues <- sum(eigenvalues < -tolerance)
		results$n_zero_eigenvalues <- sum(abs(eigenvalues) < tolerance)
	}, error = function(e) {
		warning("Eigenvalue computation failed: ", e$message)
		results$eigenvalues <- NA
		results$min_eigenvalue <- NA
		results$max_eigenvalue <- NA
		results$n_negative_eigenvalues <- NA
		results$n_zero_eigenvalues <- NA
	})
	
	# Matrix definiteness (check for NA first)
	if (!is.na(results$min_eigenvalue)) {
		results$is_positive_definite <- results$min_eigenvalue > tolerance
		results$is_positive_semidefinite <- results$min_eigenvalue > -tolerance
		results$is_invertible <- results$min_eigenvalue > tolerance
	} else {
		results$is_positive_definite <- NA
		results$is_positive_semidefinite <- NA
		results$is_invertible <- NA
	}
	
	# Condition number analysis
	if (!is.na(results$is_invertible) && results$is_invertible) {
		results$condition_number <- results$max_eigenvalue / results$min_eigenvalue
		results$condition_category <- classify_condition_number(results$condition_number)
		results$log_condition_number <- log10(results$condition_number)
	} else {
		results$condition_number <- Inf
		results$condition_category <- "Singular"
		results$log_condition_number <- Inf
	}
	
	# Matrix determinant and rank
	if (!is.na(results$min_eigenvalue)) {
		results$determinant <- prod(results$eigenvalues)
		results$log_determinant <- sum(log(results$eigenvalues[results$eigenvalues > tolerance]))
		results$rank <- sum(results$eigenvalues > tolerance)
		results$rank_deficiency <- nrow(matrix) - results$rank
	} else {
		results$determinant <- NA
		results$log_determinant <- NA
		results$rank <- NA
		results$rank_deficiency <- NA
	}
	
	# Numerical stability measures
	if (!is.na(results$min_eigenvalue)) {
		results$spectral_radius <- max(abs(results$eigenvalues))
		results$nuclear_norm <- sum(abs(results$eigenvalues))
	} else {
		results$spectral_radius <- NA
		results$nuclear_norm <- NA
	}
	
	results$frobenius_norm <- sqrt(sum(matrix^2))
	
	# Correlation-specific properties
	if (nrow(matrix) > 1) {
		upper_tri_vals <- matrix[upper.tri(matrix)]
		results$max_off_diagonal <- max(abs(upper_tri_vals))
		results$min_off_diagonal <- min(abs(upper_tri_vals))
		results$mean_correlation <- mean(upper_tri_vals)
		results$n_high_correlations <- sum(abs(upper_tri_vals) > 0.8)
	} else {
		results$max_off_diagonal <- 0
		results$min_off_diagonal <- 0
		results$mean_correlation <- 0
		results$n_high_correlations <- 0
	}
	
	# Overall validity assessment
	results$is_valid_correlation_matrix <- 
		!is.na(results$is_symmetric) && results$is_symmetric && 
		!is.na(results$has_unit_diagonal) && results$has_unit_diagonal && 
		!is.na(results$is_positive_semidefinite) && results$is_positive_semidefinite &&
		(!is.na(results$max_off_diagonal) && results$max_off_diagonal <= 1 + tolerance)
	
	return(results)
}

#' Classify condition number for interpretation (FIXED VERSION)
#'
#' @param cond_num Condition number
#' @return Character classification
classify_condition_number <- function(cond_num) {
	if (is.na(cond_num)) return("Unknown")
	if (is.infinite(cond_num)) return("Singular")
	if (cond_num <= 10) return("Well-conditioned")
	if (cond_num <= 100) return("Moderately ill-conditioned")
	if (cond_num <= 1000) return("Ill-conditioned")
	return("Severely ill-conditioned")
}
