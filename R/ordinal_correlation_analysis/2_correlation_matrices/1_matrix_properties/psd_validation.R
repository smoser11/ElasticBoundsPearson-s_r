# psd_validation.R
# Positive semidefinite matrix validation and correction

#' Test if matrix is positive semidefinite
is_positive_semidefinite <- function(matrix, tolerance = 1e-10) {
	eigen_values <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
	return(min(eigen_values) >= -tolerance)
}

#' Make matrix positive semidefinite
make_matrix_psd_robust <- function(matrix, min_eigenvalue = 1e-8) {
	# Extract this logic from the original enforce_psd section
	# in construct_correlation_matrix()
	#' Make matrix positive semidefinite (robust version)
		# Input validation
		if (!is.matrix(matrix) || nrow(matrix) != ncol(matrix)) {
			stop("Input must be a square matrix")
		}
		
		n <- nrow(matrix)
		
		# Ensure symmetry
		matrix <- (matrix + t(matrix)) / 2
		
		# Eigendecomposition with error handling
		tryCatch({
			eigen_decomp <- eigen(matrix, symmetric = TRUE)
			eigenvalues <- eigen_decomp$values
			eigenvectors <- eigen_decomp$vectors
			
			# Check if adjustment is needed
			if (min(eigenvalues) < -1e-12) {  # Use small tolerance
				# Adjust eigenvalues
				eigenvalues[eigenvalues < min_eigenvalue] <- min_eigenvalue
				
				# Reconstruct
				matrix_new <- eigenvectors %*% diag(eigenvalues) %*% t(eigenvectors)
				
				# Rescale diagonal to 1 (for correlation matrices)
				diag_elements <- diag(matrix_new)
				if (all(diag_elements > 0)) {
					scaling <- diag(1/sqrt(diag_elements))
					matrix_new <- scaling %*% matrix_new %*% scaling
				}
				
				# Final cleanup
				matrix_new <- (matrix_new + t(matrix_new)) / 2
				diag(matrix_new) <- 1
				
				return(matrix_new)
			} else {
				# No adjustment needed
				return(matrix)
			}
		}, error = function(e) {
			warning("Eigendecomposition failed: ", e$message)
			return(matrix)  # Return original matrix if decomposition fails
		})
	}

#' Make matrix positive semidefinite using eigenvalue adjustment
#'
#' @param matrix Input correlation matrix
#' @param min_eigenvalue Minimum allowed eigenvalue
#' @return Adjusted positive semidefinite matrix
make_matrix_psd <- function(matrix, min_eigenvalue = 1e-8) {
	# Input validation
	if (!is.matrix(matrix)) {
		stop("Input must be a matrix")
	}
	
	if (nrow(matrix) != ncol(matrix)) {
		stop("Input must be a square matrix")
	}
	
	# Store original dimensions and names
	orig_dim <- dim(matrix)
	orig_rownames <- rownames(matrix)
	orig_colnames <- colnames(matrix)
	
	# Ensure matrix is symmetric
	matrix <- (matrix + t(matrix)) / 2
	
	# Check if already positive semidefinite
	eigenvalues <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
	
	if (min(eigenvalues) >= -1e-12) {
		# Already PSD, return as-is
		return(matrix)
	}
	
	# Need to adjust - do full eigendecomposition
	eigen_decomp <- eigen(matrix, symmetric = TRUE)
	eigenvalues <- eigen_decomp$values
	eigenvectors <- eigen_decomp$vectors
	
	# Adjust negative eigenvalues
	eigenvalues[eigenvalues < min_eigenvalue] <- min_eigenvalue
	
	# Reconstruct matrix
	matrix_new <- eigenvectors %*% diag(eigenvalues) %*% t(eigenvectors)
	
	# For correlation matrices: rescale to ensure unit diagonal
	diag_elements <- diag(matrix_new)
	if (all(diag_elements > 0)) {
		scaling_factors <- 1/sqrt(diag_elements)
		scaling_matrix <- diag(scaling_factors)
		matrix_new <- scaling_matrix %*% matrix_new %*% scaling_matrix
	}
	
	# Final cleanup
	matrix_new <- (matrix_new + t(matrix_new)) / 2
	diag(matrix_new) <- 1
	
	# Restore names if they existed
	if (!is.null(orig_rownames)) {
		rownames(matrix_new) <- orig_rownames
	}
	if (!is.null(orig_colnames)) {
		colnames(matrix_new) <- orig_colnames
	}
	
	# Ensure result is a matrix
	if (!is.matrix(matrix_new)) {
		stop("Function failed to return a matrix")
	}
	
	return(matrix_new)
}

# psd_validation.R
# Positive semidefinite matrix validation and correction functions
# Extracted from correlation_matrix_test.R

#' Test if matrix is positive semidefinite
#'
#' @param matrix Matrix to test
#' @param tolerance Numerical tolerance for eigenvalue tests
#' @return TRUE if matrix is positive semidefinite
is_positive_semidefinite <- function(matrix, tolerance = 1e-10) {
	eigen_values <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
	return(min(eigen_values) >= -tolerance)
}



#' Check matrix properties related to positive semidefiniteness
#'
#' @param matrix Matrix to analyze
#' @return List with PSD-related properties
analyze_psd_properties <- function(matrix) {
	eigen_values <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
	
	return(list(
		is_psd = min(eigen_values) >= -1e-10,
		min_eigenvalue = min(eigen_values),
		max_eigenvalue = max(eigen_values),
		n_negative_eigenvalues = sum(eigen_values < -1e-10),
		eigenvalue_range = max(eigen_values) - min(eigen_values),
		needs_adjustment = min(eigen_values) < -1e-10
	))
}