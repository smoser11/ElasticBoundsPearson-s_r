# 2_correlation_matrices/1_matrix_properties/condition_number_analysis.R
# Condition number analysis functions

#' Classify condition number for interpretation
classify_condition_number <- function(cond_num) {
	if (is.infinite(cond_num)) return("Singular")
	if (cond_num <= 10) return("Well-conditioned")
	if (cond_num <= 100) return("Moderately ill-conditioned")
	if (cond_num <= 1000) return("Ill-conditioned")
	return("Severely ill-conditioned")
}