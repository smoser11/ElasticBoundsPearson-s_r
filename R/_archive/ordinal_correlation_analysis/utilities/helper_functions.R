# utilities/helper_functions.R
# Common utility functions

#' Safe division avoiding division by zero
safe_divide <- function(numerator, denominator, default_value = NA) {
	ifelse(abs(denominator) < 1e-15, default_value, numerator / denominator)
}

#' Format numbers for reporting
format_number <- function(x, digits = 3) {
	if (is.na(x)) return("NA")
	return(formatC(round(x, digits), format = "f", digits = digits))
}