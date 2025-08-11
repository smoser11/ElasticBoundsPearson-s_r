# utilities/data_loading.R
# Functions for loading and preprocessing BES data

#' Load BES data with proper formatting and validation
#'
#' @param file_path Path to BES data file
#' @param validate_data Perform data validation
#' @return Cleaned BES data frame
load_bes_data <- function(file_path, validate_data = TRUE) {
	if (!file.exists(file_path)) {
		stop("BES data file not found: ", file_path)
	}
	
	cat("Loading BES data from:", file_path, "\n")
	
	# Load data (adjust based on your file format)
	if (grepl("\\.csv$", file_path)) {
		bes_data <- read.csv(file_path, stringsAsFactors = FALSE)
	} else if (grepl("\\.rds$", file_path)) {
		bes_data <- readRDS(file_path)
	} else if (grepl("\\.dta$", file_path)) {
		bes_data <- haven::read_dta(file_path)
	} else {
		stop("Unsupported file format. Use .csv, .rds, or .dta")
	}
	
	cat("Loaded", nrow(bes_data), "rows and", ncol(bes_data), "columns\n")
	
	# Basic data cleaning
	bes_data <- clean_bes_data(bes_data)
	
	# Validation if requested
	if (validate_data) {
		validation_results <- validate_bes_dataset(bes_data)
		if (!validation_results$valid) {
			warning("Data validation issues found:\n", paste(validation_results$errors, collapse = "\n"))
		} else {
			cat("✓ Data validation passed\n")
		}
	}
	
	return(bes_data)
}

#' Clean and standardize BES data
#'
#' @param bes_data Raw BES data frame
#' @return Cleaned data frame
clean_bes_data <- function(bes_data) {
	# Ensure required columns exist
	required_cols <- c("var1", "var2", "corr", "nobs", "var1cats", "var2cats")
	missing_cols <- setdiff(required_cols, names(bes_data))
	
	if (length(missing_cols) > 0) {
		stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
	}
	
	# Convert variable identifiers to character
	bes_data$var1 <- as.character(bes_data$var1)
	bes_data$var2 <- as.character(bes_data$var2)
	
	# Ensure numeric columns are properly formatted
	numeric_cols <- c("corr", "nobs", "var1cats", "var2cats")
	for (col in numeric_cols) {
		if (col %in% names(bes_data)) {
			bes_data[[col]] <- as.numeric(bes_data[[col]])
		}
	}
	
	# Remove rows with invalid data
	valid_rows <- !is.na(bes_data$corr) & 
		!is.na(bes_data$nobs) & 
		bes_data$nobs > 0 &
		abs(bes_data$corr) <= 1
	
	if (sum(!valid_rows) > 0) {
		cat("Removing", sum(!valid_rows), "rows with invalid data\n")
		bes_data <- bes_data[valid_rows, ]
	}
	
	# Sort by variable pairs for consistent ordering
	bes_data <- bes_data[order(bes_data$var1, bes_data$var2), ]
	
	return(bes_data)
}

#' Validate entire BES dataset
#'
#' @param bes_data BES data frame
#' @return Validation results
validate_bes_dataset <- function(bes_data) {
	errors <- c()
	warnings <- c()
	
	# Check data dimensions
	if (nrow(bes_data) == 0) {
		errors <- c(errors, "Dataset is empty")
	}
	
	# Check for required columns
	required_cols <- c("var1", "var2", "corr", "nobs", "var1cats", "var2cats")
	missing_cols <- setdiff(required_cols, names(bes_data))
	if (length(missing_cols) > 0) {
		errors <- c(errors, paste("Missing columns:", paste(missing_cols, collapse = ", ")))
	}
	
	# Validate correlations
	invalid_corr <- sum(abs(bes_data$corr) > 1, na.rm = TRUE)
	if (invalid_corr > 0) {
		errors <- c(errors, paste(invalid_corr, "correlations outside [-1, 1] range"))
	}
	
	# Check sample sizes
	small_samples <- sum(bes_data$nobs < 30, na.rm = TRUE)
	if (small_samples > 0) {
		warnings <- c(warnings, paste(small_samples, "pairs have small sample sizes (< 30)"))
	}
	
	# Check frequency/proportion columns
	freq_cols <- grep("freq[0-9]+", names(bes_data), value = TRUE)
	prop_cols <- grep("prop[0-9]+", names(bes_data), value = TRUE)
	
	if (length(freq_cols) == 0 && length(prop_cols) == 0) {
		warnings <- c(warnings, "No frequency or proportion columns found")
	}
	
	# Check for duplicate variable pairs
	if (anyDuplicated(paste(bes_data$var1, bes_data$var2))) {
		errors <- c(errors, "Duplicate variable pairs found")
	}
	
	return(list(
		valid = length(errors) == 0,
		errors = errors,
		warnings = warnings,
		n_rows = nrow(bes_data),
		n_valid_pairs = sum(!is.na(bes_data$corr))
	))
}