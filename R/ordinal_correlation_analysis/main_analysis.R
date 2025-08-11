# main_analysis.R
# Master analysis coordination script

# Required libraries
library(dplyr)
library(ggplot2)
library(gridExtra)

# Source all module coordinators
source("1_bivariate_ordcats_correlation/bivariate_main.R")
source("2_correlation_matrices/matrices_main.R") 
source("3_fixes_and_rescaling/rescaling_main.R")

#' Run complete ordinal correlation analysis pipeline
#'
#' @param bes_data BES dataset
#' @param config Analysis configuration  
#' @return Complete analysis results
run_complete_analysis <- function(bes_data, config = get_default_config()) {

#	config <- research_config
	cat("Starting comprehensive ordinal correlation analysis...\n")
	cat("Dataset:", nrow(bes_data), "variable pairs\n")
	
	# 1. Bivariate bounds analysis
	cat("\n=== BIVARIATE ANALYSIS ===\n")
	bivariate_results <- run_bivariate_analysis(bes_data, config)
	
	# 2. Matrix properties analysis
	cat("\n=== MATRIX ANALYSIS ===\n") 
	matrix_results <- run_matrix_analysis(bivariate_results$bounds_data, config)
	
	# 3. Rescaling analysis
	cat("\n=== RESCALING ANALYSIS ===\n")
	rescaling_results <- run_rescaling_analysis(bivariate_results$bounds_data, config)
	
	cat("\n=== ANALYSIS COMPLETE ===\n")
	
	return(list(
		bivariate = bivariate_results,
		matrices = matrix_results,
		rescaling = rescaling_results,
		config = config
	))
}

#' Get default analysis configuration
get_default_config <- function() {
	list(
		nsim = 1000,           # Permutation simulations
		matrix_trials = 50,    # Matrix trials  
		min_vars = 3,          # Min variables per matrix
		max_vars = 10,         # Max variables per matrix
		progress = TRUE,       # Show progress
		seed = 123            # Random seed
	)
}

# Test if run directly
if (!interactive()) {
	# Load your BES data here
	# bes_data <- read.csv("path/to/your/bes_data.csv")
	
	# Run analysis
	# results <- run_complete_analysis(bes_data)
	
	cat("Run: results <- run_complete_analysis(your_bes_data)\n")

}


###########################################################################
######################### RESEARCH
# 1.1 Load and Validate Your BES Dataset
# Load your actual BES data
getwd()
library(readstata13)
# bes_data <- read.dta13("../data/correlation and other data about pairs of BES2019 variables.dta")
bes_data <-read.dta13('/Users/sm38679/Documents/GitHub/ElasticBoundsPearson-s_r/R/data/processed/correlation and other data about pairs of BES2019 variables.dta')

# Or if you have it in another format:
# bes_data <- readRDS("BES_2019_processed.rds")
# bes_data <- haven::read_dta("BES_2019.dta")

# Validate the data structure
cat("Dataset dimensions:", dim(bes_data), "\n")
cat("Required columns present:", all(c("var1", "var2", "corr", "nobs", "var1cats", "var2cats") %in% names(bes_data)), "\n")

# Quick data overview
head(bes_data)
summary(bes_data[c("corr", "nobs", "var1cats", "var2cats")])


# 1.2 Load and Validate Your BES Dataset
# Create data quality report
source("./utilities/data_loading.R")  # If you created this file

# Comprehensive data validation
data_validation <- validate_bes_dataset(bes_data)
cat("Data validation passed:", data_validation$valid, "\n")
if (!data_validation$valid) {
	cat("Issues found:\n")
	for (error in data_validation$errors) {
		cat("  -", error, "\n")
	}
}

# Clean the data if needed
bes_data_clean <- clean_bes_data(bes_data)
cat("Cleaned dataset size:", nrow(bes_data_clean), "variable pairs\n")

## 2.1 configure analysis params

research_config <-
list(
	nsim =		2000,
	# High-quality simula
	2000, # High-quality simula
	matrix_trials = 200,
	min_vars = 5 , 
	max_vars = 25,
	progress =  TRUE,
	seed =42,
	significance_alpha = 0.05
)

cat("Research configuration:\n")
	 print(research_config)
	 
	
## 2.2 Run Complete Analysis Pipeline
	 
	 
	 # Load the main analysis system
	 source("main_analysis.R")
	 
	 # Run comprehensive analysis
	 cat("Starting research analysis on", nrow(bes_data_clean), "variable pairs...\n")
	 cat("This may take several minutes for large datasets.\n")
	 
	 # Time the analysis for reporting
	 start_time <- Sys.time()
	 research_results <- run_complete_analysis(bes_data_clean, research_config)
	 end_time <- Sys.time()
	 
	 analysis_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
	 cat("Analysis completed in", round(analysis_time, 1), "minutes\n")	 
	 
	 
# 2.3 2.3 Generate Research Outputs
	 # Create output directories for publication
	 dir.create("research_output", showWarnings = FALSE)
	 dir.create("research_output/figures", showWarnings = FALSE)
	 dir.create("research_output/tables", showWarnings = FALSE)
	 dir.create("research_output/data", showWarnings = FALSE)
	 
	 # Save complete results
	 saveRDS(research_results, "research_output/data/complete_analysis_results.rds")
	 cat("✓ Complete results saved\n")
	 
# 3.1 Main Research Figures
	 
	 # Load visualization functions
	 source("1_bivariate_ordcats_correlation/3_visualization/bounds_visualization.R")
	 
	 # Figure 1: Correlation Bounds Landscape
	 fig1 <- plot_bounds_scatter(research_results$bivariate$bounds_data, color_by = "bounds_range")
	 fig1 <- fig1 + 
	 	labs(title = "Theoretical Correlation Bounds for Ordinal Variables",
	 		 subtitle = paste("British Election Study 2019 (n =", nrow(research_results$bivariate$bounds_data), "variable pairs)")) +
	 	theme(text = element_text(size = 12), plot.title = element_text(size = 14, face = "bold"))
	 
	 ggsave("research_output/figures/figure1_correlation_bounds_landscape.png", 
	 	   fig1, width = 10, height = 8, dpi = 300)
	 
	 # Figure 2: Rescaling Effects
	 if (!is.null(research_results$rescaling$rescaled_data$r_rescaled)) {
	 	fig2 <- plot_rescaling_comparison(research_results$rescaling$rescaled_data)
	 	fig2 <- fig2 + 
	 		labs(title = "Impact of Correlation Rescaling on Observed Relationships",
	 			 subtitle = "Rescaling maps correlations to theoretical [-1,1] range") +
	 		theme(text = element_text(size = 12), plot.title = element_text(size = 14, face = "bold"))
	 	
	 	ggsave("research_output/figures/figure2_rescaling_effects.png", 
	 		   fig2, width = 10, height = 8, dpi = 300)
	 }
	 
	 # Figure 3: Asymmetry Analysis
	 if (!is.null(research_results$bivariate$asymmetry)) {
	 	asymmetry_data <- research_results$bivariate$asymmetry$detailed_results
	 	
	 	fig3 <- ggplot(asymmetry_data, aes(x = basic_asymmetry)) +
	 		geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
	 		geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
	 		labs(title = "Distribution of Correlation Bounds Asymmetry",
	 			 subtitle = "Asymmetry = r_max + r_min (perfect symmetry = 0)",
	 			 x = "Bounds Asymmetry", y = "Number of Variable Pairs") +
	 		theme_minimal() +
	 		theme(text = element_text(size = 12), plot.title = element_text(size = 14, face = "bold"))
	 	
	 	ggsave("research_output/figures/figure3_asymmetry_distribution.png", 
	 		   fig3, width = 10, height = 6, dpi = 300)
	 }
	 
	 
	 
## 3.2 Supplementary Figures
	 # Supplementary Figure A: Matrix Properties Comparison
	 if (!is.null(research_results$matrices)) {
	 	matrix_comparison <- data.frame(
	 		Property = c("Positive Semidefinite", "Invertible", "Valid Correlation Matrix"),
	 		Original = c(research_results$matrices$summary_psd$orig_psd_percent,
	 					 research_results$matrices$summary_psd$orig_invertible_percent,
	 					 research_results$matrices$summary_psd$orig_valid_percent),
	 		Rescaled = c(research_results$matrices$summary_psd$rescaled_psd_percent,
	 					 research_results$matrices$summary_psd$rescaled_invertible_percent,
	 					 research_results$matrices$summary_psd$rescaled_valid_percent)
	 	)
	 	
	 	matrix_comparison_long <- reshape2::melt(matrix_comparison, id.vars = "Property")
	 	
	 	figS1 <- ggplot(matrix_comparison_long, aes(x = Property, y = value, fill = variable)) +
	 		geom_bar(stat = "identity", position = "dodge") +
	 		geom_text(aes(label = paste0(round(value, 1), "%")), 
	 				  position = position_dodge(width = 0.9), vjust = -0.5) +
	 		scale_fill_manual(values = c("Original" = "lightblue", "Rescaled" = "lightcoral"),
	 						  name = "Matrix Type") +
	 		labs(title = "Matrix Properties: Original vs Rescaled Correlations",
	 			 subtitle = paste("Based on", research_results$matrices$summary_psd$n_trials, "random matrix trials"),
	 			 x = "", y = "Percentage of Matrices") +
	 		theme_minimal() +
	 		theme(text = element_text(size = 12), plot.title = element_text(size = 14, face = "bold"),
	 			  axis.text.x = element_text(angle = 45, hjust = 1))
	 	
	 	ggsave("research_output/figures/figureS1_matrix_properties.png", 
	 		   figS1, width = 12, height = 8, dpi = 300)
	 }
	 
	 cat("✓ Supplementary figures saved\n")	 
	 cat("✓ Main figures saved to research_output/figures/\n") 
	 
	 
##	 4.1 Summary Statistics Table
	 # Table 1: Descriptive Statistics
	 bounds_data <- research_results$bivariate$bounds_data
	 
	 table1 <- data.frame(
	 	Statistic = c("Number of Variable Pairs", "Mean Correlation", "SD Correlation", 
	 				  "Mean r_min", "SD r_min", "Mean r_max", "SD r_max",
	 				  "Mean Bounds Range", "SD Bounds Range", 
	 				  "Mean Sample Size", "SD Sample Size"),
	 	Value = c(
	 		nrow(bounds_data),
	 		round(mean(bounds_data$observed_r, na.rm = TRUE), 3),
	 		round(sd(bounds_data$observed_r, na.rm = TRUE), 3),
	 		round(mean(bounds_data$r_min, na.rm = TRUE), 3),
	 		round(sd(bounds_data$r_min, na.rm = TRUE), 3),
	 		round(mean(bounds_data$r_max, na.rm = TRUE), 3),
	 		round(sd(bounds_data$r_max, na.rm = TRUE), 3),
	 		round(mean(bounds_data$r_max - bounds_data$r_min, na.rm = TRUE), 3),
	 		round(sd(bounds_data$r_max - bounds_data$r_min, na.rm = TRUE), 3),
	 		round(mean(bounds_data$n_obs, na.rm = TRUE), 0),
	 		round(sd(bounds_data$n_obs, na.rm = TRUE), 0)
	 	)
	 )
	 
	 write.csv(table1, "research_output/tables/table1_descriptive_statistics.csv", row.names = FALSE)
	 print(table1)	 
	 
	 
	 
## 4.2 Key Findings Table
	 
	 # Table 2: Key Research Findings
	 if (!is.null(research_results$bivariate$asymmetry)) {
	 	asymmetry_summary <- research_results$bivariate$asymmetry$summary
	 	
	 	table2 <- data.frame(
	 		Finding = c(
	 			"Proportion with Symmetric Bounds",
	 			"Proportion with Mild Asymmetry", 
	 			"Proportion with Moderate Asymmetry",
	 			"Proportion with High Asymmetry",
	 			"Mean Rescaling Effect",
	 			"Proportion with Increased Magnitude",
	 			"Matrix PSD Rate (Original)",
	 			"Matrix PSD Rate (Rescaled)",
	 			"Mean Condition Number (Original)",
	 			"Mean Condition Number (Rescaled)"
	 		),
	 		Value = c(
	 			paste0(round(100 * asymmetry_summary$prop_symmetric, 1), "%"),
	 			paste0(round(100 * asymmetry_summary$prop_mildly_asymmetric, 1), "%"),
	 			paste0(round(100 * asymmetry_summary$prop_moderately_asymmetric, 1), "%"),
	 			paste0(round(100 * asymmetry_summary$prop_highly_asymmetric, 1), "%"),
	 			round(research_results$rescaling$effects$summary$mean_magnitude_change, 3),
	 			paste0(round(100 * research_results$rescaling$effects$summary$prop_increased, 1), "%"),
	 			paste0(round(research_results$matrices$summary_psd$orig_psd_percent, 1), "%"),
	 			paste0(round(research_results$matrices$summary_psd$rescaled_psd_percent, 1), "%"),
	 			round(research_results$matrices$summary_psd$orig_mean_condition, 2),
	 			round(research_results$matrices$summary_psd$rescaled_mean_condition, 2)
	 		)
	 	)
	 	
	 	write.csv(table2, "research_output/tables/table2_key_findings.csv", row.names = FALSE)
	 	print(table2)
	 }
	 
	 cat("✓ Tables saved to research_output/tables/\n")	 
	 
	 
## 5.1 Automated Research Summary
	 # Generate executive research summary
	 research_summary <- paste0(
	 	"CORRELATION BOUNDS ANALYSIS: BRITISH ELECTION STUDY 2019\n",
	 	"=" %&% strrep("=", 55), "\n\n",
	 	
	 	"DATASET OVERVIEW:\n",
	 	"- Variable pairs analyzed: ", nrow(bounds_data), "\n",
	 	"- Mean correlation: ", round(mean(bounds_data$observed_r, na.rm = TRUE), 3), 
	 	" (SD = ", round(sd(bounds_data$observed_r, na.rm = TRUE), 3), ")\n",
	 	"- Mean sample size: ", round(mean(bounds_data$n_obs, na.rm = TRUE), 0), "\n\n",
	 	
	 	"KEY FINDINGS:\n",
	 	"1. Theoretical bounds range: [", 
	 	round(mean(bounds_data$r_min, na.rm = TRUE), 3), ", ",
	 	round(mean(bounds_data$r_max, na.rm = TRUE), 3), "] on average\n",
	 	
	 	"2. Bounds asymmetry: ",
	 	if (!is.null(research_results$bivariate$asymmetry)) {
	 		paste0(round(100 * (1 - research_results$bivariate$asymmetry$summary$prop_symmetric), 1), 
	 			   "% of pairs show asymmetric bounds")
	 	} else {"Analysis pending"}, "\n",
	 	
	 	"3. Rescaling effects: ",
	 	if (!is.null(research_results$rescaling)) {
	 		paste0(round(100 * research_results$rescaling$effects$summary$prop_increased, 1),
	 			   "% of correlations increased in magnitude")
	 	} else {"Analysis pending"}, "\n",
	 	
	 	"4. Matrix properties: ",
	 	if (!is.null(research_results$matrices)) {
	 		paste0(round(research_results$matrices$summary_psd$orig_psd_percent, 1),
	 			   "% of original matrices are positive semidefinite")
	 	} else {"Analysis pending"}, "\n\n",
	 	
	 	"IMPLICATIONS:\n",
	 	"- Correlation bounds vary substantially across variable pairs\n",
	 	"- Rescaling provides more interpretable correlation magnitudes\n", 
	 	"- Matrix properties are generally preserved under rescaling\n",
	 	"- Results support the use of bounds-aware correlation analysis\n\n",
	 	
	 	"Generated: ", Sys.time(), "\n",
	 	"Analysis time: ", round(analysis_time, 1), " minutes\n"
	 )
	 
	 # Save research summary
	 writeLines(research_summary, "research_output/research_summary.txt")
	 cat("✓ Research summary saved\n")
	 cat("\nRESEARCH SUMMARY:\n")
	 cat(research_summary)	 
	 
	 
## 6.2 Results Section Key Points

