
## take cross-tabes -- simluated or BES or whatever, and bootstrap the 2-D uncertainty   of $(\hat{r_{min}}, \hat{r_{max}}$)






# Build FH upper (comonotonic) joint from marginals pX, pY (both sum to 1)
fh_upper <- function(pX, pY){
	cx <- cumsum(pX); cy <- cumsum(pY)
	i <- j <- 1; px <- py <- 0
	J <- matrix(0, length(pX), length(pY))
	while(i <= length(pX) && j <= length(pY)){
		nx <- cx[i]; ny <- cy[j]
		J[i,j] <- max(0, min(nx, ny) - max(px, py))
		if(nx < ny){ px <- nx; i <- i + 1 } 
		else if(nx > ny){ py <- ny; j <- j + 1 }
		else { px <- nx; py <- ny; i <- i + 1; j <- j + 1 }
	}
	J
}

# FH "lower" (antitone) by pairing X with reversed Y, then un-reversing columns
fh_lower <- function(pX, pY){
	Jrev <- fh_upper(pX, rev(pY))
	Jrev[, ncol(Jrev):1, drop=FALSE]
}

# Pearson r from a joint table and score vectors
r_from_joint <- function(J, x_scores, y_scores){
	pX <- rowSums(J); pY <- colSums(J)
	EX <- sum(pX * x_scores); EY <- sum(pY * y_scores)
	EXY <- sum(J * outer(x_scores, y_scores))
	cov <- EXY - EX * EY
	sdX <- sqrt(sum(pX * x_scores^2) - EX^2)
	sdY <- sqrt(sum(pY * y_scores^2) - EY^2)
	cov / (sdX * sdY)
}

# One-shot estimator from data vectors x,y (ordinal categories 1..K etc.)
extreme_rs <- function(x, y, x_scores = NULL, y_scores = NULL){
	if(is.null(x_scores)) x_scores <- sort(unique(x))
	if(is.null(y_scores)) y_scores <- sort(unique(y))
	# align scores to contiguous indices
	x_levels <- x_scores; y_levels <- y_scores
	pX <- prop.table(table(factor(x, levels = x_levels)))
	pY <- prop.table(table(factor(y, levels = y_levels)))
	Jmax <- fh_upper(as.numeric(pX), as.numeric(pY))
	Jmin <- fh_lower(as.numeric(pX), as.numeric(pY))
	r_max <- r_from_joint(Jmax, x_levels, y_levels)
	r_min <- r_from_joint(Jmin, x_levels, y_levels)
	c(r_min = r_min, r_max = r_max)
}

# Bootstrap wrapper
boot_extremes <- function(x, y, B = 1000, seed = 1){
	set.seed(seed)
	est0 <- extreme_rs(x, y)
	n <- length(x); out <- matrix(NA_real_, B, 2, dimnames=list(NULL, names(est0)))
	for(b in 1:B){
		idx <- sample.int(n, n, replace = TRUE)
		out[b,] <- extreme_rs(x[idx], y[idx])
	}
	list(
		est = est0,
		boot = out,
		ci = apply(out, 2, quantile, probs = c(.025, .975))
	)
}


## WRAPING THE ABOVE:
# Simplified bootstrap analysis using your existing functions

# Simplified bootstrap analysis using your existing functions

# Convert contingency table matrix to individual observation vectors
# Simplified bootstrap analysis using your existing functions

# Simplified bootstrap analysis using your existing functions

# Simplified bootstrap analysis using your existing functions

# Simplified bootstrap analysis using your existing functions
# Simplified bootstrap analysis using your existing functions

# Simplified bootstrap analysis using your existing functions

# Convert contingency table matrix to individual observation vectors
contingency_to_vectors <- function(table_matrix) {
	K1 <- nrow(table_matrix)
	K2 <- ncol(table_matrix)
	
	# Create vectors to hold individual observations
	x_vec <- c()
	y_vec <- c()
	
	# Convert table back to individual observations
	for(i in 1:K1) {
		for(j in 1:K2) {
			count <- table_matrix[i, j]
			if(count > 0) {
				x_vec <- c(x_vec, rep(i, count))
				y_vec <- c(y_vec, rep(j, count))
			}
		}
	}
	
	return(list(x = x_vec, y = y_vec))
}

# Check if a table is suitable for bootstrap analysis
is_table_valid_for_bootstrap <- function(table_matrix) {
	# Check 1: Must have observations in at least 2 rows and 2 columns
	row_sums <- rowSums(table_matrix)
	col_sums <- colSums(table_matrix)
	
	non_empty_rows <- sum(row_sums > 0)
	non_empty_cols <- sum(col_sums > 0)
	
	if(non_empty_rows < 2 || non_empty_cols < 2) {
		return(list(valid = FALSE, 
					reason = paste("Only", non_empty_rows, "non-empty rows and", 
								   non_empty_cols, "non-empty columns")))
	}
	
	# Check 2: Must have at least some minimum number of observations
	total_obs <- sum(table_matrix)
	if(total_obs < 10) {
		return(list(valid = FALSE, 
					reason = paste("Too few observations:", total_obs)))
	}
	
	# Check 3: Convert to vectors and check for variance
	vectors <- contingency_to_vectors(table_matrix)
	
	if(length(unique(vectors$x)) < 2 || length(unique(vectors$y)) < 2) {
		return(list(valid = FALSE, 
					reason = "No variance in x or y"))
	}
	
	# Check 4: Check for extreme sparsity - need sufficient observations in non-empty cells
	non_empty_cells <- sum(table_matrix > 0)
	if(non_empty_cells < 3) {
		return(list(valid = FALSE, 
					reason = paste("Too sparse: only", non_empty_cells, "non-empty cells")))
	}
	
	# Check 5: Check if largest cell dominates too much (>95% of observations)
	max_cell_prop <- max(table_matrix) / total_obs
	if(max_cell_prop > 0.95) {
		return(list(valid = FALSE, 
					reason = paste("One cell dominates:", round(max_cell_prop * 100, 1), "% of observations")))
	}
	
	# Check 6: Try a test run of extreme_rs to see if it works
	test_result <- tryCatch({
		extreme_rs(vectors$x, vectors$y)
	}, error = function(e) {
		return(list(error = e$message))
	})
	
	if(is.list(test_result) && !is.null(test_result$error)) {
		return(list(valid = FALSE, 
					reason = paste("extreme_rs failed:", test_result$error)))
	}
	
	return(list(valid = TRUE, reason = "OK"))
}
# Apply bootstrap to a single configuration
# B = number of bootstrap samples (matches your boot_extremes function)
bootstrap_configuration <- function(sim_results, config_id, B = 1000, 
									seed_offset = 1000, store_full_samples = FALSE) {
	
	if(config_id > length(sim_results$tables)) {
		stop("config_id out of range")
	}
	
	config_tables <- sim_results$tables[[config_id]]
	config_info <- sim_results$config_info[config_id, ]
	
	cat(sprintf("Bootstrapping config %d: K1=%d, K2=%d, N=%d\n", 
				config_id, config_info$K1, config_info$K2, config_info$N))
	cat("Number of tables:", length(config_tables), "\n")
	
	
	# Storage for summary results
	summary_stats <- data.frame()
	full_bootstrap_samples <- list()  # NEW: Store full samples if requested
	skipped_tables <- 0
	
	# Process each table
	for(sim_id in 1:length(config_tables)) {
		if(sim_id %% 50 == 0) {
			cat("Processing table", sim_id, "of", length(config_tables), "\n")
		}
		
		# Check if table is suitable for bootstrap
		table_check <- is_table_valid_for_bootstrap(config_tables[[sim_id]])
		
		if(!table_check$valid) {
			cat("Skipping table", sim_id, ":", table_check$reason, "\n")
			skipped_tables <- skipped_tables + 1
			next
		}
		
		# Convert table to vectors
		vectors <- contingency_to_vectors(config_tables[[sim_id]])
		
		# Set unique seed for reproducibility
		table_seed <- seed_offset + config_id * 10000 + sim_id
		
		# Use your existing boot_extremes function!
		boot_result <- tryCatch({
			# Suppress warnings and messages during bootstrap
			suppressWarnings(suppressMessages(
				boot_extremes(vectors$x, vectors$y, B = B, seed = table_seed)
			))
		}, error = function(e) {
			# Silently skip problematic tables - just increment skip counter
			return(NULL)
		})
		
		# Skip this table if bootstrap failed
		if(is.null(boot_result)) {
			skipped_tables <- skipped_tables + 1
			next
		}
		
		# NEW: Store full bootstrap samples if requested
		if(store_full_samples) {
			full_bootstrap_samples[[length(full_bootstrap_samples) + 1]] <- list(
				config_id = config_id,
				sim_id = sim_id,
				bootstrap_samples = boot_result$boot,  # The full B x 2 matrix
				original_estimates = boot_result$est
			)
		}
		
		# Extract summary statistics INCLUDING joint relationship
		
		# Extract summary statistics INCLUDING joint relationship
		boot_cor <- cor(boot_result$boot[,"r_min"], boot_result$boot[,"r_max"])
		boot_cov <- cov(boot_result$boot[,"r_min"], boot_result$boot[,"r_max"])
		
		summary_stats <- rbind(summary_stats, data.frame(
			config_id = config_id,
			sim_id = sim_id,
			K1 = config_info$K1,
			K2 = config_info$K2,
			N = config_info$N,
			r_min_est = boot_result$est["r_min"],
			r_max_est = boot_result$est["r_max"],
			r_range = boot_result$est["r_max"] - boot_result$est["r_min"],
			r_min_ci_lower = boot_result$ci["2.5%", "r_min"],
			r_min_ci_upper = boot_result$ci["97.5%", "r_min"],
			r_max_ci_lower = boot_result$ci["2.5%", "r_max"],
			r_max_ci_upper = boot_result$ci["97.5%", "r_max"],
			r_min_ci_width = boot_result$ci["97.5%", "r_min"] - boot_result$ci["2.5%", "r_min"],
			r_max_ci_width = boot_result$ci["97.5%", "r_max"] - boot_result$ci["2.5%", "r_max"],
			# NEW: Joint relationship statistics
			boot_cor_rmin_rmax = boot_cor,
			boot_cov_rmin_rmax = boot_cov,
			boot_sd_rmin = sd(boot_result$boot[,"r_min"]),
			boot_sd_rmax = sd(boot_result$boot[,"r_max"]),
			boot_mean_rmin = mean(boot_result$boot[,"r_min"]),
			boot_mean_rmax = mean(boot_result$boot[,"r_max"])
		))
	}
	
	cat("Bootstrap complete for configuration", config_id, "\n")
	cat("Successfully processed:", nrow(summary_stats), "tables\n")
	cat("Skipped:", skipped_tables, "tables\n\n")
	
	result <- list(
		config_info = config_info,
		summary_stats = summary_stats,
		skipped_count = skipped_tables,
		success_count = nrow(summary_stats)
	)
	
	# NEW: Add full samples if stored
	if(store_full_samples) {
		result$full_bootstrap_samples <- full_bootstrap_samples
		cat("Stored full bootstrap samples for", length(full_bootstrap_samples), "tables\n")
	}
	
	return(result)
}

# Apply bootstrap to ALL configurations
bootstrap_all_configurations <- function(sim_results, B = 1000, 
										 seed_offset = 1000, store_full_samples = FALSE) {
	
	num_configs <- length(sim_results$tables)
	all_summaries <- data.frame()
	
	cat("Starting bootstrap analysis for", num_configs, "configurations\n")
	cat("Bootstrap samples per table:", B, "\n\n")
	cat("Bootstrap samples per table:", B, "\n")
	if(store_full_samples) {
		cat("WARNING: Storing full samples will create very large files!\n")
	}
	cat("\n")
	
	start_time <- Sys.time()
	
	for(config_id in 1:num_configs) {
		config_result <- bootstrap_configuration(sim_results, config_id, 
												 B = B, seed_offset = seed_offset,
												 store_full_samples = store_full_samples)
		# Combine summary statistics
		all_summaries <- rbind(all_summaries, config_result$summary_stats)
		
		# Progress update
		elapsed <- as.numeric(Sys.time() - start_time, units = "mins")
		remaining <- elapsed * (num_configs - config_id) / config_id
		cat(sprintf("Completed %d/%d configs. Elapsed: %.1f min, Est. remaining: %.1f min\n", 
					config_id, num_configs, elapsed, remaining))
	}
	
	end_time <- Sys.time()
	total_time <- as.numeric(end_time - start_time, units = "mins")
	
	cat(sprintf("\nBootstrap analysis complete! Total time: %.1f minutes\n", total_time))
	cat("Total bootstrap samples:", nrow(all_summaries) * B, "\n")
	
	return(list(
		config_info = sim_results$config_info,
		summary_stats = all_summaries,
		analysis_info = list(
			B = B,
			seed_offset = seed_offset,
			total_time_mins = total_time,
			store_full_samples = store_full_samples  # ← This line should be there
			
		)
	))
}

# Analyze a single table WITH bootstrap uncertainty
analyze_single_table_with_bootstrap <- function(table_matrix, B = 1000, seed = 123) {
	vectors <- contingency_to_vectors(table_matrix)
	
	# Get both point estimates and bootstrap uncertainty
	boot_result <- boot_extremes(vectors$x, vectors$y, B = B, seed = seed)
	
	cat("Table dimensions:", nrow(table_matrix), "x", ncol(table_matrix), "\n")
	cat("Sample size:", sum(table_matrix), "\n")
	cat("r_min estimate:", round(boot_result$est["r_min"], 4), "\n")
	cat("r_max estimate:", round(boot_result$est["r_max"], 4), "\n")
	cat("Range:", round(boot_result$est["r_max"] - boot_result$est["r_min"], 4), "\n")
	cat("\nBootstrap Confidence Intervals (95%):\n")
	cat("r_min CI: [", round(boot_result$ci["2.5%", "r_min"], 4), ",", 
		round(boot_result$ci["97.5%", "r_min"], 4), "]\n")
	cat("r_max CI: [", round(boot_result$ci["2.5%", "r_max"], 4), ",", 
		round(boot_result$ci["97.5%", "r_max"], 4), "]\n")
	
	return(boot_result)
}

# Quick analysis of a single table (no bootstrap, just extremes)
analyze_single_table <- function(table_matrix) {
	vectors <- contingency_to_vectors(table_matrix)
	extremes <- extreme_rs(vectors$x, vectors$y)
	
	cat("Table dimensions:", nrow(table_matrix), "x", ncol(table_matrix), "\n")
	cat("Sample size:", sum(table_matrix), "\n")
	cat("r_min:", round(extremes["r_min"], 4), "\n")
	cat("r_max:", round(extremes["r_max"], 4), "\n")
	cat("Range:", round(extremes["r_max"] - extremes["r_min"], 4), "\n")
	
	return(extremes)
}

# Save bootstrap results
save_bootstrap_results <- function(boot_analysis, 
								   filename = "bootstrap_analysis.rds") {
	saveRDS(boot_analysis, file = filename)
	cat("Bootstrap analysis saved to:", filename, "\n")
	cat("File size:", round(file.size(filename) / 1024^2, 2), "MB\n")
}

# Function to extract and analyze joint bootstrap samples
extract_joint_bootstrap_samples <- function(boot_analysis, config_id, sim_id) {
	# Get the specific bootstrap result
	if(config_id > length(boot_analysis$boot_results)) {
		stop("config_id out of range")
	}
	if(sim_id > length(boot_analysis$boot_results[[config_id]])) {
		stop("sim_id out of range")
	}
	
	# This would require storing the full boot_results, which we're not doing currently
	# But we can add this functionality if needed
	cat("To extract full joint samples, need to modify bootstrap functions to store boot_results\n")
	cat("Currently only storing summary statistics\n")
}

# Function to analyze joint relationships across configurations
analyze_joint_relationships <- function(boot_analysis) {
	stats <- boot_analysis$summary_stats
	
	cat("Joint Relationship Analysis\n")
	cat("==========================\n")
	
	# Overall statistics
	cat("Overall correlation between r_min and r_max bootstrap samples:\n")
	cat("Mean correlation:", round(mean(stats$boot_cor_rmin_rmax, na.rm = TRUE), 4), "\n")
	cat("SD correlation:", round(sd(stats$boot_cor_rmin_rmax, na.rm = TRUE), 4), "\n")
	cat("Range:", round(range(stats$boot_cor_rmin_rmax, na.rm = TRUE), 4), "\n\n")
	
	# By configuration
	config_joint_summary <- aggregate(
		cbind(boot_cor_rmin_rmax, boot_cov_rmin_rmax, r_range) ~ K1 + K2 + N, 
		data = stats, 
		FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
							sd = sd(x, na.rm = TRUE),
							min = min(x, na.rm = TRUE),
							max = max(x, na.rm = TRUE))
	)
	
	cat("Configuration summary (mean correlation, covariance, range):\n")
	print(config_joint_summary)
	
	# Find configurations with most/least dependence
	mean_cors <- aggregate(boot_cor_rmin_rmax ~ K1 + K2 + N, data = stats, mean, na.rm = TRUE)
	
	most_dependent <- mean_cors[which.max(mean_cors$boot_cor_rmin_rmax), ]
	least_dependent <- mean_cors[which.min(mean_cors$boot_cor_rmin_rmax), ]
	
	cat("\nMost dependent r_min/r_max (highest correlation):\n")
	cat("K1 =", most_dependent$K1, ", K2 =", most_dependent$K2, 
		", N =", most_dependent$N, ", Cor =", round(most_dependent$boot_cor_rmin_rmax, 4), "\n")
	
	cat("Least dependent r_min/r_max (lowest correlation):\n") 
	cat("K1 =", least_dependent$K1, ", K2 =", least_dependent$K2,
		", N =", least_dependent$N, ", Cor =", round(least_dependent$boot_cor_rmin_rmax, 4), "\n")
	
	return(list(
		overall_stats = list(
			mean_cor = mean(stats$boot_cor_rmin_rmax, na.rm = TRUE),
			sd_cor = sd(stats$boot_cor_rmin_rmax, na.rm = TRUE),
			range_cor = range(stats$boot_cor_rmin_rmax, na.rm = TRUE)
		),
		config_summary = config_joint_summary,
		most_dependent = most_dependent,
		least_dependent = least_dependent
	))
}
# stats <- boot_analysis$summary_stats
# 
# cat("Bootstrap Analysis Summary\n")
# cat("========================\n")
# cat("Total configurations:", length(unique(stats$config_id)), "\n")
# cat("Total tables analyzed:", nrow(stats), "\n")
# cat("Bootstrap samples per table:", boot_analysis$analysis_info$B, "\n")
# cat("Analysis time:", round(boot_analysis$analysis_info$total_time_mins, 1), "minutes\n\n")
# 
# # Summary by configuration
# config_summary <- aggregate(cbind(r_min_est, r_max_est, r_range) ~ K1 + K2 + N, 
# 							data = stats, FUN = function(x) c(mean = mean(x), sd = sd(x)))
# 
# print(config_summary)
# 
# return(config_summary)
# }



# Example usage:
# 

# 
# Quick check of a single table (no bootstrap)
single_table <- sim_data$tables[[1]][[1]]
analyze_single_table(single_table)

# Single table WITH bootstrap uncertainty (B = number of bootstrap samples)
analyze_single_table_with_bootstrap(single_table, B = 500)

# Bootstrap analysis for ALL configurations (now includes joint relationships!)
full_analysis <- bootstrap_all_configurations(sim_data, B = 10, store_full_samples = TRUE)

# Analyze the joint r_min/r_max relationships
joint_analysis <- analyze_joint_relationships(full_analysis)

# Check which configurations have most/least correlated r_min and r_max
print(joint_analysis$most_dependent)
print(joint_analysis$least_dependent)

# Save and summarize
library(here)
save_bootstrap_results(full_analysis, here("R", "ordinal_correlation_analysis", "data", "processed", "bootstrap_results_sim.rds"))
summary <- summarize_bootstrap_results(full_analysis)
summary


#####################################
# VISUALIZE dependence


# Visualize joint uncertainty distributions from bootstrap results

# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)

# Visualize joint uncertainty distributions from bootstrap results

# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)

# Function to visualize joint distributions using summary statistics
visualize_joint_uncertainty_summary <- function(boot_results_file) {
	
	# Load the bootstrap results
	if(is.character(boot_results_file)) {
		boot_analysis <- readRDS(boot_results_file)
	} else {
		boot_analysis <- boot_results_file
	}
	
	stats <- boot_analysis$summary_stats
	
	cat("Visualizing joint uncertainty for", nrow(stats), "tables across", 
		length(unique(stats$config_id)), "configurations\n")
	
	# 1. Correlation heatmap by configuration
	config_cors <- aggregate(boot_cor_rmin_rmax ~ K1 + K2 + N, 
							 data = stats, mean, na.rm = TRUE)
	
	p1 <- ggplot(config_cors, aes(x = factor(K1), y = factor(K2), fill = boot_cor_rmin_rmax)) +
		geom_tile() +
		facet_wrap(~paste("N =", N), ncol = 3) +
		scale_fill_viridis_c(name = "Mean\nCorrelation") +
		labs(title = "Mean r_min/r_max Correlation by Configuration",
			 x = "K1 (categories var1)", y = "K2 (categories var2)") +
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1))
	
	# 2. Scatter plot of r_min vs r_max estimates colored by correlation
	p2 <- ggplot(stats, aes(x = r_min_est, y = r_max_est, color = boot_cor_rmin_rmax)) +
		geom_point(alpha = 0.6, size = 1) +
		geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
		scale_color_viridis_c(name = "Bootstrap\nCorrelation") +
		labs(title = "r_min vs r_max Estimates",
			 subtitle = "Points colored by bootstrap correlation",
			 x = "r_min estimate", y = "r_max estimate") +
		theme_minimal()
	
	# 3. Distribution of correlations
	p3 <- ggplot(stats, aes(x = boot_cor_rmin_rmax)) +
		geom_histogram(bins = 30, alpha = 0.7, fill = "steelblue") +
		geom_vline(xintercept = mean(stats$boot_cor_rmin_rmax, na.rm = TRUE), 
				   color = "red", linetype = "dashed") +
		labs(title = "Distribution of r_min/r_max Bootstrap Correlations",
			 x = "Bootstrap Correlation", y = "Count") +
		theme_minimal()
	
	# 4. Range vs correlation
	p4 <- ggplot(stats, aes(x = r_range, y = boot_cor_rmin_rmax, color = factor(N))) +
		geom_point(alpha = 0.6) +
		scale_color_viridis_d(name = "Sample\nSize") +
		labs(title = "Range vs Bootstrap Correlation",
			 x = "r_max - r_min", y = "Bootstrap Correlation") +
		theme_minimal()
	
	# Combine plots
	combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2)
	
	return(list(
		correlation_heatmap = p1,
		scatter_plot = p2,
		correlation_dist = p3,
		range_vs_correlation = p4,
		combined = combined_plot
	))
}

# Unified function to create confidence regions with method choice
visualize_joint_confidence_regions <- function(boot_results_file, 
											   method = c("empirical", "parametric"),
											   selected_tables = NULL, 
											   confidence_level = 0.95,
											   max_tables = 20) {
	
	method <- match.arg(method)
	
	# Load the bootstrap results
	if(is.character(boot_results_file)) {
		boot_analysis <- readRDS(boot_results_file)
	} else {
		boot_analysis <- boot_results_file
	}
	
	if(method == "empirical") {
		cat("Using EMPIRICAL method (actual bootstrap samples)\n")
		
		# Check if full samples are available
		if(is.null(boot_analysis$full_bootstrap_samples)) {
			cat("Full bootstrap samples not available. Falling back to parametric method.\n")
			method <- "parametric"
		} else {
			return(visualize_empirical_method(boot_analysis, selected_tables, 
											  confidence_level, max_tables))
		}
	}
	
	if(method == "parametric") {
		cat("Using PARAMETRIC method (summary statistics + ellipses)\n")
		return(visualize_parametric_method(boot_analysis, selected_tables, 
										   confidence_level, max_tables))
	}
}

# Internal function for empirical method
visualize_empirical_method <- function(boot_analysis, selected_tables, 
									   confidence_level, max_tables) {
	
	full_samples <- boot_analysis$full_bootstrap_samples
	cat("Found full bootstrap samples for", length(full_samples), "tables\n")
	
	# Select subset if specified
	if(!is.null(selected_tables)) {
		if(max(selected_tables) > length(full_samples)) {
			selected_tables <- selected_tables[selected_tables <= length(full_samples)]
		}
		full_samples <- full_samples[selected_tables]
	} else if(length(full_samples) > max_tables) {
		selected_indices <- sample(length(full_samples), max_tables)
		full_samples <- full_samples[selected_indices]
		cat("Sampling", max_tables, "tables from", length(boot_analysis$full_bootstrap_samples), "for visualization\n")
	}
	
	# Create the empirical data
	all_points <- data.frame()
	region_data <- data.frame()
	
	for(i in 1:length(full_samples)) {
		table_data <- full_samples[[i]]
		boot_samples <- table_data$bootstrap_samples
		
		# Add all bootstrap points
		points_df <- data.frame(
			r_min = boot_samples[, "r_min"],
			r_max = boot_samples[, "r_max"],
			table_id = paste("Config", table_data$config_id, "Table", table_data$sim_id),
			config_id = table_data$config_id,
			sim_id = table_data$sim_id
		)
		all_points <- rbind(all_points, points_df)
		
		# Create empirical confidence region using Mahalanobis distance
		center_x <- mean(boot_samples[, "r_min"])
		center_y <- mean(boot_samples[, "r_max"])
		cov_matrix <- cov(boot_samples)
		
		# Handle potential singular covariance matrices
		if(det(cov_matrix) > 1e-10) {
			mahal_dist <- mahalanobis(boot_samples, c(center_x, center_y), cov_matrix)
			cutoff <- quantile(mahal_dist, confidence_level)
			in_region <- mahal_dist <= cutoff
			boundary_points <- boot_samples[in_region, ]
			
			# Create convex hull
			if(nrow(boundary_points) > 3) {
				hull_indices <- chull(boundary_points[, "r_min"], boundary_points[, "r_max"])
				hull_points <- boundary_points[hull_indices, ]
				
				region_df <- data.frame(
					r_min = hull_points[, "r_min"],
					r_max = hull_points[, "r_max"],
					table_id = paste("Config", table_data$config_id, "Table", table_data$sim_id),
					config_id = table_data$config_id,
					sim_id = table_data$sim_id
				)
				region_data <- rbind(region_data, region_df)
			}
		}
	}
	
	# Create plots
	p1 <- create_confidence_plot(all_points, region_data, confidence_level, 
								 "Empirical", "convex hull regions")
	
	p2 <- ggplot(all_points, aes(x = r_min, y = r_max)) +
		geom_density_2d_filled(alpha = 0.7) +
		geom_point(alpha = 0.5, size = 0.5) +
		facet_wrap(~table_id, scales = "free") +
		geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
		labs(title = "Empirical Joint Density by Table",
			 x = "r_min", y = "r_max") +
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1),
			  strip.text = element_text(size = 8))
	
	return(list(
		method = "empirical",
		confidence_plot = p1,
		density_plot = p2,
		sample_data = all_points,
		region_data = region_data
	))
}

# Internal function for parametric method
visualize_parametric_method <- function(boot_analysis, selected_tables, 
										confidence_level, max_tables) {
	
	stats <- boot_analysis$summary_stats
	
	# Select subset if specified
	if(!is.null(selected_tables)) {
		if(max(selected_tables) > nrow(stats)) {
			selected_tables <- selected_tables[selected_tables <= nrow(stats)]
		}
		stats <- stats[selected_tables, ]
	} else if(nrow(stats) > max_tables) {
		stats <- stats[sample(nrow(stats), max_tables), ]
		cat("Sampling", max_tables, "tables from", nrow(boot_analysis$summary_stats), "for visualization\n")
	}
	
	# Create parametric ellipses using summary statistics
	create_parametric_ellipse <- function(row) {
		rho <- row$boot_cor_rmin_rmax
		sx <- row$boot_sd_rmin
		sy <- row$boot_sd_rmax
		mx <- row$boot_mean_rmin
		my <- row$boot_mean_rmax
		
		# Handle edge cases
		if(is.na(rho) || sx <= 0 || sy <= 0) return(NULL)
		
		# Generate ellipse points
		theta <- seq(0, 2*pi, length.out = 100)
		chi2_val <- qchisq(confidence_level, 2)
		
		# Create correlation matrix and eigendecomposition
		corr_matrix <- matrix(c(1, rho, rho, 1), 2, 2)
		eigen_decomp <- eigen(corr_matrix)
		
		# Scale by standard deviations
		scale_matrix <- diag(c(sx, sy))
		transform_matrix <- scale_matrix %*% eigen_decomp$vectors %*% 
			diag(sqrt(pmax(eigen_decomp$values, 0)))
		
		# Generate ellipse
		ellipse_points <- sqrt(chi2_val) * cbind(cos(theta), sin(theta))
		transformed_points <- t(transform_matrix %*% t(ellipse_points))
		
		data.frame(
			r_min = transformed_points[, 1] + mx,
			r_max = transformed_points[, 2] + my,
			table_id = paste("Config", row$config_id, "Table", row$sim_id),
			config_id = row$config_id,
			sim_id = row$sim_id
		)
	}
	
	# Create ellipse data
	ellipse_list <- lapply(1:nrow(stats), function(i) create_parametric_ellipse(stats[i, ]))
	ellipse_data <- do.call(rbind, ellipse_list[!sapply(ellipse_list, is.null)])
	
	# Create point data for centers
	point_data <- data.frame(
		r_min = stats$boot_mean_rmin,
		r_max = stats$boot_mean_rmax,
		table_id = paste("Config", stats$config_id, "Table", stats$sim_id),
		config_id = stats$config_id,
		sim_id = stats$sim_id
	)
	
	# Create plots
	p1 <- create_confidence_plot(point_data, ellipse_data, confidence_level, 
								 "Parametric", "elliptical regions")
	
	# Correlation vs range plot
	p2 <- ggplot(stats, aes(x = r_range, y = boot_cor_rmin_rmax, 
							color = factor(K1), shape = factor(K2))) +
		geom_point(size = 3, alpha = 0.7) +
		scale_color_viridis_d(name = "K1") +
		labs(title = "Bootstrap Correlation vs Range by Configuration",
			 x = "r_max - r_min", y = "Bootstrap Correlation",
			 shape = "K2") +
		theme_minimal()
	
	return(list(
		method = "parametric", 
		confidence_plot = p1,
		correlation_plot = p2,
		sample_data = point_data,
		region_data = ellipse_data,
		summary_stats = stats
	))
}

# Helper function to create consistent plots
create_confidence_plot <- function(point_data, region_data, confidence_level, 
								   method_name, region_type) {
	
	ggplot() +
		geom_polygon(data = region_data, 
					 aes(x = r_min, y = r_max, group = table_id, 
					 	fill = factor(config_id)), 
					 alpha = 0.3) +
		geom_point(data = point_data, 
				   aes(x = r_min, y = r_max, color = factor(config_id)), 
				   alpha = 0.8, size = 1.5) +
		geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
		scale_fill_viridis_d(name = "Config ID") +
		scale_color_viridis_d(name = "Config ID") +
		labs(title = paste0(method_name, " Joint Distribution of (r_min, r_max)"),
			 subtitle = paste0(confidence_level * 100, "% confidence ", region_type),
			 x = "r_min", y = "r_max") +
		theme_minimal() +
		theme(legend.position = "bottom")
}

# Function to analyze empirical correlations and patterns
analyze_empirical_patterns <- function(boot_results_file) {
	
	# Load the bootstrap results
	if(is.character(boot_results_file)) {
		boot_analysis <- readRDS(boot_results_file)
	} else {
		boot_analysis <- boot_results_file
	}
	
	# Check if full samples are available
	if(is.null(boot_analysis$full_bootstrap_samples)) {
		stop("Full bootstrap samples not available. Rerun bootstrap with store_full_samples = TRUE")
	}
	
	full_samples <- boot_analysis$full_bootstrap_samples
	
	# Extract patterns from each table
	pattern_summary <- data.frame()
	
	for(i in 1:length(full_samples)) {
		table_data <- full_samples[[i]]
		boot_samples <- table_data$bootstrap_samples
		
		# Calculate empirical statistics
		empirical_cor <- cor(boot_samples[, "r_min"], boot_samples[, "r_max"])
		empirical_cov <- cov(boot_samples[, "r_min"], boot_samples[, "r_max"])
		
		# Shape of the distribution (skewness, kurtosis if needed)
		range_r_min <- diff(range(boot_samples[, "r_min"]))
		range_r_max <- diff(range(boot_samples[, "r_max"]))
		
		pattern_summary <- rbind(pattern_summary, data.frame(
			config_id = table_data$config_id,
			sim_id = table_data$sim_id,
			empirical_cor = empirical_cor,
			empirical_cov = empirical_cov,
			range_r_min = range_r_min,
			range_r_max = range_r_max,
			mean_r_min = mean(boot_samples[, "r_min"]),
			mean_r_max = mean(boot_samples[, "r_max"]),
			sd_r_min = sd(boot_samples[, "r_min"]),
			sd_r_max = sd(boot_samples[, "r_max"])
		))
	}
	
	# Merge with config info
	config_patterns <- merge(pattern_summary, boot_analysis$config_info, by = "config_id")
	
	return(config_patterns)
}

# Function to modify bootstrap analysis to store full joint samples
# (For future use - requires modifying the bootstrap_configuration function)
store_full_bootstrap_samples <- function(sim_results, config_ids, B = 1000) {
	cat("This function shows how to modify the analysis to store full samples\n")
	cat("Current code only stores summary statistics\n")
	cat("To get full joint distributions, modify bootstrap_configuration to store:\n")
	cat("boot_result$boot (the full B x 2 matrix of bootstrap samples)\n")
	
	# Example of what the modification would look like:
	cat("\nExample modification to bootstrap_configuration:\n")
	cat("# Add this after computing boot_result:\n")
	cat("full_samples[[sim_id]] <- boot_result$boot\n")
	cat("# Then return: list(summary_stats = summary_stats, full_samples = full_samples)\n")
}

# Quick visualization function
quick_joint_viz <- function(boot_results_file) {
	plots <- visualize_joint_uncertainty_summary(boot_results_file)
	return(plots$combined)
}

# Example usage:
# 
# # Step 1: Run bootstrap WITH full samples (WARNING: Large files!)
# full_analysis <- bootstrap_all_configurations(sim_data, B = 1000, store_full_samples = TRUE)
# save_bootstrap_results(full_analysis, "bootstrap_with_samples.rds")
# 
# # Step 2: Visualize empirical joint distributions
# empirical_plots <- visualize_empirical_joint_distribution("bootstrap_with_samples.rds")
# print(empirical_plots$joint_plot)
# print(empirical_plots$individual_densities)
# 
# # Step 3: Analyze empirical patterns
# patterns <- analyze_empirical_patterns("bootstrap_with_samples.rds")
# print(head(patterns))
# 
# # Step 4: Basic visualization using summary statistics (smaller files)
# basic_plots <- visualize_joint_uncertainty_summary("bootstrap_results.rds")
# print(basic_plots$combined)

# Example usage:
# 
# Basic visualization using summary statistics
plots <- visualize_joint_uncertainty_summary("bootstrap_results.rds")

# Individual plots
print(plots$correlation_heatmap)
print(plots$scatter_plot)

# Confidence ellipses (approximation)
ellipse_plot <- visualize_confidence_ellipses("bootstrap_results.rds")
print(ellipse_plot)

# Empirical method (uses actual bootstrap samples → convex hulls)
empirical_viz <- visualize_joint_confidence_regions(here("R", "ordinal_correlation_analysis", "data", "processed", "bootstrap_with_samples.rds"),
													method = "empirical")
empirical_viz

# Parametric method (uses summary stats → ellipses) 
parametric_viz <- visualize_joint_confidence_regions(here("R", "ordinal_correlation_analysis", "data", "processed", "bootstrap_with_samples.rds"), 
													 method = "parametric")
parametric_viz

# Auto-fallback: tries empirical, falls back to parametric if no full samples
auto_viz <- visualize_joint_confidence_regions("bootstrap_results.rds", 
											   method = "empirical")


# Load and visualize your results
plots <- visualize_joint_uncertainty_summary("bootstrap_results.rds")

# See individual plots
print(plots$correlation_heatmap)
print(plots$scatter_plot)

# Get confidence ellipses
ellipse_plot <- visualize_confidence_ellipses("bootstrap_results.rds")
print(ellipse_plot)

# Quick combined view
quick_joint_viz("bootstrap_results.rds")

#### USE EMPERICAL JOINT SE --- this will take a long time

# # Step 1: Run bootstrap WITH full samples (WARNING: Large files!)
full_analysis <- bootstrap_all_configurations(sim_data, B = 1000, store_full_samples = TRUE)
save_bootstrap_results(full_analysis, here("R", "ordinal_correlation_analysis", "data", "processed", "bootstrap_with_samples.rds"))
# 
# # Step 2: Visualize empirical joint distributions
empirical_plots <- visualize_empirical_joint_distribution(here("R", "ordinal_correlation_analysis", "data", "processed", "bootstrap_with_samples.rds"))
# print(empirical_plots$joint_plot)
# print(empirical_plots$individual_densities)
# 
# # Step 3: Analyze empirical patterns
# patterns <- analyze_empirical_patterns("bootstrap_with_samples.rds")
# print(head(patterns))
# 
# # Step 4: Basic visualization using summary statistics (smaller files)
# basic_plots <- visualize_joint_uncertainty_summary("bootstrap_results.rds")
# print(basic_plots$combined)

