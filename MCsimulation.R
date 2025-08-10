#--------------------------------------
# Cross-Tab Space Simulation
#--------------------------------------
library(ggplot2)
library(gridExtra)
library(dplyr)

# Set parameters for the simulation
K_values <- c(4, 5, 6, 7, 10, 11)
sample_sizes <- c(100, 500, 1000, 2000)
n_simulations <- 10  # Number of random cross-tabs per combination

# Functions for computing correlation bounds (from previous work)
max_corr_analytical <- function(pX, pY) {
	x <- 0:(length(pX)-1)
	y <- 0:(length(pY)-1)
	
	muX <- sum(x*pX)
	muY <- sum(y*pY)
	sX <- sqrt(sum(x^2*pX) - muX^2)
	sY <- sqrt(sum(y^2*pY) - muY^2)
	
	if(sX == 0 || sY == 0) return(0)  # Handle degenerate cases
	
	FX <- cumsum(pX)
	FY <- cumsum(pY)
	joint <- matrix(0, length(x), length(y))
	i <- j <- 1
	lastX <- lastY <- 0
	
	while(i <= length(x) && j <= length(y)) {
		m  <- min(FX[i], FY[j])
		M <- max(lastX, lastY)
		joint[i,j] <- m - M
		
		if(FX[i] < FY[j]) {
			lastX <- FX[i]
			i <- i + 1
		} else if(FX[i] > FY[j]) {
			lastY <- FY[j]
			j <- j + 1
		} else {
			lastX <- FX[i]
			lastY <- FY[j]
			i <- i+1
			j <- j+1
		}
	}
	
	EXY <- sum(outer(x, y, FUN="*") * joint)
	(EXY - muX*muY) / (sX * sY)
}

min_corr_analytical <- function(pX, pY) {
	x <- 0:(length(pX)-1)
	y <- 0:(length(pY)-1)
	
	# Create sample vectors
	N <- 10000  # Large sample for approximation
	countsX <- round(pX * N)
	countsY <- round(pY * N)
	
	# Adjust for rounding errors
	countsX[1] <- countsX[1] + (N - sum(countsX))
	countsY[1] <- countsY[1] + (N - sum(countsY))
	
	x_vec <- rep(x, countsX)
	y_vec <- rep(y, countsY)
	
	# Countermonotonic pairing
	x_sorted <- sort(x_vec, decreasing = TRUE)
	y_sorted <- sort(y_vec, decreasing = FALSE)
	
	cor(x_sorted, y_sorted)
}

# Function to generate random cross-tab
generate_random_crosstab <- function(K1, K2, N) {
	# Generate random counts using multinomial distribution
	# Each cell gets a random probability, then we sample N observations
	cell_probs <- runif(K1 * K2)
	cell_probs <- cell_probs / sum(cell_probs)
	
	# Sample N observations into the K1*K2 cells
	cell_counts <- as.vector(rmultinom(1, N, cell_probs))
	
	# Convert to matrix form
	crosstab <- matrix(cell_counts, nrow = K1, ncol = K2)
	
	return(crosstab)
}

# Function to extract marginals from cross-tab
get_marginals <- function(crosstab) {
	marginal_X <- rowSums(crosstab) / sum(crosstab)
	marginal_Y <- colSums(crosstab) / sum(crosstab)
	
	return(list(pX = marginal_X, pY = marginal_Y))
}

# Function to calculate empirical correlation from cross-tab
empirical_correlation <- function(crosstab) {
	K1 <- nrow(crosstab)
	K2 <- ncol(crosstab)
	
	# Create vectors of observations
	x_vals <- c()
	y_vals <- c()
	
	for(i in 1:K1) {
		for(j in 1:K2) {
			if(crosstab[i,j] > 0) {
				x_vals <- c(x_vals, rep(i-1, crosstab[i,j]))  # 0-indexed
				y_vals <- c(y_vals, rep(j-1, crosstab[i,j]))  # 0-indexed
			}
		}
	}
	
	if(length(x_vals) < 2) return(0)
	return(cor(x_vals, y_vals))
}

# Initialize results storage
results <- data.frame()

# Main simulation loop
set.seed(42)  # For reproducibility
cat("Starting cross-tab simulation...\n")

sim_count <- 0
total_sims <- length(K_values) * length(K_values) * length(sample_sizes) * n_simulations

for(K1 in K_values) {
	for(K2 in K_values) {
		for(N in sample_sizes) {
			for(sim in 1:n_simulations) {
				sim_count <- sim_count + 1
				if(sim_count %% 50 == 0) {
					cat(sprintf("Progress: %d/%d (%.1f%%)\n", 
								sim_count, total_sims, 100*sim_count/total_sims))
				}
				
				# Generate random cross-tab
				crosstab <- generate_random_crosstab(K1, K2, N)
				
				# Extract marginals
				marginals <- get_marginals(crosstab)
				pX <- marginals$pX
				pY <- marginals$pY
				
				# Calculate theoretical bounds
				r_max <- max_corr_analytical(pX, pY)
				r_min <- min_corr_analytical(pX, pY)
				
				# Calculate empirical correlation
				r_empirical <- empirical_correlation(crosstab)
				
				# Calculate asymmetry measures
				delta <- r_max + r_min
				range_r <- r_max - r_min
				
				# Calculate marginal asymmetries (TV distance from symmetry)
				tv_X <- 0.5 * sum(abs(pX - rev(pX)))
				tv_Y <- 0.5 * sum(abs(pY - rev(pY)))
				
				# Store results
				results <- rbind(results, data.frame(
					K1 = K1,
					K2 = K2,
					N = N,
					simulation = sim,
					r_max = r_max,
					r_min = r_min,
					r_empirical = r_empirical,
					delta = delta,
					range_r = range_r,
					tv_X = tv_X,
					tv_Y = tv_Y,
					symmetric_dims = (K1 == K2),
					stringsAsFactors = FALSE
				))
			}
		}
	}
}

cat("Simulation complete!\n")

#--------------------------------------
# Analysis and Visualization
#--------------------------------------

# Summary statistics
cat("\nSummary Statistics:\n")
cat("===================\n")
print(summary(results[, c("r_max", "r_min", "delta", "range_r", "tv_X", "tv_Y")]))

# When is delta significantly different from 0?
results$delta_significant <- abs(results$delta) > 0.01

cat("\nProportion of cases where |r_max + r_min| > 0.01:\n")
prop_asymmetric <- mean(results$delta_significant)
cat(sprintf("%.3f (%.1f%%)\n", prop_asymmetric, 100*prop_asymmetric))

# Analyze by dimension symmetry
cat("\nAsymmetry by dimension matching:\n")
asymmetry_by_dims <- results %>%
	group_by(symmetric_dims) %>%
	summarise(
		mean_abs_delta = mean(abs(delta)),
		prop_asymmetric = mean(delta_significant),
		.groups = 'drop'
	)
print(asymmetry_by_dims)

#--------------------------------------
# Visualizations - r_min vs r_max plots
#--------------------------------------

# Add derived variables for plotting
results$max_tv <- pmax(results$tv_X, results$tv_Y)
results$total_categories <- results$K1 + results$K2
results$asymmetry_level <- cut(abs(results$delta), 
							   breaks = c(0, 0.01, 0.05, 0.1, Inf),
							   labels = c("Symmetric", "Low", "Medium", "High"))

# Plot 1: r_min vs r_max colored by asymmetry level
p1 <- ggplot(results, aes(x = r_max, y = r_min, color = asymmetry_level)) +
	geom_point(alpha = 0.6, size = 1.5) +
	geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "black") +
	scale_color_manual(values = c("Symmetric" = "blue", 
								  "Low" = "green", 
								  "Medium" = "orange", 
								  "High" = "red")) +
	labs(title = "r_min vs r_max Colored by Bound Asymmetry",
		 subtitle = "Dashed line shows perfect symmetry (r_min = -r_max)",
		 x = "r_max", y = "r_min", color = "Asymmetry\nLevel") +
	theme_minimal() +
	theme(legend.position = "right")

# Plot 2: r_min vs r_max colored by total number of categories
p2 <- ggplot(results, aes(x = r_max, y = r_min, color = total_categories)) +
	geom_point(alpha = 0.6, size = 1.5) +
	geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "black") +
	scale_color_gradient(low = "lightblue", high = "darkred") +
	labs(title = "r_min vs r_max Colored by Total Categories",
		 subtitle = "Darker colors indicate more categories (K1 + K2)",
		 x = "r_max", y = "r_min", color = "Total\nCategories") +
	theme_minimal() +
	theme(legend.position = "right")

# Plot 3: r_min vs r_max colored by dimension symmetry
results$dimension_type <- ifelse(results$K1 == results$K2, "Square", "Rectangular")
p3 <- ggplot(results, aes(x = r_max, y = r_min, color = dimension_type)) +
	geom_point(alpha = 0.6, size = 1.5) +
	geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "black") +
	scale_color_manual(values = c("Square" = "purple", "Rectangular" = "orange")) +
	labs(title = "r_min vs r_max Colored by Table Shape",
		 subtitle = "Square tables (K1=K2) vs Rectangular tables (K1≠K2)",
		 x = "r_max", y = "r_min", color = "Table\nShape") +
	theme_minimal() +
	theme(legend.position = "right")

p3
# Plot 4: r_min vs r_max colored by sample size
p4 <- ggplot(results, aes(x = r_max, y = r_min, color = factor(N))) +
	geom_point(alpha = 0.6, size = 1.5) +
	geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "black") +
	scale_color_brewer(type = "qual", palette = "Set1") +
	labs(title = "r_min vs r_max Colored by Sample Size",
		 subtitle = "Different sample sizes show similar correlation bounds",
		 x = "r_max", y = "r_min", color = "Sample\nSize") +
	theme_minimal() +
	theme(legend.position = "right")

p4
# Plot 5: Distribution of delta values
p5 <- ggplot(results, aes(x = delta)) +
	geom_histogram(bins = 50, alpha = 0.7, fill = "skyblue") +
	geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
	labs(title = "Distribution of Δ = r_max + r_min",
		 subtitle = "Values near 0 indicate symmetric bounds",
		 x = "Δ (r_max + r_min)", y = "Frequency") +
	theme_minimal()
p5
# Plot 6: Delta vs marginal asymmetry
p6 <- ggplot(results, aes(x = max_tv, y = abs(delta), color = total_categories)) +
	geom_point(alpha = 0.6) +
	geom_smooth(method = "loess", se = FALSE, color = "black") +
	scale_color_gradient(low = "lightblue", high = "darkred") +
	labs(title = "Bound Asymmetry vs Marginal Asymmetry",
		 subtitle = "Relationship between marginal and bound asymmetries",
		 x = "Max TV Distance (X or Y)", y = "|Δ| (|r_max + r_min|)",
		 color = "Total\nCategories") +
	theme_minimal()

# Arrange the main plots
grid.arrange(p1, p2, p3, p4, ncol = 2)

# Show additional analysis plots
grid.arrange(p5, p6, ncol = 2)

#--------------------------------------
# Summary table for key findings
#--------------------------------------

# Create summary by dimension combinations
summary_by_dims <- results %>%
	group_by(K1, K2) %>%
	summarise(
		n_sims = n(),
		mean_r_max = mean(r_max),
		mean_r_min = mean(r_min),
		mean_delta = mean(delta),
		mean_range = mean(range_r),
		prop_asymmetric = mean(delta_significant),
		.groups = 'drop'
	) %>%
	arrange(K1, K2)

cat("\nSummary by Table Dimensions:\n")
cat("============================\n")
print(summary_by_dims)

# Identify most extreme cases
cat("\nMost asymmetric cases (top 10 by |delta|):\n")
names(results)
top_asymmetric <- results %>% 
	select(K1, K2, N, r_max, r_min, delta, tv_X, tv_Y) %>%
	arrange(desc(abs(delta))) %>%
	head(10) 
print(top_asymmetric)

