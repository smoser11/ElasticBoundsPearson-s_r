# Monte Carlo simulation for contingency tables
# Sweeps over different configurations of categories and sample sizes


sample_uniform_table <- function(K1, K2, N) {
	K <- K1 * K2
	if (K == 1L) return(matrix(N, nrow = K1, ncol = K2))
	
	# Use the standard "stars and bars" method for generating uniform random counts
	# This ensures all counts are non-negative and sum to N
	if (K - 1L > 0) {
		cuts <- sort(sample.int(N + K - 1L, K - 1L))
		# Add boundary points
		cuts <- c(0L, cuts, N + K)
		
		# Calculate differences to get the counts
		x <- integer(K)
		for (i in seq_len(K)) {
			x[i] <- cuts[i + 1L] - cuts[i] - 1L
		}
	} else {
		x <- N
	}
	
	# Ensure non-negative counts and correct sum
	x <- pmax(x, 0L)  # Ensure no negative values
	
	# Adjust if sum doesn't match N due to rounding
	diff <- N - sum(x)
	if (diff != 0) {
		# Add/subtract the difference to a random cell
		cell_to_adjust <- sample(K, 1)
		x[cell_to_adjust] <- x[cell_to_adjust] + diff
		# Ensure no negative values after adjustment
		x[cell_to_adjust] <- max(x[cell_to_adjust], 0L)
	}
	
	matrix(x, nrow = K1, ncol = K2, byrow = TRUE)
}

# Function to generate a single random contingency table
generate_random_table <- function(K1, K2, N) {
	# Initialize empty table
	table_matrix <- matrix(0, nrow = K1, ncol = K2)
	
	# Distribute all N observations
	remaining <- N
	
	while(remaining > 0) {
		# Randomly decide how much to allocate (between 1 and remaining)
		# Using a geometric-like distribution to allow various allocation sizes
		allocation <- sample(1:remaining, 1)
		
		# Randomly choose a cell
		cell_row <- sample(1:K1, 1)
		cell_col <- sample(1:K2, 1)
		
		# Add allocation to chosen cell
		table_matrix[cell_row, cell_col] <- table_matrix[cell_row, cell_col] + allocation
		
		# Update remaining count
		remaining <- remaining - allocation
	}
	
	return(table_matrix)
}

# Main simulation function
run_mc_simulation <- function(numsims = 1000, seed = 123, 
							  categories_var1 = c(4, 5, 6, 7, 10, 11,13,14),
							  categories_var2 = c(4, 5, 6, 7, 10, 11,13,14),
							  sample_sizes = c(100, 200, 500, 1000, 2000, 5000)) {
	set.seed(seed)
	
	# Use user-provided parameter ranges
	categories1 <- categories_var1
	categories2 <- categories_var2
	
	# Create storage structure
	# Using a nested list: results[[config_id]][[sim_id]] = table_matrix
	results <- list()
	config_info <- data.frame()
	
	config_id <- 1
	
	cat("Starting Monte Carlo simulation...\n")
	cat("Categories for var1:", paste(categories1, collapse = ", "), "\n")
	cat("Categories for var2:", paste(categories2, collapse = ", "), "\n")
	cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
	cat("Total configurations:", length(categories1) * length(categories2) * length(sample_sizes), "\n\n")
	
	# Sweep over all combinations
	for(K1 in categories1) {
		for(K2 in categories2) {
			for(N in sample_sizes) {
				
				cat(sprintf("Configuration %d: K1=%d, K2=%d, N=%d\n", 
							config_id, K1, K2, N))
				
				# Store configuration info
				config_info <- rbind(config_info, 
									 data.frame(config_id = config_id,
									 		   K1 = K1, K2 = K2, N = N))
				
				# Generate numsims tables for this configuration
				config_tables <- list()
				
				for(sim in 1:numsims) {
					config_tables[[sim]] <- sample_uniform_table(K1, K2, N)
				}
				
				# Store results for this configuration
				results[[config_id]] <- config_tables
				
				config_id <- config_id + 1
			}
		}
	}
	
	cat("\nSimulation complete!\n")
	cat("Generated", length(results), "configurations with", numsims, "tables each\n")
	cat("Total tables generated:", length(results) * numsims, "\n")
	
	return(list(
		tables = results,
		config_info = config_info,
		numsims = numsims,
		seed = seed
	))
}

# Function to save simulation results
save_simulation <- function(sim_results, filename = "mc_contingency_simulation.rds") {
	saveRDS(sim_results, file = filename)
	cat("Simulation saved to:", filename, "\n")
	cat("File size:", round(file.size(filename) / 1024^2, 2), "MB\n")
}

# Function to load simulation results
load_simulation <- function(filename = "mc_contingency_simulation.rds") {
	if(!file.exists(filename)) {
		stop("File does not exist: ", filename)
	}
	sim_results <- readRDS(filename)
	cat("Simulation loaded from:", filename, "\n")
	cat("Configurations:", nrow(sim_results$config_info), "\n")
	cat("Simulations per config:", sim_results$numsims, "\n")
	return(sim_results)
}

# Function to examine a specific table
examine_table <- function(sim_results, config_id, sim_id) {
	if(config_id > length(sim_results$tables)) {
		stop("config_id out of range")
	}
	if(sim_id > length(sim_results$tables[[config_id]])) {
		stop("sim_id out of range")
	}
	
	table_matrix <- sim_results$tables[[config_id]][[sim_id]]
	config <- sim_results$config_info[config_id, ]
	
	cat("Configuration:", config_id, "\n")
	cat("K1 =", config$K1, ", K2 =", config$K2, ", N =", config$N, "\n")
	cat("Simulation:", sim_id, "\n")
	cat("Table sum:", sum(table_matrix), "\n\n")
	
	print(table_matrix)
	return(table_matrix)
}

# Example usage:
# Run simulation with default categories (this will take some time!)
# sim_data <- run_mc_simulation(numsims = 100)  
#
# Run with custom categories
# sim_data <- run_mc_simulation(numsims = 100, 
#                              categories_var1 = c(3, 4, 5), 
#                              categories_var2 = c(4, 6, 8),
#                              sample_sizes = c(100, 500))
# 
# Save the results
# save_simulation(sim_data, "my_contingency_sim.rds")
#
# Load results later
# loaded_data <- load_simulation("my_contingency_sim.rds")
#
# Examine a specific table
# examine_table(loaded_data, config_id = 1, sim_id = 1)


# Example usage (commented out to prevent execution when sourcing):
# sim_data <- run_mc_simulation(numsims = 10) 
# library(here)
# save_simulation(sim_data, here("R", "ordinal_correlation_analysis", "data", "raw", "MCsim.rds"))
# examine_table(sim_data, config_id = 1, sim_id = 1)
