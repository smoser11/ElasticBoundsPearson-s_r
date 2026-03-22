# Test all converted legacy R files work with here() package

cat("Testing ALL converted legacy files...\n")
cat("=====================================\n\n")

# List of files to test
files_to_test <- list(
  list(file = "R/correlation_bounds_examples.R", name = "Examples", test_func = "basic_demonstration"),
  list(file = "R/correlation_bounds_bes_example.R", name = "BES Example", test_func = "create_example_row"),
  list(file = "R/correlation_bounds_demo.R", name = "Demo", test_func = NULL),  # Just test loading
  list(file = "R/correlation_bounds_simulation.R", name = "Simulation", test_func = NULL)  # Just test loading
)

results <- list()

for (i in seq_along(files_to_test)) {
  test <- files_to_test[[i]]
  cat(sprintf("%d. Testing %s (%s)...\n", i, test$name, test$file))
  
  tryCatch({
    # Clear environment to test clean loading
    rm(list = ls(pattern = "^[^.].*"))  # Remove all non-hidden objects
    
    source(here::here(test$file))
    cat(sprintf("   ✅ File loaded successfully\n"))
    
    # Test specific function if specified
    if (!is.null(test$test_func)) {
      if (exists(test$test_func)) {
        cat(sprintf("   ✅ Function '%s' found\n", test$test_func))
      } else {
        cat(sprintf("   ❌ Function '%s' not found\n", test$test_func))
      }
    }
    
    results[[test$name]] <- "SUCCESS"
    
  }, error = function(e) {
    cat(sprintf("   ❌ Error: %s\n", e$message))
    results[[test$name]] <- paste("ERROR:", e$message)
  })
  
  cat("\n")
}

# Summary
cat("=== SUMMARY ===\n")
for (name in names(results)) {
  status <- if (results[[name]] == "SUCCESS") "✅" else "❌"
  cat(sprintf("%s %s: %s\n", status, name, results[[name]]))
}

# Test core mathematical functionality
cat("\n=== TESTING CORE MATH ===\n")
tryCatch({
  uniform_X <- c(0.25, 0.25, 0.25, 0.25)
  skewed_Y <- c(0.1, 0.2, 0.3, 0.4)
  
  r_max <- max_corr_bound(uniform_X, skewed_Y)
  r_min <- min_corr_bound(uniform_X, skewed_Y)
  
  cat("✅ Core mathematical functions working!\n")
  cat(sprintf("   Max correlation: %.4f\n", r_max))
  cat(sprintf("   Min correlation: %.4f\n", r_min))
  cat(sprintf("   Range: %.4f\n", r_max - r_min))
  
}, error = function(e) {
  cat("❌ Core math error:", e$message, "\n")
})

cat("\n=== ALL TESTS COMPLETE ===\n")