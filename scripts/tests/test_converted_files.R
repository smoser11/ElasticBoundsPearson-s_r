# Test that the converted files work correctly with here()

cat("Testing converted files with here() package...\n")

# Test 1: Test examples file
cat("\n1. Testing correlation_bounds_examples.R...\n")
tryCatch({
  source(here::here("R", "correlation_bounds_examples.R"))
  cat("✅ Examples file loaded successfully\n")
  
  # Test a function exists
  if (exists("basic_demonstration")) {
    cat("✅ basic_demonstration() function found\n")
  } else {
    cat("❌ basic_demonstration() function not found\n")
  }
}, error = function(e) {
  cat("❌ Error loading examples:", e$message, "\n")
})

# Test 2: Test BES example file  
cat("\n2. Testing correlation_bounds_bes_example.R...\n")
tryCatch({
  source(here::here("R", "correlation_bounds_bes_example.R"))
  cat("✅ BES example file loaded successfully\n")
  
  # Test a function exists
  if (exists("create_example_row")) {
    cat("✅ create_example_row() function found\n")
  } else {
    cat("❌ create_example_row() function not found\n")
  }
}, error = function(e) {
  cat("❌ Error loading BES example:", e$message, "\n")
})

# Test 3: Test BES functions file
cat("\n3. Testing correlation_bounds_bes.R...\n") 
tryCatch({
  source(here::here("R", "correlation_bounds_bes.R"))
  cat("✅ BES functions file loaded successfully\n")
  
  # Test a function exists
  if (exists("calculate_corr_bounds")) {
    cat("✅ calculate_corr_bounds() function found\n")
  } else {
    cat("❌ calculate_corr_bounds() function not found\n")
  }
}, error = function(e) {
  cat("❌ Error loading BES functions:", e$message, "\n")
})

cat("\n=== CONVERSION TEST COMPLETE ===\n")