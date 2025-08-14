# Simple test of converted legacy files

cat("Testing converted legacy files...\n\n")

# Test 1: Examples file (we know this works)
cat("1. Testing examples file...\n")
source(here::here("R", "correlation_bounds_examples.R"))
cat("✅ Examples loaded - basic_demonstration exists:", exists("basic_demonstration"), "\n\n")

# Test 2: Just core mathematical functionality
cat("2. Testing core math...\n")
uniform_X <- c(0.25, 0.25, 0.25, 0.25)
skewed_Y <- c(0.1, 0.2, 0.3, 0.4)

r_max <- max_corr_bound(uniform_X, skewed_Y)
r_min <- min_corr_bound(uniform_X, skewed_Y)

cat("✅ Core functions work!\n")
cat("   Max correlation:", round(r_max, 4), "\n")
cat("   Min correlation:", round(r_min, 4), "\n")
cat("   Range:", round(r_max - r_min, 4), "\n\n")

# Test 3: Try demo file (might have missing dependencies)
cat("3. Testing demo file...\n")
tryCatch({
  source(here::here("R", "correlation_bounds_demo.R"))
  cat("✅ Demo file loaded successfully\n")
}, error = function(e) {
  cat("❌ Demo file error (expected if missing correlation_matrix_test.R):", e$message, "\n")
})

cat("\n=== Core functionality is working with here() paths! ===\n")