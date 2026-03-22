# Test the converted modular structure files

cat("Testing converted modular structure files...\n")
cat("==========================================\n\n")

# Test if main analysis loads (expect data errors but path resolution should work)
cat("1. Testing main_analysis.R path resolution...\n")
tryCatch({
  # Don't run the full analysis, just test if paths resolve
  source(here::here("R", "ordinal_correlation_analysis", "main_analysis.R"))
  cat("✅ Main analysis file sourced successfully\n")
  
  # Check if key functions exist
  if (exists("run_complete_analysis")) {
    cat("✅ run_complete_analysis() function found\n")
  } else {
    cat("❌ run_complete_analysis() function not found\n")
  }
  
}, error = function(e) {
  cat("❌ Error in main_analysis.R:", e$message, "\n")
})

cat("\n2. Testing module coordinators...\n")

# Test bivariate main
tryCatch({
  source(here::here("R", "ordinal_correlation_analysis", "1_bivariate_ordcats_correlation", "bivariate_main.R"))
  cat("✅ Bivariate coordinator loaded\n")
}, error = function(e) {
  cat("❌ Bivariate coordinator error:", e$message, "\n")
})

# Test matrices main  
tryCatch({
  source(here::here("R", "ordinal_correlation_analysis", "2_correlation_matrices", "matrices_main.R"))
  cat("✅ Matrices coordinator loaded\n")
}, error = function(e) {
  cat("❌ Matrices coordinator error:", e$message, "\n")
})

# Test rescaling main
tryCatch({
  source(here::here("R", "ordinal_correlation_analysis", "3_fixes_and_rescaling", "rescaling_main.R"))
  cat("✅ Rescaling coordinator loaded\n")
}, error = function(e) {
  cat("❌ Rescaling coordinator error:", e$message, "\n")
})

cat("\n=== MODULAR PATH CONVERSION TEST COMPLETE ===\n")
cat("Note: Some functions may be missing but here() paths should work\n")