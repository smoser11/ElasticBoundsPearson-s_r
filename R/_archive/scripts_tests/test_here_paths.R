# Test that here package works correctly
if (!require(here, quietly = TRUE)) {
  install.packages("here")
}
library(here)

cat("Project root detected by here():", here(), "\n")
cat("This should be: /Users/sm38679/Documents/GitHub/ElasticBoundsPearson-s_r\n")

# Test some key paths
cat("Core file path:", here("R", "correlation_bounds_core.R"), "\n")
cat("BES data path:", here("R", "ordinal_correlation_analysis", "data", "processed"), "\n")
cat("Figures path:", here("R", "ordinal_correlation_analysis", "output", "figures"), "\n")

# Verify core file exists
if (file.exists(here("R", "correlation_bounds_core.R"))) {
  cat("✅ Core file found!\n")
} else {
  cat("❌ Core file not found!\n")
}