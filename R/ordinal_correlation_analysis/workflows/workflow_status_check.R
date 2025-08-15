# ================================================================
# WORKFLOW STATUS CHECK
# Quick overview of all workflow results and cache status
# ================================================================

library(here)
source(here("R", "ordinal_correlation_analysis", "utilities", "cache_management.R"))

cat("📊 ORDINAL CORRELATION ANALYSIS WORKFLOW STATUS\n")
cat("================================================\n\n")

# ================================================================
# CHECK CACHE STATUS ACROSS ALL WORKFLOWS
# ================================================================

cat("🗂️  CACHE STATUS CHECK\n")
cat("=====================\n")

# Check cache directories
cache_dirs <- list(
  "Raw Data" = here("R", "ordinal_correlation_analysis", "data", "raw"),
  "Processed Data" = here("R", "ordinal_correlation_analysis", "data", "processed"), 
  "Analysis Reports" = here("R", "ordinal_correlation_analysis", "output", "reports"),
  "Figures" = here("R", "ordinal_correlation_analysis", "output", "figures"),
  "Tables" = here("R", "ordinal_correlation_analysis", "output", "tables")
)

total_files <- 0
total_size_mb <- 0

for (dir_name in names(cache_dirs)) {
  dir_path <- cache_dirs[[dir_name]]
  
  if (dir.exists(dir_path)) {
    files <- list.files(dir_path, recursive = TRUE)
    if (length(files) > 0) {
      file_paths <- file.path(dir_path, files)
      dir_size_mb <- round(sum(file.size(file_paths), na.rm = TRUE) / 1024^2, 1)
      
      cat("📁", dir_name, ":", length(files), "files (", dir_size_mb, "MB )\n")
      
      # Show most recent files in each directory
      if (length(files) > 0) {
        file_info <- file.info(file_paths)
        recent_files <- files[order(file_info$mtime, decreasing = TRUE)[1:min(3, length(files))]]
        for (f in recent_files) {
          cat("   └─", f, "\n")
        }
      }
      
      total_files <- total_files + length(files)
      total_size_mb <- total_size_mb + dir_size_mb
    } else {
      cat("📁", dir_name, ": empty\n")
    }
  } else {
    cat("📁", dir_name, ": directory not found\n")
  }
  cat("\n")
}

cat("💾 Total cache: ", total_files, " files, ", total_size_mb, " MB\n\n")

# ================================================================
# CHECK WORKFLOW COMPLETION STATUS
# ================================================================

cat("✅ WORKFLOW COMPLETION STATUS\n")
cat("=============================\n")

workflows <- list(
  "Monte Carlo Simulation" = list(
    data_file = "MCsim_numsims.*\\.rds",
    analysis_file = "mc_bounds_analysis.*\\.rds",
    summary_file = "mc_summary_data.*\\.rds"
  ),
  "BES Data Analysis" = list(
    bounds_file = "bes_bounds.*\\.rds", 
    analysis_file = "bes_analysis.*\\.rds"
  ),
  "Bootstrap Uncertainty" = list(
    bootstrap_file = "bootstrap_results.*\\.rds",
    ci_file = "confidence_intervals.*\\.rds"
  )
)

for (workflow_name in names(workflows)) {
  cat("🔬", workflow_name, ":\n")
  workflow_files <- workflows[[workflow_name]]
  
  completed_steps <- 0
  total_steps <- length(workflow_files)
  
  for (file_type in names(workflow_files)) {
    pattern <- workflow_files[[file_type]]
    
    # Check in relevant directories
    found_files <- c()
    
    for (check_dir in c("data/raw", "data/processed", "output/reports")) {
      full_dir <- here("R", "ordinal_correlation_analysis", check_dir)
      if (dir.exists(full_dir)) {
        files <- list.files(full_dir, pattern = pattern, full.names = TRUE)
        found_files <- c(found_files, files)
      }
    }
    
    if (length(found_files) > 0) {
      # Get most recent file
      most_recent <- found_files[which.max(file.mtime(found_files))]
      file_age <- round(as.numeric(Sys.time() - file.mtime(most_recent), units = "hours"), 1)
      file_size <- round(file.size(most_recent) / 1024^2, 1)
      
      cat("   ✅", file_type, ":", basename(most_recent), "(", file_age, "hrs old,", file_size, "MB )\n")
      completed_steps <- completed_steps + 1
    } else {
      cat("   ❌", file_type, ": not found\n")
    }
  }
  
  completion_pct <- round(100 * completed_steps / total_steps, 0)
  cat("   📊 Completion:", completion_pct, "% (", completed_steps, "/", total_steps, " steps)\n\n")
}

# ================================================================
# CHECK VISUALIZATION AND TABLE OUTPUTS
# ================================================================

cat("📈 OUTPUT STATUS\n")
cat("================\n")

# Check figures
figures_dir <- here("R", "ordinal_correlation_analysis", "output", "figures")
if (dir.exists(figures_dir)) {
  figures <- list.files(figures_dir, pattern = "\\.png$")
  cat("📊 Figures generated:", length(figures), "\n")
  
  if (length(figures) > 0) {
    # Group by workflow type
    mc_figures <- grep("mc_", figures, value = TRUE)
    bes_figures <- grep("bes_", figures, value = TRUE) 
    bootstrap_figures <- grep("bootstrap_", figures, value = TRUE)
    
    if (length(mc_figures) > 0) cat("   └─ Monte Carlo:", length(mc_figures), "figures\n")
    if (length(bes_figures) > 0) cat("   └─ BES Analysis:", length(bes_figures), "figures\n")
    if (length(bootstrap_figures) > 0) cat("   └─ Bootstrap:", length(bootstrap_figures), "figures\n")
  }
} else {
  cat("📊 Figures: directory not found\n")
}

# Check tables  
tables_dir <- here("R", "ordinal_correlation_analysis", "output", "tables")
if (dir.exists(tables_dir)) {
  tables <- list.files(tables_dir, pattern = "\\.(csv|html|md)$")
  cat("📋 Tables generated:", length(tables), "\n")
  
  if (length(tables) > 0) {
    csv_tables <- grep("\\.csv$", tables, value = TRUE)
    html_tables <- grep("\\.html$", tables, value = TRUE)
    md_tables <- grep("\\.md$", tables, value = TRUE)
    
    if (length(csv_tables) > 0) cat("   └─ CSV format:", length(csv_tables), "files\n")
    if (length(html_tables) > 0) cat("   └─ HTML format:", length(html_tables), "files\n") 
    if (length(md_tables) > 0) cat("   └─ Markdown format:", length(md_tables), "files\n")
  }
} else {
  cat("📋 Tables: directory not found\n")
}

cat("\n")

# ================================================================
# WORKFLOW RECOMMENDATIONS
# ================================================================

cat("💡 WORKFLOW RECOMMENDATIONS\n")
cat("===========================\n")

# Check what should be run next
mc_data_exists <- length(list.files(here("R", "ordinal_correlation_analysis", "data", "raw"), 
                                   pattern = "MCsim.*\\.rds")) > 0

bes_data_exists <- file.exists(here("R", "ordinal_correlation_analysis", "data", "processed", 
                                   "correlation and other data about pairs of BES2019 variables.dta"))

if (!mc_data_exists) {
  cat("🎲 Run workflow_bivariate_mc.R to generate Monte Carlo simulation data\n")
}

if (!bes_data_exists) {
  cat("🇬🇧 Ensure BES data file is in data/processed/ directory\n")
} else {
  bes_analysis_exists <- length(list.files(here("R", "ordinal_correlation_analysis", "output", "reports"), 
                                          pattern = "bes_analysis.*\\.rds")) > 0
  if (!bes_analysis_exists) {
    cat("🇬🇧 Run workflow_bivariate_bes.R to analyze BES data\n")
  }
}

bootstrap_exists <- length(list.files(here("R", "ordinal_correlation_analysis", "data", "processed"), 
                                     pattern = "bootstrap_results.*\\.rds")) > 0

if (mc_data_exists && !bootstrap_exists) {
  cat("🔄 Run workflow_bivariate_bootstrap.R for uncertainty analysis\n")
}

# Check if ready for matrix analysis
bounds_data_exists <- length(list.files(here("R", "ordinal_correlation_analysis", "output", "reports"), 
                                        pattern = "bes_bounds.*\\.rds")) > 0

if (bounds_data_exists) {
  cat("🔢 Ready to run matrix analysis workflows\n")
  cat("⚖️  Ready to run rescaling analysis workflows\n")
}

cat("\n")

# ================================================================
# QUICK CACHE MANAGEMENT OPTIONS
# ================================================================

cat("🧹 CACHE MANAGEMENT\n")
cat("===================\n")
cat("To clean old cache files (older than 30 days):\n")
cat("   clean_old_cache(here('R', 'ordinal_correlation_analysis'), days_old = 30, dry_run = TRUE)\n\n")

cat("To force regenerate all cached results:\n") 
cat("   # Set force_regenerate = TRUE in workflow parameter sections\n\n")

cat("To check detailed cache statistics:\n")
cat("   report_cache_stats()\n\n")

cat("📝 For detailed workflow documentation, see:\n")
cat("   - workflow_bivariate_mc.R for Monte Carlo simulation\n")
cat("   - workflow_bivariate_bes.R for BES data analysis\n") 
cat("   - workflow_bivariate_bootstrap.R for uncertainty analysis\n\n")