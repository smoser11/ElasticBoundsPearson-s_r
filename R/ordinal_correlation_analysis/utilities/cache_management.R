# cache_management.R
# Smart caching utilities for workflow orchestration

library(here)

#' Cache-or-compute pattern for expensive operations
#' 
#' @param cache_file Full path to cache file
#' @param compute_func Function that performs the computation
#' @param force_regenerate Logical, whether to ignore existing cache
#' @return Result of computation (from cache or fresh)
cache_or_compute <- function(cache_file, compute_func, force_regenerate = FALSE) {
  
  if (file.exists(cache_file) && !force_regenerate) {
    cat("📂 Loading cached result:", basename(cache_file), "\n")
    cat("   File size:", format(file.size(cache_file), units = "MB"), "\n")
    cat("   Modified:", format(file.mtime(cache_file), "%Y-%m-%d %H:%M"), "\n")
    return(readRDS(cache_file))
    
  } else {
    cat("⚡ Computing new result for:", basename(cache_file), "\n")
    
    # Ensure directory exists
    dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
    
    # Time the computation
    start_time <- Sys.time()
    result <- compute_func()
    end_time <- Sys.time()
    
    # Save result
    saveRDS(result, cache_file)
    
    # Report
    computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    cat("💾 Cached result saved:", basename(cache_file), "\n")
    cat("   Computation time:", round(computation_time, 1), "seconds\n")
    cat("   File size:", format(file.size(cache_file), units = "MB"), "\n")
    
    return(result)
  }
}

#' Generate parameter-aware cache filename
#' 
#' @param base_name Base name for the cache file
#' @param params List of parameters that affect the computation
#' @param extension File extension (default: "rds")
#' @return Full cache filename with parameters encoded
generate_cache_filename <- function(base_name, params, extension = "rds") {
  # Create parameter string
  param_string <- paste(names(params), params, sep = "", collapse = "_")
  param_string <- gsub("[^A-Za-z0-9_]", "", param_string)  # Remove special characters
  
  filename <- paste0(base_name, "_", param_string, ".", extension)
  return(filename)
}

#' Check cache status for multiple files
#' 
#' @param cache_files Vector of cache file paths
#' @return Data frame with cache status information
check_cache_status <- function(cache_files) {
  status_df <- data.frame(
    file = basename(cache_files),
    exists = file.exists(cache_files),
    size_mb = ifelse(file.exists(cache_files), 
                     round(file.size(cache_files) / 1024^2, 2), 
                     NA),
    modified = ifelse(file.exists(cache_files),
                      as.character(file.mtime(cache_files)),
                      NA),
    stringsAsFactors = FALSE
  )
  
  return(status_df)
}

#' Clean old cache files based on age
#' 
#' @param cache_dir Directory containing cache files
#' @param days_old Remove files older than this many days
#' @param dry_run If TRUE, just report what would be deleted
clean_old_cache <- function(cache_dir, days_old = 30, dry_run = TRUE) {
  
  if (!dir.exists(cache_dir)) {
    cat("Cache directory doesn't exist:", cache_dir, "\n")
    return(invisible())
  }
  
  cache_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
  
  if (length(cache_files) == 0) {
    cat("No cache files found in:", cache_dir, "\n")
    return(invisible())
  }
  
  # Find old files
  file_ages <- as.numeric(Sys.time() - file.mtime(cache_files), units = "days")
  old_files <- cache_files[file_ages > days_old]
  
  if (length(old_files) == 0) {
    cat("No cache files older than", days_old, "days found\n")
    return(invisible())
  }
  
  if (dry_run) {
    cat("Would delete", length(old_files), "cache files older than", days_old, "days:\n")
    for (f in old_files) {
      age <- round(file_ages[cache_files == f], 1)
      size <- round(file.size(f) / 1024^2, 1)
      cat("  -", basename(f), "(", age, "days old,", size, "MB )\n")
    }
    cat("Run with dry_run = FALSE to actually delete\n")
  } else {
    deleted_size <- sum(file.size(old_files))
    unlink(old_files)
    cat("Deleted", length(old_files), "cache files, freed", 
        round(deleted_size / 1024^2, 1), "MB\n")
  }
}

#' Report cache directory statistics
#' 
#' @param cache_dir Directory to analyze
report_cache_stats <- function(cache_dir = here("R", "ordinal_correlation_analysis")) {
  
  cat("📊 CACHE STATISTICS\n")
  cat("===================\n")
  
  # Find all cache directories
  cache_dirs <- c(
    file.path(cache_dir, "data", "raw"),
    file.path(cache_dir, "data", "processed"), 
    file.path(cache_dir, "output", "reports"),
    file.path(cache_dir, "output", "figures"),
    file.path(cache_dir, "output", "tables")
  )
  
  total_files <- 0
  total_size <- 0
  
  for (dir in cache_dirs) {
    if (dir.exists(dir)) {
      files <- list.files(dir, pattern = "\\.(rds|csv|png)$", recursive = TRUE)
      if (length(files) > 0) {
        file_paths <- file.path(dir, files)
        dir_size <- sum(file.size(file_paths), na.rm = TRUE)
        
        cat("📁", basename(dir), ":", length(files), "files,", 
            round(dir_size / 1024^2, 1), "MB\n")
        
        total_files <- total_files + length(files)
        total_size <- total_size + dir_size
      }
    }
  }
  
  cat("-------------------\n")
  cat("📊 Total:", total_files, "files,", round(total_size / 1024^2, 1), "MB\n")
}