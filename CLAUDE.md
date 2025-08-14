# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a statistical research project exploring **theoretical bounds on Pearson's correlation coefficient for categorical variables** with fixed marginal distributions. The core insight is that correlation bounds are "elastic" - they stretch and contract based on marginal distributions, affecting both the range of possible correlations and statistical inference.

## Common Commands

### New Modular Architecture (Preferred)
```r
# Run complete analysis pipeline using new structure
source("R/ordinal_correlation_analysis/main_analysis.R")

# Load BES data and run comprehensive analysis
bes_data <- load_bes_data()  # From utilities/data_loading.R
results <- run_complete_analysis(bes_data)

# Run specific modules
source("R/ordinal_correlation_analysis/1_bivariate_ordcats_correlation/bivariate_main.R")
bivariate_results <- run_bivariate_analysis(bes_data, config)
```

### Legacy Functions (Still Functional)
```r
# Source core functions from legacy structure
source("R/correlation_bounds_core.R")
source("R/correlation-bounds-visualization.R") 
source("R/correlation_bounds_examples.R")

# Run comprehensive examples
results <- run_all_examples()

# Run specific demonstrations
basic_results <- basic_demonstration()
range_results <- range_comparison_demonstration()

# Analyze real data (BES2019)
source("R/correlation_bounds_bes_example.R")
```

### Document Generation
```bash
# Render latest manuscript version (from project root)
quarto render paper/ElasticBounds-r_v4b.qmd
quarto render paper/ElasticBounds-r_v4b.qmd --to pdf
quarto render paper/ElasticBounds-r_v4b.qmd --to docx

# Render specific versions
quarto render paper/ElasticBounds-r_v3b.qmd --to pdf
```

### Required R Packages
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "gridExtra", 
                   "Matrix", "matrixcalc", "MASS", "readstata13"))
```

## Code Architecture

### New Modular Structure (Primary)

The codebase has been reorganized into a sophisticated modular architecture under `R/ordinal_correlation_analysis/`:

#### 1. Bivariate Analysis (`1_bivariate_ordcats_correlation/`)
- **Core bounds computation**: `1_rmin_rmax_rhat/1_monte_carlo_simulation/core_bounds_functions.R`
- **Real-world validation**: `1_rmin_rmax_rhat/2_bes_illustrative_example/bes_data_analysis.R`
- **Bootstrap uncertainty**: `1_rmin_rmax_rhat/bootstrapJoint_r_minmax.R`
- **Asymmetry analysis**: `2_asymmetry_analysis/asymmetry_measures.R`
- **Visualization**: `3_visualization/bounds_visualization.R`
- **Module coordinator**: `bivariate_main.R`

#### 2. Correlation Matrices (`2_correlation_matrices/`)
- **Matrix diagnostics**: `1_matrix_properties/` (condition numbers, PSD validation, invertibility)
- **Matrix construction**: `2_matrix_construction/random_trials.R`
- **Module coordinator**: `matrices_main.R`

#### 3. Rescaling Methods (`3_fixes_and_rescaling/`)
- **Simple rescaling**: `1_simple_rescaling/linear_rescaling.R`
- **Advanced methods**: `2_advanced_rescaling/` (future extensions)
- **Module coordinator**: `rescaling_main.R`

#### 4. Utilities and Data Management
- **Data loading**: `utilities/data_loading.R`
- **Helper functions**: `utilities/helper_functions.R`
- **Processed datasets**: `data/processed/` (bootstrap results, BES correlations)
- **Generated outputs**: `output/` (figures, reports, tables)

### Legacy Functions (Still Functional)

#### Core Mathematical Functions (`R/correlation_bounds_core.R`)
- `max_corr_bound()` - Computes maximum correlation using comonotonic coupling (Fréchet–Hoeffding upper bound)
- `min_corr_bound()` - Computes minimum correlation using anti-comonotonic coupling  
- `simulate_permutation_r()` - Generates permutation distributions under null hypothesis
- `analyze_all_corr_bounds()` - Comprehensive analysis function for datasets

#### Visualization System (`R/correlation-bounds-visualization.R`)
- `plot_permutation_distribution()` - Shows empirical distributions with theoretical bounds
- `plot_bounds_summary()` - Summary plots comparing bounds and confidence intervals
- `plot_range_comparison()` - Comparison plots across scenarios
- `plot_significance_comparison()` - Compares different significance testing approaches

#### Example Framework (`R/correlation_bounds_examples.R`)
- `basic_demonstration()` - Explores uniform vs extreme distributions
- `range_comparison_demonstration()` - Compares theoretical vs empirical ranges
- `run_all_examples()` - Executes comprehensive example suite

#### Applied Analysis (`R/correlation_bounds_bes.R`)
- Functions for analyzing British Election Study (BES) 2019 data
- Real-world validation of theoretical bounds with survey data

## Development Patterns

### File Organization
- **New modular structure**: `R/ordinal_correlation_analysis/` contains reorganized production code
- **Legacy functions**: Root of `R/` directory contains original working implementations
- **Experimental work**: `R/scratchWork/` for development versions and prototypes
- **Academic output**: `paper/` directory for Quarto manuscripts (latest: ElasticBounds-r_v4b.qmd)
- **Research documentation**: `Notes/` and `Meetings/` for mathematical proofs and research discussions

### Mathematical Framework
The project implements Fréchet–Hoeffding bounds adapted for categorical variables:
- **Comonotonic coupling**: Both variables sorted in same direction (maximum correlation)
- **Anti-comonotonic coupling**: Variables sorted in opposite directions (minimum correlation)
- **Permutation inference**: Random permutation preserves marginals for null hypothesis testing

### Data Structures
Functions typically work with:
- Marginal distributions as probability vectors or count vectors
- Results returned as lists with `$summary` and `$simulations` components
- Visualization functions expect data frames with specific column naming conventions

### Significance Testing Innovation
The project explores how theoretical bounds affect statistical inference by comparing:
- Traditional t-tests for correlation significance
- Permutation-based randomization tests that respect marginal constraints
- How bounds influence interpretation of "significant" correlations

### Architecture Transition
The project is transitioning from legacy functions to a new modular architecture:
- **Current state**: Both architectures coexist and are functional
- **Legacy code**: Fully working, well-tested, documented in research papers
- **New modular code**: Enhanced organization, improved maintainability, extended functionality
- **Data flow**: New structure processes data through `main_analysis.R` → module coordinators → specific analyses
- **Recommendation**: Use new modular structure for new development; legacy functions remain available for validation

## Manuscript Generation

The project uses Quarto for academic manuscript generation with:
- PDF output via pdflatex with mathematical theorem environments
- Bibliography management with APA style (`apa.csl`)
- Multiple output formats (PDF, DOCX, HTML) configured in `quarto.yml`
- Version control of manuscripts with descriptive commit messages

## Working with Both Architectures

### For New Development
```r
# Use the modular architecture for new features
source("R/ordinal_correlation_analysis/main_analysis.R")

# Access specific modules as needed
source("R/ordinal_correlation_analysis/utilities/helper_functions.R")
```

### For Validation and Comparison
```r
# Use legacy functions to validate new implementations
source("R/correlation_bounds_core.R")
legacy_result <- max_corr_bound(marginals)

# Compare with new modular implementation
source("R/ordinal_correlation_analysis/1_bivariate_ordcats_correlation/1_rmin_rmax_rhat/1_monte_carlo_simulation/core_bounds_functions.R")
new_result <- compute_max_bound(marginals)

# Verify consistency
all.equal(legacy_result, new_result)
```

### Data Locations
- **Raw data**: `R/ordinal_correlation_analysis/data/raw/MCsim.rds`
- **Processed BES data**: `R/ordinal_correlation_analysis/data/processed/`
- **Bootstrap results**: `R/ordinal_correlation_analysis/data/processed/bootstrap_results_sim.rds`

## Key Research Questions

1. How do marginal distributions constrain possible correlation values?
2. When do empirical bounds closely approximate theoretical bounds?
3. How do these constraints affect statistical significance testing?
4. What are the implications for survey research and categorical data analysis?

This codebase represents sophisticated statistical research combining rigorous mathematical theory with practical applications in survey data analysis. The dual architecture provides both stability (legacy code) and extensibility (modular structure) for ongoing research development.