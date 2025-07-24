# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a statistical research project exploring **theoretical bounds on Pearson's correlation coefficient for categorical variables** with fixed marginal distributions. The core insight is that correlation bounds are "elastic" - they stretch and contract based on marginal distributions, affecting both the range of possible correlations and statistical inference.

## Common Commands

### Running Analysis
```r
# Source core functions
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
# Render academic manuscripts (from project root)
quarto render paper/ElasticBounds-r_v3b.qmd
quarto render paper/ElasticBounds-r_v3b.qmd --to pdf
quarto render paper/ElasticBounds-r_v3b.qmd --to docx
```

### Required R Packages
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "gridExtra", 
                   "Matrix", "matrixcalc", "MASS", "readstata13"))
```

## Code Architecture

### Core Mathematical Functions (`R/correlation_bounds_core.R`)
- `max_corr_bound()` - Computes maximum correlation using comonotonic coupling (Fréchet–Hoeffding upper bound)
- `min_corr_bound()` - Computes minimum correlation using anti-comonotonic coupling  
- `simulate_permutation_r()` - Generates permutation distributions under null hypothesis
- `analyze_all_corr_bounds()` - Comprehensive analysis function for datasets

### Visualization System (`R/correlation-bounds-visualization.R`)
- `plot_permutation_distribution()` - Shows empirical distributions with theoretical bounds
- `plot_bounds_summary()` - Summary plots comparing bounds and confidence intervals
- `plot_range_comparison()` - Comparison plots across scenarios
- `plot_significance_comparison()` - Compares different significance testing approaches

### Example Framework (`R/correlation_bounds_examples.R`)
- `basic_demonstration()` - Explores uniform vs extreme distributions
- `range_comparison_demonstration()` - Compares theoretical vs empirical ranges
- `run_all_examples()` - Executes comprehensive example suite

### Applied Analysis (`R/correlation_bounds_bes.R`)
- Functions for analyzing British Election Study (BES) 2019 data
- Real-world validation of theoretical bounds with survey data

## Development Patterns

### File Organization
- **Main implementations**: `R/` directory contains production code
- **Experimental work**: `R/scratchWork/` for development versions
- **Academic output**: `paper/` directory for Quarto manuscripts
- **Research notes**: `Notes/` and `Meetings/` for documentation

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

## Manuscript Generation

The project uses Quarto for academic manuscript generation with:
- PDF output via pdflatex with mathematical theorem environments
- Bibliography management with APA style (`apa.csl`)
- Multiple output formats (PDF, DOCX, HTML) configured in `quarto.yml`
- Version control of manuscripts with descriptive commit messages

## Key Research Questions

1. How do marginal distributions constrain possible correlation values?
2. When do empirical bounds closely approximate theoretical bounds?
3. How do these constraints affect statistical significance testing?
4. What are the implications for survey research and categorical data analysis?

This codebase represents sophisticated statistical research combining rigorous mathematical theory with practical applications in survey data analysis.