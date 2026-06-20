# What Each R File Actually Does (Plain English)

**Ignore the `ordinal_correlation_analysis/` folder for now - it's broken. Use these files instead:**

## Core Functions (The Mathematical Engine)

### `R/correlation_bounds_core.R` ⭐ **MOST IMPORTANT**
**What it does:** Calculates the theoretical minimum and maximum correlation possible between two categorical variables when you know their marginal distributions.

**Key functions:**
- `max_corr_bound()` - "What's the highest correlation I could possibly get?"
- `min_corr_bound()` - "What's the lowest correlation I could possibly get?"
- `simulate_permutation_r()` - "What correlations would I get by random chance?"

**Why it matters:** This is your core research contribution - showing that correlation "bounds" are elastic and depend on marginal distributions.

### `R/correlation-bounds-visualization.R` ⭐ **FOR MAKING PLOTS**
**What it does:** Creates all the pretty plots and figures for your papers.

**Key functions:**
- `plot_permutation_distribution()` - Shows histogram of permutation correlations with bounds
- `plot_bounds_summary()` - Summary plots comparing different scenarios
- `plot_range_comparison()` - Compares correlation ranges across conditions

**Use when:** You want to visualize your results for papers or presentations.

## Examples and Demonstrations

### `R/correlation_bounds_examples.R` ⭐ **START HERE FOR LEARNING**
**What it does:** Runs comprehensive examples that demonstrate the theory with concrete numbers.

**Key functions:**
- `run_all_examples()` - Runs everything, shows you how the theory works
- `basic_demonstration()` - Simple example comparing uniform vs skewed distributions
- `range_comparison_demonstration()` - Shows how different marginals affect correlation ranges

**Use when:** You want to understand what your research actually shows, or create examples for teaching.

### `R/correlation_bounds_demo.R`
**What it does:** Interactive demonstrations, probably similar to examples but more focused.

## Real Data Analysis

### `R/correlation_bounds_bes.R` ⭐ **FOR REAL DATA**
**What it does:** Functions specifically designed to analyze the British Election Study (BES) 2019 dataset.

**Key functions:**
- Functions to load BES data
- Apply correlation bounds analysis to real survey data
- Validate your theory with actual data

### `R/correlation_bounds_bes_example.R` ⭐ **BES WALKTHROUGH**  
**What it does:** Step-by-step example of analyzing BES data, probably the script you run to generate your paper's empirical results.

**Use when:** You want to replicate your BES analysis or understand how the theory applies to real survey data.

## Advanced Features

### `R/correlation_bounds_simulation.R`
**What it does:** Monte Carlo simulations to test your theory under different conditions.

**Use when:** You want to test robustness or generate simulation studies for your papers.

### `R/correlation_matrix_example.R` & `R/correlation_matrix_test.R`
**What it does:** Tests what happens when you try to build full correlation matrices using your bounds.

**The problem:** Individual pairwise bounds might not be compatible when building a full matrix.
**Use when:** You're working on the matrix-level implications of your theory.

## Your Research in Simple Terms

**What you discovered:** 
- Correlation coefficients between categorical variables aren't just affected by the relationship between the variables
- They're also constrained by the marginal distributions (how many people fall into each category)
- These constraints are "elastic" - they stretch and contract based on the marginals
- This affects statistical significance testing and interpretation

**Why it matters:**
- Survey researchers need to account for these bounds when interpreting correlations
- Traditional significance tests might be misleading
- Your bounds provide the "reference frame" for understanding what correlation values actually mean

## Recommended Workflow

1. **Start with:** `SIMPLE_WORKING_EXAMPLE.R` (I just created this for you)
2. **Learn the theory:** Run `R/correlation_bounds_examples.R` functions
3. **See real applications:** Use `R/correlation_bounds_bes_example.R` 
4. **Make plots:** Use functions from `R/correlation-bounds-visualization.R`
5. **Write papers:** Generate figures and results from the above

**Ignore for now:** The entire `ordinal_correlation_analysis/` folder - it's causing more problems than it solves!