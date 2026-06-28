# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a statistical research project exploring **theoretical bounds on Pearson's correlation coefficient for categorical variables** with fixed marginal distributions. The core insight is that correlation bounds are "elastic" - they stretch and contract based on marginal distributions, affecting both the range of possible correlations and statistical inference.

## Common Commands

### Run the full pipeline
```r
# From the project root (where ElasticBoundsPearson-s_r.Rproj lives):
source("R/00_run_pipeline.R")
```
```bash
# Equivalently, from a shell at the project root:
Rscript R/00_run_pipeline.R
```
This runs `R/01_bes_compute_bounds.R` through `R/07_permutation_zscore.R` in order and regenerates everything in `output/data/` and `output/figures/` that the manuscript pulls in.

### Run a single pipeline stage
Each numbered script in `R/` is self-contained (each redefines whatever helper functions it needs rather than sourcing a shared library), so any stage can be run on its own:
```r
source("R/02_mc_simulation.R")       # just the Monte Carlo simulation
source("R/05_hypothesis_testing.R")  # just the permutation-test figure
```

### Render the manuscript
```bash
# From the project root
quarto render paper/ElasticBounds-r_v6b.qmd --to pdf
quarto render paper/ElasticBounds-r_v6b.qmd --to docx
```
The manuscript does not execute R code inline — it pulls in pre-generated figures from `output/figures/*.pdf` via `knitr::include_graphics()`. Run `R/00_run_pipeline.R` first (or after changing any analysis script) so the figures are current before rendering.

### Required R Packages
```r
install.packages(c("haven", "dplyr", "tidyr", "ggplot2", "patchwork", "scales"))
```
(`scales` ships as a transitive dependency of `ggplot2` but is called directly —
`scales::pseudo_log_trans()` — in `03_exploratory_plots.R`/`04_additional_analyses.R`.)

## Code Architecture

### The pipeline (`R/`)

The live, current analysis pipeline is a flat set of numbered scripts directly under `R/`. Numbering reflects execution order; each stage reads/writes the top-level `output/` directory.

| Script | Purpose | Reads | Writes |
|---|---|---|---|
| `00_run_pipeline.R` | Master script — sources 01 through 05 in order | — | — |
| `01_bes_compute_bounds.R` | Computes r_min, r_max, and constraint metrics (C1, C2, entropy, TV, Bhattacharyya, overlap, asymmetry, skew) for all 7,503 BES 2019 variable pairs | `R/data/raw/bes2019_pairs.dta` | `output/data/bes_bounds.rds` |
| `02_mc_simulation.R` | Monte Carlo simulation over random Dirichlet marginals for K ∈ {4,5,6,7,10,11}, computing the same bound/constraint metrics | — | `output/data/mc_bounds.rds` |
| `03_exploratory_plots.R` | Core exploratory figures (r_min/r_max cloud, C1 distributions, rescaling comparisons) | `output/data/*.rds` | `output/figures/01_*.pdf`–`09_*.pdf` |
| `04_additional_analyses.R` | "Spur" structure analysis + asymmetry-magnitude-vs-K analysis | `output/data/*.rds` | `output/figures/10_*.pdf`–`14_*.pdf` |
| `05_hypothesis_testing.R` | Permutation-test vs. t-test comparison figure for the manuscript's Hypothesis Testing section | — | `output/figures/15_*.pdf` |
| `06_correlation_adjustment.R` | Kendall's tau-b bounds (tau_min/tau_max) for all BES pairs; five candidate adjustments to r_obs (raw, U1_conservative, tight_linear, split_warrens, U3) per Warrens (2013)'s correction-for-chance/correction-for-maximum-value framework; synthetic + real-BES-data PSD (invertibility) sweeps for Pearson r and Kendall's tau matrices under each adjustment | `R/data/raw/bes2019_pairs.dta` | `output/data/bes_tau_bounds.rds`, `output/data/adjustment_psd_synthetic.rds`, `output/data/adjustment_psd_bes_real.rds` |
| `07_permutation_zscore.R` | Permutation-null z-scores for all BES pairs: `z_raw` (closed-form, r_obs standardized by the exact 1/sqrt(N-1) permutation-null sd — shape-invariant, depends only on N) and `z_split` (simulated, mean-centered z-score on split_warrens — folds the bound-tightness correction and the sampling-noise correction into one statistic) | `R/data/raw/bes2019_pairs.dta` | `output/data/bes_perm_zscores.rds` |

Note: prose in the "Reads" / "Writes" columns omits the implicit project-root working directory; all paths in these scripts are relative to the repo root (run via `Rscript` from the root, or in an R session opened at the `.Rproj`).

### In-development modules (not yet wired into the pipeline)

| Script | Purpose |
|---|---|
| `02a_mc_sampling_uniform_crosstab.R` | Uniform sampling of J×K cross-tab tables with fixed total N (stars-and-bars method) |
| `02b_mc_sampling_crosstab_given_r.R` | Builds on 02a: samples J×K tables conditioned on a target Pearson's r via simulated annealing + Metropolis swap moves |
| `02c_mc_sampling_margin_fixed_polytope.R` | Margin-fixed (not just N-fixed) uniform sampling over the feasible polytope, via a 2x2 corner-swap Metropolis chain with a variable-magnitude ("bigger swap step") proposal; provides `percentile_in_polytope()` and a polytope-mean/sd-based `polytope_zscore()` — the "Version B" companion to `07_permutation_zscore.R`'s closed-form `z_raw`. Includes built-in validation against exact enumeration and a real-BES-pairs demo. |

These implement the "K×J cell-allocation-with-replacement" algorithm for the expanded Monte Carlo design (square + rectangular pairs, K = 3–11). They are standalone and not yet called by `02_mc_simulation.R` or the master pipeline — run them directly with `source("R/02a_mc_sampling_uniform_crosstab.R")` etc. once ready to integrate.

### Data (`R/data/raw/`)
- `bes2019_pairs.dta` — British Election Study 2019 pairwise variable data, the only data dependency the live pipeline needs.

### Archived / superseded code (`R/_archive/`)
Everything below was explored at some point but is **not** part of the live pipeline. Kept for reference/history, not for use in new work:
- `core/` — an earlier reusable function library (`correlation_bounds_core.R`, `correlation_bounds_visualization.R`, `correlation_bounds_bes.R`, `correlation_matrix_test.R`). The numbered pipeline scripts redefine the bound functions they need inline rather than sourcing this.
- `ordinal_correlation_analysis/` — a deep nested modular architecture (begun, never finished — large parts of it reference files that were never written). Superseded by the flat numbered pipeline above.
- `examples_legacy/` — older worked-example/demo scripts, plus one early draft of the crosstab-given-r sampler (the current version lives at `R/02b_mc_sampling_crosstab_given_r.R`).
- `scratch/` — informal dev/prototype scripts.
- `scripts_demos/`, `scripts_tests/` — leftovers from an earlier `here()`-path-conversion effort.
- `WHAT_EACH_FILE_DOES.md`, `r_code_reorganization_notes.md`, `correlation-bounds-readme.md` — superseded planning/orientation docs from earlier reorganizations.

### Output (`output/`)
- `output/data/*.rds` — computed bounds data frames (BES + Monte Carlo)
- `output/figures/*.pdf` — all manuscript figures, numbered to match the section/script that produces them

### Visualization conventions
Several quantities in this project are bounded (r_min ∈ [-1,0), r_max ∈ (0,1], abs_asym ≥ 0) and pile up
hard against those bounds, so a plain linear axis often overplots the most interesting region. When adding
or editing a figure, check the actual distribution first and pick a transform that fits how the variable
is bounded, rather than defaulting to a raw linear scale:
- **Variables strictly bounded away from a hard limit** (e.g. r_min is never exactly -1, r_max is never
  exactly 1 in the BES/MC data) → log-transform the *distance to the bound* (`1 + r_min`, `1 - r_max`), as
  in Figures 1-2 of `03_exploratory_plots.R`. A plain `log10()` is safe here because the distance is always
  strictly positive.
- **Variables that can legitimately equal their bound** (e.g. abs_asym can be exactly 0) → do not use a
  plain log (undefined/-Inf at 0); use `scales::pseudo_log_trans()` instead, as in Figure 13 of
  `04_additional_analyses.R`. It behaves linearly near 0 and logarithmically beyond a chosen `sigma`, so it
  compresses a long right tail without blowing up at the boundary.
- Always update the subtitle/axis labels to say what's plotted (e.g. "distance from bound, log scale") —
  don't silently change a figure's axes without explaining the new geometry, since e.g. a symmetry
  reference line's "above/below" interpretation can flip under one of these transforms.

## Mathematical Framework
The project implements Fréchet–Hoeffding bounds adapted for categorical variables:
- **Comonotonic coupling**: Both variables sorted in same direction (maximum correlation)
- **Anti-comonotonic coupling**: Variables sorted in opposite directions (minimum correlation)
- **Permutation inference**: Random permutation preserves marginals for null hypothesis testing

## Data Structures
Functions typically work with:
- Marginal distributions as probability vectors or count vectors
- Results returned as data frames with one row per variable pair (BES) or per simulated pair (Monte Carlo)
- Constraint-severity metrics: C1 = (|r_min| + r_max)/2, C2 = (r_min² + r_max²)/2
- Marginal-shape measures: Shannon entropy, total variation distance vs. reversed marginal, Bhattacharyya coefficient, overlap coefficient

## Significance Testing Innovation
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
- The current/live manuscript is `paper/ElasticBounds-r_v6b.qmd`; the `paper/` directory contains many earlier versioned drafts (v1–v6a and variants) kept for history

## Key Research Questions

1. How do marginal distributions constrain possible correlation values?
2. When do empirical bounds closely approximate theoretical bounds?
3. How do these constraints affect statistical significance testing?
4. What are the implications for survey research and categorical data analysis?

This codebase represents statistical research combining rigorous mathematical theory with practical applications in survey data analysis. The live pipeline (`R/00_run_pipeline.R` plus the numbered scripts in `R/`) is intentionally flat and self-contained so each stage can be read, run, and edited independently.
