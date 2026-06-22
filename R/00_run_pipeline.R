# =============================================================================
# 00_run_pipeline.R  --  MASTER SCRIPT
#
# Runs the full ElasticBounds analysis pipeline end-to-end, in order:
#
#   01_bes_compute_bounds.R    Loads BES 2019 data, computes r_min/r_max and
#                              constraint metrics for all 7,503 variable pairs.
#                              -> output/data/bes_bounds.rds
#
#   02_mc_simulation.R         Monte Carlo simulation over random marginal
#                              distributions (K = 3,4,5,6,7,10,11), computing the
#                              same bound/constraint metrics.
#                              -> output/data/mc_bounds.rds
#
#   03_exploratory_plots.R     Core exploratory figures from the BES + MC data
#                              (cloud plots, C1 distributions, rescaling, etc).
#                              -> output/figures/01_*.pdf ... 09_*.pdf
#
#   04_additional_analyses.R   "Spur" structure analysis + asymmetry-vs-K
#                              analysis.
#                              -> output/figures/10_*.pdf ... 14_*.pdf
#
#   05_hypothesis_testing.R    Permutation-test vs t-test comparison figure
#                              for the manuscript's Hypothesis Testing section.
#                              -> output/figures/15_*.pdf
#
# All five scripts are also self-contained and can be run individually
# (each one redefines whatever helper functions it needs, by design, so you
# can read/run/edit any single stage without chasing dependencies across
# files). Running this master script just runs all five in the right order
# and regenerates everything the manuscript (paper/ElasticBounds-r_v6b.qmd)
# pulls in via knitr::include_graphics().
#
# NOT included in this auto-run sequence:
#   02a_mc_sampling_uniform_crosstab.R   } Standalone, in-development modules
#   02b_mc_sampling_crosstab_given_r.R   } for the K x J cell-allocation MC
#                                         } sampler (uniform-table and
#                                         } target-r sampling). Not yet wired
#                                         } into 02_mc_simulation.R's output.
#                                         } Run these directly when ready:
#                                         }   source("R/02a_mc_sampling_uniform_crosstab.R")
#                                         }   source("R/02b_mc_sampling_crosstab_given_r.R")
#
# USAGE:
#   From the project root (where ElasticBoundsPearson-s_r.Rproj lives):
#     Rscript R/00_run_pipeline.R
#   ...or from an R/RStudio session already opened at the project root:
#     source("R/00_run_pipeline.R")
# =============================================================================

if (!file.exists("ElasticBoundsPearson-s_r.Rproj")) {
  stop(
    "00_run_pipeline.R must be run with the working directory set to the ",
    "project root (where ElasticBoundsPearson-s_r.Rproj lives), e.g.:\n",
    "  Rscript R/00_run_pipeline.R\n",
    "Current working directory is: ", getwd()
  )
}

pipeline_steps <- c(
  "R/01_bes_compute_bounds.R",
  "R/02_mc_simulation.R",
  "R/03_exploratory_plots.R",
  "R/04_additional_analyses.R",
  "R/05_hypothesis_testing.R"
)

dir.create("output/data",    recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

cat("=============================================================\n")
cat(" ElasticBounds pipeline -- running", length(pipeline_steps), "steps\n")
cat("=============================================================\n\n")

pipeline_start <- Sys.time()

for (step in pipeline_steps) {
  cat("-------------------------------------------------------------\n")
  cat("STEP:", step, "\n")
  cat("-------------------------------------------------------------\n")
  step_start <- Sys.time()
  source(step, echo = FALSE)
  step_time <- round(difftime(Sys.time(), step_start, units = "secs"), 1)
  cat("\n>> Finished", step, "in", step_time, "sec\n\n")
}

total_time <- round(difftime(Sys.time(), pipeline_start, units = "secs"), 1)
cat("=============================================================\n")
cat(" Pipeline complete in", total_time, "sec\n")
cat(" Data:    output/data/bes_bounds.rds, output/data/mc_bounds.rds\n")
cat(" Figures: output/figures/*.pdf (", length(list.files("output/figures", pattern = "\\.pdf$")), "files)\n")
cat("=============================================================\n")
