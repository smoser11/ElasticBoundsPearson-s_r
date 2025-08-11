# Complete R Project Reorganization for Ordinal Correlation Analysis

## Directory Structure with File Mapping

```
ordinal_correlation_analysis/
├── 1_bivariate_ordcats_correlation/
│   ├── 1_rmin_rmax_rhat/
│   │   ├── 1_monte_carlo_simulation/
│   │   │   ├── core_bounds_functions.R          # From: correlation_bounds_core.R
│   │   │   ├── simulation_studies.R             # From: correlation_bounds_simulation.R
│   │   │   └── validation_tests.R               # From: correlation_bounds_examples.R
│   │   ├── 2_bes_illustrative_example/
│   │   │   ├── bes_data_analysis.R              # From: correlation_bounds_bes.R
│   │   │   ├── example_demonstrations.R         # From: correlation_bounds_bes_example.R
│   │   │   └── case_studies.R                   # From: correlation_bounds_demo.R
│   │   └── bounds_analysis_main.R               # Integration module
│   ├── 2_asymmetry_analysis/
│   │   ├── asymmetry_measures.R                 # New synthesis
│   │   ├── total_variation_functions.R          # New functionality
│   │   └── asymmetry_visualization.R            # From: correlation-bounds-visualization.R
│   ├── 3_visualization/
│   │   ├── bounds_visualization.R               # From: correlation-bounds-visualization.R
│   │   ├── distribution_plots.R                 # From: correlation_bounds_bes.R (plotting functions)
│   │   └── interactive_plots.R                  # Enhanced from examples
│   └── bivariate_main.R                         # Module coordinator
├── 2_correlation_matrices/
│   ├── 1_matrix_properties/
│   │   ├── invertibility_tests.R                # From: correlation_matrix_test.R
│   │   ├── condition_number_analysis.R          # From: correlation_matrix_test.R
│   │   ├── psd_validation.R                     # From: correlation_matrix_test.R
│   │   └── matrix_diagnostics.R                 # Synthesis of matrix functions
│   ├── 2_matrix_construction/
│   │   ├── matrix_builder.R                     # From: correlation_matrix_test.R
│   │   ├── random_trials.R                      # From: correlation_matrix_test.R
│   │   └── validation_framework.R               # From: correlation_matrix_example.R
│   └── matrices_main.R                          # Module coordinator
├── 3_fixes_and_rescaling/
│   ├── 1_simple_rescaling/
│   │   ├── linear_rescaling.R                   # From: correlation_bounds_bes.R
│   │   └── rescaling_validation.R               # From: correlation_bounds_demo.R
│   ├── 2_advanced_rescaling/
│   │   ├── sophisticated_methods.R              # New advanced approaches
│   │   ├── preserving_properties.R              # Matrix property preservation
│   │   └── alternative_approaches.R             # Future extensions
│   └── rescaling_main.R                         # Module coordinator
├── utilities/
│   ├── helper_functions.R                       # Common utilities
│   ├── data_loading.R                           # Data handling functions
│   └── plotting_themes.R                        # Consistent plot styling
├── data/
│   ├── raw/                                     # Original BES data
│   └── processed/                               # Cleaned/transformed data
├── output/
│   ├── figures/                                 # Generated plots
│   ├── tables/                                  # Analysis results
│   └── reports/                                 # Full analysis reports
└── main_analysis.R                              # Master analysis script
```

