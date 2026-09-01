# Simulation Study Code

This folder reproduces the primary simulation study for the manuscript. It implements the final design-based adaptive composite estimator, the design-based SRS comparator, metric-specific logit-Wald confidence intervals, and the outcome-dependent verification sensitivity analysis.

## Requirements

Use R 4.5.1 or later with these packages installed:

```r
install.packages(c(
  "dplyr", "tidyr", "readr", "ggplot2", "patchwork",
  "forcats", "scales", "svglite", "ragg"
))
```

## Run the simulation

From this folder, run the complete simulation study with 2,000 Monte Carlo replicates per scenario:

```sh
Rscript run_simulation.R 2000 8 outputs
Rscript summarize_results.R outputs
Rscript make_manuscript_figures.R outputs figures
```

The arguments to `run_simulation.R` are, in order: the number of replicates, the number of worker processes, and the output directory. The manuscript uses 2,000 replicates, eight workers, and 324 total settings: 162 primary probability-sampling scenarios plus 162 matched outcome-dependent-verification stress-test settings. Runtime depends on the computing environment. A small implementation check can be run first:

```sh
Rscript run_simulation.R 20 2 smoke_outputs
Rscript summarize_results.R smoke_outputs
```

## Outputs

`run_simulation.R` writes the scenario grid, scenario-level performance summaries, implementation checks, and one raw `.rds` file per scenario under `outputs/raw/`. `summarize_results.R` produces collapsed summaries, efficiency ratios, coverage diagnostics, and the matched sequential-versus-composite comparison. `make_manuscript_figures.R` writes publication-format PDF, SVG, TIFF, and PNG figures plus source-data CSV files.

The output directories are excluded from version control because the complete simulation results are large. The simulation is fully data-generating and does not require patient-level data.

## Reproducibility note

Set an R random seed before running if an exact rerun from a fixed seed is desired. Parallel execution is used on Unix-like systems; on Windows, the simulation runs serially to preserve compatibility.
