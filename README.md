# Design-based validation of EHR algorithms under partial verification

This repository contains the analysis code accompanying the manuscript
*Design-Based Validation Method of EHR Algorithms to Address Partial
Verification Bias*. It provides two self-contained components:

- `simulation_study/` reproduces the finite-population simulation study,
  including the adaptive composite estimator, design-based SRS comparator,
  confidence intervals, and outcome-dependent-verification stress test.
- `real_data_application_template/` is a portable template for applying the
  analysis to a new, de-identified cohort. It includes synthetic inputs for a
  structural smoke test but does **not** include Duke University Health System
  EHR data, chart-review data, patient identifiers, or site-specific
  preprocessing code.

## Repository layout

| Directory | Contents |
| --- | --- |
| `simulation_study/` | Simulation, summaries, and manuscript-figure code. |
| `real_data_application_template/` | Reusable two-stage validation workflow, input contract, and synthetic example. |

Each directory has its own README with required R packages, commands, inputs,
outputs, and implementation details.

## Quick start

Install R (version 4.5.1 or later recommended), then install the packages
listed in the component README that you intend to run.

To run a short simulation smoke test:

```sh
cd simulation_study
Rscript run_simulation.R 20 2 smoke_outputs
Rscript summarize_results.R smoke_outputs
Rscript make_manuscript_figures.R smoke_outputs smoke_figures
```

To test the reusable application template with non-identifiable synthetic
data:

```sh
cd real_data_application_template
Rscript generate_synthetic_example.R
VALIDATION_ENDPOINT_YEARS=3 \
VALIDATION_M1_VALUES=75 \
VALIDATION_REPLICATES=2 \
VALIDATION_WORKERS=1 \
MULTI_ENDPOINT_OUTPUT_ROOT=outputs_smoke_test \
Rscript run_multi_endpoint_validation.R
```

The repeated-sampling command is a tutorial audit requiring complete
gold-standard outcomes for all eligible patients. It is not an operational
workflow for a usual partially verified cohort.

## Data availability and privacy

The Duke University Health System data used in the application are restricted
and are not included. To use the real-data template, provide two de-identified
CSV files following the contract in
[`real_data_application_template/input/README.md`](real_data_application_template/input/README.md).
The repository ignores the expected input filenames and all derived outputs to
help prevent accidental disclosure.

## Reproducibility

The simulation study is fully reproducible without external data. For an exact
rerun, record the R version, package versions, random seed, worker count, and
simulation settings. The real-data application requires a user-supplied cohort
and a prespecified endpoint, algorithm threshold, risk-stratification rule,
and chart-review budget.

## Citation

Please cite the accompanying manuscript when using this code. A persistent
software citation and license will be added with the public release.

## License

License terms are pending approval from the copyright holders. Do not reuse or
redistribute this code until a license is added to this repository.
