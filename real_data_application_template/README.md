# Reproducible template for the real-data validation application

This folder contains the final composite-estimator implementation used for the
real-data repeated-sampling evaluation. It intentionally contains no EHR data,
chart-review records, patient identifiers, selected-chart lists, or site-specific
preprocessing code. To reproduce the tutorial-style repeated-sampling analysis,
the analyst must supply a complete-gold-standard cohort from their own setting.

## What this template does

The workflow evaluates a binary algorithm in a fixed cohort when all
algorithm-positive patients are reviewed and algorithm-negative patients are
sampled. It implements:

- SRSWOR of algorithm-negative patients;
- an optional Phase 2 risk-stratified Neyman allocation;
- the adaptive composite estimator when the prespecified variance-ratio trigger
  is met, and pooled SRS HT estimation otherwise; and
- design-based logit-Wald CIs for sensitivity, specificity, NPV, and accuracy,
  together with a Wilson CI for PPV.

When the adaptive path is evaluated, the allocation reserves one Phase 2
review for each risk stratum that still contains unsampled algorithm-negative
patients. This maintains the positive conditional inclusion probabilities
needed by the sequential component of the composite estimator.

`run_multi_endpoint_validation.R` can draw one validation sample or multiple
repeated samples from the supplied complete-gold-standard cohort. The latter is
an additional methodological audit that is only possible when gold-standard
outcomes are known for every cohort member, as in the tutorial application.
This release is not an operational workflow for an ordinary partially verified
cohort in which outcomes are unavailable until after chart selection.

## Private input files

Place two **de-identified** CSV files in `input/`:

1. `analytic_cohort.csv`: chart-review truth and endpoint evaluability;
2. `algorithm_predictions.csv`: the algorithm binary output and predicted risk.

Their required schema and accepted formats are described in
[`input/README.md`](input/README.md). The repository ignores these two filenames
so that private data are not accidentally committed. Use environment variables
`TRUTH_CSV` and `ALGORITHM_CSV` if the files must remain outside this folder.

For a structural smoke test, create synthetic example inputs:

```bash
Rscript generate_synthetic_example.R
```

The generated rows are simulated and must never be used for scientific results.

## Run one sampling draw

From this folder, with a complete-gold-standard analytic cohort in `input/`:

```bash
VALIDATION_ENDPOINT_YEARS=3 \
VALIDATION_M1_VALUES=75 \
VALIDATION_REPLICATES=1 \
VALIDATION_WORKERS=1 \
MULTI_ENDPOINT_OUTPUT_ROOT=outputs_single_run \
Rscript run_multi_endpoint_validation.R
```

Before running, set the user-modifiable settings near the top of
`stage1_validation_analysis.R`, especially the endpoint definition, the total
negative-review precision target, the number of risk strata, and the adaptive
trigger threshold. The algorithm-positive/negative label must correspond to the
prespecified decision rule being validated.

## Repeated-sampling tutorial audit

The supplied analytic cohort must have gold-standard outcomes for every eligible
patient. The following is a small smoke test, not a target analysis:

```bash
VALIDATION_ENDPOINT_YEARS=3,4,5 \
VALIDATION_M1_VALUES=50,75,100 \
VALIDATION_REPLICATES=2 \
VALIDATION_WORKERS=1 \
MULTI_ENDPOINT_OUTPUT_ROOT=outputs_smoke_test \
Rscript run_multi_endpoint_validation.R
```

To make the formatted repeated-sampling table, set the same output location:

```bash
MULTI_ENDPOINT_OUTPUT_ROOT=outputs_smoke_test \
RESULTS_TABLE_M1=75 \
Rscript make_repeated_sampling_results_table.R
```

## Main files

| File | Purpose |
| --- | --- |
| `stage1_validation_analysis.R` | Cohort checks, sample-size planning, Phase 1 sampling, risk stratification, and optional Phase 2 allocation. |
| `stage2_compare_designs.R` | SRS and adaptive composite estimation, variance calculation, and CIs. |
| `run_multi_endpoint_validation.R` | Reproducible endpoint and replicate orchestrator. |
| `make_repeated_sampling_results_table.R` | Summary table for a complete-cohort repeated-sampling audit. |
| `generate_synthetic_example.R` | Creates non-identifiable synthetic inputs for a smoke test. |

## Software

R (version 4.2 or later) with `tidyverse`, `readr`, `janitor`, and `lubridate`.
Record `sessionInfo()` with any substantive analysis.
