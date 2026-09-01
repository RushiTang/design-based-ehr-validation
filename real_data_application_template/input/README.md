# Input data contract

Do not place identifiable data in a public repository. The two CSV files below
must contain one row per patient and a unique shared patient key. Column names
are case-insensitive and are standardized to snake case by the scripts.

## `analytic_cohort.csv`

This file provides chart-review gold-standard status. The supplied reproduction
workflow requires the endpoint to be available for the complete analytic cohort,
solely to permit repeated sampling and comparison with the known finite-
population truth. It is therefore intended for tutorial-style validation audits,
not direct deployment in an ordinarily partially verified cohort.

| Column | Required | Format | Meaning |
| --- | --- | --- | --- |
| `patient_id` | Yes | String, unique, nonmissing | De-identified patient key used to merge the two files. |
| `diagnosis_date` | Yes | `YYYY-MM-DD` | Baseline date used to define fixed endpoints. |
| `true_recur_3yr`, `true_recur_4yr`, `true_recur_5yr` | For each analyzed endpoint | `0` or `1` | Gold-standard event by the specified endpoint. Use the naming pattern `true_recur_<N>yr` for other endpoints. |
| `evaluable_3yr_primary`, `evaluable_4yr_primary`, `evaluable_5yr_primary` | Recommended | `TRUE`/`FALSE` or `1`/`0` | Indicates adequate follow-up or an observed event to classify the endpoint. Use the naming pattern `evaluable_<N>yr_primary`. |
| `algorithm_eligible` | No | `TRUE`/`FALSE` or `1`/`0` | Optional cohort eligibility flag; missing means all rows are eligible. |
| `exclude` | No | `TRUE`/`FALSE` or `1`/`0` | Optional exclusion flag; missing means no rows are excluded. |

The scripts also support deriving endpoint labels from `recurrence_by_pull`,
`recurrence_date_by_pull`, and `followup_years_primary`, but supplying the
prespecified binary endpoint columns is clearer and is recommended for reuse.

## `algorithm_predictions.csv`

| Column | Required | Format | Meaning |
| --- | --- | --- | --- |
| `patient_id` | Yes | String, unique, nonmissing | Same de-identified key as `analytic_cohort.csv`. |
| `algo_positive` | Yes | `0` or `1` | Binary algorithm result used to define the all-reviewed positive group and the sampled negative group. |
| `predicted_risk` | Yes | Numeric, preferably in `[0,1]` | Continuous predicted probability or a monotone risk score used only to form the prespecified risk strata. |
| `predicted_recurrence_time` | Recommended for multiple fixed endpoints | Numeric | Predicted time from baseline. The unit is set by `algo_pred_time_unit` in `stage1_validation_analysis.R`; the supplied template uses months. |
| `predicted_recurrence_date` | Alternative | `YYYY-MM-DD` | Optional predicted event date. If supplied, it takes precedence over predicted time. |

For endpoint-specific binary classifications, an algorithm-positive patient is
considered positive at endpoint `N` only when the predicted time/date lies from
baseline through that endpoint. Patients originally labeled algorithm-negative
remain negative at every endpoint.

## Quality checks before analysis

1. Ensure one unique `patient_id` per file and resolve duplicates before merging.
2. Confirm that patient IDs are de-identified and carry no direct identifiers.
3. Verify the endpoint definition, administrative censoring, and evaluability
   rule with the clinical review protocol.
4. Confirm that `algo_positive` implements the prespecified decision threshold.
5. Confirm that the risk score was produced before chart review and was not
   modified using gold-standard outcomes.
