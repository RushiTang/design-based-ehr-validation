# Creates non-identifiable synthetic files solely to test the input schema and
# workflow. The values do not represent an EHR cohort and must not be used in
# a manuscript, presentation, or clinical decision.

set.seed(20260829)

dir.create("input", showWarnings = FALSE, recursive = TRUE)
n <- 600L
patient_id <- sprintf("SYN%04d", seq_len(n))
diagnosis_date <- as.Date("2017-01-01") + sample.int(365 * 4, n, replace = TRUE)
predicted_risk <- pmin(pmax(rbeta(n, 1.4, 8), 0.001), 0.999)
algo_positive <- as.integer(predicted_risk >= 0.18)
predicted_recurrence_time <- ifelse(
  algo_positive == 1L,
  pmax(0.1, rexp(n, rate = 1 / 34)),
  NA_real_
)

event_3yr <- rbinom(n, 1, plogis(-3.0 + 3.2 * predicted_risk))
event_4yr <- pmax(event_3yr, rbinom(n, 1, plogis(-3.4 + 3.1 * predicted_risk)))
event_5yr <- pmax(event_4yr, rbinom(n, 1, plogis(-3.7 + 3.0 * predicted_risk)))

analytic_cohort <- data.frame(
  patient_id = patient_id,
  diagnosis_date = format(diagnosis_date, "%Y-%m-%d"),
  true_recur_3yr = event_3yr,
  true_recur_4yr = event_4yr,
  true_recur_5yr = event_5yr,
  evaluable_3yr_primary = TRUE,
  evaluable_4yr_primary = TRUE,
  evaluable_5yr_primary = TRUE,
  algorithm_eligible = TRUE,
  exclude = FALSE
)

algorithm_predictions <- data.frame(
  patient_id = patient_id,
  algo_positive = algo_positive,
  predicted_risk = predicted_risk,
  predicted_recurrence_time = predicted_recurrence_time
)

write.csv(analytic_cohort, "input/analytic_cohort.csv", row.names = FALSE)
write.csv(algorithm_predictions, "input/algorithm_predictions.csv", row.names = FALSE)
message("Wrote synthetic example inputs to input/. These files are ignored by Git.")
