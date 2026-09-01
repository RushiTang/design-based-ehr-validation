############################################################
# Phase 1 validation analysis for a user-supplied EHR algorithm
#
# Required private inputs (not distributed with this repository):
#   1) input/analytic_cohort.csv: one row per eligible patient, with the
#      chart-review gold-standard endpoint(s) and evaluability indicator(s).
#   2) input/algorithm_predictions.csv: one row per patient, with the binary
#      algorithm label and continuous predicted risk score.
#
# Read input/README.md before replacing the synthetic example files. The
# expected column names are set in the "User-modifiable settings" block below.
# This script never uploads, shares, or writes the supplied patient-level data
# outside the local output directory selected by the analyst.
#
# Outputs:
#   stage1_validation_outputs/*.csv
#
# Design summary:
# - All algorithm-positive patients are verified, so their inclusion
#   probability is 1 and their design weight is 1.
# - Algorithm-negative patients are sampled by SRSWOR in Stage 1.
#   Their inclusion probability is m1 / N_neg and their design weight is
#   N_neg / m1.
# - Naive complete-case metrics from the reviewed charts are printed for
#   transparency, but they are biased under partial verification because
#   algorithm positives and algorithm negatives have different sampling
#   fractions.
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(janitor)
  library(readr)
})

# -------------------------------------------------------------------------
# User-modifiable settings
# -------------------------------------------------------------------------
script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- "--file="
script_file_match <- script_args[startsWith(script_args, script_file_arg)]
script_path <- if (length(script_file_match) > 0) {
  sub(script_file_arg, "", script_file_match[[1]])
} else {
  NA_character_
}
script_root <- if (!is.na(script_path)) {
  dirname(normalizePath(script_path))
} else {
  normalizePath(getwd())
}
analysis_root <- normalizePath(
  Sys.getenv("VALIDATION_PROJECT_ROOT", unset = script_root),
  winslash = "/",
  mustWork = TRUE
)
setwd(analysis_root)

# Analysis endpoint. Setting VALIDATION_ENDPOINT_YEAR lets the same script be
# rerun for any prespecified binary endpoint without editing its analysis body.
endpoint_year <- as.integer(Sys.getenv("VALIDATION_ENDPOINT_YEAR", unset = "3"))
if (is.na(endpoint_year) || endpoint_year <= 0) {
  stop("VALIDATION_ENDPOINT_YEAR must be a positive integer.")
}
endpoint_label <- paste0(endpoint_year, "yr")
# Each invocation analyzes one prespecified endpoint. The multi-endpoint runner
# calls this script separately for every requested value, so do not require
# endpoint columns that the analyst did not ask to evaluate.
endpoint_years_for_reporting <- endpoint_year

validation_replicate <- as.integer(Sys.getenv("VALIDATION_REPLICATE", unset = "1"))
if (is.na(validation_replicate) || validation_replicate <= 0) {
  stop("VALIDATION_REPLICATE must be a positive integer.")
}

# Place private, de-identified files in input/ or pass their locations through
# TRUTH_CSV and ALGORITHM_CSV. Do not commit patient-level files to a public
# repository.
truth_csv <- Sys.getenv("TRUTH_CSV", unset = file.path("input", "analytic_cohort.csv"))
algorithm_output_csv <- Sys.getenv(
  "ALGORITHM_CSV",
  unset = file.path("input", "algorithm_predictions.csv")
)

output_dir <- "stage1_validation_outputs"

output_dir <- Sys.getenv("STAGE1_OUTPUT_DIR", unset = output_dir)

# Key variable names. These are cleaned with janitor::clean_names() after
# import, so use snake_case names here. Edit only this block when your input
# uses different names; see input/README.md for definitions and valid formats.
truth_id_var <- "patient_id"
truth_diagnosis_date_var <- "diagnosis_date"
truth_outcome_var <- paste0("true_recur_", endpoint_year, "yr")
truth_evaluable_var <- paste0("evaluable_", endpoint_year, "yr_primary")
truth_algorithm_eligible_var <- "algorithm_eligible"
truth_exclude_var <- "exclude"

algo_id_var <- "patient_id"
algo_binary_var <- "algo_positive"
algo_risk_var <- "predicted_risk"

# Optional timing fields. If a predicted recurrence date is present, it takes
# precedence. Otherwise, if predicted time is present, it is interpreted using
# algo_pred_time_unit. If neither is present, the script uses the binary
# algorithm flag directly as the current endpoint algorithm label and emits a
# warning.
algo_pred_date_var <- NA_character_
algo_pred_time_var <- "predicted_recurrence_time"
algo_pred_time_unit <- "months"  # one of: "months", "days", "years"

# Set this only if the endpoint labels must be derived from dated source data.
# It is not used when true_recur_<endpoint>yr columns are supplied directly.
admin_censor_date <- as.Date("2026-01-01")

seed <- as.integer(Sys.getenv("VALIDATION_SEED", unset = "20260523"))
if (is.na(seed)) {
  stop("VALIDATION_SEED must be an integer.")
}
# Stage 1 negative pilot size. The multi-run wrapper varies this to assess how
# pilot information affects the Neyman allocation and downstream performance.
m1 <- as.integer(Sys.getenv("VALIDATION_M1", unset = "100"))
if (is.na(m1) || m1 < 2) {
  stop("VALIDATION_M1 must be an integer of at least 2.")
}

# Let the script calculate total planned negative sample size m from the
# current analytic cohort. Set FALSE if you want to force m_manual.
auto_calculate_m <- TRUE
m_manual <- 250L

# Precision targets for choosing total negative sample size m.
# m is chosen among algorithm-negative patients only.
sens_halfwidth_target <- 0.10
spec_halfwidth_target <- 0.05
npv_halfwidth_target  <- 0.05
acc_halfwidth_target  <- 0.05
conf_level_for_m <- 0.95

# Use 25-chart grid to match the original planning convention.
# Set to 1L if you want the exact smallest integer m.
m_round_to <- 25L

# Placeholder; overwritten after analytic cohort is created.
m <- m_manual
m2 <- m - m1

H <- 5L
trigger_threshold <- 1.1

stop_on_duplicate_ids <- TRUE

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

clean_requested_name <- function(x) {
  if (length(x) == 0 || is.na(x) || identical(x, "")) {
    return(NA_character_)
  }

  janitor::make_clean_names(x)
}

resolve_col <- function(data, requested, candidates = character(),
                        label, required = TRUE) {
  requested_clean <- clean_requested_name(requested)
  candidate_clean <- janitor::make_clean_names(candidates)

  if (!is.na(requested_clean) && requested_clean %in% names(data)) {
    return(requested_clean)
  }

  found <- candidate_clean[candidate_clean %in% names(data)]
  if (length(found) > 0) {
    message("Using inferred column for ", label, ": ", found[[1]])
    return(found[[1]])
  }

  if (required) {
    stop("Could not find required column for ", label,
         ". Requested: ", requested,
         ". Available columns: ", paste(names(data), collapse = ", "))
  }

  NA_character_
}

standardize_patient_id <- function(x) {
  str_to_upper(str_squish(as.character(x)))
}

is_blank_text <- function(x) {
  y <- str_to_lower(str_squish(as.character(x)))
  is.na(y) | y %in% c("", "na", "n/a", "missing", "unknown", "unk",
                      "not available", ".", "-", "--")
}

parse_date_flexible <- function(x) {
  x_chr <- str_squish(as.character(x))
  x_chr[is_blank_text(x_chr)] <- NA_character_

  out <- rep(as.Date(NA), length(x_chr))

  numeric_value <- suppressWarnings(as.numeric(x_chr))
  numeric_like <- str_detect(x_chr, "^\\d+(\\.\\d+)?$")
  year_only_like <- numeric_like &
    str_detect(x_chr, "^\\d{4}$") &
    !is.na(numeric_value) &
    numeric_value >= 1900 &
    numeric_value <= 2100

  excel_serial_like <- numeric_like &
    !is.na(numeric_value) &
    !year_only_like &
    numeric_value >= 1 &
    numeric_value <= 80000

  out[excel_serial_like] <- as.Date(floor(numeric_value[excel_serial_like]),
                                    origin = "1899-12-30")

  remaining <- is.na(out) & !is.na(x_chr)
  if (any(remaining)) {
    parsed <- suppressWarnings(lubridate::parse_date_time(
      x_chr[remaining],
      orders = c("ymd", "mdy", "dmy", "Ymd", "mdY", "dmY",
                 "ymd HMS", "mdy HMS", "dmy HMS"),
      tz = "UTC"
    ))
    out[remaining] <- as.Date(parsed)
  }

  out
}

standardize_binary <- function(x) {
  if (is.logical(x)) {
    return(case_when(
      is.na(x) ~ NA_integer_,
      x ~ 1L,
      !x ~ 0L,
      TRUE ~ NA_integer_
    ))
  }

  y <- str_to_lower(str_squish(as.character(x)))
  y <- str_replace_all(y, "[^a-z0-9.\\-]+", " ")

  numeric_value <- suppressWarnings(as.numeric(y))

  case_when(
    y %in% c("yes", "y", "true", "t", "recurrent", "recurrence",
             "positive", "pos", "case") ~ 1L,
    y %in% c("no", "n", "false", "f", "non recurrent", "nonrecurrent",
             "none", "negative", "neg", "control") ~ 0L,
    !is.na(numeric_value) & numeric_value == 1 ~ 1L,
    !is.na(numeric_value) & numeric_value == 0 ~ 0L,
    TRUE ~ NA_integer_
  )
}

standardize_logical <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  y <- str_to_lower(str_squish(as.character(x)))
  case_when(
    y %in% c("true", "t", "yes", "y", "1") ~ TRUE,
    y %in% c("false", "f", "no", "n", "0") ~ FALSE,
    TRUE ~ NA
  )
}

parse_numeric_flexible <- function(x) {
  suppressWarnings(readr::parse_number(as.character(x)))
}

add_truth_endpoint_columns <- function(data, endpoint_years) {
  required_cols <- c("diagnosis_date", "recurrence_by_pull",
                     "recurrence_date_by_pull", "followup_years_primary")
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    warning(
      "Cannot derive all endpoint-specific truth columns because these columns ",
      "are missing: ", paste(missing_cols, collapse = ", "),
      ". Existing truth columns will be used when available."
    )
    return(data)
  }

  diagnosis_date <- parse_date_flexible(data$diagnosis_date)
  recurrence_by_pull <- standardize_binary(data$recurrence_by_pull)
  recurrence_date_by_pull <- parse_date_flexible(data$recurrence_date_by_pull)
  followup_years_primary <- parse_numeric_flexible(data$followup_years_primary)

  for (yr in endpoint_years) {
    endpoint_date <- diagnosis_date %m+% years(yr)
    true_col <- paste0("true_recur_", yr, "yr")
    evaluable_col <- paste0("evaluable_", yr, "yr_primary")

    true_recur <- case_when(
      recurrence_by_pull == 1L &
        !is.na(recurrence_date_by_pull) &
        !is.na(endpoint_date) &
        recurrence_date_by_pull <= endpoint_date ~ 1L,
      recurrence_by_pull == 1L &
        !is.na(recurrence_date_by_pull) &
        !is.na(endpoint_date) &
        recurrence_date_by_pull > endpoint_date ~ 0L,
      recurrence_by_pull == 0L &
        !is.na(followup_years_primary) &
        followup_years_primary >= yr ~ 0L,
      TRUE ~ NA_integer_
    )

    data[[true_col]] <- true_recur
    data[[evaluable_col]] <- !is.na(true_recur)
  }

  data
}

safe_divide <- function(numerator, denominator) {
  ifelse(is.na(denominator) | denominator == 0, NA_real_,
         numerator / denominator)
}

sample_variance <- function(x) {
  x <- as.numeric(x)
  if (sum(!is.na(x)) <= 1) {
    return(NA_real_)
  }
  stats::var(x, na.rm = TRUE)
}

sample_covariance <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) <= 1) {
    return(NA_real_)
  }
  stats::cov(as.numeric(x[ok]), as.numeric(y[ok]))
}

metric_from_totals <- function(tp, fp, fn, tn, denominator_n) {
  tibble(
    metric = c("sensitivity", "specificity", "ppv", "npv", "accuracy"),
    estimate = c(
      safe_divide(tp, tp + fn),
      safe_divide(tn, tn + fp),
      safe_divide(tp, tp + fp),
      safe_divide(tn, tn + fn),
      safe_divide(tp + tn, denominator_n)
    )
  )
}

logit_ci <- function(p_hat, var_p, level = 0.95, eps = 1e-6) {
  if (is.na(p_hat) || is.na(var_p) || var_p < 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  z <- qnorm(1 - (1 - level) / 2)
  p_for_logit <- min(max(p_hat, eps), 1 - eps)
  se_logit <- sqrt(var_p / (p_for_logit^2 * (1 - p_for_logit)^2))
  ci <- plogis(qlogis(p_for_logit) + c(-1, 1) * z * se_logit)

  c(lower = ci[[1]], upper = ci[[2]])
}

wilson_ci <- function(x, n, level = 0.95) {
  if (is.na(x) || is.na(n) || n <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  z <- qnorm(1 - (1 - level) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half_width <- z * sqrt((p * (1 - p) / n) + z^2 / (4 * n^2)) / denom

  c(lower = max(0, center - half_width),
    upper = min(1, center + half_width))
}

print_duplicate_ids <- function(data, id_col, label) {
  dupes <- data %>%
    count(.data[[id_col]], name = "n") %>%
    filter(!is.na(.data[[id_col]]), n > 1) %>%
    arrange(desc(n), .data[[id_col]])

  if (nrow(dupes) > 0) {
    cat("\nDuplicate patient IDs in ", label, "\n", sep = "")
    print(dupes, n = Inf)
  } else {
    cat("\nNo duplicate patient IDs in ", label, ".\n", sep = "")
  }

  dupes
}

allocate_neyman_constrained <- function(stratum_summary, target_m) {
  alloc <- stratum_summary %>%
    mutate(
      allocation_basis = N_h * S_h,
      allocation_basis = if_else(is.na(allocation_basis), 0, allocation_basis),
      # Conditional Stage 2 estimation requires every stratum with unsampled
      # patients to receive at least one Stage 2 review.
      m_h_min_total = pmin(N_h, m1_h + 1L)
    )

  if (sum(alloc$N_h, na.rm = TRUE) < target_m) {
    stop("Target total negative sample m exceeds N_neg.")
  }

  if (sum(alloc$m1_h, na.rm = TRUE) > target_m) {
    stop("Stage 1 negative sample already exceeds target total negative sample.")
  }

  if (sum(alloc$m_h_min_total, na.rm = TRUE) > target_m) {
    stop(
      "The Stage 2 budget is too small to give each partially sampled risk ",
      "stratum one additional review. Increase the total negative-review ",
      "budget or reduce the number of risk strata."
    )
  }

  if (sum(alloc$allocation_basis, na.rm = TRUE) <= 0) {
    alloc <- alloc %>% mutate(allocation_basis = N_h)
  }

  alloc <- alloc %>%
    mutate(
      m_h_continuous = target_m * allocation_basis / sum(allocation_basis),
      m_h_neyman_total = floor(m_h_continuous),
      m_h_neyman_total = pmax(m_h_neyman_total, m_h_min_total),
      m_h_neyman_total = pmin(m_h_neyman_total, N_h)
    )

  iteration <- 0L
  max_iterations <- target_m * 20L + 100L

  while (sum(alloc$m_h_neyman_total) < target_m) {
    iteration <- iteration + 1L
    if (iteration > max_iterations) {
      stop("Neyman allocation did not converge while adding samples.")
    }

    candidates <- which(alloc$m_h_neyman_total < alloc$N_h)
    if (length(candidates) == 0) {
      stop("No stratum has remaining capacity for Neyman allocation.")
    }

    priority <- alloc$m_h_continuous - alloc$m_h_neyman_total
    priority[-candidates] <- -Inf
    idx <- which.max(priority)
    alloc$m_h_neyman_total[[idx]] <- alloc$m_h_neyman_total[[idx]] + 1L
  }

  while (sum(alloc$m_h_neyman_total) > target_m) {
    iteration <- iteration + 1L
    if (iteration > max_iterations) {
      stop("Neyman allocation did not converge while removing samples.")
    }

    candidates <- which(alloc$m_h_neyman_total > alloc$m_h_min_total)
    if (length(candidates) == 0) {
      stop("No stratum can be reduced without dropping Stage 1 charts.")
    }

    priority <- alloc$m_h_neyman_total - alloc$m_h_continuous
    priority[-candidates] <- -Inf
    idx <- which.max(priority)
    alloc$m_h_neyman_total[[idx]] <- alloc$m_h_neyman_total[[idx]] - 1L
  }

  alloc %>%
    mutate(
      m_h_neyman_total = as.integer(m_h_neyman_total),
      m2_h = as.integer(m_h_neyman_total - m1_h)
    )
}

sample_within_strata <- function(data, allocation, n_col, seed_value) {
  set.seed(seed_value)

  data %>%
    left_join(allocation %>% select(risk_stratum, n_to_sample = all_of(n_col)),
              by = "risk_stratum") %>%
    group_by(risk_stratum) %>%
    group_modify(~ {
      n_take <- unique(.x$n_to_sample)
      n_take <- n_take[!is.na(n_take)]
      if (length(n_take) == 0 || n_take[[1]] <= 0) {
        return(.x[0, , drop = FALSE])
      }
      dplyr::slice_sample(.x, n = n_take[[1]])
    }) %>%
    ungroup() %>%
    select(-n_to_sample)
}

# -------------------------------------------------------------------------
# Step 1. Read, clean, standardize, and merge data
# -------------------------------------------------------------------------

truth_raw <- readr::read_csv(truth_csv, show_col_types = FALSE) %>%
  janitor::clean_names()

algo_raw <- readr::read_csv(algorithm_output_csv, show_col_types = FALSE) %>%
  janitor::clean_names()

# Reuse analyst-supplied binary endpoint labels whenever available. Derivation
# from recurrence dates and follow-up is only a fallback for source datasets
# that have not yet been converted to the required endpoint-specific format.
requested_truth_col <- paste0("true_recur_", endpoint_year, "yr")
if (!requested_truth_col %in% names(truth_raw)) {
  truth_raw <- add_truth_endpoint_columns(truth_raw, endpoint_years_for_reporting)
}

truth_id_col <- resolve_col(
  truth_raw,
  truth_id_var,
  c("mrn", "medical_record_number", "patient_id", "patid"),
  "chart-review patient ID"
)
truth_dx_col <- resolve_col(
  truth_raw,
  truth_diagnosis_date_var,
  c("diagnosis_date", "date_diagnosis", "date_of_diagnosis"),
  "diagnosis date"
)
truth_outcome_col <- resolve_col(
  truth_raw,
  truth_outcome_var,
  c(truth_outcome_var, "true_recur_3yr", "true_recurrence_3yr",
    "three_year_recurrence"),
  paste0("true ", endpoint_year, "-year recurrence")
)
truth_evaluable_col <- resolve_col(
  truth_raw,
  truth_evaluable_var,
  c(truth_evaluable_var, "evaluable_3yr_primary", "evaluable_3yr"),
  paste0(endpoint_year, "-year evaluability flag"),
  required = FALSE
)
truth_algorithm_eligible_col <- resolve_col(
  truth_raw,
  truth_algorithm_eligible_var,
  c("algorithm_eligible_1yr_primary", "algorithm_eligible_1yr"),
  "1-year algorithm eligibility flag",
  required = FALSE
)
truth_exclude_col <- resolve_col(
  truth_raw,
  truth_exclude_var,
  c("exclude_primary", "exclude"),
  "primary exclusion flag",
  required = FALSE
)

algo_id_col <- resolve_col(
  algo_raw,
  algo_id_var,
  c("mrn", "medical_record_number", "patient_id", "patid"),
  "algorithm patient ID"
)
algo_binary_col <- resolve_col(
  algo_raw,
  algo_binary_var,
  c("algo_positive", "algorithm_positive", "recur_pred",
    "recurrence", "predicted_recurrence", "algorithm_recurrence"),
  "algorithm binary recurrence output"
)
algo_risk_col <- resolve_col(
  algo_raw,
  algo_risk_var,
  c("pred_prob", "predicted_risk", "recur_probability", "prob",
    "risk_score", "probability"),
  "algorithm predicted risk score"
)
algo_pred_date_col <- resolve_col(
  algo_raw,
  algo_pred_date_var,
  c("predicted_recurrence_date", "pred_recurrence_date",
    "algorithm_recurrence_date", "recurrence_date_predicted"),
  "optional predicted recurrence date",
  required = FALSE
)
algo_pred_time_col <- resolve_col(
  algo_raw,
  algo_pred_time_var,
  c("predtime", "pred_time", "predicted_time_to_recurrence",
    "predicted_recurrence_time", "months_to_recurrence"),
  "optional predicted recurrence time",
  required = FALSE
)

truth <- truth_raw %>%
  mutate(
    patient_id = standardize_patient_id(.data[[truth_id_col]]),
    diagnosis_date = parse_date_flexible(.data[[truth_dx_col]]),
    true_recur_3yr = standardize_binary(.data[[truth_outcome_col]]),
    evaluable_3yr_primary = if (!is.na(truth_evaluable_col)) {
      standardize_logical(.data[[truth_evaluable_col]])
    } else {
      warning(endpoint_year, "-year evaluability column not found; ",
              "assuming all truth rows are evaluable.")
      TRUE
    },
    algorithm_eligible_1yr_primary = if (!is.na(truth_algorithm_eligible_col)) {
      standardize_logical(.data[[truth_algorithm_eligible_col]])
    } else {
      warning("1-year algorithm eligibility column not found; assuming all truth rows are eligible.")
      TRUE
    },
    exclude_primary = if (!is.na(truth_exclude_col)) {
      standardize_logical(.data[[truth_exclude_col]])
    } else {
      warning("Primary exclusion column not found; assuming no truth rows are excluded.")
      FALSE
    }
  )

algo <- algo_raw %>%
  mutate(
    patient_id = standardize_patient_id(.data[[algo_id_col]]),
    algo_binary_raw = .data[[algo_binary_col]],
    algo_positive_original = standardize_binary(algo_binary_raw),
    predicted_risk = parse_numeric_flexible(.data[[algo_risk_col]]),
    predicted_recurrence_date = if (!is.na(algo_pred_date_col)) {
      parse_date_flexible(.data[[algo_pred_date_col]])
    } else {
      as.Date(NA)
    },
    predicted_recurrence_time = if (!is.na(algo_pred_time_col)) {
      parse_numeric_flexible(.data[[algo_pred_time_col]])
    } else {
      NA_real_
    }
  )

truth_dupes <- print_duplicate_ids(truth, "patient_id", "chart-review truth file")
algo_dupes <- print_duplicate_ids(algo, "patient_id", "algorithm output file")

if (stop_on_duplicate_ids && (nrow(truth_dupes) > 0 || nrow(algo_dupes) > 0)) {
  stop("Duplicate patient IDs found. Resolve duplicates or set stop_on_duplicate_ids = FALSE.")
}

truth_ids <- truth %>% distinct(patient_id)
algo_ids <- algo %>% distinct(patient_id)

merge_counts <- tibble(
  truth_rows = nrow(truth),
  algorithm_rows = nrow(algo),
  matched_patient_ids = nrow(inner_join(truth_ids, algo_ids, by = "patient_id")),
  truth_without_algorithm = nrow(anti_join(truth_ids, algo_ids, by = "patient_id")),
  algorithm_without_truth = nrow(anti_join(algo_ids, truth_ids, by = "patient_id"))
)

cat("\nMerge counts before cohort filtering\n")
print(merge_counts)

readr::write_csv(anti_join(truth_ids, algo_ids, by = "patient_id"),
                 file.path(output_dir, "unmatched_truth_patient_ids.csv"))
readr::write_csv(anti_join(algo_ids, truth_ids, by = "patient_id"),
                 file.path(output_dir, "unmatched_algorithm_patient_ids.csv"))

merged_pre_filter <- truth %>%
  inner_join(algo, by = "patient_id", suffix = c("_truth", "_algo"))

cohort_filter_counts <- tibble(
  step = c("truth_rows", "algorithm_rows", "matched_rows_before_filter"),
  n = c(nrow(truth), nrow(algo), nrow(merged_pre_filter))
)

analytic <- merged_pre_filter %>%
  filter(
    exclude_primary == FALSE,
    algorithm_eligible_1yr_primary == TRUE,
    evaluable_3yr_primary == TRUE
  )

cohort_filter_counts <- bind_rows(
  cohort_filter_counts,
  tibble(step = "matched_rows_after_primary_analytic_filter",
         n = nrow(analytic))
)

cat("\nCohort counts before and after primary analytic filtering\n")
print(cohort_filter_counts, n = Inf)

missing_required <- analytic %>%
  filter(
    is.na(patient_id) |
      is.na(diagnosis_date) |
      is.na(true_recur_3yr) |
      is.na(algo_positive_original)
  ) %>%
  select(patient_id, diagnosis_date, true_recur_3yr, algo_positive_original)

if (nrow(missing_required) > 0) {
  print(missing_required, n = min(25, nrow(missing_required)))
  stop("Analytic rows have missing required ID/date/outcome/algorithm fields.")
}

# -------------------------------------------------------------------------
# Step 2. Define algorithm-predicted recurrence for the analysis endpoint
# -------------------------------------------------------------------------

analytic <- analytic %>%
  mutate(
    analysis_endpoint_year = endpoint_year,
    analysis_endpoint_label = endpoint_label,
    analysis_endpoint_date = diagnosis_date %m+% years(endpoint_year),
    algo_prediction_timing_source = case_when(
      !is.na(algo_pred_date_col) ~ "predicted_recurrence_date",
      !is.na(algo_pred_time_col) ~ paste0("predicted_recurrence_time_", algo_pred_time_unit),
      TRUE ~ "binary_algorithm_flag_only"
    )
  )

if (!is.na(algo_pred_date_col)) {
  analytic <- analytic %>%
    mutate(
      alg_recur_3yr = case_when(
        algo_positive_original == 1L &
          !is.na(predicted_recurrence_date) &
          predicted_recurrence_date >= diagnosis_date &
          predicted_recurrence_date <= analysis_endpoint_date ~ 1L,
        !is.na(algo_positive_original) ~ 0L,
        TRUE ~ NA_integer_
      )
    )
} else if (!is.na(algo_pred_time_col)) {
  time_threshold <- case_when(
    algo_pred_time_unit == "months" ~ 12 * endpoint_year,
    algo_pred_time_unit == "days" ~ 365.25 * endpoint_year,
    algo_pred_time_unit == "years" ~ endpoint_year,
    TRUE ~ NA_real_
  )

  if (is.na(time_threshold)) {
    stop("algo_pred_time_unit must be one of: months, days, years.")
  }

  n_positive_missing_time <- analytic %>%
    filter(algo_positive_original == 1L, is.na(predicted_recurrence_time)) %>%
    nrow()

  if (n_positive_missing_time > 0) {
    warning(n_positive_missing_time,
            " algorithm-positive rows have missing predicted recurrence time; ",
            "they will be classified as alg_recur_3yr = 0.")
  }

  analytic <- analytic %>%
    mutate(
      alg_recur_3yr = case_when(
        algo_positive_original == 1L &
          !is.na(predicted_recurrence_time) &
          predicted_recurrence_time >= 0 &
          predicted_recurrence_time <= time_threshold ~ 1L,
        !is.na(algo_positive_original) ~ 0L,
        TRUE ~ NA_integer_
      )
    )
} else {
  warning("No predicted recurrence time/date column was found. ",
          "Using the binary algorithm recurrence flag directly as the ",
          endpoint_label, " algorithm label.")
  analytic <- analytic %>%
    mutate(alg_recur_3yr = algo_positive_original)
}

if (any(is.na(analytic$alg_recur_3yr))) {
  stop("Analysis endpoint algorithm label contains missing values after standardization.")
}

analytic[[paste0("alg_recur_", endpoint_year, "yr")]] <- analytic$alg_recur_3yr

alg_difference_n <- analytic %>%
  summarise(n = sum(algo_positive_original != alg_recur_3yr, na.rm = TRUE)) %>%
  pull(n)

cat("\nRows where algo_positive_original differs from the ",
    endpoint_label, " algorithm label: ",
    alg_difference_n, "\n", sep = "")

algorithm_prediction_crosstab <- analytic %>%
  count(algo_positive_original, alg_recur_3yr, true_recur_3yr,
        name = "n") %>%
  mutate(endpoint_year = endpoint_year,
         endpoint_label = endpoint_label,
         validation_replicate = validation_replicate,
         .before = 1) %>%
  arrange(algo_positive_original, alg_recur_3yr, true_recur_3yr)

cat("\nCross-tab of original algorithm label, ",
    endpoint_label, " algorithm label, and truth\n", sep = "")
print(algorithm_prediction_crosstab, n = Inf)

# -------------------------------------------------------------------------
# Step 2b. Calculate planned total negative sample size m
# -------------------------------------------------------------------------
# m is the total number of algorithm-negative charts to review across
# Stage 1 + Stage 2. It should be calculated in the same analytic cohort used
# for the validation estimand, not in the broader raw algorithm-output cohort.

calculate_precision_table_for_m <- function(data,
                                            m1,
                                            sens_target,
                                            spec_target,
                                            npv_target,
                                            acc_target,
                                            conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  
  plan_df <- data %>%
    mutate(
      planning_prob = pmin(pmax(predicted_risk, 1e-6), 1 - 1e-6),
      sampling_positive = algo_positive_original == 1L,
      sampling_negative = algo_positive_original == 0L,
      predicted_positive_endpoint = alg_recur_3yr == 1L
    )
  
  N_plan <- nrow(plan_df)
  N_pos_plan <- sum(plan_df$sampling_positive, na.rm = TRUE)
  N_neg_plan <- sum(plan_df$sampling_negative, na.rm = TRUE)
  
  if (N_neg_plan < m1) {
    stop("N_neg_plan is smaller than m1; cannot calculate m.")
  }
  
  if (any(is.na(plan_df$planning_prob))) {
    stop("predicted_risk has missing values in analytic cohort; resolve before planning m.")
  }
  
  # Fully verified original algorithm-positive patients contribute fixed
  # expected counts. Among these patients, the endpoint-specific
  # predicted-positive label is alg_recur_3yr.
  TP_fixed_guess <- sum(
    plan_df$planning_prob[
      plan_df$sampling_positive & plan_df$predicted_positive_endpoint
    ],
    na.rm = TRUE
  )
  
  FP_fixed_guess <- sum(
    1 - plan_df$planning_prob[
      plan_df$sampling_positive & plan_df$predicted_positive_endpoint
    ],
    na.rm = TRUE
  )
  
  FN_fixed_guess <- sum(
    plan_df$planning_prob[
      plan_df$sampling_positive & !plan_df$predicted_positive_endpoint
    ],
    na.rm = TRUE
  )
  
  TN_fixed_guess <- sum(
    1 - plan_df$planning_prob[
      plan_df$sampling_positive & !plan_df$predicted_positive_endpoint
    ],
    na.rm = TRUE
  )
  
  # Uncertainty comes from sampled original algorithm-negatives.
  p_fn_neg_guess <- mean(
    plan_df$planning_prob[plan_df$sampling_negative],
    na.rm = TRUE
  )
  p_tn_neg_guess <- 1 - p_fn_neg_guess
  
  m_seq <- seq(m1, N_neg_plan, by = 1L)
  
  precision_table <- map_dfr(m_seq, function(m_candidate) {
    f <- m_candidate / N_neg_plan
    
    var_p_fn <- (1 - f) * p_fn_neg_guess * (1 - p_fn_neg_guess) /
      (m_candidate - 1)
    
    var_FN_neg <- N_neg_plan^2 * var_p_fn
    var_TN_neg <- var_FN_neg
    cov_TN_FN_neg <- -var_FN_neg
    
    FN_neg_hat <- N_neg_plan * p_fn_neg_guess
    TN_neg_hat <- N_neg_plan * p_tn_neg_guess
    
    TP_hat <- TP_fixed_guess
    FP_hat <- FP_fixed_guess
    FN_hat <- FN_fixed_guess + FN_neg_hat
    TN_hat <- TN_fixed_guess + TN_neg_hat
    
    sensitivity <- TP_hat / (TP_hat + FN_hat)
    specificity <- TN_hat / (TN_hat + FP_hat)
    npv <- TN_hat / (TN_hat + FN_hat)
    accuracy <- (TP_hat + TN_hat) / N_plan
    
    d_sens_d_fn <- TP_hat / (TP_hat + FN_hat)^2
    var_sens <- d_sens_d_fn^2 * var_FN_neg
    
    d_spec_d_tn <- FP_hat / (TN_hat + FP_hat)^2
    var_spec <- d_spec_d_tn^2 * var_TN_neg
    
    d_npv_d_tn <- FN_hat / (TN_hat + FN_hat)^2
    d_npv_d_fn <- -TN_hat / (TN_hat + FN_hat)^2
    var_npv <- d_npv_d_tn^2 * var_TN_neg +
      d_npv_d_fn^2 * var_FN_neg +
      2 * d_npv_d_tn * d_npv_d_fn * cov_TN_FN_neg
    
    var_acc <- var_TN_neg / N_plan^2
    
    tibble(
      m = m_candidate,
      N_plan = N_plan,
      N_pos_plan = N_pos_plan,
      N_neg_plan = N_neg_plan,
      TP_fixed_guess = TP_fixed_guess,
      FP_fixed_guess = FP_fixed_guess,
      FN_fixed_guess = FN_fixed_guess,
      TN_fixed_guess = TN_fixed_guess,
      p_fn_neg_guess = p_fn_neg_guess,
      p_tn_neg_guess = p_tn_neg_guess,
      sensitivity = sensitivity,
      specificity = specificity,
      npv = npv,
      accuracy = accuracy,
      sens_halfwidth = z * sqrt(var_sens),
      spec_halfwidth = z * sqrt(var_spec),
      npv_halfwidth = z * sqrt(var_npv),
      acc_halfwidth = z * sqrt(var_acc)
    )
  })
  
  precision_table
}

sample_size_precision_table <- calculate_precision_table_for_m(
  data = analytic,
  m1 = m1,
  sens_target = sens_halfwidth_target,
  spec_target = spec_halfwidth_target,
  npv_target = npv_halfwidth_target,
  acc_target = acc_halfwidth_target,
  conf = conf_level_for_m
)

first_m_for_metric <- function(table, halfwidth_col, target) {
  out <- table %>%
    filter(.data[[halfwidth_col]] <= target) %>%
    slice_head(n = 1) %>%
    pull(m)
  
  if (length(out) == 0) NA_integer_ else as.integer(out[[1]])
}

m_for_sens <- first_m_for_metric(
  sample_size_precision_table,
  "sens_halfwidth",
  sens_halfwidth_target
)

m_for_spec <- first_m_for_metric(
  sample_size_precision_table,
  "spec_halfwidth",
  spec_halfwidth_target
)

m_for_npv <- first_m_for_metric(
  sample_size_precision_table,
  "npv_halfwidth",
  npv_halfwidth_target
)

m_for_acc <- first_m_for_metric(
  sample_size_precision_table,
  "acc_halfwidth",
  acc_halfwidth_target
)

m_precision_rule_exact_row <- sample_size_precision_table %>%
  filter(
    sens_halfwidth <= sens_halfwidth_target,
    spec_halfwidth <= spec_halfwidth_target,
    npv_halfwidth <= npv_halfwidth_target,
    acc_halfwidth <= acc_halfwidth_target
  ) %>%
  slice_head(n = 1)

m_precision_rule_grid_row <- sample_size_precision_table %>%
  filter(m %% m_round_to == 0L) %>%
  filter(
    sens_halfwidth <= sens_halfwidth_target,
    spec_halfwidth <= spec_halfwidth_target,
    npv_halfwidth <= npv_halfwidth_target,
    acc_halfwidth <= acc_halfwidth_target
  ) %>%
  slice_head(n = 1)

if (nrow(m_precision_rule_exact_row) == 0) {
  stop("No m satisfies all precision targets up to N_neg.")
}

if (nrow(m_precision_rule_grid_row) == 0) {
  stop("No grid-rounded m satisfies all precision targets up to N_neg.")
}

sample_size_m_summary <- tibble(
  endpoint_year = endpoint_year,
  endpoint_label = endpoint_label,
  validation_replicate = validation_replicate,
  N_plan = m_precision_rule_exact_row$N_plan,
  N_pos_plan = m_precision_rule_exact_row$N_pos_plan,
  N_neg_plan = m_precision_rule_exact_row$N_neg_plan,
  m1 = m1,
  m_for_sens = m_for_sens,
  m_for_spec = m_for_spec,
  m_for_npv = m_for_npv,
  m_for_acc = m_for_acc,
  m_recommend_exact_by_max = max(
    c(m_for_sens, m_for_spec, m_for_npv, m_for_acc),
    na.rm = TRUE
  ),
  m_precision_rule_exact = m_precision_rule_exact_row$m,
  sens_halfwidth_at_exact_m = m_precision_rule_exact_row$sens_halfwidth,
  spec_halfwidth_at_exact_m = m_precision_rule_exact_row$spec_halfwidth,
  npv_halfwidth_at_exact_m = m_precision_rule_exact_row$npv_halfwidth,
  acc_halfwidth_at_exact_m = m_precision_rule_exact_row$acc_halfwidth,
  m_precision_rule_grid = m_precision_rule_grid_row$m,
  sens_halfwidth_at_grid_m = m_precision_rule_grid_row$sens_halfwidth,
  spec_halfwidth_at_grid_m = m_precision_rule_grid_row$spec_halfwidth,
  npv_halfwidth_at_grid_m = m_precision_rule_grid_row$npv_halfwidth,
  acc_halfwidth_at_grid_m = m_precision_rule_grid_row$acc_halfwidth,
  m_round_to = m_round_to,
  auto_calculate_m = auto_calculate_m
)

if (auto_calculate_m) {
  m <- as.integer(sample_size_m_summary$m_precision_rule_grid)
} else {
  m <- as.integer(m_manual)
}

m2 <- m - m1

if (m < m1) {
  stop("Calculated m = ", m, " is smaller than m1 = ", m1, ".")
}

sample_size_m_summary <- sample_size_m_summary %>%
  mutate(
    planned_total_negative_sample_m = m,
    planned_stage2_negative_sample_m2 = m2,
    .after = m1
  )

cat("\nNegative sample-size planning summary\n")
print(sample_size_m_summary, width = Inf)

cat("\nUsing total planned negative sample size m = ", m,
    " and Stage 2 negative sample size m2 = ", m2, ".\n", sep = "")

readr::write_csv(
  sample_size_m_summary,
  file.path(output_dir, "sample_size_m_summary.csv")
)

readr::write_csv(
  sample_size_precision_table,
  file.path(output_dir, "sample_size_precision_table.csv")
)

# -------------------------------------------------------------------------
# Step 3. Define full analytic algorithm-positive and negative frames
# -------------------------------------------------------------------------

N <- nrow(analytic)
N_pos <- sum(analytic$algo_positive_original == 1L)
N_neg <- sum(analytic$algo_positive_original == 0L)

if (N_neg < m1) {
  stop("N_neg = ", N_neg, " is smaller than m1 = ", m1,
       "; cannot draw the Stage 1 negative SRS sample.")
}

if (m > N_neg) {
  stop("Total planned negative sample m = ", m,
       " exceeds N_neg = ", N_neg, ".")
}

cat("\nAnalytic cohort size\n")
cat("N = ", N, "\nN_pos = ", N_pos, "\nN_neg = ", N_neg, "\n", sep = "")

positive_frame <- analytic %>% filter(algo_positive_original == 1L)
negative_frame <- analytic %>% filter(algo_positive_original == 0L)

positive_truth_counts <- positive_frame %>%
  summarise(
    TP_pos = sum(true_recur_3yr == 1L & alg_recur_3yr == 1L),
    FP_pos = sum(true_recur_3yr == 0L & alg_recur_3yr == 1L)
  )

cat("\nObserved TP/FP among fully verified algorithm-positive patients\n")
print(positive_truth_counts)

# -------------------------------------------------------------------------
# Step 4. Draw Stage 1 validation sample
# -------------------------------------------------------------------------

set.seed(seed)
stage1_negative_ids <- negative_frame %>%
  slice_sample(n = m1, replace = FALSE) %>%
  pull(patient_id)

analytic <- analytic %>%
  mutate(
    stage1_sample =
      algo_positive_original == 1L | patient_id %in% stage1_negative_ids,
    sample_component = case_when(
      algo_positive_original == 1L ~ "all_algorithm_positive",
      patient_id %in% stage1_negative_ids ~ "stage1_srs_negative",
      TRUE ~ "not_sampled_stage1"
    )
  )

# -------------------------------------------------------------------------
# Step 5. Create Stage 1 design weights
# -------------------------------------------------------------------------

pi_neg_stage1 <- m1 / N_neg
w_neg_stage1 <- N_neg / m1

analytic <- analytic %>%
  mutate(
    pi_stage1 = case_when(
      algo_positive_original == 1L ~ 1,
      sample_component == "stage1_srs_negative" ~ pi_neg_stage1,
      TRUE ~ NA_real_
    ),
    weight_stage1 = case_when(
      algo_positive_original == 1L ~ 1,
      sample_component == "stage1_srs_negative" ~ w_neg_stage1,
      TRUE ~ NA_real_
    )
  )

stage1_data <- analytic %>% filter(stage1_sample)
sampled_neg_stage1 <- stage1_data %>% filter(sample_component == "stage1_srs_negative")

stage1_sample_size_summary <- tibble(
  endpoint_year = endpoint_year,
  endpoint_label = endpoint_label,
  validation_replicate = validation_replicate,
  N_analytic = N,
  N_algorithm_positive = N_pos,
  N_algorithm_negative = N_neg,
  planned_stage1_negative_sample_m1 = m1,
  planned_total_negative_sample_m = m,
  planned_stage2_negative_sample_m2 = m2,
  stage1_sample_size_total = nrow(stage1_data),
  positives_reviewed = sum(stage1_data$sample_component == "all_algorithm_positive"),
  sampled_negatives_reviewed = sum(stage1_data$sample_component == "stage1_srs_negative"),
  pi_neg_stage1 = pi_neg_stage1,
  weight_neg_stage1 = w_neg_stage1
)

cat("\nStage 1 sample size summary\n")
print(stage1_sample_size_summary)

# -------------------------------------------------------------------------
# Step 6. Stage 1 design-based and naive performance metrics
# -------------------------------------------------------------------------

weighted_totals <- stage1_data %>%
  summarise(
    TP_hat = sum(weight_stage1 * (alg_recur_3yr == 1L & true_recur_3yr == 1L)),
    FP_hat = sum(weight_stage1 * (alg_recur_3yr == 1L & true_recur_3yr == 0L)),
    FN_hat = sum(weight_stage1 * (alg_recur_3yr == 0L & true_recur_3yr == 1L)),
    TN_hat = sum(weight_stage1 * (alg_recur_3yr == 0L & true_recur_3yr == 0L))
  )

TP_hat <- weighted_totals$TP_hat
FP_hat <- weighted_totals$FP_hat
FN_hat <- weighted_totals$FN_hat
TN_hat <- weighted_totals$TN_hat

stage1_metrics_weighted <- metric_from_totals(
  TP_hat, FP_hat, FN_hat, TN_hat, denominator_n = N
) %>%
  mutate(endpoint_year = endpoint_year,
         endpoint_label = endpoint_label,
         validation_replicate = validation_replicate,
         .before = 1) %>%
  rename(weighted_estimate_stage1 = estimate)

# These naive metrics are complete-case summaries of the reviewed Stage 1
# sample. They are biased for population performance under partial
# verification because all positives are reviewed but only a small fraction of
# negatives is reviewed.
naive_totals <- stage1_data %>%
  summarise(
    TP = sum(alg_recur_3yr == 1L & true_recur_3yr == 1L),
    FP = sum(alg_recur_3yr == 1L & true_recur_3yr == 0L),
    FN = sum(alg_recur_3yr == 0L & true_recur_3yr == 1L),
    TN = sum(alg_recur_3yr == 0L & true_recur_3yr == 0L)
  )

stage1_metrics_naive <- metric_from_totals(
  naive_totals$TP,
  naive_totals$FP,
  naive_totals$FN,
  naive_totals$TN,
  denominator_n = nrow(stage1_data)
) %>%
  mutate(endpoint_year = endpoint_year,
         endpoint_label = endpoint_label,
         validation_replicate = validation_replicate,
         .before = 1) %>%
  rename(naive_complete_case_estimate_stage1 = estimate)

stage1_metrics_comparison <- stage1_metrics_weighted %>%
  full_join(stage1_metrics_naive,
            by = c("endpoint_year", "endpoint_label",
                   "validation_replicate", "metric"))

cat("\nStage 1 weighted and naive metrics\n")
print(stage1_metrics_comparison, n = Inf)

# -------------------------------------------------------------------------
# Step 7. Design-based variance and confidence intervals
# -------------------------------------------------------------------------

sampled_neg_stage1 <- sampled_neg_stage1 %>%
  mutate(
    fn_indicator = as.integer(alg_recur_3yr == 0L & true_recur_3yr == 1L),
    tn_indicator = as.integer(alg_recur_3yr == 0L & true_recur_3yr == 0L)
  )

m1_actual <- nrow(sampled_neg_stage1)
f_stage1 <- m1_actual / N_neg

# The pilot FN rate among sampled negatives drives sensitivity uncertainty
# because all algorithm positives are fully verified and most residual
# uncertainty is about how many recurrences were missed among algorithm
# negatives.
Var_FN_hat <- N_neg^2 * (1 - f_stage1) *
  sample_variance(sampled_neg_stage1$fn_indicator) / m1_actual
Var_TN_hat <- N_neg^2 * (1 - f_stage1) *
  sample_variance(sampled_neg_stage1$tn_indicator) / m1_actual
Cov_TN_FN_hat <- N_neg^2 * (1 - f_stage1) *
  sample_covariance(sampled_neg_stage1$tn_indicator,
                    sampled_neg_stage1$fn_indicator) / m1_actual

Var_sensitivity <- (safe_divide(TP_hat, (TP_hat + FN_hat)^2))^2 * Var_FN_hat
Var_specificity <- (safe_divide(FP_hat, (TN_hat + FP_hat)^2))^2 * Var_TN_hat
Var_accuracy <- Var_TN_hat / N^2

den_npv <- TN_hat + FN_hat
grad_tn_npv <- safe_divide(FN_hat, den_npv^2)
grad_fn_npv <- -safe_divide(TN_hat, den_npv^2)
Var_NPV <- grad_tn_npv^2 * Var_TN_hat +
  grad_fn_npv^2 * Var_FN_hat +
  2 * grad_tn_npv * grad_fn_npv * Cov_TN_FN_hat

weighted_metric_lookup <- stage1_metrics_weighted %>%
  select(metric, weighted_estimate_stage1) %>%
  deframe()

ci_sensitivity <- logit_ci(weighted_metric_lookup[["sensitivity"]],
                           Var_sensitivity)
ci_specificity <- logit_ci(weighted_metric_lookup[["specificity"]],
                           Var_specificity)
ci_npv <- logit_ci(weighted_metric_lookup[["npv"]], Var_NPV)
ci_accuracy <- logit_ci(weighted_metric_lookup[["accuracy"]], Var_accuracy)

ppv_n <- sum(stage1_data$alg_recur_3yr == 1L)
ppv_x <- sum(stage1_data$alg_recur_3yr == 1L & stage1_data$true_recur_3yr == 1L)
ci_ppv <- wilson_ci(ppv_x, ppv_n)
se_ppv <- sqrt(
  weighted_metric_lookup[["ppv"]] * (1 - weighted_metric_lookup[["ppv"]]) /
    ppv_n
)

if (any(stage1_data$sample_component == "stage1_srs_negative" &
        stage1_data$alg_recur_3yr == 1L)) {
  warning("At least one sampled algorithm-negative row has alg_recur_3yr = 1. ",
          "The PPV Wilson interval assumes predicted-positive rows are fully verified.")
}

stage1_metrics_ci <- tibble(
  endpoint_year = endpoint_year,
  endpoint_label = endpoint_label,
  validation_replicate = validation_replicate,
  metric = c("sensitivity", "specificity", "ppv", "npv", "accuracy"),
  estimate = c(
    weighted_metric_lookup[["sensitivity"]],
    weighted_metric_lookup[["specificity"]],
    weighted_metric_lookup[["ppv"]],
    weighted_metric_lookup[["npv"]],
    weighted_metric_lookup[["accuracy"]]
  ),
  standard_error = c(
    sqrt(Var_sensitivity),
    sqrt(Var_specificity),
    se_ppv,
    sqrt(Var_NPV),
    sqrt(Var_accuracy)
  ),
  ci_lower = c(
    ci_sensitivity[["lower"]],
    ci_specificity[["lower"]],
    ci_ppv[["lower"]],
    ci_npv[["lower"]],
    ci_accuracy[["lower"]]
  ),
  ci_upper = c(
    ci_sensitivity[["upper"]],
    ci_specificity[["upper"]],
    ci_ppv[["upper"]],
    ci_npv[["upper"]],
    ci_accuracy[["upper"]]
  ),
  variance_method = c(
    "Delta method using SRSWOR variance of FN total among sampled negatives",
    "Delta method using SRSWOR variance of TN total among sampled negatives",
    "Wilson/binomial among verified predicted-positive patients",
    "Delta method ratio estimator using TN/FN covariance among sampled negatives",
    "Delta method using SRSWOR variance of TN total among sampled negatives"
  )
)

cat("\nStage 1 design-based CI table\n")
print(stage1_metrics_ci, n = Inf)

# -------------------------------------------------------------------------
# Step 8. Create risk strata among algorithm-negative patients
# -------------------------------------------------------------------------

if (all(is.na(negative_frame$predicted_risk))) {
  stop("All algorithm-negative predicted_risk values are missing; cannot create risk strata.")
}

median_negative_risk <- median(negative_frame$predicted_risk, na.rm = TRUE)
n_missing_negative_risk <- sum(is.na(negative_frame$predicted_risk))

if (n_missing_negative_risk > 0) {
  warning(n_missing_negative_risk,
          " algorithm-negative patients have missing predicted risk. ",
          "They will be assigned strata using the median negative risk.")
}

negative_strata <- analytic %>%
  filter(algo_positive_original == 0L) %>%
  mutate(risk_for_strata = coalesce(predicted_risk, median_negative_risk)) %>%
  arrange(risk_for_strata, patient_id) %>%
  mutate(risk_stratum = ntile(row_number(), H)) %>%
  select(patient_id, risk_stratum, risk_for_strata)

analytic <- analytic %>%
  left_join(negative_strata, by = "patient_id")

stage1_data <- analytic %>% filter(stage1_sample)
sampled_neg_stage1 <- stage1_data %>%
  filter(sample_component == "stage1_srs_negative") %>%
  mutate(
    fn_indicator = as.integer(alg_recur_3yr == 0L & true_recur_3yr == 1L),
    tn_indicator = as.integer(alg_recur_3yr == 0L & true_recur_3yr == 0L)
  )

cat("\nAlgorithm-negative risk stratum counts\n")
print(analytic %>%
        filter(algo_positive_original == 0L) %>%
        count(risk_stratum, name = "N_h") %>%
        arrange(risk_stratum),
      n = Inf)

# -------------------------------------------------------------------------
# Step 9. Pilot FN analysis by risk stratum
# -------------------------------------------------------------------------

stratum_index <- tibble(risk_stratum = seq_len(H))

stage1_pilot_fn_by_stratum <- stratum_index %>%
  left_join(
    analytic %>%
      filter(algo_positive_original == 0L) %>%
      count(risk_stratum, name = "N_h"),
    by = "risk_stratum"
  ) %>%
  left_join(
    sampled_neg_stage1 %>%
      group_by(risk_stratum) %>%
      summarise(
        m1_h = n(),
        c1_h = sum(fn_indicator),
        d1_h = sum(tn_indicator),
        .groups = "drop"
      ),
    by = "risk_stratum"
  ) %>%
  mutate(
    endpoint_year = endpoint_year,
    endpoint_label = endpoint_label,
    validation_replicate = validation_replicate,
    across(c(N_h, m1_h, c1_h, d1_h), ~ replace_na(.x, 0L)),
    p_FN_h_raw = if_else(m1_h > 0, c1_h / m1_h, NA_real_),
    p_FN_h_jeffreys = (c1_h + 0.5) / (m1_h + 1),
    p_TN_h_jeffreys = 1 - p_FN_h_jeffreys,
    S_h = sqrt(p_FN_h_jeffreys * (1 - p_FN_h_jeffreys))
  )

c1_total <- sum(stage1_pilot_fn_by_stratum$c1_h)
p_FN_overall_raw <- c1_total / m1_actual
p_FN_overall_jeffreys <- (c1_total + 0.5) / (m1_actual + 1)

cat("\nStage 1 pilot false-negative summary by risk stratum\n")
print(stage1_pilot_fn_by_stratum, n = Inf)

# -------------------------------------------------------------------------
# Step 10. Projected variance ratio Rhat(m) and Neyman allocation
# -------------------------------------------------------------------------

f_srs <- m / N_neg

# For the projected Rhat comparison, use a single internally consistent pilot
# model for both SRS and Neyman. The descriptive overall Jeffreys rate above
# uses one 0.5/0.5 prior for the entire pilot sample. Neyman uses stratum-level
# Jeffreys rates, which effectively apply that weak prior within each stratum.
# With rare false negatives, comparing overall-smoothed SRS to stratum-smoothed
# Neyman can make the optimized allocation look artificially worse. The
# stratum-calibrated rate below keeps the SRS variance on the same smoothed
# stratum scale as Neyman.
p_FN_srs_projected <- weighted.mean(
  stage1_pilot_fn_by_stratum$p_FN_h_jeffreys,
  stage1_pilot_fn_by_stratum$N_h
)

Var_SRS_FN <- N_neg^2 * (1 - f_srs) *
  p_FN_srs_projected * (1 - p_FN_srs_projected) / (m - 1)

neyman_allocation <- allocate_neyman_constrained(
  stage1_pilot_fn_by_stratum,
  target_m = m
)

if (sum(neyman_allocation$m2_h) != m2) {
  stop("Neyman Stage 2 allocation sums to ",
       sum(neyman_allocation$m2_h), " but expected m2 = ", m2, ".")
}

Var_Ney_FN <- neyman_allocation %>%
  mutate(
    f_h = m_h_neyman_total / N_h,
    variance_component = case_when(
      N_h <= 0 ~ 0,
      m_h_neyman_total >= N_h ~ 0,
      m_h_neyman_total > 1 ~
        N_h^2 * (1 - f_h) *
        p_FN_h_jeffreys * (1 - p_FN_h_jeffreys) /
        (m_h_neyman_total - 1),
      TRUE ~ NA_real_
    )
  ) %>%
  summarise(value = sum(variance_component, na.rm = TRUE)) %>%
  pull(value)

sens_derivative <- safe_divide(TP_hat, (TP_hat + FN_hat)^2)
Var_SRS_Sens <- sens_derivative^2 * Var_SRS_FN
Var_Ney_Sens <- sens_derivative^2 * Var_Ney_FN
Rhat_FN <- safe_divide(Var_SRS_FN, Var_Ney_FN)
Rhat_Sens <- safe_divide(Var_SRS_Sens, Var_Ney_Sens)

# Rhat > 1.1 means the projected sensitivity variance under SRS is at least
# 10% larger than under Neyman allocation. That margin is used as a practical
# trigger so Stage 2 only switches designs when the expected gain is meaningful.
Rhat <- Rhat_FN
trigger_neyman <- !is.na(Rhat) && Rhat > trigger_threshold

projected_variance_ratio_rhat <- tibble(
  endpoint_year = endpoint_year,
  endpoint_label = endpoint_label,
  validation_replicate = validation_replicate,
  m = m,
  m1 = m1,
  m2 = m2,
  N_neg = N_neg,
  p_FN_overall_raw = p_FN_overall_raw,
  p_FN_overall_jeffreys = p_FN_overall_jeffreys,
  p_FN_srs_projected = p_FN_srs_projected,
  Var_SRS_FN = Var_SRS_FN,
  Var_Ney_FN = Var_Ney_FN,
  Var_SRS_Sens = Var_SRS_Sens,
  Var_Ney_Sens = Var_Ney_Sens,
  Rhat_FN = Rhat_FN,
  Rhat_Sens = Rhat_Sens,
  Rhat = Rhat,
  trigger_threshold = trigger_threshold,
  trigger_neyman = trigger_neyman
)

cat("\nProjected variance ratio Rhat(m)\n")
print(as.data.frame(projected_variance_ratio_rhat))

stage2_allocation_recommendation <- neyman_allocation %>%
  select(
    risk_stratum, N_h, m1_h, c1_h, p_FN_h_jeffreys, S_h,
    m_h_min_total, m_h_neyman_total, m2_h
  ) %>%
  mutate(
    endpoint_year = endpoint_year,
    endpoint_label = endpoint_label,
    validation_replicate = validation_replicate,
    trigger_neyman = trigger_neyman,
    stage2_method_recommended = if_else(trigger_neyman,
                                        "neyman_by_risk_stratum",
                                        "srs_among_remaining_negatives")
  )

cat("\nRecommended Stage 2 allocation table\n")
print(stage2_allocation_recommendation, n = Inf)

# -------------------------------------------------------------------------
# Step 11. Generate Stage 2 and final validation sample lists
#         under BOTH scenarios:
#         1) Neyman allocation
#         2) Continued SRS among remaining algorithm-negatives
# -------------------------------------------------------------------------

remaining_negative_frame <- analytic %>%
  filter(algo_positive_original == 0L, stage1_sample == FALSE)

if (nrow(remaining_negative_frame) < m2) {
  stop("Fewer than m2 remaining algorithm-negative patients are available.")
}

# ---- Scenario A: Stage 2 Neyman allocation ----
# Uses the recommended m2_h by risk stratum, regardless of trigger_neyman.
stage2_neyman_sample <- sample_within_strata(
  remaining_negative_frame,
  stage2_allocation_recommendation,
  n_col = "m2_h",
  seed_value = seed + 1L
) %>%
  mutate(
    stage2_component = "stage2_neyman_negative",
    stage2_design = "neyman"
  )

# ---- Scenario B: Stage 2 SRS allocation ----
# Samples m2 additional algorithm-negative patients from all remaining negatives.
set.seed(seed + 2L)

stage2_srs_sample <- remaining_negative_frame %>%
  slice_sample(n = m2, replace = FALSE) %>%
  mutate(
    stage2_component = "stage2_srs_negative",
    stage2_design = "srs"
  )

# ---- Compare realized Stage 2 stratum counts ----
stage2_neyman_counts_by_stratum <- stage2_neyman_sample %>%
  count(risk_stratum, name = "m2_h_used_neyman")

stage2_srs_counts_by_stratum <- stage2_srs_sample %>%
  count(risk_stratum, name = "m2_h_used_srs")

stage2_allocation_recommendation <- stage2_allocation_recommendation %>%
  left_join(stage2_neyman_counts_by_stratum, by = "risk_stratum") %>%
  left_join(stage2_srs_counts_by_stratum, by = "risk_stratum") %>%
  mutate(
    m2_h_used_neyman = replace_na(m2_h_used_neyman, 0L),
    m2_h_used_srs = replace_na(m2_h_used_srs, 0L),
    trigger_neyman = trigger_neyman,
    trigger_threshold = trigger_threshold,
    stage2_method_if_following_trigger = if_else(
      trigger_neyman,
      "neyman_by_risk_stratum",
      "srs_among_remaining_negatives"
    )
  )

# ---- Stage 1 negative sample, common to both designs ----
stage1_negative_sample <- analytic %>%
  filter(sample_component == "stage1_srs_negative") %>%
  mutate(
    validation_phase = "stage1",
    final_sample_component = "stage1_srs_negative"
  )

# ---- Final negative sample under Neyman ----
final_negative_sample_neyman <- bind_rows(
  stage1_negative_sample %>%
    mutate(stage2_design = "neyman"),
  stage2_neyman_sample %>%
    mutate(
      validation_phase = "stage2",
      final_sample_component = stage2_component
    )
) %>%
  arrange(validation_phase, risk_stratum, patient_id)

# ---- Final negative sample under SRS ----
final_negative_sample_srs <- bind_rows(
  stage1_negative_sample %>%
    mutate(stage2_design = "srs"),
  stage2_srs_sample %>%
    mutate(
      validation_phase = "stage2",
      final_sample_component = stage2_component
    )
) %>%
  arrange(validation_phase, risk_stratum, patient_id)

# ---- Final all-patient validation sample under Neyman ----
final_validation_sample_neyman <- bind_rows(
  analytic %>%
    filter(algo_positive_original == 1L) %>%
    mutate(
      validation_phase = "stage1",
      final_sample_component = "all_algorithm_positive",
      stage2_design = "neyman"
    ),
  final_negative_sample_neyman
) %>%
  arrange(final_sample_component, risk_stratum, patient_id)

# ---- Final all-patient validation sample under SRS ----
final_validation_sample_srs <- bind_rows(
  analytic %>%
    filter(algo_positive_original == 1L) %>%
    mutate(
      validation_phase = "stage1",
      final_sample_component = "all_algorithm_positive",
      stage2_design = "srs"
    ),
  final_negative_sample_srs
) %>%
  arrange(final_sample_component, risk_stratum, patient_id)

# ---- Add indicators back to analytic dataset ----
stage2_neyman_ids <- stage2_neyman_sample$patient_id
stage2_srs_ids <- stage2_srs_sample$patient_id

analytic <- analytic %>%
  mutate(
    stage2_sample_neyman = patient_id %in% stage2_neyman_ids,
    stage2_sample_srs = patient_id %in% stage2_srs_ids,
    
    final_validation_sample_neyman =
      patient_id %in% final_validation_sample_neyman$patient_id,
    final_validation_sample_srs =
      patient_id %in% final_validation_sample_srs$patient_id,
    
    final_sample_component_neyman = case_when(
      sample_component == "all_algorithm_positive" ~ "all_algorithm_positive",
      sample_component == "stage1_srs_negative" ~ "stage1_srs_negative",
      patient_id %in% stage2_neyman_ids ~ "stage2_neyman_negative",
      TRUE ~ "not_in_final_validation_sample"
    ),
    
    final_sample_component_srs = case_when(
      sample_component == "all_algorithm_positive" ~ "all_algorithm_positive",
      sample_component == "stage1_srs_negative" ~ "stage1_srs_negative",
      patient_id %in% stage2_srs_ids ~ "stage2_srs_negative",
      TRUE ~ "not_in_final_validation_sample"
    )
  )

# ---- Sanity checks ----
if (nrow(stage2_neyman_sample) != m2) {
  stop("Neyman Stage 2 sample size is ", nrow(stage2_neyman_sample),
       " but expected m2 = ", m2, ".")
}

if (nrow(stage2_srs_sample) != m2) {
  stop("SRS Stage 2 sample size is ", nrow(stage2_srs_sample),
       " but expected m2 = ", m2, ".")
}

if (nrow(final_negative_sample_neyman) != m) {
  stop("Final Neyman negative sample size is ",
       nrow(final_negative_sample_neyman),
       " but expected m = ", m, ".")
}

if (nrow(final_negative_sample_srs) != m) {
  stop("Final SRS negative sample size is ",
       nrow(final_negative_sample_srs),
       " but expected m = ", m, ".")
}

# -------------------------------------------------------------------------
# Step 12. Save and print QC outputs
# -------------------------------------------------------------------------

cohort_size_summary <- tibble(
  endpoint_year = endpoint_year,
  endpoint_label = endpoint_label,
  validation_replicate = validation_replicate,
  N_analytic = N,
  N_algorithm_positive = N_pos,
  N_algorithm_negative = N_neg,
  planned_stage1_negative_sample_m1 = m1,
  Stage_1_sample_size_total = nrow(stage1_data),
  positives_reviewed = sum(stage1_data$sample_component == "all_algorithm_positive"),
  sampled_negatives_reviewed = sum(stage1_data$sample_component == "stage1_srs_negative"),
  planned_total_negative_sample_m = m,
  planned_stage2_negative_sample_m2 = m2,
  
  stage2_neyman_sample_size = nrow(stage2_neyman_sample),
  final_negative_sample_size_neyman = nrow(final_negative_sample_neyman),
  final_validation_sample_size_neyman = nrow(final_validation_sample_neyman),
  
  stage2_srs_sample_size = nrow(stage2_srs_sample),
  final_negative_sample_size_srs = nrow(final_negative_sample_srs),
  final_validation_sample_size_srs = nrow(final_validation_sample_srs),
  
  trigger_neyman = trigger_neyman,
  trigger_threshold = trigger_threshold,
  stage2_method_if_following_trigger = if_else(
    trigger_neyman,
    "neyman_by_risk_stratum",
    "srs_among_remaining_negatives"
  )
)

readr::write_csv(
  analytic,
  file.path(output_dir, "merged_clean_validation_data.csv")
)

readr::write_csv(
  stage1_data,
  file.path(output_dir, "stage1_validation_sample.csv")
)

readr::write_csv(
  stage1_metrics_weighted,
  file.path(output_dir, "stage1_metrics_weighted.csv")
)

readr::write_csv(
  stage1_metrics_naive,
  file.path(output_dir, "stage1_metrics_naive.csv")
)

readr::write_csv(
  stage1_metrics_comparison,
  file.path(output_dir, "stage1_metrics_comparison.csv")
)

readr::write_csv(
  stage1_metrics_ci,
  file.path(output_dir, "stage1_metrics_ci.csv")
)

readr::write_csv(
  stage1_pilot_fn_by_stratum,
  file.path(output_dir, "stage1_pilot_fn_by_stratum.csv")
)

readr::write_csv(
  projected_variance_ratio_rhat,
  file.path(output_dir, "projected_variance_ratio_rhat.csv")
)

readr::write_csv(
  stage2_allocation_recommendation,
  file.path(output_dir, "stage2_allocation_recommendation.csv")
)

# New: save BOTH Stage 2 samples
readr::write_csv(
  stage2_neyman_sample,
  file.path(output_dir, "stage2_neyman_sample_list.csv")
)

readr::write_csv(
  stage2_srs_sample,
  file.path(output_dir, "stage2_srs_sample_list.csv")
)

# New: save BOTH final negative samples
readr::write_csv(
  final_negative_sample_neyman,
  file.path(output_dir, "final_negative_sample_list_neyman.csv")
)

readr::write_csv(
  final_negative_sample_srs,
  file.path(output_dir, "final_negative_sample_list_srs.csv")
)

# New: save BOTH final validation samples
readr::write_csv(
  final_validation_sample_neyman,
  file.path(output_dir, "final_validation_sample_list_neyman.csv")
)

readr::write_csv(
  final_validation_sample_srs,
  file.path(output_dir, "final_validation_sample_list_srs.csv")
)

# Optional backward-compatible files based on trigger decision
if (trigger_neyman) {
  readr::write_csv(
    stage2_neyman_sample,
    file.path(output_dir, "stage2_sample_list.csv")
  )
  readr::write_csv(
    final_negative_sample_neyman,
    file.path(output_dir, "final_negative_sample_list.csv")
  )
  readr::write_csv(
    final_validation_sample_neyman,
    file.path(output_dir, "final_validation_sample_list.csv")
  )
} else {
  readr::write_csv(
    stage2_srs_sample,
    file.path(output_dir, "stage2_sample_list.csv")
  )
  readr::write_csv(
    final_negative_sample_srs,
    file.path(output_dir, "final_negative_sample_list.csv")
  )
  readr::write_csv(
    final_validation_sample_srs,
    file.path(output_dir, "final_validation_sample_list.csv")
  )
}

readr::write_csv(
  cohort_size_summary,
  file.path(output_dir, "cohort_size_summary.csv")
)

readr::write_csv(
  weighted_totals,
  file.path(output_dir, "stage1_weighted_confusion_totals.csv")
)

readr::write_csv(
  naive_totals,
  file.path(output_dir, "stage1_naive_confusion_totals.csv")
)

readr::write_csv(
  algorithm_prediction_crosstab,
  file.path(output_dir, "algorithm_prediction_crosstab.csv")
)

cat("\nCohort size summary\n")
print(cohort_size_summary, width = Inf)

cat("\nStage 1 weighted metrics table\n")
print(stage1_metrics_weighted, n = Inf)

cat("\nStage 1 naive complete-case metrics table\n")
print(stage1_metrics_naive, n = Inf)

cat("\nStage 1 variance/CI table\n")
print(stage1_metrics_ci, n = Inf)

cat("\nPilot negative FN summary by risk stratum\n")
print(stage1_pilot_fn_by_stratum, n = Inf)

cat("\nProjected Rhat table\n")
print(projected_variance_ratio_rhat, width = Inf)

cat("\nStage 2 allocation table with both realized samples\n")
print(stage2_allocation_recommendation, n = Inf, width = Inf)

cat("\nSaved Stage 1/Stage 2 validation outputs to: ",
    normalizePath(output_dir), "\n", sep = "")
