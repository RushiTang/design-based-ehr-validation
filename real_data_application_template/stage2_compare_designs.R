############################################################
# Stage 2 final validation analysis:
# Compare the adaptive composite design with final SRS
#
# Inputs from Stage 1 script:
#   stage1_validation_outputs/merged_clean_validation_data.csv
#   stage1_validation_outputs/final_validation_sample_list_neyman.csv
#   stage1_validation_outputs/final_validation_sample_list_srs.csv
#   stage1_validation_outputs/stage2_allocation_recommendation.csv
#
# Outputs:
#   stage2_validation_outputs/*.csv
#
# Design logic:
# - All original algorithm-positive patients are a census and have weight 1.
# - Under final SRS, all sampled original algorithm-negative patients share
#   weight N_neg / m.
# - If the adaptive rule triggers, estimation uses the fixed-weight lambda-pool
#   composite of the Stage 1 HT estimator and the sequential conditional HT
#   estimator based on the Stage 2 Neyman sample.
# - If the adaptive rule does not trigger, its final sample and estimator are
#   the same pooled SRS procedure used by the SRS comparator.
# - Metric-specific logit-Wald intervals use the corresponding design variance.
#   A finite-population hypergeometric safeguard is used when a final SRS sample
#   contains zero false negatives. PPV retains its Wilson interval.
# - Performance metrics use algorithm recurrence versus true recurrence at the
#   selected x-year endpoint.
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
  library(lubridate)
})

# -------------------------------------------------------------------------
# User settings
# -------------------------------------------------------------------------

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- "--file="
script_file_match <- script_args[startsWith(script_args, script_file_arg)]
script_path <- if (length(script_file_match) > 0) {
  sub(script_file_arg, "", script_file_match[[1]])
} else {
  NA_character_
}
if (!is.na(script_path)) {
  setwd(dirname(normalizePath(script_path)))
}

stage1_output_dir <- Sys.getenv("STAGE1_OUTPUT_DIR",
                                unset = "stage1_validation_outputs")
output_dir <- Sys.getenv("STAGE2_OUTPUT_DIR",
                         unset = "stage2_validation_outputs")

endpoint_year <- as.integer(Sys.getenv("VALIDATION_ENDPOINT_YEAR", unset = "3"))
if (is.na(endpoint_year) || endpoint_year <= 0) {
  stop("VALIDATION_ENDPOINT_YEAR must be a positive integer.")
}
endpoint_label <- paste0(endpoint_year, "yr")

validation_replicate <- as.integer(Sys.getenv("VALIDATION_REPLICATE", unset = "1"))
if (is.na(validation_replicate) || validation_replicate <= 0) {
  stop("VALIDATION_REPLICATE must be a positive integer.")
}

full_analytic_csv <- file.path(stage1_output_dir, "merged_clean_validation_data.csv")
final_neyman_csv <- file.path(stage1_output_dir, "final_validation_sample_list_neyman.csv")
final_srs_csv <- file.path(stage1_output_dir, "final_validation_sample_list_srs.csv")
allocation_csv <- file.path(stage1_output_dir, "stage2_allocation_recommendation.csv")
rhat_csv <- file.path(stage1_output_dir, "projected_variance_ratio_rhat.csv")
sample_size_m_summary_csv <- file.path(stage1_output_dir, "sample_size_m_summary.csv")

id_var <- "patient_id"
truth_var <- paste0("true_recur_", endpoint_year, "yr")
algorithm_endpoint_var <- paste0("alg_recur_", endpoint_year, "yr")
original_algorithm_positive_var <- "algo_positive_original"
stratum_var <- "risk_stratum"
algo_pred_time_unit <- "months"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

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

parse_numeric_flexible <- function(x) {
  suppressWarnings(readr::parse_number(as.character(x)))
}

add_endpoint_columns <- function(data, endpoint_year) {
  truth_col <- paste0("true_recur_", endpoint_year, "yr")
  evaluable_col <- paste0("evaluable_", endpoint_year, "yr_primary")
  alg_col <- paste0("alg_recur_", endpoint_year, "yr")

  if (!truth_col %in% names(data)) {
    required_truth_cols <- c("diagnosis_date", "recurrence_by_pull",
                             "recurrence_date_by_pull",
                             "followup_years_primary")
    missing_truth_cols <- setdiff(required_truth_cols, names(data))

    if (length(missing_truth_cols) > 0) {
      stop("Cannot derive ", truth_col, " because merged analytic data are missing: ",
           paste(missing_truth_cols, collapse = ", "))
    }

    diagnosis_date <- parse_date_flexible(data$diagnosis_date)
    recurrence_by_pull <- standardize_binary(data$recurrence_by_pull)
    recurrence_date_by_pull <- parse_date_flexible(data$recurrence_date_by_pull)
    followup_years_primary <- parse_numeric_flexible(data$followup_years_primary)
    endpoint_date <- diagnosis_date %m+% years(endpoint_year)

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
        followup_years_primary >= endpoint_year ~ 0L,
      TRUE ~ NA_integer_
    )

    data[[truth_col]] <- true_recur
    data[[evaluable_col]] <- !is.na(true_recur)
  }

  if (!alg_col %in% names(data)) {
    if (!original_algorithm_positive_var %in% names(data)) {
      stop("Cannot derive ", alg_col, " because ",
           original_algorithm_positive_var, " is missing.")
    }

    algo_positive_original <- standardize_binary(data[[original_algorithm_positive_var]])

    if ("predicted_recurrence_date" %in% names(data) &&
        any(!is.na(data$predicted_recurrence_date))) {
      if (!"diagnosis_date" %in% names(data)) {
        stop("Cannot derive ", alg_col,
             " from predicted_recurrence_date without diagnosis_date.")
      }

      diagnosis_date <- parse_date_flexible(data$diagnosis_date)
      predicted_recurrence_date <- parse_date_flexible(data$predicted_recurrence_date)
      endpoint_date <- diagnosis_date %m+% years(endpoint_year)

      data[[alg_col]] <- case_when(
        algo_positive_original == 1L &
          !is.na(predicted_recurrence_date) &
          predicted_recurrence_date >= diagnosis_date &
          predicted_recurrence_date <= endpoint_date ~ 1L,
        !is.na(algo_positive_original) ~ 0L,
        TRUE ~ NA_integer_
      )
    } else if ("predicted_recurrence_time" %in% names(data) &&
               any(!is.na(data$predicted_recurrence_time))) {
      predicted_recurrence_time <- parse_numeric_flexible(data$predicted_recurrence_time)
      time_threshold <- case_when(
        algo_pred_time_unit == "months" ~ 12 * endpoint_year,
        algo_pred_time_unit == "days" ~ 365.25 * endpoint_year,
        algo_pred_time_unit == "years" ~ endpoint_year,
        TRUE ~ NA_real_
      )

      if (is.na(time_threshold)) {
        stop("algo_pred_time_unit must be one of: months, days, years.")
      }

      data[[alg_col]] <- case_when(
        algo_positive_original == 1L &
          !is.na(predicted_recurrence_time) &
          predicted_recurrence_time >= 0 &
          predicted_recurrence_time <= time_threshold ~ 1L,
        !is.na(algo_positive_original) ~ 0L,
        TRUE ~ NA_integer_
      )
    } else {
      warning("No predicted recurrence time/date column was found. ",
              "Using the binary algorithm flag directly as ", alg_col, ".")
      data[[alg_col]] <- algo_positive_original
    }
  }

  data
}

safe_divide <- function(numerator, denominator) {
  ifelse(is.na(denominator) | denominator == 0, NA_real_,
         numerator / denominator)
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

hypergeom_zero_upper <- function(population_size, sample_size,
                                 level = 0.95) {
  if (is.na(population_size) || is.na(sample_size) || sample_size <= 0L ||
      population_size <= 0L || sample_size > population_size) {
    return(NA_integer_)
  }

  # Central two-sided exact interval: P(X = 0 | K) >= alpha / 2.
  log_threshold <- log((1 - level) / 2)
  log_prob_zero <- function(total_cases) {
    dhyper(0, total_cases, population_size - total_cases,
           sample_size, log = TRUE)
  }

  low <- 0L
  high <- as.integer(population_size)
  while (low < high) {
    middle <- as.integer(ceiling((low + high) / 2))
    if (log_prob_zero(middle) >= log_threshold) {
      low <- middle
    } else {
      high <- middle - 1L
    }
  }
  low
}

metric_interval_from_frame_fn <- function(metric, fixed_totals,
                                          frame_size, fn_upper,
                                          denominator_n) {
  totals_at <- function(frame_fn) {
    c(
      TP_hat = fixed_totals[["TP_hat"]],
      FP_hat = fixed_totals[["FP_hat"]],
      FN_hat = fixed_totals[["FN_hat"]] + frame_fn,
      TN_hat = fixed_totals[["TN_hat"]] + frame_size - frame_fn
    )
  }

  metric_value <- function(frame_fn) {
    totals <- totals_at(frame_fn)
    metric_from_totals(
      totals[["TP_hat"]], totals[["FP_hat"]],
      totals[["FN_hat"]], totals[["TN_hat"]], denominator_n
    ) %>%
      filter(.data$metric == .env$metric) %>%
      pull(estimate)
  }

  c(lower = metric_value(fn_upper), upper = metric_value(0))
}

weighted_auc <- function(truth, score, weights = NULL) {
  if (is.null(weights)) {
    weights <- rep(1, length(truth))
  }

  data <- tibble(
    truth = standardize_binary(truth),
    score = as.numeric(score),
    weight = as.numeric(weights)
  ) %>%
    filter(!is.na(truth), !is.na(score), !is.na(weight), weight > 0)

  if (nrow(data) == 0 ||
      sum(data$truth == 1L) == 0 ||
      sum(data$truth == 0L) == 0) {
    return(NA_real_)
  }

  score_groups <- data %>%
    group_by(score) %>%
    summarise(
      w_pos = sum(weight[truth == 1L]),
      w_neg = sum(weight[truth == 0L]),
      .groups = "drop"
    ) %>%
    arrange(score) %>%
    mutate(cum_neg_before = lag(cumsum(w_neg), default = 0))

  total_pos_weight <- sum(score_groups$w_pos)
  total_neg_weight <- sum(score_groups$w_neg)

  sum(score_groups$w_pos *
        (score_groups$cum_neg_before + 0.5 * score_groups$w_neg)) /
    (total_pos_weight * total_neg_weight)
}

add_cell_indicators <- function(data) {
  data %>%
    mutate(
      cell_tp = as.integer(alg_recur_3yr == 1L & true_recur_3yr == 1L),
      cell_fp = as.integer(alg_recur_3yr == 1L & true_recur_3yr == 0L),
      cell_fn = as.integer(alg_recur_3yr == 0L & true_recur_3yr == 1L),
      cell_tn = as.integer(alg_recur_3yr == 0L & true_recur_3yr == 0L)
    )
}

safe_cov_matrix <- function(data, indicator_cols) {
  x <- data %>%
    select(all_of(indicator_cols)) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()

  if (nrow(x) <= 1) {
    warning("A variance stratum has <=1 sampled row; returning zero covariance.")
    out <- matrix(0, nrow = length(indicator_cols), ncol = length(indicator_cols))
    dimnames(out) <- list(indicator_cols, indicator_cols)
    return(out)
  }

  out <- stats::cov(x, use = "pairwise.complete.obs")
  out[is.na(out)] <- 0
  out
}

metric_gradients <- function(totals, denominator_n) {
  tp <- totals[["TP_hat"]]
  fp <- totals[["FP_hat"]]
  fn <- totals[["FN_hat"]]
  tn <- totals[["TN_hat"]]

  # Gradient order must match: TP, FP, FN, TN.
  list(
    sensitivity = c(
      safe_divide(fn, (tp + fn)^2),
      0,
      -safe_divide(tp, (tp + fn)^2),
      0
    ),
    specificity = c(
      0,
      -safe_divide(tn, (tn + fp)^2),
      0,
      safe_divide(fp, (tn + fp)^2)
    ),
    ppv = c(
      safe_divide(fp, (tp + fp)^2),
      -safe_divide(tp, (tp + fp)^2),
      0,
      0
    ),
    npv = c(
      0,
      0,
      -safe_divide(tn, (tn + fn)^2),
      safe_divide(fn, (tn + fn)^2)
    ),
    accuracy = c(1 / denominator_n, 0, 0, 1 / denominator_n)
  )
}

quadratic_variance <- function(gradient, variance_matrix) {
  gradient <- matrix(gradient, ncol = 1)
  as.numeric(t(gradient) %*% variance_matrix %*% gradient)
}

read_clean_csv <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  readr::read_csv(path, show_col_types = FALSE) %>%
    janitor::clean_names()
}

build_final_sample <- function(sample_path, full_analytic, design_label) {
  sample_ids <- read_clean_csv(sample_path)

  if (!id_var %in% names(sample_ids)) {
    stop("Final sample file does not contain ", id_var, ": ", sample_path)
  }

  metadata_cols <- intersect(
    c(id_var, "validation_phase", "final_sample_component",
      "stage2_design", "stage2_component"),
    names(sample_ids)
  )

  sample_ids <- sample_ids %>%
    distinct(.data[[id_var]], .keep_all = TRUE) %>%
    select(all_of(metadata_cols))

  out <- full_analytic %>%
    semi_join(sample_ids %>% select(all_of(id_var)), by = id_var) %>%
    left_join(sample_ids, by = id_var, suffix = c("", "_sample_file")) %>%
    mutate(stage2_design = design_label)

  if (nrow(out) != nrow(sample_ids)) {
    stop(
      "Could not recover every sampled patient for ", design_label,
      ". sample file rows = ", nrow(sample_ids),
      ", joined rows = ", nrow(out), "."
    )
  }

  out
}

add_design_weights <- function(sample_data, full_analytic, design_label) {
  N_neg <- full_analytic %>%
    filter(algo_positive_original == 0L) %>%
    nrow()

  if (design_label == "srs") {
    m_neg <- sample_data %>%
      filter(algo_positive_original == 0L) %>%
      nrow()

    sample_data %>%
      mutate(
        design_weight = case_when(
          algo_positive_original == 1L ~ 1,
          algo_positive_original == 0L ~ N_neg / m_neg,
          TRUE ~ NA_real_
        ),
        design_inclusion_probability = case_when(
          algo_positive_original == 1L ~ 1,
          algo_positive_original == 0L ~ m_neg / N_neg,
          TRUE ~ NA_real_
        )
      )
  } else if (design_label == "neyman") {
    stratum_weights <- full_analytic %>%
      filter(algo_positive_original == 0L) %>%
      count(risk_stratum, name = "N_h") %>%
      left_join(
        sample_data %>%
          filter(algo_positive_original == 0L) %>%
          count(risk_stratum, name = "m_h"),
        by = "risk_stratum"
      ) %>%
      mutate(
        m_h = replace_na(m_h, 0L),
        stratum_weight = if_else(m_h > 0, N_h / m_h, NA_real_),
        stratum_pi = if_else(N_h > 0, m_h / N_h, NA_real_)
      )

    if (any(is.na(stratum_weights$stratum_weight))) {
      print(stratum_weights)
      stop("At least one negative risk stratum has m_h = 0 under Neyman.")
    }

    sample_data %>%
      left_join(
        stratum_weights %>%
          select(risk_stratum, N_h, m_h, stratum_weight, stratum_pi),
        by = "risk_stratum"
      ) %>%
      mutate(
        design_weight = case_when(
          algo_positive_original == 1L ~ 1,
          algo_positive_original == 0L ~ stratum_weight,
          TRUE ~ NA_real_
        ),
        design_inclusion_probability = case_when(
          algo_positive_original == 1L ~ 1,
          algo_positive_original == 0L ~ stratum_pi,
          TRUE ~ NA_real_
        )
      )
  } else {
    stop("Unknown design label: ", design_label)
  }
}

estimate_design <- function(sample_data, full_analytic, design_label,
                            sampling_design, triggered) {
  N_total <- nrow(full_analytic)
  N_neg <- full_analytic %>% filter(algo_positive_original == 0L) %>% nrow()
  indicator_cols <- c("cell_tp", "cell_fp", "cell_fn", "cell_tn")

  sample_weighted <- sample_data %>%
    add_design_weights(full_analytic, sampling_design) %>%
    add_cell_indicators()

  fixed_census <- full_analytic %>%
    filter(algo_positive_original == 1L) %>%
    add_cell_indicators()
  fixed_totals <- fixed_census %>%
    summarise(
      TP_hat = sum(cell_tp, na.rm = TRUE),
      FP_hat = sum(cell_fp, na.rm = TRUE),
      FN_hat = sum(cell_fn, na.rm = TRUE),
      TN_hat = sum(cell_tn, na.rm = TRUE)
    ) %>%
    unlist(use.names = TRUE)

  neg_sample <- sample_weighted %>%
    filter(algo_positive_original == 0L)
  m_neg <- nrow(neg_sample)
  if (m_neg <= 1L) {
    stop("At least two sampled algorithm-negative patients are required.")
  }

  if (sampling_design == "srs") {
    neg_matrix <- neg_sample %>%
      select(all_of(indicator_cols)) %>%
      as.matrix()
    neg_totals <- N_neg * colMeans(neg_matrix)
    variance_matrix <- N_neg^2 * (1 - m_neg / N_neg) *
      safe_cov_matrix(neg_sample, indicator_cols) / m_neg
    lambda_pool <- NA_real_
    estimator_method <- "Pooled SRS HT"
    variance_method <- "SRSWOR among final sampled algorithm-negatives"
  } else if (sampling_design == "neyman" && triggered) {
    stage1_neg <- neg_sample %>%
      filter(final_sample_component == "stage1_srs_negative")
    stage2_neg <- neg_sample %>%
      filter(str_detect(final_sample_component, "^stage2_neyman"))
    m1_actual <- nrow(stage1_neg)

    if (m1_actual <= 1L || nrow(stage2_neg) <= 1L) {
      stop("Triggered composite estimation requires nontrivial Stage 1 and Stage 2 samples.")
    }

    stratum_frame <- full_analytic %>%
      filter(algo_positive_original == 0L) %>%
      count(risk_stratum, name = "N_h") %>%
      left_join(stage1_neg %>% count(risk_stratum, name = "n1_h"),
                by = "risk_stratum") %>%
      left_join(stage2_neg %>% count(risk_stratum, name = "n2_h"),
                by = "risk_stratum") %>%
      mutate(
        n1_h = replace_na(n1_h, 0L),
        n2_h = replace_na(n2_h, 0L),
        remaining_h = N_h - n1_h
      )

    if (any(stratum_frame$remaining_h > 0L & stratum_frame$n2_h <= 0L)) {
      stop("Triggered design left a nonempty remainder stratum without a Stage 2 sample.")
    }

    stage1_matrix <- stage1_neg %>%
      select(all_of(indicator_cols)) %>%
      as.matrix()
    stage1_ht <- N_neg * colMeans(stage1_matrix)
    stage1_variance <- N_neg^2 * (1 - m1_actual / N_neg) *
      safe_cov_matrix(stage1_neg, indicator_cols) / m1_actual

    sequential_totals <- colSums(stage1_matrix)
    sequential_variance <- matrix(
      0, nrow = length(indicator_cols), ncol = length(indicator_cols),
      dimnames = list(indicator_cols, indicator_cols)
    )

    for (index in seq_len(nrow(stratum_frame))) {
      h <- stratum_frame$risk_stratum[[index]]
      remaining_h <- stratum_frame$remaining_h[[index]]
      n2_h <- stratum_frame$n2_h[[index]]
      if (remaining_h <= 0L) next

      stage2_h <- stage2_neg %>% filter(risk_stratum == h)
      stage2_matrix <- stage2_h %>%
        select(all_of(indicator_cols)) %>%
        as.matrix()
      sequential_totals <- sequential_totals +
        remaining_h * colMeans(stage2_matrix)
      sequential_variance <- sequential_variance +
        remaining_h^2 * (1 - n2_h / remaining_h) *
        safe_cov_matrix(stage2_h, indicator_cols) / n2_h
    }

    lambda_pool <- m1_actual * (N_neg - m_neg) /
      (m_neg * (N_neg - m1_actual))
    neg_totals <- lambda_pool * stage1_ht +
      (1 - lambda_pool) * sequential_totals
    variance_matrix <- lambda_pool^2 * stage1_variance +
      (1 - lambda_pool)^2 * sequential_variance
    estimator_method <- "Lambda-pool composite HT"
    variance_method <- paste0(
      "Fixed-weight composite variance: Stage 1 SRSWOR + ",
      "sequential conditional Stage 2 stratified SRSWOR"
    )
  } else {
    stop("Unsupported design state: sampling_design = ", sampling_design,
         ", triggered = ", triggered, ".")
  }

  names(neg_totals) <- indicator_cols
  names_for_matrix <- c("TP_hat", "FP_hat", "FN_hat", "TN_hat")
  totals <- tibble(
    TP_hat = fixed_totals[["TP_hat"]] + neg_totals[["cell_tp"]],
    FP_hat = fixed_totals[["FP_hat"]] + neg_totals[["cell_fp"]],
    FN_hat = fixed_totals[["FN_hat"]] + neg_totals[["cell_fn"]],
    TN_hat = fixed_totals[["TN_hat"]] + neg_totals[["cell_tn"]]
  )
  if (abs(sum(unlist(totals, use.names = FALSE)) - N_total) > 1e-8) {
    stop("Estimated confusion-matrix totals do not sum to the analytic cohort size.")
  }
  dimnames(variance_matrix) <- list(names_for_matrix, names_for_matrix)

  metrics <- metric_from_totals(
    totals$TP_hat,
    totals$FP_hat,
    totals$FN_hat,
    totals$TN_hat,
    denominator_n = N_total
  ) %>%
    mutate(endpoint_year = endpoint_year,
           endpoint_label = endpoint_label,
           validation_replicate = validation_replicate,
           design = design_label,
           trigger_neyman = triggered,
           stage2_design_used = sampling_design,
           estimator_method = estimator_method,
           lambda_pool = lambda_pool,
           .before = 1)

  weighted_auroc <- if ("predicted_risk" %in% names(sample_weighted)) {
    weighted_auc(
      sample_weighted$true_recur_3yr,
      sample_weighted$predicted_risk,
      sample_weighted$design_weight
    )
  } else {
    NA_real_
  }

  naive_totals <- sample_weighted %>%
    summarise(
      TP = sum(cell_tp, na.rm = TRUE),
      FP = sum(cell_fp, na.rm = TRUE),
      FN = sum(cell_fn, na.rm = TRUE),
      TN = sum(cell_tn, na.rm = TRUE)
    )

  naive_metrics <- metric_from_totals(
    naive_totals$TP,
    naive_totals$FP,
    naive_totals$FN,
    naive_totals$TN,
    denominator_n = nrow(sample_weighted)
  ) %>%
    mutate(endpoint_year = endpoint_year,
           endpoint_label = endpoint_label,
           validation_replicate = validation_replicate,
           design = design_label,
           .before = 1)

  gradients <- metric_gradients(totals, N_total)
  metric_variances <- tibble(
    metric = names(gradients),
    variance = map_dbl(gradients, quadratic_variance, variance_matrix),
    standard_error = sqrt(pmax(variance, 0))
  )

  # PPV is based on a census of predicted-positive patients if every
  # alg_recur_3yr-positive row is among original algorithm positives.
  ppv_den <- full_analytic %>% filter(alg_recur_3yr == 1L) %>% nrow()
  ppv_num <- full_analytic %>%
    filter(alg_recur_3yr == 1L, true_recur_3yr == 1L) %>%
    nrow()
  ppv_ci <- wilson_ci(ppv_num, ppv_den)

  frame_is_algorithm_negative <- full_analytic %>%
    filter(algo_positive_original == 0L) %>%
    summarise(all_negative = all(alg_recur_3yr == 0L)) %>%
    pull(all_negative)
  observed_frame_fn <- sum(neg_sample$cell_fn, na.rm = TRUE)
  zero_fn_fallback_eligible <- sampling_design == "srs" &&
    frame_is_algorithm_negative && observed_frame_fn == 0L
  fn_upper_exact <- if (zero_fn_fallback_eligible) {
    hypergeom_zero_upper(N_neg, m_neg)
  } else {
    NA_integer_
  }

  ci_table <- metrics %>%
    left_join(metric_variances, by = "metric") %>%
    rowwise() %>%
    mutate(
      fallback_used = metric != "ppv" && zero_fn_fallback_eligible,
      fallback_lower = if (fallback_used) {
        metric_interval_from_frame_fn(
          metric, fixed_totals, N_neg, fn_upper_exact, N_total
        )[["lower"]]
      } else {
        NA_real_
      },
      fallback_upper = if (fallback_used) {
        metric_interval_from_frame_fn(
          metric, fixed_totals, N_neg, fn_upper_exact, N_total
        )[["upper"]]
      } else {
        NA_real_
      },
      ci_lower = if (metric == "ppv") {
        ppv_ci[["lower"]]
      } else if (fallback_used) {
        fallback_lower
      } else {
        logit_ci(estimate, variance)[["lower"]]
      },
      ci_upper = if (metric == "ppv") {
        ppv_ci[["upper"]]
      } else if (fallback_used) {
        fallback_upper
      } else {
        logit_ci(estimate, variance)[["upper"]]
      },
      ci_method = case_when(
        metric == "ppv" ~ "Wilson/binomial",
        fallback_used ~ "Finite-population hypergeometric zero-FN safeguard",
        TRUE ~ "Metric-specific logit-Wald"
      ),
      variance_method = if_else(
        metric == "ppv",
        "Wilson/binomial among census of predicted-positive patients",
        variance_method
      )
    ) %>%
    select(-fallback_lower, -fallback_upper) %>%
    ungroup()

  sample_summary <- tibble(
    endpoint_year = endpoint_year,
    endpoint_label = endpoint_label,
    validation_replicate = validation_replicate,
    design = design_label,
    trigger_neyman = triggered,
    stage2_design_used = sampling_design,
    estimator_method = estimator_method,
    lambda_pool = lambda_pool,
    zero_fn_fallback_eligible = zero_fn_fallback_eligible,
    exact_fn_upper = fn_upper_exact,
    N_total = N_total,
    N_algorithm_positive = nrow(full_analytic %>% filter(algo_positive_original == 1L)),
    N_algorithm_negative = N_neg,
    final_sample_n = nrow(sample_weighted),
    final_positive_census_n = sum(sample_weighted$algo_positive_original == 1L),
    final_negative_sample_n = sum(sample_weighted$algo_positive_original == 0L),
    stage1_negative_n = sum(sample_weighted$final_sample_component == "stage1_srs_negative", na.rm = TRUE),
    stage2_negative_n = sum(str_detect(sample_weighted$final_sample_component, "^stage2_"), na.rm = TRUE),
    min_negative_weight = min(sample_weighted$design_weight[sample_weighted$algo_positive_original == 0L]),
    max_negative_weight = max(sample_weighted$design_weight[sample_weighted$algo_positive_original == 0L])
  )

  stratum_summary <- sample_weighted %>%
    filter(algo_positive_original == 0L) %>%
    group_by(risk_stratum) %>%
    summarise(
      endpoint_year = endpoint_year,
      endpoint_label = endpoint_label,
      validation_replicate = validation_replicate,
      design = design_label,
      m_h_final = n(),
      m_h_stage1 = sum(final_sample_component == "stage1_srs_negative", na.rm = TRUE),
      m_h_stage2 = sum(str_detect(final_sample_component, "^stage2_"), na.rm = TRUE),
      observed_fn = sum(cell_fn, na.rm = TRUE),
      observed_tn = sum(cell_tn, na.rm = TRUE),
      observed_fn_rate = observed_fn / m_h_final,
      mean_weight = mean(design_weight, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      full_analytic %>%
        filter(algo_positive_original == 0L) %>%
        count(risk_stratum, name = "N_h"),
      by = "risk_stratum"
    ) %>%
    relocate(endpoint_year, endpoint_label, validation_replicate, design,
             risk_stratum, N_h)

  list(
    sample_weighted = sample_weighted,
    totals = totals %>%
      mutate(endpoint_year = endpoint_year,
             endpoint_label = endpoint_label,
             validation_replicate = validation_replicate,
             design = design_label,
             trigger_neyman = triggered,
             stage2_design_used = sampling_design,
             estimator_method = estimator_method,
             lambda_pool = lambda_pool,
             .before = 1),
    metrics = metrics,
    naive_metrics = naive_metrics,
    ci_table = ci_table,
    auroc = tibble(endpoint_year = endpoint_year,
                   endpoint_label = endpoint_label,
                   validation_replicate = validation_replicate,
                   design = design_label,
                   metric = "auroc",
                   estimate = weighted_auroc),
    sample_summary = sample_summary,
    stratum_summary = stratum_summary,
    variance_matrix = as_tibble(variance_matrix, rownames = "cell")
  )
}

# -------------------------------------------------------------------------
# Read inputs and standardize core variables
# -------------------------------------------------------------------------

endpoint_evaluable_var <- paste0("evaluable_", endpoint_year, "yr_primary")

full_analytic <- read_clean_csv(full_analytic_csv) %>%
  add_endpoint_columns(endpoint_year)

if (endpoint_evaluable_var %in% names(full_analytic)) {
  n_before_endpoint_filter <- nrow(full_analytic)
  full_analytic <- full_analytic %>%
    filter(.data[[endpoint_evaluable_var]] == TRUE)
  n_after_endpoint_filter <- nrow(full_analytic)

  if (n_after_endpoint_filter < n_before_endpoint_filter) {
    message("Filtered ", n_before_endpoint_filter - n_after_endpoint_filter,
            " rows not evaluable for the ", endpoint_label, " endpoint.")
  }
}

full_analytic <- full_analytic %>%
  mutate(
    true_recur_3yr = standardize_binary(.data[[truth_var]]),
    alg_recur_3yr = standardize_binary(.data[[algorithm_endpoint_var]]),
    algo_positive_original =
      standardize_binary(.data[[original_algorithm_positive_var]]),
    risk_stratum = as.integer(.data[[stratum_var]])
  )

required_cols <- c(id_var, "true_recur_3yr", "alg_recur_3yr",
                   "algo_positive_original", "risk_stratum")
missing_cols <- setdiff(required_cols, names(full_analytic))
if (length(missing_cols) > 0) {
  stop("merged_clean_validation_data.csv is missing: ",
       paste(missing_cols, collapse = ", "))
}

if (any(is.na(full_analytic$true_recur_3yr)) ||
    any(is.na(full_analytic$alg_recur_3yr)) ||
    any(is.na(full_analytic$algo_positive_original))) {
  stop("Core binary validation variables contain missing values.")
}

sample_size_m_summary <- if (file.exists(sample_size_m_summary_csv)) {
  read_clean_csv(sample_size_m_summary_csv)
} else {
  tibble()
}

planned_m1 <- if (nrow(sample_size_m_summary) > 0 &&
                  "m1" %in% names(sample_size_m_summary)) {
  as.integer(sample_size_m_summary$m1[[1]])
} else {
  NA_integer_
}
planned_m <- if (nrow(sample_size_m_summary) > 0 &&
                 "planned_total_negative_sample_m" %in%
                 names(sample_size_m_summary)) {
  as.integer(sample_size_m_summary$planned_total_negative_sample_m[[1]])
} else if (nrow(sample_size_m_summary) > 0 &&
           "m_precision_rule_grid" %in% names(sample_size_m_summary)) {
  as.integer(sample_size_m_summary$m_precision_rule_grid[[1]])
} else {
  NA_integer_
}
planned_m2 <- if (nrow(sample_size_m_summary) > 0 &&
                  "planned_stage2_negative_sample_m2" %in%
                  names(sample_size_m_summary)) {
  as.integer(sample_size_m_summary$planned_stage2_negative_sample_m2[[1]])
} else if (!is.na(planned_m) && !is.na(planned_m1)) {
  as.integer(planned_m - planned_m1)
} else {
  NA_integer_
}

rhat <- if (file.exists(rhat_csv)) {
  read_clean_csv(rhat_csv)
} else {
  tibble()
}
rhat_value <- if (nrow(rhat) > 0 && "rhat" %in% names(rhat)) {
  rhat$rhat[[1]]
} else {
  NA_real_
}
trigger_neyman <- if (nrow(rhat) > 0 &&
                       "trigger_neyman" %in% names(rhat)) {
  isTRUE(rhat$trigger_neyman[[1]])
} else {
  !is.na(rhat_value) && rhat_value > 1.1
}

add_design_parameter_columns <- function(data) {
  data %>%
    mutate(
      m1 = .env$planned_m1,
      planned_m = .env$planned_m,
      planned_m2 = .env$planned_m2,
      .after = validation_replicate
    )
}

neyman_sample <- build_final_sample(final_neyman_csv, full_analytic, "neyman")
srs_sample <- build_final_sample(final_srs_csv, full_analytic, "srs")
adaptive_sample <- if (trigger_neyman) neyman_sample else srs_sample
adaptive_sampling_design <- if (trigger_neyman) "neyman" else "srs"

# -------------------------------------------------------------------------
# Estimate final validation performance under each design
# -------------------------------------------------------------------------

neyman_results <- estimate_design(
  adaptive_sample, full_analytic, "neyman",
  sampling_design = adaptive_sampling_design,
  triggered = trigger_neyman
)
srs_results <- estimate_design(
  srs_sample, full_analytic, "srs",
  sampling_design = "srs",
  triggered = FALSE
)

if (!trigger_neyman) {
  adaptive_totals <- unlist(neyman_results$totals[1, c(
    "TP_hat", "FP_hat", "FN_hat", "TN_hat"
  )], use.names = FALSE)
  srs_totals <- unlist(srs_results$totals[1, c(
    "TP_hat", "FP_hat", "FN_hat", "TN_hat"
  )], use.names = FALSE)
  if (max(abs(adaptive_totals - srs_totals)) > 1e-10) {
    stop("Non-triggered adaptive estimates must equal pooled SRS estimates.")
  }
}

if (trigger_neyman &&
    (!is.finite(neyman_results$sample_summary$lambda_pool[[1]]) ||
     neyman_results$sample_summary$lambda_pool[[1]] < 0 ||
     neyman_results$sample_summary$lambda_pool[[1]] > 1)) {
  stop("Triggered lambda-pool weight is outside [0, 1].")
}

all_ci <- bind_rows(neyman_results$ci_table, srs_results$ci_table)
invalid_ci <- all_ci %>%
  filter(!is.na(ci_lower), !is.na(ci_upper)) %>%
  filter(ci_lower < 0 | ci_upper > 1 | ci_lower > ci_upper)
if (nrow(invalid_ci) > 0L) {
  stop("At least one confidence interval is outside [0, 1] or reversed.")
}

weighted_metrics <- bind_rows(neyman_results$metrics, srs_results$metrics) %>%
  add_design_parameter_columns()
naive_metrics <- bind_rows(neyman_results$naive_metrics, srs_results$naive_metrics) %>%
  add_design_parameter_columns()
metrics_ci <- bind_rows(neyman_results$ci_table, srs_results$ci_table) %>%
  add_design_parameter_columns()
weighted_auroc <- bind_rows(neyman_results$auroc, srs_results$auroc) %>%
  add_design_parameter_columns()
confusion_totals <- bind_rows(neyman_results$totals, srs_results$totals) %>%
  add_design_parameter_columns()
sample_summary <- bind_rows(neyman_results$sample_summary, srs_results$sample_summary) %>%
  add_design_parameter_columns()
stratum_summary <- bind_rows(neyman_results$stratum_summary, srs_results$stratum_summary) %>%
  add_design_parameter_columns()

oracle_recurrence_prevalence <- mean(full_analytic$true_recur_3yr == 1L,
                                     na.rm = TRUE)

# Oracle/full-cohort metrics are available only because the working analytic
# file already contains chart-review truth for every analytic row. In a real
# partially verified workflow, this is a benchmark/sanity check, not something
# available before all charts are reviewed.
oracle_totals <- full_analytic %>%
  add_cell_indicators() %>%
  summarise(
    TP = sum(cell_tp),
    FP = sum(cell_fp),
    FN = sum(cell_fn),
    TN = sum(cell_tn)
  )

oracle_metrics <- metric_from_totals(
  oracle_totals$TP,
  oracle_totals$FP,
  oracle_totals$FN,
  oracle_totals$TN,
  denominator_n = nrow(full_analytic)
) %>%
  mutate(endpoint_year = endpoint_year,
         endpoint_label = endpoint_label,
         validation_replicate = validation_replicate,
         design = "oracle_full_analytic",
         .before = 1)

oracle_auroc <- tibble(
  endpoint_year = endpoint_year,
  endpoint_label = endpoint_label,
  validation_replicate = validation_replicate,
  design = "oracle_full_analytic",
  metric = "auroc",
  estimate = if ("predicted_risk" %in% names(full_analytic)) {
    weighted_auc(full_analytic$true_recur_3yr, full_analytic$predicted_risk)
  } else {
    NA_real_
  }
)

design_vs_oracle <- weighted_metrics %>%
  left_join(
    oracle_metrics %>%
      select(endpoint_year, endpoint_label, validation_replicate,
             metric, oracle_estimate = estimate),
    by = c("endpoint_year", "endpoint_label",
           "validation_replicate", "metric")
  ) %>%
  mutate(
    difference_from_oracle = estimate - oracle_estimate,
    absolute_difference_from_oracle = abs(difference_from_oracle)
  ) %>%
  arrange(metric, design)

comparison_wide <- metrics_ci %>%
  select(endpoint_year, endpoint_label, validation_replicate,
         m1, planned_m, planned_m2,
         design, metric, estimate,
         standard_error, ci_lower, ci_upper) %>%
  pivot_wider(
    id_cols = c(endpoint_year, endpoint_label, validation_replicate,
                m1, planned_m, planned_m2, metric),
    names_from = design,
    values_from = c(estimate, standard_error, ci_lower, ci_upper),
    names_glue = "{.value}_{design}"
  ) %>%
  left_join(
    oracle_metrics %>%
      select(endpoint_year, endpoint_label, validation_replicate,
             metric, oracle_estimate = estimate),
    by = c("endpoint_year", "endpoint_label",
           "validation_replicate", "metric")
  ) %>%
  mutate(
    neyman_minus_srs = estimate_neyman - estimate_srs,
    neyman_abs_error = abs(estimate_neyman - oracle_estimate),
    srs_abs_error = abs(estimate_srs - oracle_estimate)
  )

auroc_comparison_wide <- weighted_auroc %>%
  select(endpoint_year, endpoint_label, validation_replicate,
         m1, planned_m, planned_m2,
         design, metric, estimate) %>%
  mutate(
    standard_error = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_
  ) %>%
  pivot_wider(
    id_cols = c(endpoint_year, endpoint_label, validation_replicate,
                m1, planned_m, planned_m2, metric),
    names_from = design,
    values_from = c(estimate, standard_error, ci_lower, ci_upper),
    names_glue = "{.value}_{design}"
  ) %>%
  left_join(
    oracle_auroc %>%
      select(endpoint_year, endpoint_label, validation_replicate,
             metric, oracle_estimate = estimate),
    by = c("endpoint_year", "endpoint_label",
           "validation_replicate", "metric")
  ) %>%
  mutate(
    neyman_minus_srs = estimate_neyman - estimate_srs,
    neyman_abs_error = abs(estimate_neyman - oracle_estimate),
    srs_abs_error = abs(estimate_srs - oracle_estimate)
  )

comparison_wide <- bind_rows(comparison_wide, auroc_comparison_wide) %>%
  mutate(rhat = rhat_value, .after = endpoint_label) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("sensitivity", "specificity", "ppv", "npv",
                 "accuracy", "auroc")
    )
  ) %>%
  arrange(metric) %>%
  mutate(metric = as.character(metric))

metric_summary_long <- metrics_ci %>%
  left_join(
    oracle_metrics %>%
      select(endpoint_year, endpoint_label, validation_replicate,
             metric, oracle_estimate = estimate),
    by = c("endpoint_year", "endpoint_label",
           "validation_replicate", "metric")
  ) %>%
  transmute(
    endpoint_year,
    endpoint_label,
    validation_replicate,
    m1,
    planned_m,
    planned_m2,
    rhat = rhat_value,
    trigger_neyman,
    stage2_design_used,
    estimator_method,
    lambda_pool,
    ci_method,
    zero_fn_fallback_used = fallback_used,
    oracle_recurrence_prevalence = oracle_recurrence_prevalence,
    design,
    metric,
    estimate,
    standard_error,
    ci_lower,
    ci_upper,
    oracle_estimate,
    bias = estimate - oracle_estimate,
    sd = standard_error,
    mse = (estimate - oracle_estimate)^2 + standard_error^2,
    rmse = sqrt(mse),
    ci_width = ci_upper - ci_lower,
    covered_oracle = ci_lower <= oracle_estimate & oracle_estimate <= ci_upper
  )

auroc_summary_long <- weighted_auroc %>%
  left_join(
    oracle_auroc %>%
      select(endpoint_year, endpoint_label, validation_replicate,
             metric, oracle_estimate = estimate),
    by = c("endpoint_year", "endpoint_label",
           "validation_replicate", "metric")
  ) %>%
  transmute(
    endpoint_year,
    endpoint_label,
    validation_replicate,
    m1,
    planned_m,
    planned_m2,
    rhat = rhat_value,
    trigger_neyman = if_else(design == "neyman", .env$trigger_neyman, FALSE),
    stage2_design_used = if_else(
      design == "neyman", .env$adaptive_sampling_design, "srs"
    ),
    estimator_method = "Descriptive pooled weighted AUROC",
    lambda_pool = if_else(design == "neyman" & .env$trigger_neyman,
                          neyman_results$sample_summary$lambda_pool[[1]],
                          NA_real_),
    ci_method = NA_character_,
    zero_fn_fallback_used = FALSE,
    oracle_recurrence_prevalence = oracle_recurrence_prevalence,
    design,
    metric,
    estimate,
    standard_error = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    oracle_estimate,
    bias = estimate - oracle_estimate,
    sd = NA_real_,
    mse = (estimate - oracle_estimate)^2,
    rmse = sqrt(mse),
    ci_width = NA_real_,
    covered_oracle = NA
  )

endpoint_summary_long <- bind_rows(metric_summary_long, auroc_summary_long) %>%
  arrange(endpoint_year, m1, design, metric)

endpoint_summary_wide <- endpoint_summary_long %>%
  pivot_wider(
    id_cols = c(endpoint_year, endpoint_label, validation_replicate,
                m1, planned_m, planned_m2, rhat, trigger_neyman,
                stage2_design_used, estimator_method, lambda_pool,
                oracle_recurrence_prevalence, design),
    names_from = metric,
    values_from = c(estimate, standard_error, ci_lower, ci_upper,
                    oracle_estimate, bias, sd, mse, rmse, ci_width,
                    covered_oracle),
    names_glue = "{metric}_{.value}"
  ) %>%
  arrange(endpoint_year, m1, design)

stage2_overlap <- tibble(
  endpoint_year = endpoint_year,
  endpoint_label = endpoint_label,
  validation_replicate = validation_replicate,
  m1 = planned_m1,
  planned_m = planned_m,
  planned_m2 = planned_m2,
  neyman_final_negative_n =
    neyman_results$sample_weighted %>% filter(algo_positive_original == 0L) %>% nrow(),
  srs_final_negative_n =
    srs_results$sample_weighted %>% filter(algo_positive_original == 0L) %>% nrow(),
  overlapping_final_negative_n =
    length(intersect(
      neyman_results$sample_weighted %>%
        filter(algo_positive_original == 0L) %>%
        pull(patient_id),
      srs_results$sample_weighted %>%
        filter(algo_positive_original == 0L) %>%
        pull(patient_id)
    )),
  neyman_stage2_negative_n =
    neyman_results$sample_weighted %>%
    filter(final_sample_component == "stage2_neyman_negative") %>%
    nrow(),
  srs_stage2_negative_n =
    srs_results$sample_weighted %>%
    filter(final_sample_component == "stage2_srs_negative") %>%
    nrow(),
  overlapping_stage2_negative_n =
    length(intersect(
      neyman_results$sample_weighted %>%
        filter(final_sample_component == "stage2_neyman_negative") %>%
        pull(patient_id),
      srs_results$sample_weighted %>%
        filter(final_sample_component == "stage2_srs_negative") %>%
        pull(patient_id)
    ))
)

allocation <- if (file.exists(allocation_csv)) {
  read_clean_csv(allocation_csv)
} else {
  tibble()
}

# -------------------------------------------------------------------------
# Save outputs
# -------------------------------------------------------------------------

write_csv(neyman_results$sample_weighted,
          file.path(output_dir, "final_validation_sample_neyman_weighted.csv"))
write_csv(srs_results$sample_weighted,
          file.path(output_dir, "final_validation_sample_srs_weighted.csv"))
write_csv(weighted_metrics,
          file.path(output_dir, "stage2_weighted_metrics_by_design.csv"))
write_csv(naive_metrics,
          file.path(output_dir, "stage2_naive_metrics_by_design.csv"))
write_csv(metrics_ci,
          file.path(output_dir, "stage2_metrics_ci_by_design.csv"))
write_csv(weighted_auroc,
          file.path(output_dir, "stage2_weighted_auroc_by_design.csv"))
write_csv(confusion_totals,
          file.path(output_dir, "stage2_weighted_confusion_totals_by_design.csv"))
write_csv(sample_summary,
          file.path(output_dir, "stage2_sample_summary_by_design.csv"))
write_csv(stratum_summary,
          file.path(output_dir, "stage2_negative_stratum_summary_by_design.csv"))
write_csv(oracle_metrics,
          file.path(output_dir, "oracle_full_analytic_metrics.csv"))
write_csv(oracle_auroc,
          file.path(output_dir, "oracle_full_analytic_auroc.csv"))
write_csv(design_vs_oracle,
          file.path(output_dir, "stage2_design_vs_oracle.csv"))
write_csv(comparison_wide,
          file.path(output_dir, "stage2_neyman_vs_srs_comparison_wide.csv"))
write_csv(endpoint_summary_long,
          file.path(output_dir, "stage2_endpoint_summary_long.csv"))
write_csv(endpoint_summary_wide,
          file.path(output_dir, "stage2_endpoint_summary_wide.csv"))
write_csv(stage2_overlap,
          file.path(output_dir, "stage2_sample_overlap_summary.csv"))
write_csv(neyman_results$variance_matrix,
          file.path(output_dir, "variance_matrix_neyman.csv"))
write_csv(srs_results$variance_matrix,
          file.path(output_dir, "variance_matrix_srs.csv"))

if (nrow(allocation) > 0) {
  write_csv(allocation,
            file.path(output_dir, "stage2_allocation_recommendation_input.csv"))
}
if (nrow(rhat) > 0) {
  write_csv(rhat,
            file.path(output_dir, "projected_rhat_input.csv"))
}

# -------------------------------------------------------------------------
# Print compact report
# -------------------------------------------------------------------------

cat("\nStage 2 sample summary by design\n")
print(sample_summary, width = Inf)

cat("\nWeighted final validation metrics with design-based CIs\n")
print(metrics_ci, n = Inf, width = Inf)

cat("\nNeyman vs SRS comparison, wide table\n")
print(comparison_wide, n = Inf, width = Inf)

cat("\nOracle full analytic metrics, available only because truth exists in this working file\n")
print(oracle_metrics, n = Inf)

cat("\nEndpoint summary long table\n")
print(endpoint_summary_long, n = Inf, width = Inf)

cat("\nNegative sample overlap summary\n")
print(stage2_overlap, width = Inf)

cat("\nNegative stratum summary by design\n")
print(stratum_summary, n = Inf, width = Inf)

cat("\nSaved Stage 2 comparison outputs to: ",
    normalizePath(output_dir), "\n", sep = "")
