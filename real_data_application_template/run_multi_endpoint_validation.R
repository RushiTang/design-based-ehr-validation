############################################################
# Run Phase 1 -> Phase 2 validation for multiple prespecified binary endpoints
# and optional repeated-sampling replicates.
#
# A single replicate implements the proposed validation design. More than one
# replicate is appropriate only for a methodological audit when gold-standard
# outcomes are known for the complete cohort; ordinary partially verified EHR
# applications should use VALIDATION_REPLICATES=1.
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- "--file="
script_file_match <- script_args[startsWith(script_args, script_file_arg)]
script_path <- if (length(script_file_match) > 0) {
  sub(script_file_arg, "", script_file_match[[1]])
} else {
  NA_character_
}
script_root <- if (!is.na(script_path)) {
  dirname(normalizePath(script_path, winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}
if (!is.na(script_path)) {
  setwd(script_root)
}

analysis_root <- normalizePath(
  Sys.getenv("VALIDATION_PROJECT_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)
setwd(analysis_root)
message("Using analysis folder: ", analysis_root)

parse_endpoint_years <- function(x) {
  out <- strsplit(x, ",", fixed = TRUE)[[1]] %>%
    str_trim() %>%
    discard(~ .x == "") %>%
    as.integer()

  if (length(out) == 0 || any(is.na(out)) || any(out <= 0)) {
    stop("VALIDATION_ENDPOINT_YEARS must be a comma-separated list of positive integers.")
  }

  sort(unique(out))
}

parse_positive_integer_values <- function(x, env_name) {
  out <- strsplit(x, ",", fixed = TRUE)[[1]] %>%
    str_trim() %>%
    discard(~ .x == "") %>%
    as.integer()

  if (length(out) == 0 || any(is.na(out)) || any(out <= 0)) {
    stop(env_name, " must be a comma-separated list of positive integers.")
  }

  sort(unique(out))
}

parse_positive_integer <- function(x, env_name) {
  out <- as.integer(x)
  if (is.na(out) || out <= 0) {
    stop(env_name, " must be a positive integer.")
  }
  out
}

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

sd_or_na <- function(x) {
  if (sum(!is.na(x)) <= 1) NA_real_ else sd(x, na.rm = TRUE)
}

quantile_or_na <- function(x, prob) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  unname(quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 1))
}

first_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA)
  }

  x[[1]]
}

normalize_output_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

resolve_script_path <- function(script) {
  candidates <- c(file.path(script_root, script), script)
  candidates <- candidates[file.exists(candidates)]

  if (length(candidates) == 0) {
    stop(
      "Could not find ", script, " in script folder: ", script_root,
      "\nPut run_multi_endpoint_validation.R, stage1_validation_analysis.R, ",
      "and stage2_compare_designs.R in the same folder, or set ",
      "VALIDATION_PROJECT_ROOT to that folder before sourcing this script."
    )
  }

  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

check_analysis_folder <- function() {
  missing_paths <- character()

  truth_csv <- Sys.getenv("TRUTH_CSV", unset = file.path("input", "analytic_cohort.csv"))
  algorithm_csv <- Sys.getenv(
    "ALGORITHM_CSV",
    unset = file.path("input", "algorithm_predictions.csv")
  )

  if (!file.exists(truth_csv)) {
    missing_paths <- c(missing_paths, truth_csv)
  }
  if (!file.exists(algorithm_csv)) {
    missing_paths <- c(missing_paths, algorithm_csv)
  }

  if (length(missing_paths) > 0) {
    stop(
      "The analysis folder is missing required input(s):\n- ",
      paste(missing_paths, collapse = "\n- "),
      "\nCurrent analysis folder: ",
      getwd(),
      "\nPlace de-identified input files in input/, or set TRUTH_CSV and ",
      "ALGORITHM_CSV to their local locations. See input/README.md."
    )
  }
}

endpoint_years <- parse_endpoint_years(
  Sys.getenv("VALIDATION_ENDPOINT_YEARS", unset = "3,4,5")
)
m1_values <- parse_positive_integer_values(
  Sys.getenv("VALIDATION_M1_VALUES", unset = "50,75,100,150"),
  "VALIDATION_M1_VALUES"
)
validation_replicates <- parse_positive_integer(
  Sys.getenv("VALIDATION_REPLICATES", unset = "200"),
  "VALIDATION_REPLICATES"
)
validation_workers <- parse_positive_integer(
  Sys.getenv("VALIDATION_WORKERS", unset = "1"),
  "VALIDATION_WORKERS"
)
base_seed <- parse_positive_integer(
  Sys.getenv("VALIDATION_BASE_SEED", unset = "20260523"),
  "VALIDATION_BASE_SEED"
)

output_root <- Sys.getenv("MULTI_ENDPOINT_OUTPUT_ROOT",
                          unset = "multi_endpoint_validation_outputs")
check_analysis_folder()
stage1_script <- resolve_script_path("stage1_validation_analysis.R")
stage2_script <- resolve_script_path("stage2_compare_designs.R")
plot_output_dir <- file.path(output_root, "endpoint_replicate_summary_plots")
log_output_dir <- file.path(output_root, "run_logs")

dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_output_dir, showWarnings = FALSE, recursive = TRUE)

parse_env_vector <- function(env) {
  pieces <- strsplit(env, "=", fixed = TRUE)
  keys <- map_chr(pieces, ~ .x[[1]])
  values <- map_chr(pieces, ~ paste(.x[-1], collapse = "="))

  stats::setNames(as.list(values), keys)
}

with_env_vars <- function(env, code) {
  env_list <- parse_env_vector(env)
  env_names <- names(env_list)
  old_values <- Sys.getenv(env_names, unset = NA_character_)

  on.exit({
    for (env_name in env_names) {
      old_value <- old_values[[env_name]]
      if (is.na(old_value)) {
        Sys.unsetenv(env_name)
      } else {
        do.call(Sys.setenv, stats::setNames(list(old_value), env_name))
      }
    }
  }, add = TRUE)

  for (env_name in env_names) {
    do.call(Sys.setenv, stats::setNames(list(env_list[[env_name]]), env_name))
  }

  force(code)
}

with_working_dir <- function(path, code) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(path)
  force(code)
}

tail_log <- function(log_path, n = 40L) {
  if (!file.exists(log_path)) {
    return("Log file was not created.")
  }

  lines <- readLines(log_path, warn = FALSE)
  if (length(lines) == 0) {
    return("Log file is empty.")
  }

  paste(tail(lines, n), collapse = "\n")
}

source_script_to_log <- function(script, env, log_path) {
  dir.create(dirname(log_path), showWarnings = FALSE, recursive = TRUE)

  if (file.exists(log_path)) {
    file.remove(log_path)
  }

  log_con <- file(log_path, open = "wt")
  output_sunk <- FALSE
  message_sunk <- FALSE
  closed <- FALSE

  close_log <- function() {
    if (closed) {
      return(invisible(NULL))
    }
    if (message_sunk) {
      sink(type = "message")
      message_sunk <<- FALSE
    }
    if (output_sunk) {
      sink(type = "output")
      output_sunk <<- FALSE
    }
    close(log_con)
    closed <<- TRUE
    invisible(NULL)
  }

  on.exit(close_log(), add = TRUE)

  tryCatch({
    sink(log_con, type = "output")
    output_sunk <- TRUE
    sink(log_con, type = "message")
    message_sunk <- TRUE

    with_env_vars(env, {
      with_working_dir(dirname(script), {
        sys.source(
          basename(script),
          envir = new.env(parent = globalenv()),
          keep.source = FALSE
        )
      })
    })

    close_log()
    0L
  }, error = function(e) {
    cat("\nError while sourcing ", script, ":\n", sep = "")
    cat(conditionMessage(e), "\n", sep = "")
    close_log()
    1L
  })
}

run_rscript <- function(script, env, log_path) {
  status <- source_script_to_log(script, env, log_path)

  if (!identical(status, 0L)) {
    stop(script, " failed with exit status ", status,
         ". See log: ", normalizePath(log_path, mustWork = FALSE),
         "\n\nLast lines from log:\n", tail_log(log_path))
  }
}

read_replicate_summary <- function(summary_path,
                                   endpoint_year,
                                   endpoint_label,
                                   m1_value,
                                   replicate_index,
                                   replicate_seed) {
  if (!file.exists(summary_path)) {
    stop("Expected summary file was not created: ", summary_path)
  }

  dat <- readr::read_csv(summary_path, show_col_types = FALSE) %>%
    janitor::clean_names()

  if (!"validation_replicate" %in% names(dat)) {
    dat$validation_replicate <- replicate_index
  }
  if (!"m1" %in% names(dat)) {
    dat$m1 <- m1_value
  }
  if (!"planned_m" %in% names(dat)) {
    dat$planned_m <- NA_integer_
  }
  if (!"planned_m2" %in% names(dat)) {
    dat$planned_m2 <- NA_integer_
  }
  if (!"rmse" %in% names(dat) && "mse" %in% names(dat)) {
    dat$rmse <- sqrt(dat$mse)
  }
  if (!"oracle_recurrence_prevalence" %in% names(dat)) {
    dat$oracle_recurrence_prevalence <- NA_real_
  }
  if (!"standard_error" %in% names(dat) && "sd" %in% names(dat)) {
    dat$standard_error <- dat$sd
  }
  if (!"ci_width" %in% names(dat) &&
      all(c("ci_lower", "ci_upper") %in% names(dat))) {
    dat$ci_width <- dat$ci_upper - dat$ci_lower
  }
  if (!"covered_oracle" %in% names(dat) &&
      all(c("ci_lower", "ci_upper", "oracle_estimate") %in% names(dat))) {
    dat$covered_oracle <- with(
      dat,
      ifelse(
        is.na(ci_lower) | is.na(ci_upper) | is.na(oracle_estimate),
        NA,
        ci_lower <= oracle_estimate & oracle_estimate <= ci_upper
      )
    )
  }
  if (!"trigger_neyman" %in% names(dat)) {
    dat$trigger_neyman <- NA
  }
  if (!"stage2_design_used" %in% names(dat)) {
    dat$stage2_design_used <- NA_character_
  }
  if (!"estimator_method" %in% names(dat)) {
    dat$estimator_method <- NA_character_
  }
  if (!"ci_method" %in% names(dat)) {
    dat$ci_method <- NA_character_
  }
  if (!"zero_fn_fallback_used" %in% names(dat)) {
    dat$zero_fn_fallback_used <- FALSE
  }
  if (!"lambda_pool" %in% names(dat)) {
    dat$lambda_pool <- NA_real_
  }

  dat %>%
    mutate(
      endpoint_year = endpoint_year,
      endpoint_label = endpoint_label,
      m1 = m1_value,
      validation_replicate = replicate_index,
      validation_seed = replicate_seed,
      .before = 1
    )
}

run_one_replicate <- function(endpoint_year, m1_value, replicate_index) {
  endpoint_label <- paste0(endpoint_year, "yr")
  m1_label <- paste0("m1_", stringr::str_pad(m1_value, width = 3, pad = "0"))
  replicate_label <- sprintf("replicate_%03d", replicate_index)
  endpoint_root <- file.path(output_root, m1_label, endpoint_label, replicate_label)
  stage1_output_dir <- file.path(endpoint_root, "stage1_validation_outputs")
  stage2_output_dir <- file.path(endpoint_root, "stage2_validation_outputs")
  replicate_seed <- base_seed + endpoint_year * 100000L +
    m1_value * 1000L + replicate_index

  dir.create(stage1_output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(stage2_output_dir, showWarnings = FALSE, recursive = TRUE)

  stage1_output_dir <- normalize_output_path(stage1_output_dir)
  stage2_output_dir <- normalize_output_path(stage2_output_dir)

  env <- c(
    paste0("VALIDATION_ENDPOINT_YEAR=", endpoint_year),
    paste0("VALIDATION_M1=", m1_value),
    paste0("VALIDATION_REPLICATE=", replicate_index),
    paste0("VALIDATION_SEED=", replicate_seed),
    paste0("STAGE1_OUTPUT_DIR=", stage1_output_dir),
    paste0("STAGE2_OUTPUT_DIR=", stage2_output_dir)
  )

  log_prefix <- paste(m1_label, endpoint_label, replicate_label, sep = "_")
  stage1_log <- file.path(log_output_dir, paste0(log_prefix, "_stage1.log"))
  stage2_log <- file.path(log_output_dir, paste0(log_prefix, "_stage2.log"))

  cat("Running ", endpoint_label, " ", m1_label, " ", replicate_label,
      " with seed ", replicate_seed, "\n", sep = "")

  run_rscript(stage1_script, env, stage1_log)
  run_rscript(stage2_script, env, stage2_log)

  read_replicate_summary(
    file.path(stage2_output_dir, "stage2_endpoint_summary_long.csv"),
    endpoint_year = endpoint_year,
    endpoint_label = endpoint_label,
    m1_value = m1_value,
    replicate_index = replicate_index,
    replicate_seed = replicate_seed
  )
}

replicate_grid <- tidyr::expand_grid(
  endpoint_year = endpoint_years,
  m1 = m1_values,
  validation_replicate = seq_len(validation_replicates)
)

cat("\nRunning endpoints ",
    paste(paste0(endpoint_years, "yr"), collapse = ", "),
    " with m1 values ",
    paste(m1_values, collapse = ", "),
    " and ", validation_replicates,
    " replicates per endpoint/m1 setting using ", validation_workers,
    " worker(s).\n", sep = "")

run_grid_row <- function(index) {
  tryCatch(
    list(
      result = run_one_replicate(
        replicate_grid$endpoint_year[[index]],
        replicate_grid$m1[[index]],
        replicate_grid$validation_replicate[[index]]
      ),
      error = NULL
    ),
    error = function(e) list(result = NULL, error = conditionMessage(e))
  )
}

if (validation_workers > 1L && .Platform$OS.type != "unix") {
  warning("Parallel workers require a Unix-like system; using one worker.")
  validation_workers <- 1L
}

grid_results <- if (validation_workers > 1L) {
  parallel::mclapply(
    seq_len(nrow(replicate_grid)),
    run_grid_row,
    mc.cores = validation_workers,
    mc.preschedule = FALSE
  )
} else {
  lapply(seq_len(nrow(replicate_grid)), run_grid_row)
}

failed_indices <- which(vapply(
  grid_results, function(x) !is.null(x$error), logical(1)
))
if (length(failed_indices) > 0L) {
  failure_text <- vapply(failed_indices, function(index) {
    row <- replicate_grid[index, ]
    paste0(
      row$endpoint_year, "yr, m1=", row$m1,
      ", replicate=", row$validation_replicate,
      ": ", grid_results[[index]]$error
    )
  }, character(1))
  stop("One or more validation runs failed:\n- ",
       paste(failure_text, collapse = "\n- "))
}

multi_endpoint_replicate_estimates_long <- bind_rows(
  lapply(grid_results, `[[`, "result")
) %>%
  arrange(endpoint_year, m1, validation_replicate, design, metric)

if (!"rhat" %in% names(multi_endpoint_replicate_estimates_long)) {
  multi_endpoint_replicate_estimates_long$rhat <- NA_real_
}

monte_carlo_summary_long <- multi_endpoint_replicate_estimates_long %>%
  group_by(endpoint_year, endpoint_label, m1, design, metric) %>%
  summarise(
    planned_m = first_or_na(planned_m),
    planned_m2 = first_or_na(planned_m2),
    n_replicates = n_distinct(validation_replicate),
    n_estimates = sum(!is.na(estimate)),
    oracle_recurrence_prevalence =
      mean_or_na(oracle_recurrence_prevalence),
    oracle_estimate = mean_or_na(oracle_estimate),
    mean_estimate = mean_or_na(estimate),
    empirical_bias = mean_or_na(estimate - oracle_estimate),
    empirical_sd = sd_or_na(estimate),
    empirical_mse = mean_or_na((estimate - oracle_estimate)^2),
    empirical_rmse = sqrt(empirical_mse),
    mean_single_replicate_rmse = mean_or_na(rmse),
    mc_ci_lower = quantile_or_na(estimate, 0.025),
    mc_ci_upper = quantile_or_na(estimate, 0.975),
    mean_single_replicate_se = mean_or_na(standard_error),
    mean_single_replicate_ci_width = mean_or_na(ci_width),
    coverage_n = sum(!is.na(covered_oracle)),
    empirical_coverage = mean_or_na(as.numeric(covered_oracle)),
    mean_rhat = mean_or_na(rhat),
    adaptive_trigger_rate = mean_or_na(as.numeric(trigger_neyman)),
    zero_fn_fallback_rate = mean_or_na(as.numeric(zero_fn_fallback_used)),
    mean_lambda_pool_when_triggered = mean_or_na(
      ifelse(trigger_neyman, lambda_pool, NA_real_)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    normal_ci_lower = mean_estimate - 1.96 * empirical_sd,
    normal_ci_upper = mean_estimate + 1.96 * empirical_sd,
    monte_carlo_se_mean = empirical_sd / sqrt(n_estimates),
    mean_estimate_ci_lower = mean_estimate - 1.96 * monte_carlo_se_mean,
    mean_estimate_ci_upper = mean_estimate + 1.96 * monte_carlo_se_mean,
    mc_bias_ci_lower = mc_ci_lower - oracle_estimate,
    mc_bias_ci_upper = mc_ci_upper - oracle_estimate
  ) %>%
  arrange(endpoint_year, m1, design, metric)

monte_carlo_summary_wide <- monte_carlo_summary_long %>%
  pivot_wider(
    id_cols = c(endpoint_year, endpoint_label, m1, planned_m, planned_m2,
                oracle_recurrence_prevalence),
    names_from = c(design, metric),
    values_from = c(
      oracle_estimate,
      mean_estimate,
      empirical_bias,
      empirical_sd,
      empirical_mse,
      empirical_rmse,
      mean_single_replicate_rmse,
      mc_ci_lower,
      mc_ci_upper,
      monte_carlo_se_mean,
      mean_estimate_ci_lower,
      mean_estimate_ci_upper,
      mean_single_replicate_se,
      mean_single_replicate_ci_width,
      empirical_coverage,
      mean_rhat
    ),
    names_glue = "{design}_{metric}_{.value}"
  ) %>%
  arrange(endpoint_year, m1)

replicate_long_path <- file.path(
  output_root,
  "multi_endpoint_replicate_estimates_long.csv"
)
summary_long_path <- file.path(
  output_root,
  "multi_endpoint_monte_carlo_summary_long.csv"
)
summary_wide_path <- file.path(
  output_root,
  "multi_endpoint_monte_carlo_summary_wide.csv"
)

readr::write_csv(multi_endpoint_replicate_estimates_long, replicate_long_path)
readr::write_csv(monte_carlo_summary_long, summary_long_path)
readr::write_csv(monte_carlo_summary_wide, summary_wide_path)

metric_labels <- c(
  sensitivity = "Sensitivity",
  specificity = "Specificity",
  ppv = "PPV",
  npv = "NPV",
  accuracy = "Accuracy",
  auroc = "AUROC"
)
design_labels <- c(
  neyman = "Adaptive Neyman",
  srs = "SRS"
)

reference_m1 <- if (100L %in% m1_values) 100L else m1_values[[1]]

plot_df <- monte_carlo_summary_long %>%
  filter(metric %in% names(metric_labels)) %>%
  mutate(
    endpoint_label = factor(endpoint_label, levels = paste0(endpoint_years, "yr")),
    m1_label = factor(paste0("m1=", m1),
                      levels = paste0("m1=", m1_values)),
    metric_label = factor(metric_labels[metric],
                          levels = unname(metric_labels)),
    design_label = recode(design, !!!design_labels)
  )

endpoint_line_df <- plot_df %>%
  filter(m1 == reference_m1)

pd <- position_dodge(width = 0.35)
has_multiple_endpoints <- length(endpoint_years) > 1
estimate_design_line <- if (has_multiple_endpoints) {
  geom_line(linewidth = 0.8, position = pd)
} else {
  NULL
}
plain_design_line <- if (has_multiple_endpoints) {
  geom_line(linewidth = 0.8)
} else {
  NULL
}
estimate_oracle_line <- if (has_multiple_endpoints) {
  geom_line(
    aes(x = endpoint_label, y = oracle_estimate, group = 1),
    color = "gray30",
    linewidth = 0.7,
    linetype = "dashed",
    inherit.aes = FALSE,
    data = endpoint_line_df %>%
      distinct(endpoint_label, metric_label, oracle_estimate)
  )
} else {
  NULL
}

estimate_plot <- ggplot(
  endpoint_line_df,
  aes(x = endpoint_label, y = mean_estimate,
      color = design_label, group = design_label)
) +
  geom_errorbar(aes(ymin = mc_ci_lower, ymax = mc_ci_upper),
                width = 0.12, linewidth = 0.5, position = pd) +
  estimate_design_line +
  geom_point(size = 2.2, position = pd) +
  estimate_oracle_line +
  geom_point(
    aes(x = endpoint_label, y = oracle_estimate),
    color = "gray30",
    size = 1.8,
    inherit.aes = FALSE,
    data = endpoint_line_df %>%
      distinct(endpoint_label, metric_label, oracle_estimate)
  ) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Adaptive Neyman" = "#1B7F79", "SRS" = "#B85C38")) +
  labs(
    x = "Endpoint",
    y = "Mean estimate across replicates",
    color = "Design",
    title = paste0("Mean Validation Performance Across ",
                   validation_replicates, " Real-Data Replicates"),
    subtitle = paste0("Reference m1 = ", reference_m1,
                      "; error bars show empirical 2.5th to 97.5th percentiles; dashed gray line is the full-cohort oracle")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot"
  )

bias_plot <- ggplot(
  endpoint_line_df,
  aes(x = endpoint_label, y = empirical_bias,
      color = design_label, group = design_label)
) +
  geom_hline(yintercept = 0, color = "gray35", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mc_bias_ci_lower, ymax = mc_bias_ci_upper),
                width = 0.12, linewidth = 0.5, position = pd) +
  estimate_design_line +
  geom_point(size = 2.2, position = pd) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Adaptive Neyman" = "#1B7F79", "SRS" = "#B85C38")) +
  labs(
    x = "Endpoint",
    y = "Mean estimate minus oracle",
    color = "Design",
    title = "Empirical Bias Across Real-Data Replicates",
    subtitle = paste0("Reference m1 = ", reference_m1,
                      "; error bars show empirical replicate percentiles after subtracting the oracle")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot"
  )

sd_plot <- ggplot(
  endpoint_line_df,
  aes(x = endpoint_label, y = empirical_sd,
      color = design_label, group = design_label)
) +
  plain_design_line +
  geom_point(size = 2.2) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Adaptive Neyman" = "#1B7F79", "SRS" = "#B85C38")) +
  labs(
    x = "Endpoint",
    y = "Empirical SD across replicates",
    color = "Design",
    title = "Empirical Sampling SD Across Real-Data Replicates",
    subtitle = paste0("Reference m1 = ", reference_m1)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot"
  )

rmse_plot <- ggplot(
  endpoint_line_df,
  aes(x = endpoint_label, y = empirical_rmse,
      color = design_label, group = design_label)
) +
  plain_design_line +
  geom_point(size = 2.2) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Adaptive Neyman" = "#1B7F79", "SRS" = "#B85C38")) +
  labs(
    x = "Endpoint",
    y = "Empirical RMSE across replicates",
    color = "Design",
    title = "Empirical RMSE Across Real-Data Replicates",
    subtitle = paste0("Reference m1 = ", reference_m1)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot"
  )

ci_width_plot <- endpoint_line_df %>%
  filter(!is.na(mean_single_replicate_ci_width), metric != "auroc") %>%
  ggplot(aes(x = endpoint_label, y = mean_single_replicate_ci_width,
             color = design_label, group = design_label)) +
  plain_design_line +
  geom_point(size = 2.2) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Adaptive Neyman" = "#1B7F79", "SRS" = "#B85C38")) +
  labs(
    x = "Endpoint",
    y = "Average single-replicate 95% CI width",
    color = "Design",
    title = "Average Design-Based CI Width Across Replicates",
    subtitle = paste0("Reference m1 = ", reference_m1)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot"
  )

coverage_plot <- endpoint_line_df %>%
  filter(!is.na(empirical_coverage), metric != "auroc") %>%
  ggplot(aes(x = endpoint_label, y = empirical_coverage,
             color = design_label, group = design_label)) +
  geom_hline(yintercept = 0.95, color = "gray35", linewidth = 0.4,
             linetype = "dashed") +
  plain_design_line +
  geom_point(size = 2.2) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Adaptive Neyman" = "#1B7F79", "SRS" = "#B85C38")) +
  labs(
    x = "Endpoint",
    y = "Empirical coverage",
    color = "Design",
    title = "Single-Replicate CI Coverage Across Replicates",
    subtitle = paste0("Reference m1 = ", reference_m1,
                      "; dashed gray line marks 95% nominal coverage")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot"
  )

heatmap_theme <- theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 9),
    plot.title.position = "plot"
  )

bias_heatmap <- plot_df %>%
  ggplot(aes(x = m1_label, y = endpoint_label, fill = empirical_bias)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = if_else(is.na(empirical_bias), "",
                                sprintf("%.3f", empirical_bias))),
            size = 2.4) +
  facet_grid(metric_label ~ design_label) +
  scale_x_discrete(labels = function(x) str_replace(x, "^m1=", "")) +
  scale_fill_gradient2(low = "#B85C38", mid = "white", high = "#1B7F79",
                       midpoint = 0, na.value = "grey90") +
  labs(
    x = "Stage 1 negative pilot size m1",
    y = "Endpoint",
    fill = "Bias",
    title = "Empirical Bias by Endpoint and m1",
    subtitle = "Bias is mean estimate minus full-cohort oracle across replicates"
  ) +
  heatmap_theme

rmse_heatmap <- plot_df %>%
  ggplot(aes(x = m1_label, y = endpoint_label, fill = empirical_rmse)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = if_else(is.na(empirical_rmse), "",
                                sprintf("%.3f", empirical_rmse))),
            size = 2.4) +
  facet_grid(metric_label ~ design_label) +
  scale_x_discrete(labels = function(x) str_replace(x, "^m1=", "")) +
  scale_fill_gradient(low = "white", high = "#08519C", na.value = "grey90") +
  labs(
    x = "Stage 1 negative pilot size m1",
    y = "Endpoint",
    fill = "RMSE",
    title = "Empirical RMSE by Endpoint and m1",
    subtitle = "RMSE is sqrt(mean((estimate - oracle)^2)) across replicates"
  ) +
  heatmap_theme

coverage_heatmap <- plot_df %>%
  filter(metric != "auroc", !is.na(empirical_coverage)) %>%
  ggplot(aes(x = m1_label, y = endpoint_label, fill = empirical_coverage)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", empirical_coverage)), size = 2.4) +
  facet_grid(metric_label ~ design_label) +
  scale_x_discrete(labels = function(x) str_replace(x, "^m1=", "")) +
  scale_fill_gradient2(low = "#B85C38", mid = "white", high = "#1B7F79",
                       midpoint = 0.95, na.value = "grey90") +
  labs(
    x = "Stage 1 negative pilot size m1",
    y = "Endpoint",
    fill = "Coverage",
    title = "Design-Based CI Coverage by Endpoint and m1",
    subtitle = "Coverage is the fraction of replicate CIs covering the full-cohort oracle"
  ) +
  heatmap_theme

design_difference_df <- monte_carlo_summary_long %>%
  filter(metric %in% names(metric_labels)) %>%
  select(endpoint_year, endpoint_label, m1, metric, design,
         empirical_bias, empirical_rmse, empirical_coverage) %>%
  pivot_wider(
    names_from = design,
    values_from = c(empirical_bias, empirical_rmse, empirical_coverage),
    names_glue = "{.value}_{design}"
  ) %>%
  mutate(
    endpoint_label = factor(endpoint_label, levels = paste0(endpoint_years, "yr")),
    m1_label = factor(paste0("m1=", m1),
                      levels = paste0("m1=", m1_values)),
    metric_label = factor(metric_labels[metric],
                          levels = unname(metric_labels)),
    rmse_neyman_minus_srs = empirical_rmse_neyman - empirical_rmse_srs,
    abs_bias_neyman_minus_srs =
      abs(empirical_bias_neyman) - abs(empirical_bias_srs),
    coverage_neyman_minus_srs =
      empirical_coverage_neyman - empirical_coverage_srs
  )

rmse_difference_heatmap <- design_difference_df %>%
  ggplot(aes(x = m1_label, y = endpoint_label,
             fill = rmse_neyman_minus_srs)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = if_else(is.na(rmse_neyman_minus_srs), "",
                                sprintf("%.3f", rmse_neyman_minus_srs))),
            size = 2.6) +
  facet_wrap(~ metric_label, scales = "free", ncol = 3) +
  scale_x_discrete(labels = function(x) str_replace(x, "^m1=", "")) +
  scale_fill_gradient2(low = "#1B7F79", mid = "white", high = "#B85C38",
                       midpoint = 0, na.value = "grey90") +
  labs(
    x = "Stage 1 negative pilot size m1",
    y = "Endpoint",
    fill = "Neyman - SRS",
    title = "RMSE Difference by Endpoint and m1",
    subtitle = "Negative values mean Adaptive Neyman has lower empirical RMSE"
  ) +
  heatmap_theme

abs_bias_difference_heatmap <- design_difference_df %>%
  ggplot(aes(x = m1_label, y = endpoint_label,
             fill = abs_bias_neyman_minus_srs)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = if_else(is.na(abs_bias_neyman_minus_srs), "",
                                sprintf("%.3f", abs_bias_neyman_minus_srs))),
            size = 2.6) +
  facet_wrap(~ metric_label, scales = "free", ncol = 3) +
  scale_x_discrete(labels = function(x) str_replace(x, "^m1=", "")) +
  scale_fill_gradient2(low = "#1B7F79", mid = "white", high = "#B85C38",
                       midpoint = 0, na.value = "grey90") +
  labs(
    x = "Stage 1 negative pilot size m1",
    y = "Endpoint",
    fill = "Neyman - SRS",
    title = "Absolute Bias Difference by Endpoint and m1",
    subtitle = "Negative values mean Adaptive Neyman has lower absolute bias"
  ) +
  heatmap_theme

coverage_difference_heatmap <- design_difference_df %>%
  filter(metric != "auroc") %>%
  ggplot(aes(x = m1_label, y = endpoint_label,
             fill = coverage_neyman_minus_srs)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = if_else(is.na(coverage_neyman_minus_srs), "",
                                sprintf("%.2f", coverage_neyman_minus_srs))),
            size = 2.6) +
  facet_wrap(~ metric_label, scales = "free", ncol = 3) +
  scale_x_discrete(labels = function(x) str_replace(x, "^m1=", "")) +
  scale_fill_gradient2(low = "#B85C38", mid = "white", high = "#1B7F79",
                       midpoint = 0, na.value = "grey90") +
  labs(
    x = "Stage 1 negative pilot size m1",
    y = "Endpoint",
    fill = "Neyman - SRS",
    title = "Coverage Difference by Endpoint and m1",
    subtitle = "Positive values mean Adaptive Neyman has higher empirical coverage"
  ) +
  heatmap_theme

ggsave(file.path(plot_output_dir, "endpoint_mean_estimates_mc_error_bars.png"),
       estimate_plot, width = 11, height = 7, dpi = 300)
ggsave(file.path(plot_output_dir, "endpoint_empirical_bias_mc_error_bars.png"),
       bias_plot, width = 11, height = 7, dpi = 300)
ggsave(file.path(plot_output_dir, "endpoint_empirical_sd.png"),
       sd_plot, width = 10, height = 6.5, dpi = 300)
ggsave(file.path(plot_output_dir, "endpoint_empirical_rmse.png"),
       rmse_plot, width = 10, height = 6.5, dpi = 300)
ggsave(file.path(plot_output_dir, "endpoint_average_ci_width.png"),
       ci_width_plot, width = 10, height = 6.5, dpi = 300)
ggsave(file.path(plot_output_dir, "endpoint_empirical_coverage.png"),
       coverage_plot, width = 10, height = 6.5, dpi = 300)
ggsave(file.path(plot_output_dir, "m1_heatmap_empirical_bias.png"),
       bias_heatmap, width = 12, height = 13, dpi = 300)
ggsave(file.path(plot_output_dir, "m1_heatmap_empirical_rmse.png"),
       rmse_heatmap, width = 12, height = 13, dpi = 300)
ggsave(file.path(plot_output_dir, "m1_heatmap_empirical_coverage.png"),
       coverage_heatmap, width = 12, height = 11, dpi = 300)
ggsave(file.path(plot_output_dir, "m1_heatmap_neyman_minus_srs_rmse.png"),
       rmse_difference_heatmap, width = 11, height = 7, dpi = 300)
ggsave(file.path(plot_output_dir, "m1_heatmap_neyman_minus_srs_abs_bias.png"),
       abs_bias_difference_heatmap, width = 11, height = 7, dpi = 300)
ggsave(file.path(plot_output_dir, "m1_heatmap_neyman_minus_srs_coverage.png"),
       coverage_difference_heatmap, width = 11, height = 7, dpi = 300)

cat("\nMonte Carlo endpoint summary\n")
print(monte_carlo_summary_long, n = Inf, width = Inf)

cat("\nSaved replicate-level and Monte Carlo outputs:\n")
cat("- ", normalizePath(replicate_long_path), "\n", sep = "")
cat("- ", normalizePath(summary_long_path), "\n", sep = "")
cat("- ", normalizePath(summary_wide_path), "\n", sep = "")

cat("\nSaved Monte Carlo endpoint summary plots:\n")
cat("- ", normalizePath(file.path(plot_output_dir,
                                  "endpoint_mean_estimates_mc_error_bars.png")),
    "\n", sep = "")
cat("- ", normalizePath(file.path(plot_output_dir,
                                  "endpoint_empirical_bias_mc_error_bars.png")),
    "\n", sep = "")
cat("- ", normalizePath(file.path(plot_output_dir,
                                  "endpoint_empirical_sd.png")),
    "\n", sep = "")
cat("- ", normalizePath(file.path(plot_output_dir,
                                  "endpoint_empirical_rmse.png")),
    "\n", sep = "")
cat("- ", normalizePath(file.path(plot_output_dir,
                                  "endpoint_average_ci_width.png")),
    "\n", sep = "")
cat("- ", normalizePath(file.path(plot_output_dir,
                                  "endpoint_empirical_coverage.png")),
    "\n", sep = "")
cat("- ", normalizePath(file.path(plot_output_dir,
                                  "m1_heatmap_empirical_bias.png")),
    "\n", sep = "")
cat("- ", normalizePath(file.path(plot_output_dir,
                                  "m1_heatmap_empirical_rmse.png")),
    "\n", sep = "")
cat("- ", normalizePath(file.path(plot_output_dir,
                                  "m1_heatmap_empirical_coverage.png")),
    "\n", sep = "")
cat("- ", normalizePath(file.path(plot_output_dir,
                                  "m1_heatmap_neyman_minus_srs_rmse.png")),
    "\n", sep = "")
cat("- ", normalizePath(file.path(plot_output_dir,
                                  "m1_heatmap_neyman_minus_srs_abs_bias.png")),
    "\n", sep = "")
cat("- ", normalizePath(file.path(plot_output_dir,
                                  "m1_heatmap_neyman_minus_srs_coverage.png")),
    "\n", sep = "")
