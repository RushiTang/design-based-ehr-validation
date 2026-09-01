############################################################
# Create a manuscript-style real-data repeated-sampling results table.
#
# Input:
#   <output_root>/multi_endpoint_monte_carlo_summary_long.csv
#
# Output:
#   <output_root>/real_world_simulation_results_table.csv
#
# The table reports mean estimates across real-data sampling replicates.
# The CI column uses the empirical 2.5th and 97.5th percentiles of the
# replicate estimates.
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
})

output_root <- Sys.getenv(
  "MULTI_ENDPOINT_OUTPUT_ROOT",
  unset = "multi_endpoint_validation_outputs_m1"
)
summary_path <- file.path(output_root, "multi_endpoint_monte_carlo_summary_long.csv")

if (!file.exists(summary_path)) {
  fallback_path <- file.path(
    "multi_endpoint_validation_outputs",
    "multi_endpoint_monte_carlo_summary_long.csv"
  )
  if (file.exists(fallback_path)) {
    summary_path <- fallback_path
    output_root <- dirname(summary_path)
  } else {
    stop(
      "Could not find multi_endpoint_monte_carlo_summary_long.csv. ",
      "Set MULTI_ENDPOINT_OUTPUT_ROOT to the folder containing that file."
    )
  }
}

m1_to_report <- Sys.getenv("RESULTS_TABLE_M1", unset = "100")

format_num <- function(x, digits = 3) {
  if_else(is.na(x), NA_character_, sprintf(paste0("%.", digits, "f"), x))
}

format_ci <- function(lower, upper, digits = 3) {
  if_else(
    is.na(lower) | is.na(upper),
    NA_character_,
    paste0("(", format_num(lower, digits), ", ", format_num(upper, digits), ")")
  )
}

format_percent <- function(x, digits = 3) {
  if_else(is.na(x), NA_character_, paste0(format_num(100 * x, digits), "%"))
}

summary_long <- readr::read_csv(summary_path, show_col_types = FALSE) %>%
  janitor::clean_names()

if (!"oracle_recurrence_prevalence" %in% names(summary_long)) {
  summary_long$oracle_recurrence_prevalence <- NA_real_
}

if ("m1" %in% names(summary_long)) {
  m1_value <- as.integer(m1_to_report)
  if (is.na(m1_value)) {
    stop("RESULTS_TABLE_M1 must be an integer.")
  }
  summary_long <- summary_long %>% filter(m1 == m1_value)
}

required_cols <- c(
  "endpoint_label", "metric", "design", "mean_rhat", "planned_m",
  "mean_estimate", "mc_ci_lower", "mc_ci_upper", "oracle_estimate"
)
missing_cols <- setdiff(required_cols, names(summary_long))
if (length(missing_cols) > 0) {
  stop("Summary file is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

metric_order <- c("ppv", "npv", "sensitivity", "specificity", "accuracy")
metric_labels <- c(
  ppv = "ppv",
  npv = "npv",
  sensitivity = "sens",
  specificity = "spec",
  accuracy = "acc"
)
design_labels <- c(
  neyman = "neyman",
  srs = "SRS"
)

results_table <- summary_long %>%
  filter(metric %in% metric_order, design %in% names(design_labels)) %>%
  mutate(
    endpoint = endpoint_label,
    metric = factor(metric, levels = metric_order),
    design = factor(recode(design, !!!design_labels),
                    levels = unname(design_labels)),
    prevalence = format_percent(oracle_recurrence_prevalence),
    rhat = format_num(mean_rhat),
    m = as.integer(planned_m),
    estimate = format_num(mean_estimate),
    CI = format_ci(mc_ci_lower, mc_ci_upper),
    oracle_estimate = format_num(oracle_estimate)
  ) %>%
  arrange(endpoint_year, metric, design) %>%
  mutate(
    metric = metric_labels[as.character(metric)],
    design = as.character(design)
  ) %>%
  select(endpoint, prevalence, metric, design, rhat, m, estimate, CI,
         oracle_estimate)

output_path <- file.path(output_root, "real_world_simulation_results_table.csv")
readr::write_csv(results_table, output_path)

metric_display_labels <- c(
  ppv = "PPV",
  npv = "NPV",
  sens = "Sensitivity",
  spec = "Specificity",
  acc = "Accuracy"
)

wide_results_table <- results_table %>%
  mutate(
    metric = metric_display_labels[metric],
    design_column = recode(design,
                           neyman = "Adaptive composite", SRS = "SRS"),
    estimate_ci = paste(estimate, CI)
  ) %>%
  select(endpoint, prevalence, rhat, m, metric, design_column, estimate_ci) %>%
  pivot_wider(names_from = design_column, values_from = estimate_ci) %>%
  left_join(
    results_table %>%
      mutate(metric = metric_display_labels[metric]) %>%
      distinct(endpoint, metric, Oracle = oracle_estimate),
    by = c("endpoint", "metric")
  ) %>%
  mutate(
    endpoint_sort = as.integer(str_remove(endpoint, "yr")),
    metric = factor(metric, levels = unname(metric_display_labels))
  ) %>%
  arrange(endpoint_sort, metric) %>%
  mutate(metric = as.character(metric)) %>%
  select(endpoint, prevalence, rhat, m, metric, `Adaptive composite`, SRS, Oracle)

wide_output_path <- file.path(
  output_root,
  "real_world_simulation_results_table_wide.csv"
)
readr::write_csv(wide_results_table, wide_output_path)

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\", "\\\\textbackslash{}", x, fixed = TRUE)
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("$", "\\\\$", x, fixed = TRUE)
  x <- gsub("#", "\\\\#", x, fixed = TRUE)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x <- gsub("{", "\\\\{", x, fixed = TRUE)
  x <- gsub("}", "\\\\}", x, fixed = TRUE)
  x
}

last_endpoint <- tail(unique(wide_results_table$endpoint), 1)

latex_rows <- wide_results_table %>%
  group_by(endpoint) %>%
  mutate(
    endpoint_cell = if_else(row_number() == 1L, endpoint, ""),
    prevalence_cell = if_else(row_number() == 1L, prevalence, ""),
    rhat_cell = if_else(row_number() == 1L, rhat, ""),
    m_cell = if_else(row_number() == 1L, as.character(m), ""),
    row_text = paste(
      latex_escape(endpoint_cell),
      latex_escape(prevalence_cell),
      latex_escape(rhat_cell),
      latex_escape(m_cell),
      latex_escape(metric),
      latex_escape(`Adaptive composite`),
      latex_escape(SRS),
      latex_escape(Oracle),
      sep = " & "
    ),
    row_text = paste0(row_text, " \\\\"),
    row_text = if_else(
      row_number() == n() & endpoint != .env$last_endpoint,
      paste0(row_text, "\n\\addlinespace"),
      row_text
    )
  ) %>%
  ungroup() %>%
  pull(row_text)

latex_table <- c(
  "% Requires \\usepackage{booktabs}",
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Real-data repeated-sampling validation performance.}",
  "\\label{tab:real_data_validation}",
  "\\begin{tabular}{lccclccc}",
  "\\toprule",
  "Endpoint & Prevalence & $\\hat{R}$ & $m$ & Metric & Adaptive composite & SRS & Oracle \\\\",
  "\\midrule",
  latex_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

latex_output_path <- file.path(
  output_root,
  "real_world_simulation_results_table_wide.tex"
)
writeLines(latex_table, latex_output_path)

backslash <- intToUtf8(92)
latex_file_lines <- readLines(latex_output_path, warn = FALSE)
latex_file_lines <- gsub(
  paste0(backslash, backslash, "%"),
  paste0(backslash, "%"),
  latex_file_lines,
  fixed = TRUE
)
writeLines(latex_file_lines, latex_output_path)

cat("\nReal-world simulation results table\n")
print(results_table, n = Inf, width = Inf)

cat("\nSaved table to: ", normalizePath(output_path, mustWork = FALSE), "\n", sep = "")
cat("Saved wide table to: ", normalizePath(wide_output_path, mustWork = FALSE), "\n", sep = "")
cat("Saved LaTeX table to: ", normalizePath(latex_output_path, mustWork = FALSE), "\n", sep = "")
if ("m1" %in% names(readr::read_csv(summary_path, n_max = 0, show_col_types = FALSE))) {
  cat("Reported m1 = ", m1_to_report, "\n", sep = "")
}
