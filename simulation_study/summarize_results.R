#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
output_dir <- if (length(args)) args[1] else file.path(script_dir, "outputs")
performance <- read.csv(file.path(output_dir, "performance_by_scenario.csv"))

primary <- performance %>%
  filter(!mar_violation,
         design %in% c("SRS_Pooled", "Adaptive_Composite"))

primary_metric_summary <- primary %>%
  group_by(design, metric) %>%
  summarise(
    mean_bias = mean(bias),
    mean_rmse = mean(rmse),
    mean_coverage = mean(coverage),
    q05_coverage = quantile(coverage, 0.05),
    minimum_coverage = min(coverage),
    proportion_below_093 = mean(coverage < 0.93),
    mean_width = mean(avg_width),
    mean_fallback_rate = mean(zero_fn_fallback_rate),
    .groups = "drop"
  )

primary_prevalence_summary <- primary %>%
  group_by(design, metric, prev) %>%
  summarise(
    mean_coverage = mean(coverage),
    q05_coverage = quantile(coverage, 0.05),
    minimum_coverage = min(coverage),
    proportion_below_093 = mean(coverage < 0.93),
    mean_width = mean(avg_width),
    mean_fallback_rate = mean(zero_fn_fallback_rate),
    .groups = "drop"
  )

primary_method_ratios <- primary %>%
  filter(metric != "ppv") %>%
  select(scenario_id, metric, design, rmse, avg_width, coverage) %>%
  pivot_wider(names_from = design,
              values_from = c(rmse, avg_width, coverage)) %>%
  mutate(
    rmse_ratio = rmse_Adaptive_Composite / rmse_SRS_Pooled,
    width_ratio = avg_width_Adaptive_Composite / avg_width_SRS_Pooled
  ) %>%
  group_by(metric) %>%
  summarise(
    mean_rmse_ratio = mean(rmse_ratio),
    median_rmse_ratio = median(rmse_ratio),
    proportion_lower_rmse = mean(rmse_ratio < 1),
    mean_width_ratio = mean(width_ratio),
    median_width_ratio = median(width_ratio),
    proportion_narrower = mean(width_ratio < 1),
    mean_coverage_difference =
      mean(coverage_Adaptive_Composite - coverage_SRS_Pooled),
    .groups = "drop"
  )

summary_by_prevalence <- performance %>%
  filter(metric == "sens") %>%
  group_by(mar_violation, design, prev) %>%
  summarise(
    mean_bias = mean(bias),
    max_abs_bias = max(abs(bias)),
    mean_rmse = mean(rmse),
    mean_coverage = mean(coverage),
    median_coverage = median(coverage),
    minimum_coverage = min(coverage),
    mean_zero_width_rate = mean(zero_width_rate),
    .groups = "drop"
  )

method_comparison <- performance %>%
  filter(!mar_violation, metric != "ppv") %>%
  select(scenario_id, prev, algo_case, assoc_case, m, m1, metric,
         design, bias, rmse, coverage, zero_width_rate) %>%
  pivot_wider(names_from = design,
              values_from = c(bias, rmse, coverage, zero_width_rate))

raw_files <- list.files(file.path(output_dir, "raw"),
                        pattern = "raw_scenario_[0-9]+[.]rds$", full.names = TRUE)
raw <- bind_rows(lapply(raw_files, readRDS))

raw_output_qa <- tibble(
  raw_files = length(raw_files),
  rows_observed = nrow(raw),
  unique_keys = nrow(distinct(raw, scenario_id, rep, design, metric)),
  finite_ordered_bounded_intervals = sum(
    is.finite(raw$ci_lo) & is.finite(raw$ci_hi) &
      raw$ci_lo <= raw$ci_hi & raw$ci_lo >= 0 & raw$ci_hi <= 1
  ),
  missing_variances = sum(!is.finite(raw$var_FN)),
  invalid_fallback_rows = sum(
    raw$zero_fn_fallback_used &
      (raw$trigger | raw$mar_violation | raw$metric == "ppv")
  )
)

adaptive_branch_metrics <- raw %>%
  filter(!mar_violation, design == "Adaptive_Composite") %>%
  mutate(branch = ifelse(trigger, "Triggered", "Not triggered")) %>%
  group_by(metric, branch, prev) %>%
  summarise(
    n = n(),
    bias = mean(est - truth),
    rmse = sqrt(mean((est - truth)^2)),
    coverage = mean(cover),
    mean_width = mean(width),
    zero_width_rate = mean(width < 1e-12),
    fallback_rate = mean(zero_fn_fallback_used),
    .groups = "drop"
  )
truth_wide <- raw %>%
  select(scenario_id, rep, design, metric, truth) %>%
  pivot_wider(names_from = metric, values_from = truth)

fn_results <- raw %>%
  filter(metric == "sens") %>%
  select(-metric, -truth) %>%
  left_join(truth_wide, by = c("scenario_id", "rep", "design")) %>%
  mutate(
    TP = 3000 / (1 + (1 - ppv) / ppv + (1 - sens) / sens +
                   ((1 - sens) / sens) * npv / (1 - npv)),
    FN_true = TP * (1 - sens) / sens,
    FN_error = FN_hat - FN_true,
    branch = ifelse(trigger, "Triggered", "Not triggered")
  )

fn_validation <- fn_results %>%
  group_by(mar_violation, design) %>%
  summarise(
    n = n(),
    FN_bias = mean(FN_error),
    empirical_variance = var(FN_error),
    mean_variance_estimate = mean(var_FN),
    variance_ratio = mean_variance_estimate / empirical_variance,
    coverage = mean(cover),
    zero_width_rate = mean(width < 1e-12),
    .groups = "drop"
  )

branch_validation <- fn_results %>%
  group_by(mar_violation, design, prev, branch) %>%
  summarise(
    n = n(),
    FN_bias = mean(FN_error),
    FN_rmse = sqrt(mean(FN_error^2)),
    coverage = mean(cover),
    zero_width_rate = mean(width < 1e-12),
    coverage_given_nonzero_width = mean(cover[width >= 1e-12]),
    mean_width = mean(width),
    .groups = "drop"
  )

worst_coverage <- performance %>%
  filter(!mar_violation, metric == "sens") %>%
  arrange(coverage) %>%
  select(scenario_id, prev, algo_case, assoc_case, m, m1, design,
         bias, rmse, coverage, zero_width_rate) %>%
  slice_head(n = 100)

write.csv(summary_by_prevalence,
          file.path(output_dir, "summary_by_prevalence.csv"), row.names = FALSE)
write.csv(primary_metric_summary,
          file.path(output_dir, "primary_metric_summary.csv"), row.names = FALSE)
write.csv(primary_prevalence_summary,
          file.path(output_dir, "primary_prevalence_summary.csv"),
          row.names = FALSE)
write.csv(primary_method_ratios,
          file.path(output_dir, "primary_method_ratios.csv"), row.names = FALSE)
write.csv(raw_output_qa,
          file.path(output_dir, "raw_output_qa.csv"), row.names = FALSE)
write.csv(adaptive_branch_metrics,
          file.path(output_dir, "adaptive_branch_metrics.csv"), row.names = FALSE)
write.csv(method_comparison,
          file.path(output_dir, "method_comparison_common_samples.csv"),
          row.names = FALSE)
write.csv(fn_validation,
          file.path(output_dir, "fn_total_validation.csv"), row.names = FALSE)
write.csv(branch_validation,
          file.path(output_dir, "branch_validation.csv"), row.names = FALSE)
write.csv(worst_coverage,
          file.path(output_dir, "worst_coverage_scenarios.csv"), row.names = FALSE)

cat("Wrote composite analysis summaries to", normalizePath(output_dir), "\n")
