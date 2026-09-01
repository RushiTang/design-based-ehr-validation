#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(forcats)
  library(scales)
  library(svglite)
  library(ragg)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
script_dir <- dirname(script_path)
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1L) args[1] else file.path(script_dir, "outputs")
output_dir <- if (length(args) >= 2L) args[2] else file.path(script_dir, "figures")
input_dir <- normalizePath(input_dir, mustWork = TRUE)
source_dir <- file.path(output_dir, "source_data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

palette <- c(
  proposed = "#168B7A",
  srs = "#4F5962",
  naive = "#A8A8A8",
  ink = "#252525",
  mid = "#747474",
  light = "#D9D9D9",
  gain = "#2A8C78",
  neutral = "#FAFAF8",
  loss = "#C45C58",
  trigger_low = "#F6F7F5",
  trigger_mid = "#A9D4C7",
  trigger_high = "#176C62"
)

metric_labels <- c(
  sens = "Sensitivity",
  spec = "Specificity",
  ppv = "PPV",
  npv = "NPV",
  acc = "Accuracy"
)
metric_levels <- c("Sensitivity", "Specificity", "PPV", "NPV", "Accuracy")
design_labels <- c(
  Adaptive_Composite = "Adaptive composite",
  SRS_Pooled = "Design-based SRS"
)
association_labels <- c(
  None = "No risk signal",
  Weak = "Weak risk signal",
  Strong = "Strong risk signal"
)

theme_nature <- function(base_size = 7, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = palette[["ink"]]),
      axis.ticks = element_line(linewidth = 0.3, colour = palette[["ink"]]),
      axis.ticks.length = unit(1.6, "pt"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4, colour = palette[["ink"]]),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.height = unit(7, "pt"),
      legend.key.width = unit(10, "pt"),
      strip.background = element_rect(fill = "white", colour = palette[["light"]], linewidth = 0.3),
      strip.text = element_text(size = base_size - 0.1, face = "bold", margin = margin(2, 2, 2, 2)),
      plot.title = element_text(size = base_size + 0.6, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = palette[["mid"]], hjust = 0),
      plot.tag = element_text(size = 9, face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(5, 6, 5, 5),
      legend.position = "bottom"
    )
}

theme_set(theme_nature())

save_pub <- function(plot, filename, width_mm, height_mm, dpi = 600) {
  base <- file.path(output_dir, filename)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in, bg = "white")
  print(plot)
  dev.off()

  grDevices::pdf(
    paste0(base, ".pdf"), width = width_in, height = height_in,
    family = "sans", bg = "white", useDingbats = FALSE
  )
  print(plot)
  dev.off()

  ragg::agg_tiff(
    paste0(base, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, background = "white", compression = "lzw"
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    paste0(base, ".png"), width = width_in, height = height_in,
    units = "in", res = 240, background = "white"
  )
  print(plot)
  dev.off()
}

quantile_summary <- function(data, value) {
  data %>%
    summarise(
      q05 = quantile({{ value }}, 0.05, na.rm = TRUE),
      q25 = quantile({{ value }}, 0.25, na.rm = TRUE),
      median = median({{ value }}, na.rm = TRUE),
      q75 = quantile({{ value }}, 0.75, na.rm = TRUE),
      q95 = quantile({{ value }}, 0.95, na.rm = TRUE),
      .groups = "drop"
    )
}

scenario_grid <- read_csv(file.path(input_dir, "scenario_grid.csv"), show_col_types = FALSE)
mar_ids <- scenario_grid %>%
  filter(!mar_violation) %>%
  pull(scenario_id)

raw_files <- file.path(input_dir, "raw", sprintf("raw_scenario_%03d.rds", mar_ids))
if (!all(file.exists(raw_files))) {
  stop("One or more B=2000 MAR raw scenario files are missing.")
}

message("Aggregating MAR raw outputs for bias, variance calibration, and trigger rate...")
raw_summaries <- lapply(raw_files, function(path) {
  dat <- readRDS(path)
  meta <- dat[1, c("scenario_id", "prev", "algo_case", "assoc_case", "m", "m1")]

  adaptive <- dat %>%
    filter(design == "Adaptive_Composite")

  bias <- adaptive %>%
    group_by(metric) %>%
    summarise(
      proposed_bias = mean(est - truth, na.rm = TRUE),
      proposed_rmse = sqrt(mean((est - truth)^2, na.rm = TRUE)),
      naive_bias = mean(naive_est - truth, na.rm = TRUE),
      naive_rmse = sqrt(mean((naive_est - truth)^2, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      scenario_id = meta$scenario_id,
      prev = meta$prev,
      algo_case = meta$algo_case,
      assoc_case = meta$assoc_case,
      m = meta$m,
      m1 = meta$m1
    )

  truth_wide <- adaptive %>%
    select(rep, metric, truth) %>%
    pivot_wider(names_from = metric, values_from = truth) %>%
    mutate(
      TP = 3000 / (1 + (1 - ppv) / ppv + (1 - sens) / sens +
                     ((1 - sens) / sens) * npv / (1 - npv)),
      FN_true = TP * (1 - sens) / sens
    ) %>%
    select(rep, FN_true)

  variance <- dat %>%
    filter(
      design %in% c("Adaptive_Composite", "SRS_Pooled"),
      metric == "sens"
    ) %>%
    left_join(truth_wide, by = "rep") %>%
    mutate(FN_error = FN_hat - FN_true) %>%
    group_by(design) %>%
    summarise(
      n_replicates = n(),
      FN_bias = mean(FN_error, na.rm = TRUE),
      FN_rmse = sqrt(mean(FN_error^2, na.rm = TRUE)),
      empirical_variance = var(FN_error, na.rm = TRUE),
      mean_variance_estimate = mean(var_FN, na.rm = TRUE),
      variance_ratio = mean_variance_estimate / empirical_variance,
      bias_mcse = sqrt(empirical_variance / n_replicates),
      standardized_bias = FN_bias / bias_mcse,
      .groups = "drop"
    ) %>%
    mutate(
      scenario_id = meta$scenario_id,
      prev = meta$prev,
      algo_case = meta$algo_case,
      assoc_case = meta$assoc_case,
      m = meta$m,
      m1 = meta$m1
    )

  trigger <- adaptive %>%
    filter(metric == "sens") %>%
    summarise(
      trigger_rate = mean(trigger, na.rm = TRUE),
      fallback_rate = mean(zero_fn_fallback_used, na.rm = TRUE),
      mean_rhat = mean(Rhat, na.rm = TRUE)
    ) %>%
    mutate(
      scenario_id = meta$scenario_id,
      prev = meta$prev,
      algo_case = meta$algo_case,
      assoc_case = meta$assoc_case,
      m = meta$m,
      m1 = meta$m1
    )

  list(bias = bias, variance = variance, trigger = trigger)
})

bias_wide <- bind_rows(lapply(raw_summaries, `[[`, "bias"))
variance_scenario <- bind_rows(lapply(raw_summaries, `[[`, "variance"))
trigger_scenario <- bind_rows(lapply(raw_summaries, `[[`, "trigger"))
rm(raw_summaries)
invisible(gc())

bias_source <- bias_wide %>%
  select(scenario_id, prev, algo_case, assoc_case, m, m1, metric,
         proposed_bias, proposed_rmse, naive_bias, naive_rmse) %>%
  pivot_longer(
    cols = c(proposed_bias, proposed_rmse, naive_bias, naive_rmse),
    names_to = c("estimator_key", ".value"),
    names_pattern = "(proposed|naive)_(bias|rmse)"
  ) %>%
  mutate(
    estimator = recode(
      estimator_key,
      proposed = "Proposed design-based",
      naive = "Naive complete-case"
    ),
    estimator = factor(estimator, levels = c("Proposed design-based", "Naive complete-case")),
    metric = factor(metric_labels[metric], levels = metric_levels),
    prevalence = factor(percent(prev, accuracy = 1), levels = c("5%", "10%", "20%")),
    bias_fraction_rmse = if_else(rmse > 0, bias / rmse, 0)
  )

write_csv(bias_source, file.path(source_dir, "figure1_bias_source.csv"))
write_csv(variance_scenario, file.path(source_dir, "figure2_fn_variance_source.csv"))
write_csv(trigger_scenario, file.path(source_dir, "figure3_trigger_source.csv"))

performance <- read_csv(file.path(input_dir, "performance_by_scenario.csv"), show_col_types = FALSE) %>%
  filter(
    !mar_violation,
    design %in% c("Adaptive_Composite", "SRS_Pooled")
  ) %>%
  mutate(
    design_label = factor(design_labels[design], levels = unname(design_labels)),
    metric_label = factor(metric_labels[metric], levels = metric_levels),
    prevalence = factor(percent(prev, accuracy = 1), levels = c("5%", "10%", "20%")),
    empirical_sd = sqrt(pmax(rmse^2 - bias^2, 0)),
    bias_mcse = empirical_sd / sqrt(n_replicates),
    standardized_bias = if_else(bias_mcse > 0, bias / bias_mcse, 0)
  )

write_csv(performance, file.path(source_dir, "figure2_scenario_performance_source.csv"))

main_table <- performance %>%
  group_by(design = design_label, metric = metric_label) %>%
  summarise(
    mean_bias = mean(bias, na.rm = TRUE),
    mean_empirical_sd = mean(empirical_sd, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_coverage = mean(coverage, na.rm = TRUE),
    q05_coverage = quantile(coverage, 0.05, na.rm = TRUE),
    minimum_coverage = min(coverage, na.rm = TRUE),
    mean_ci_width = mean(avg_width, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(main_table, file.path(source_dir, "main_simulation_summary_table.csv"))

# Figure 1: the proposed estimator versus a naive complete-case analysis.
p1 <- ggplot(bias_source, aes(x = estimator, y = bias, fill = estimator)) +
  geom_hline(yintercept = 0, colour = palette[["ink"]], linewidth = 0.3) +
  geom_boxplot(
    width = 0.62, linewidth = 0.3, outlier.shape = NA,
    colour = palette[["ink"]]
  ) +
  geom_point(
    aes(colour = estimator), position = position_jitter(width = 0.12, height = 0),
    size = 0.65, alpha = 0.34, stroke = 0
  ) +
  facet_grid(metric ~ prevalence, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "Proposed design-based" = palette[["proposed"]],
      "Naive complete-case" = palette[["naive"]]
    )
  ) +
  scale_colour_manual(
    values = c(
      "Proposed design-based" = palette[["proposed"]],
      "Naive complete-case" = palette[["mid"]]
    )
  ) +
  scale_x_discrete(labels = c(
    "Proposed design-based" = "Proposed",
    "Naive complete-case" = "Naive"
  )) +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  labs(
    title = "Design-based estimation prevents complete-case bias",
    subtitle = "Scenario-level Monte Carlo bias under the design-conforming mechanism",
    x = NULL,
    y = "Signed bias",
    fill = NULL,
    colour = NULL
  ) +
  guides(colour = "none", fill = guide_legend(nrow = 1)) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.text.y = element_text(angle = 0),
    panel.spacing = unit(4, "pt")
  )

save_pub(p1, "figure1_proposed_vs_naive_bias", 183, 190)

# Figure 2: FN-total bias, interval coverage, and variance calibration.
coverage_summary <- performance %>%
  group_by(prevalence, design_label, metric_label) %>%
  quantile_summary(coverage)

variance_summary <- variance_scenario %>%
  filter(is.finite(variance_ratio), variance_ratio > 0) %>%
  mutate(
    design_label = factor(design_labels[design], levels = unname(design_labels)),
    prevalence = factor(percent(prev, accuracy = 1), levels = c("5%", "10%", "20%"))
  ) %>%
  group_by(prevalence, design_label) %>%
  quantile_summary(variance_ratio)

fn_bias_summary <- variance_scenario %>%
  filter(is.finite(standardized_bias)) %>%
  mutate(
    design_label = factor(design_labels[design], levels = unname(design_labels)),
    prevalence = factor(percent(prev, accuracy = 1), levels = c("5%", "10%", "20%"))
  ) %>%
  group_by(prevalence, design_label) %>%
  quantile_summary(standardized_bias)

write_csv(fn_bias_summary, file.path(source_dir, "figure2_fn_bias_summary.csv"))
write_csv(coverage_summary, file.path(source_dir, "figure2_coverage_summary.csv"))
write_csv(variance_summary, file.path(source_dir, "figure2_variance_ratio_summary.csv"))

design_colours <- c(
  "Adaptive composite" = palette[["proposed"]],
  "Design-based SRS" = palette[["srs"]]
)

p2a <- ggplot(
  fn_bias_summary,
  aes(x = prevalence, y = median, colour = design_label, group = design_label)
) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.96, ymax = 1.96,
           fill = palette[["light"]], alpha = 0.22) +
  geom_hline(yintercept = 0, colour = palette[["ink"]], linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = q05, ymax = q95),
    width = 0, linewidth = 0.45,
    position = position_dodge(width = 0.18)
  ) +
  geom_point(size = 1.7, position = position_dodge(width = 0.18)) +
  scale_colour_manual(values = design_colours) +
  labs(
    title = "False-negative total bias",
    subtitle = "Median and 5th-95th percentiles across scenarios; shaded region is +/-1.96 Monte Carlo SE",
    x = "Outcome prevalence",
    y = "FN-total bias / Monte Carlo SE",
    colour = NULL
  )

p2b <- ggplot(
  performance,
  aes(x = metric_label, y = coverage, fill = design_label, colour = design_label)
) +
  geom_hline(yintercept = 0.95, colour = palette[["ink"]], linewidth = 0.35, linetype = "dashed") +
  geom_boxplot(
    width = 0.62, linewidth = 0.35, outlier.shape = NA,
    position = position_dodge(width = 0.72), alpha = 0.82
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.13, dodge.width = 0.72),
    size = 0.65, alpha = 0.32, stroke = 0
  ) +
  facet_wrap(~prevalence, nrow = 1) +
  scale_colour_manual(values = design_colours) +
  scale_fill_manual(values = design_colours) +
  coord_cartesian(ylim = c(0.87, 1.005)) +
  scale_y_continuous(breaks = c(0.88, 0.92, 0.95, 1.00), labels = label_number(accuracy = 0.01)) +
  labs(
    title = "Confidence-interval coverage across simulation scenarios",
    subtitle = "Each point represents one scenario; boxes show the median and interquartile range",
    x = NULL,
    y = "Empirical coverage",
    colour = NULL,
    fill = NULL
  ) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  guides(colour = "none", fill = guide_legend(nrow = 1))

p2c <- ggplot(
  variance_summary,
  aes(x = prevalence, y = median, colour = design_label, group = design_label)
) +
  geom_hline(yintercept = 1, colour = palette[["ink"]], linewidth = 0.35, linetype = "dashed") +
  geom_errorbar(aes(ymin = q05, ymax = q95), width = 0.04, linewidth = 0.45) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = design_colours) +
  labs(
    title = "False-negative total variance calibration",
    subtitle = "Median and 5th-95th percentiles across scenarios; dashed line denotes exact calibration",
    x = "Outcome prevalence",
    y = "Mean estimated variance / empirical variance",
    colour = NULL
  ) +
  guides(colour = "none")

figure2 <- p2b

save_pub(figure2, "figure2_estimator_and_ci_calibration", 183, 96)

figure_s1 <- p2a / p2c +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a")

save_pub(figure_s1, "figureS1_fn_estimator_audit", 183, 145)

# Figure 3: efficiency ratios and adaptive trigger behavior.
ratio_source <- performance %>%
  filter(metric != "ppv") %>%
  select(scenario_id, prev, algo_case, assoc_case, m, m1, metric,
         design, rmse, avg_width, coverage) %>%
  pivot_wider(
    names_from = design,
    values_from = c(rmse, avg_width, coverage)
  ) %>%
  mutate(
    rmse_ratio = rmse_Adaptive_Composite / rmse_SRS_Pooled,
    width_ratio = avg_width_Adaptive_Composite / avg_width_SRS_Pooled,
    coverage_difference = coverage_Adaptive_Composite - coverage_SRS_Pooled
  ) %>%
  group_by(prev, assoc_case, m, m1, metric) %>%
  summarise(
    rmse_ratio = mean(rmse_ratio, na.rm = TRUE),
    width_ratio = mean(width_ratio, na.rm = TRUE),
    coverage_difference = mean(coverage_difference, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    metric_label = factor(metric_labels[metric], levels = metric_levels),
    prevalence = factor(percent(prev, accuracy = 1), levels = c("5%", "10%", "20%")),
    association = factor(association_labels[assoc_case], levels = unname(association_labels)),
    budget = factor(
      paste0(m, "/", m1),
      levels = c("150/50", "150/100", "150/150", "275/50", "275/100", "275/150",
                 "400/50", "400/100", "400/150")
    ),
    no_stage2 = m == m1
  )

trigger_source <- trigger_scenario %>%
  group_by(prev, assoc_case, m, m1) %>%
  summarise(
    trigger_rate = mean(trigger_rate, na.rm = TRUE),
    fallback_rate = mean(fallback_rate, na.rm = TRUE),
    mean_rhat = mean(mean_rhat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    prevalence = factor(percent(prev, accuracy = 1), levels = c("5%", "10%", "20%")),
    association = factor(association_labels[assoc_case], levels = unname(association_labels)),
    budget = factor(
      paste0(m, "/", m1),
      levels = c("150/50", "150/100", "150/150", "275/50", "275/100", "275/150",
                 "400/50", "400/100", "400/150")
    ),
    no_stage2 = m == m1
  )

write_csv(ratio_source, file.path(source_dir, "figure3_efficiency_ratio_source.csv"))
write_csv(trigger_source, file.path(source_dir, "figure3_trigger_grid_source.csv"))

ratio_values <- c(ratio_source$rmse_ratio, ratio_source$width_ratio)
ratio_limit <- min(0.25, max(0.08, quantile(abs(ratio_values - 1), 0.98, na.rm = TRUE)))
ratio_limits <- c(1 - ratio_limit, 1 + ratio_limit)

ratio_heatmap <- function(data, value, title, subtitle) {
  ggplot(data, aes(x = budget, y = metric_label)) +
    geom_tile(
      data = filter(data, no_stage2),
      fill = palette[["light"]], colour = "white", linewidth = 0.3
    ) +
    geom_tile(
      data = filter(data, !no_stage2),
      aes(fill = {{ value }}), colour = "white", linewidth = 0.3
    ) +
    geom_text(
      data = filter(data, no_stage2), label = "-",
      size = 2.1, colour = palette[["mid"]]
    ) +
    facet_grid(prevalence ~ association) +
    scale_fill_gradient2(
      low = palette[["gain"]], mid = palette[["neutral"]], high = palette[["loss"]],
      midpoint = 1, limits = ratio_limits, oob = squish,
      breaks = c(ratio_limits[[1]], 1, ratio_limits[[2]]),
      labels = label_number(accuracy = 0.01)
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Total/pilot algorithm-negative reviews (m/m1)",
      y = NULL,
      fill = "Ratio"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5.8),
      axis.text.y = element_text(size = 6),
      panel.spacing = unit(2.8, "pt"),
      legend.position = "right",
      legend.key.height = unit(23, "pt")
    )
}

p3a <- ratio_heatmap(
  ratio_source, rmse_ratio,
  "RMSE ratio: adaptive composite / SRS",
  "Values below 1 favour adaptive sampling; algorithm-performance profiles are averaged"
)

p3b <- ratio_heatmap(
  ratio_source, width_ratio,
  "Mean CI-width ratio: adaptive composite / SRS",
  "Values below 1 indicate narrower intervals; algorithm-performance profiles are averaged"
)

p3c <- ggplot(trigger_source, aes(x = budget, y = "Trigger", fill = trigger_rate)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(
    aes(label = if_else(no_stage2, "-", percent(trigger_rate, accuracy = 1))),
    size = 1.9,
    colour = if_else(trigger_source$trigger_rate >= 0.45, "white", palette[["ink"]])
  ) +
  facet_grid(prevalence ~ association) +
  scale_fill_gradientn(
    colours = c(palette[["trigger_low"]], palette[["trigger_mid"]], palette[["trigger_high"]]),
    limits = c(0, 1), labels = label_percent(accuracy = 1), oob = squish
  ) +
  labs(
    title = "Adaptive trigger rate",
    subtitle = "The trigger decision is shared by all operating-characteristic estimators",
    x = "Total/pilot algorithm-negative reviews (m/m1)",
    y = NULL,
    fill = "Trigger rate"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5.8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(2.8, "pt"),
    legend.position = "right",
    legend.key.height = unit(20, "pt")
  )

figure3 <- p3a / p3b / p3c +
  plot_layout(heights = c(1.15, 1.15, 0.72)) +
  plot_annotation(tag_levels = "a")

save_pub(figure3, "figure3_adaptive_efficiency_and_trigger", 183, 255)

# Combined main-text figure: bias, efficiency ratios, and trigger behavior.
main_bias_source <- performance %>%
  mutate(
    relative_bias = if_else(empirical_sd > 0, bias / empirical_sd, 0),
    metric_label = factor(metric_label, levels = metric_levels)
  )
write_csv(main_bias_source, file.path(source_dir, "figure_main_relative_bias_source.csv"))

p_main_a <- ggplot(
  main_bias_source,
  aes(x = metric_label, y = relative_bias, fill = design_label, colour = design_label)
) +
  geom_hline(yintercept = 0, colour = palette[["ink"]], linewidth = 0.35) +
  geom_boxplot(
    width = 0.62, linewidth = 0.35, outlier.shape = NA,
    position = position_dodge(width = 0.72), alpha = 0.82
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.13, dodge.width = 0.72),
    size = 0.6, alpha = 0.28, stroke = 0
  ) +
  facet_wrap(~prevalence, nrow = 1) +
  scale_colour_manual(values = design_colours) +
  scale_fill_manual(values = design_colours) +
  coord_cartesian(ylim = c(-0.10, 0.205)) +
  scale_y_continuous(
    breaks = c(-0.10, 0, 0.10, 0.20),
    labels = label_number(accuracy = 0.01)
  ) +
  labs(
    title = "Scenario-specific bias relative to empirical SD",
    subtitle = "Values near zero indicate that systematic error is small relative to sampling variability",
    x = NULL,
    y = "Bias / empirical SD",
    fill = NULL,
    colour = NULL
  ) +
  guides(colour = "none", fill = guide_legend(nrow = 1)) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

efficiency_long <- ratio_source %>%
  select(prev, prevalence, assoc_case, association, m, m1, budget,
         no_stage2, metric, metric_label, rmse_ratio, width_ratio) %>%
  pivot_longer(
    cols = c(rmse_ratio, width_ratio),
    names_to = "outcome",
    values_to = "ratio"
  ) %>%
  mutate(
    outcome = recode(outcome, rmse_ratio = "RMSE", width_ratio = "CI width"),
    facet_row = factor(
      paste(outcome, prevalence, sep = " | "),
      levels = c(
        "RMSE | 5%", "RMSE | 10%", "RMSE | 20%",
        "CI width | 5%", "CI width | 10%", "CI width | 20%"
      )
    )
  )
write_csv(efficiency_long, file.path(source_dir, "figure_main_efficiency_source.csv"))

p_main_b <- ggplot(efficiency_long, aes(x = budget, y = metric_label)) +
  geom_tile(
    data = filter(efficiency_long, no_stage2),
    fill = palette[["light"]], colour = "white", linewidth = 0.3
  ) +
  geom_tile(
    data = filter(efficiency_long, !no_stage2),
    aes(fill = ratio), colour = "white", linewidth = 0.3
  ) +
  geom_text(
    data = filter(efficiency_long, no_stage2), label = "-",
    size = 1.9, colour = palette[["mid"]]
  ) +
  facet_grid(facet_row ~ association) +
  scale_fill_gradient2(
    low = palette[["gain"]], mid = palette[["neutral"]], high = palette[["loss"]],
    midpoint = 1, limits = ratio_limits, oob = squish,
    breaks = c(ratio_limits[[1]], 1, ratio_limits[[2]]),
    labels = label_number(accuracy = 0.01)
  ) +
  labs(
    title = "Efficiency ratios: adaptive composite / SRS",
    subtitle = "Values below 1 favour adaptive sampling; algorithm-performance profiles are averaged",
    x = "Total/pilot algorithm-negative reviews (m/m1)",
    y = NULL,
    fill = "Ratio"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5.7),
    axis.text.y = element_text(size = 5.9),
    strip.text.y = element_text(angle = 0, size = 5.8),
    panel.spacing = unit(2.2, "pt"),
    legend.position = "right",
    legend.key.height = unit(24, "pt")
  )

p_main_c <- p3c +
  labs(
    title = "Adaptive trigger rate",
    subtitle = "The trigger decision is shared by all operating-characteristic estimators"
  )

figure_main <- p_main_a / p_main_b / p_main_c +
  plot_layout(heights = c(0.78, 1.72, 0.72)) +
  plot_annotation(tag_levels = "a")

save_pub(figure_main, "figure_main_bias_efficiency_trigger", 183, 245)

message("Finished. Outputs written to: ", output_dir)
