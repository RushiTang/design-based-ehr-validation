#!/usr/bin/env Rscript

# Final simulation using the unchanged manuscript data-generating mechanism.
# The primary adaptive analysis uses the fixed-weight lambda-pool composite
# estimator, metric-specific logit-Wald CIs, and an exact finite-population
# zero-FN safeguard on non-triggered SRS paths.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

expit <- function(z) 1 / (1 + exp(-z))
logit <- function(p) log(p / (1 - p))
clamp01 <- function(p, eps = 1e-10) pmin(pmax(p, eps), 1 - eps)
p_jeff <- function(x, n) (x + 0.5) / (n + 1)

tune_intercept <- function(x, slope, target_mean) {
  uniroot(function(a0) mean(expit(a0 + slope * x)) - target_mean,
          c(-30, 30))$root
}

wilson_ci <- function(x, n, alpha = 0.05) {
  if (n <= 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - alpha / 2)
  p <- x / n
  denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  half <- z / denominator * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(max(0, center - half), min(1, center + half))
}

round_alloc <- function(total, proportions) {
  if (total == 0L) return(integer(length(proportions)))
  raw <- total * proportions
  allocation <- floor(raw)
  remainder <- total - sum(allocation)
  if (remainder > 0L) {
    add <- order(raw - allocation, decreasing = TRUE)[seq_len(remainder)]
    allocation[add] <- allocation[add] + 1L
  }
  as.integer(allocation)
}

var_fn_total <- function(N_h, n_h, S_h) {
  valid <- N_h > 1 & n_h > 0
  sum(N_h[valid]^2 * (1 - n_h[valid] / N_h[valid]) *
        S_h[valid]^2 / n_h[valid])
}

binary_total_variance <- function(y, population_size) {
  sample_size <- length(y)
  if (population_size <= 1L || sample_size <= 1L) return(NA_real_)
  population_size^2 * (1 - sample_size / population_size) *
    var(y) / sample_size
}

sample_mnar <- function(candidates, take, D, kappa) {
  if (take <= 0L) return(integer(0))
  if (take >= length(candidates)) return(candidates)
  disease <- D[candidates]
  intercept <- uniroot(
    function(a) sum(expit(a + kappa * disease)) - take,
    c(-30, 30)
  )$root
  probabilities <- expit(intercept + kappa * disease)
  selected <- candidates[runif(length(candidates)) < probabilities]
  if (length(selected) > take) selected <- sample(selected, take)
  if (length(selected) < take) {
    left <- setdiff(candidates, selected)
    selected <- c(selected, sample(left, take - length(selected)))
  }
  selected
}

draw_stage2 <- function(remaining, stratum, allocation, D, mar_violation,
                        kappa) {
  selected <- integer(0)
  H <- length(allocation)
  for (h in seq_len(H)) {
    candidates <- remaining[stratum[remaining] == h]
    take <- min(allocation[h], length(candidates))
    if (take <= 0L) next
    selected_h <- if (mar_violation) {
      sample_mnar(candidates, take, D, kappa)
    } else {
      sample(candidates, take)
    }
    selected <- c(selected, selected_h)
  }
  target <- sum(allocation)
  if (length(selected) < target) {
    left <- setdiff(remaining, selected)
    selected <- c(selected, if (mar_violation) {
      sample_mnar(left, target - length(selected), D, kappa)
    } else {
      sample(left, target - length(selected))
    })
  }
  if (length(selected) > target) selected <- selected[seq_len(target)]
  selected
}

sequential_estimate <- function(D, stratum, stage1, stage2, N_h,
                                triggered, N_negative) {
  H <- length(N_h)
  m1 <- length(stage1)
  m2 <- length(stage2)

  if (m2 == 0L) {
    FN_hat <- N_negative * mean(D[stage1])
    variance <- binary_total_variance(D[stage1], N_negative)
    return(list(FN = FN_hat, variance = variance,
                n2_min = NA_integer_, zero_stage2_fn = NA))
  }

  if (!triggered) {
    remaining_size <- N_negative - m1
    FN_hat <- sum(D[stage1]) + remaining_size * mean(D[stage2])
    variance <- binary_total_variance(D[stage2], remaining_size)
    return(list(FN = FN_hat, variance = variance,
                n2_min = m2, zero_stage2_fn = sum(D[stage2]) == 0L))
  }

  n1_h <- tabulate(stratum[stage1], nbins = H)
  n2_h <- tabulate(stratum[stage2], nbins = H)
  c1_h <- tabulate(stratum[stage1][D[stage1] == 1L], nbins = H)
  c2_h <- tabulate(stratum[stage2][D[stage2] == 1L], nbins = H)
  remaining_h <- N_h - n1_h
  if (any(n2_h <= 0L & remaining_h > 0L)) {
    stop("Triggered design left a nonempty stratum without a Stage 2 sample.")
  }

  FN_hat <- sum(c1_h + remaining_h * c2_h / n2_h)
  variance_components <- vapply(seq_len(H), function(h) {
    if (remaining_h[h] <= 0L) return(0)
    binary_total_variance(D[stage2[stratum[stage2] == h]], remaining_h[h])
  }, numeric(1))
  variance <- if (anyNA(variance_components)) NA_real_ else sum(variance_components)
  list(FN = FN_hat, variance = variance, n2_min = min(n2_h),
       zero_stage2_fn = sum(c2_h) == 0L)
}

pooled_estimate <- function(D, stratum, stage1, stage2, N_h,
                            triggered, N_negative) {
  H <- length(N_h)
  final <- c(stage1, stage2)
  if (!triggered) {
    return(list(
      FN = N_negative * mean(D[final]),
      variance = binary_total_variance(D[final], N_negative)
    ))
  }

  final_h <- tabulate(stratum[final], nbins = H)
  cases_h <- tabulate(stratum[final][D[final] == 1L], nbins = H)
  if (any(final_h <= 0L)) stop("Pooled estimator has an empty final stratum.")
  components <- vapply(seq_len(H), function(h) {
    binary_total_variance(D[final[stratum[final] == h]], N_h[h])
  }, numeric(1))
  list(
    FN = sum(N_h * cases_h / final_h),
    variance = if (anyNA(components)) NA_real_ else sum(components)
  )
}

composite_estimate <- function(D, stage1, sequential_fit, N_negative,
                               total_sample_size) {
  m1 <- length(stage1)
  lambda_pool <- m1 * (N_negative - total_sample_size) /
    (total_sample_size * (N_negative - m1))
  stage1_FN <- N_negative * mean(D[stage1])
  stage1_variance <- binary_total_variance(D[stage1], N_negative)
  list(
    FN = lambda_pool * stage1_FN +
      (1 - lambda_pool) * sequential_fit$FN,
    variance = lambda_pool^2 * stage1_variance +
      (1 - lambda_pool)^2 * sequential_fit$variance,
    lambda_pool = lambda_pool,
    stage1_FN = stage1_FN,
    stage1_variance = stage1_variance
  )
}

hypergeom_zero_upper <- function(population_size, sample_size,
                                 alpha = 0.05) {
  if (sample_size <= 0L || sample_size > population_size) return(NA_integer_)
  threshold <- log(alpha / 2)
  log_prob_zero <- function(total_cases) {
    dhyper(0, total_cases, population_size - total_cases,
           sample_size, log = TRUE)
  }
  low <- 0L
  high <- as.integer(population_size)
  while (low < high) {
    middle <- as.integer(ceiling((low + high) / 2))
    if (log_prob_zero(middle) >= threshold) {
      low <- middle
    } else {
      high <- middle - 1L
    }
  }
  low
}

transform_fn_interval <- function(FN_lo, FN_hi, N_negative, TP, FP) {
  TN_lo <- N_negative - FN_hi
  TN_hi <- N_negative - FN_lo
  N_total <- N_negative + TP + FP
  rbind(
    sens = c(TP / (TP + FN_hi), TP / (TP + FN_lo)),
    spec = c(TN_lo / (TN_lo + FP), TN_hi / (TN_hi + FP)),
    npv = c(TN_lo / N_negative, TN_hi / N_negative),
    acc = c((TP + TN_lo) / N_total, (TP + TN_hi) / N_total)
  )
}

final_intervals <- function(FN_hat, variance, N_negative, TP, FP,
                            final_sample_size, zero_final_sample_fn,
                            allow_zero_fn_fallback, alpha = 0.05) {
  if (!is.finite(variance) || variance < 0) {
    intervals <- matrix(
      NA_real_, 5, 2,
      dimnames = list(c("sens", "spec", "ppv", "npv", "acc"),
                      c("lo", "hi"))
    )
    return(list(intervals = intervals,
                methods = rep("Unavailable", 5),
                fallback_used = rep(FALSE, 5)))
  }

  estimate <- metric_estimates(FN_hat, N_negative, TP, FP)
  derivatives <- c(
    sens = -TP / (TP + FN_hat)^2,
    spec = -FP / (N_negative - FN_hat + FP)^2,
    npv = -1 / N_negative,
    acc = -1 / (N_negative + TP + FP)
  )
  intervals <- matrix(
    NA_real_, 5, 2,
    dimnames = list(c("sens", "spec", "ppv", "npv", "acc"),
                    c("lo", "hi"))
  )
  for (metric in names(derivatives)) {
    variance_metric <- derivatives[[metric]]^2 * variance
    estimate_bounded <- clamp01(estimate[[metric]])
    se_logit <- sqrt(max(variance_metric, 0)) /
      (estimate_bounded * (1 - estimate_bounded))
    center <- logit(estimate_bounded)
    intervals[metric, ] <- expit(
      center + c(-1, 1) * qnorm(1 - alpha / 2) * se_logit
    )
  }

  ppv_interval <- wilson_ci(TP, TP + FP, alpha)
  intervals["ppv", ] <- ppv_interval

  use_fallback <- allow_zero_fn_fallback && zero_final_sample_fn
  methods <- c(
    sens = "Metric-specific logit-Wald",
    spec = "Metric-specific logit-Wald",
    ppv = "Wilson/binomial",
    npv = "Metric-specific logit-Wald",
    acc = "Metric-specific logit-Wald"
  )
  fallback_used <- setNames(rep(FALSE, 5), names(methods))
  if (use_fallback) {
    FN_hi <- hypergeom_zero_upper(N_negative, final_sample_size, alpha)
    fallback_metrics <- c("sens", "spec", "npv", "acc")
    intervals[fallback_metrics, ] <- transform_fn_interval(
      0, FN_hi, N_negative, TP, FP
    )[fallback_metrics, ]
    methods[fallback_metrics] <-
      "Finite-population hypergeometric zero-FN safeguard"
    fallback_used[fallback_metrics] <- TRUE
  }

  list(intervals = intervals, methods = methods,
       fallback_used = fallback_used)
}

metric_estimates <- function(FN_hat, N_negative, TP, FP) {
  TN_hat <- N_negative - FN_hat
  N_total <- N_negative + TP + FP
  c(
    sens = TP / (TP + FN_hat),
    spec = TN_hat / (TN_hat + FP),
    ppv = TP / (TP + FP),
    npv = TN_hat / N_negative,
    acc = (TP + TN_hat) / N_total
  )
}

naive_estimates <- function(D, T, verified) {
  TP <- sum(T[verified] == 1L & D[verified] == 1L)
  FP <- sum(T[verified] == 1L & D[verified] == 0L)
  FN <- sum(T[verified] == 0L & D[verified] == 1L)
  TN <- sum(T[verified] == 0L & D[verified] == 0L)
  c(
    sens = TP / (TP + FN),
    spec = TN / (TN + FP),
    ppv = TP / (TP + FP),
    npv = TN / (TN + FN),
    acc = (TP + TN) / length(verified)
  )
}

make_rows <- function(replicate, design, truth, estimate, intervals, naive,
                      FN_hat, variance, R_hat, triggered, imbalance, n2_min,
                      zero_stage2_fn, lambda_pool = NA_real_,
                      ci_methods = NULL, fallback_used = NULL) {
  metrics <- names(truth)
  if (is.null(ci_methods)) ci_methods <- rep(NA_character_, length(metrics))
  if (is.null(fallback_used)) fallback_used <- rep(FALSE, length(metrics))
  data.frame(
    rep = replicate,
    design = design,
    metric = metrics,
    truth = unname(truth),
    est = unname(estimate[metrics]),
    ci_lo = intervals[metrics, "lo"],
    ci_hi = intervals[metrics, "hi"],
    cover = as.integer(
      truth >= intervals[metrics, "lo"] & truth <= intervals[metrics, "hi"]
    ),
    width = intervals[metrics, "hi"] - intervals[metrics, "lo"],
    naive_est = unname(naive[metrics]),
    FN_hat = FN_hat,
    var_FN = variance,
    Rhat = R_hat,
    trigger = triggered,
    imbalance = if (triggered) imbalance else NA_real_,
    n2_min = n2_min,
    zero_stage2_fn = zero_stage2_fn,
    lambda_pool = lambda_pool,
    ci_method = unname(ci_methods[metrics]),
    zero_fn_fallback_used = unname(fallback_used[metrics]),
    stringsAsFactors = FALSE
  )
}

calibrate_algorithm <- function(scenario, constants) {
  X <- rnorm(constants$N)
  beta0 <- tune_intercept(X, constants$beta1, scenario$prev)
  D <- rbinom(constants$N, 1, expit(beta0 + constants$beta1 * X))
  objective <- function(parameters) {
    q <- expit(parameters[1] + parameters[2] * D +
                 scenario$gamma3 * D * X)
    (mean(q[D == 0]) - scenario$FPR_target)^2 +
      (mean(q[D == 1]) - scenario$Se_target)^2
  }
  g0 <- logit(clamp01(scenario$FPR_target))
  g1 <- logit(clamp01(scenario$Se_target)) - g0
  optim(c(g0, g1), objective, method = "Nelder-Mead")$par
}

run_scenario <- function(scenario, B, constants, output_dir) {
  set.seed(10000L + scenario$scenario_id)
  gamma <- calibrate_algorithm(scenario, constants)
  m2 <- scenario$m - scenario$m1
  kappa <- log(constants$mnar_or)
  output <- vector("list", B)
  allocation_topups <- 0L

  for (b in seq_len(B)) {
    X <- rnorm(constants$N)
    alpha0 <- tune_intercept(X, constants$alpha1, constants$risk_mean)
    risk <- expit(alpha0 + constants$alpha1 * X)
    beta0 <- tune_intercept(X, constants$beta1, scenario$prev)
    D <- rbinom(constants$N, 1, expit(beta0 + constants$beta1 * X))
    q <- expit(gamma[1] + gamma[2] * D + scenario$gamma3 * D * X)
    T <- rbinom(constants$N, 1, q)

    negative <- which(T == 0L)
    positive <- which(T == 1L)
    N_negative <- length(negative)
    cuts <- quantile(risk[negative], seq(0, 1, length.out = constants$H + 1),
                     type = 8)
    cuts[1] <- -Inf
    cuts[length(cuts)] <- Inf
    stratum <- rep(NA_integer_, constants$N)
    stratum[negative] <- as.integer(cut(
      risk[negative], cuts, labels = FALSE, include.lowest = TRUE
    ))
    N_h <- tabulate(stratum[negative], nbins = constants$H)

    TP <- sum(T == 1L & D == 1L)
    FP <- sum(T == 1L & D == 0L)
    FN <- sum(T == 0L & D == 1L)
    TN <- sum(T == 0L & D == 0L)
    truth <- c(
      sens = TP / (TP + FN),
      spec = TN / (TN + FP),
      ppv = TP / (TP + FP),
      npv = TN / (TN + FN),
      acc = (TP + TN) / constants$N
    )

    stage1 <- sample(negative, scenario$m1)
    remaining <- setdiff(negative, stage1)
    n1_h <- tabulate(stratum[stage1], nbins = constants$H)
    c1_h <- tabulate(stratum[stage1][D[stage1] == 1L], nbins = constants$H)
    p_overall <- p_jeff(sum(c1_h), scenario$m1)
    p_h <- ifelse(n1_h > 0L, p_jeff(c1_h, n1_h), p_overall)
    S_h <- sqrt(p_h * (1 - p_h) + 1e-12)
    allocation <- round_alloc(
      m2, N_h * S_h / sum(N_h * S_h)
    )

    if (m2 > 0L) {
      n_srs <- scenario$m * N_h / N_negative
      n_pilot_expected <- scenario$m1 * N_h / N_negative
      R_hat <- var_fn_total(N_h, n_srs, S_h) /
        var_fn_total(N_h, n_pilot_expected + allocation, S_h)
      triggered <- is.finite(R_hat) && R_hat > constants$R0
    } else {
      R_hat <- 1
      triggered <- FALSE
    }
    imbalance <- if (m2 > 0L) sd(allocation / m2) else NA_real_

    if (m2 > 0L) {
      stage2_srs <- if (scenario$mar_violation) {
        sample_mnar(remaining, m2, D, kappa)
      } else {
        sample(remaining, m2)
      }
      stage2_adaptive <- if (triggered) {
        sampled <- draw_stage2(
          remaining, stratum, allocation, D, scenario$mar_violation, kappa
        )
        realized <- tabulate(stratum[sampled], nbins = constants$H)
        allocation_topups <- allocation_topups + as.integer(any(realized != allocation))
        sampled
      } else if (scenario$mar_violation) {
        sample_mnar(remaining, m2, D, kappa)
      } else {
        sample(remaining, m2)
      }
    } else {
      stage2_srs <- stage2_adaptive <- integer(0)
    }

    srs_sequential_fit <- sequential_estimate(
      D, stratum, stage1, stage2_srs, N_h, FALSE, N_negative
    )
    adaptive_sequential_fit <- sequential_estimate(
      D, stratum, stage1, stage2_adaptive, N_h, triggered, N_negative
    )
    srs_pooled_fit <- pooled_estimate(
      D, stratum, stage1, stage2_srs, N_h, FALSE, N_negative
    )
    adaptive_pooled_fit <- pooled_estimate(
      D, stratum, stage1, stage2_adaptive, N_h, triggered, N_negative
    )
    adaptive_composite_fit <- composite_estimate(
      D, stage1, adaptive_sequential_fit, N_negative, scenario$m
    )
    if (!triggered && abs(adaptive_composite_fit$FN -
                          adaptive_pooled_fit$FN) > 1e-8) {
      stop("Lambda-pool composite did not reduce to pooled SRS estimation.")
    }

    srs_naive <- naive_estimates(D, T, c(positive, stage1, stage2_srs))
    adaptive_naive <- naive_estimates(
      D, T, c(positive, stage1, stage2_adaptive)
    )

    row_for_fit <- function(label, fit, naive, design_triggered,
                            design_imbalance, n2_min, zero_stage2_fn,
                            zero_final_sample_fn,
                            lambda_pool = NA_real_) {
      estimate <- metric_estimates(fit$FN, N_negative, TP, FP)
      ci_fit <- final_intervals(
        fit$FN, fit$variance, N_negative, TP, FP,
        final_sample_size = scenario$m,
        zero_final_sample_fn = zero_final_sample_fn,
        allow_zero_fn_fallback = !design_triggered && !scenario$mar_violation
      )
      make_rows(
        b, label, truth, estimate, ci_fit$intervals, naive,
        fit$FN, fit$variance, R_hat, design_triggered, design_imbalance,
        n2_min, zero_stage2_fn, lambda_pool,
        ci_methods = ci_fit$methods,
        fallback_used = ci_fit$fallback_used
      )
    }

    zero_srs_final_fn <- sum(D[c(stage1, stage2_srs)]) == 0L
    zero_adaptive_final_fn <- sum(D[c(stage1, stage2_adaptive)]) == 0L

    output[[b]] <- bind_rows(
      row_for_fit(
        "SRS_Pooled", srs_pooled_fit, srs_naive, FALSE, NA_real_,
        srs_sequential_fit$n2_min, srs_sequential_fit$zero_stage2_fn,
        zero_srs_final_fn
      ),
      row_for_fit(
        "Adaptive_Pooled", adaptive_pooled_fit, adaptive_naive,
        triggered, imbalance, adaptive_sequential_fit$n2_min,
        adaptive_sequential_fit$zero_stage2_fn, zero_adaptive_final_fn
      ),
      row_for_fit(
        "Adaptive_Sequential", adaptive_sequential_fit, adaptive_naive,
        triggered, imbalance, adaptive_sequential_fit$n2_min,
        adaptive_sequential_fit$zero_stage2_fn, zero_adaptive_final_fn
      ),
      row_for_fit(
        "Adaptive_Composite", adaptive_composite_fit, adaptive_naive,
        triggered, imbalance, adaptive_sequential_fit$n2_min,
        adaptive_sequential_fit$zero_stage2_fn, zero_adaptive_final_fn,
        adaptive_composite_fit$lambda_pool
      )
    )
  }

  raw <- bind_rows(output) %>%
    mutate(
      scenario_id = scenario$scenario_id,
      prev = scenario$prev,
      algo_case = scenario$algo_case,
      assoc_case = scenario$assoc_case,
      mar_violation = scenario$mar_violation,
      m = scenario$m,
      m1 = scenario$m1
    )
  saveRDS(raw, file.path(
    output_dir, "raw", sprintf("raw_scenario_%03d.rds", scenario$scenario_id)
  ))

  performance <- raw %>%
    group_by(scenario_id, prev, algo_case, assoc_case, mar_violation,
             m, m1, metric, design) %>%
    summarise(
      n_replicates = n(),
      bias = mean(est - truth, na.rm = TRUE),
      rmse = sqrt(mean((est - truth)^2, na.rm = TRUE)),
      coverage = mean(cover, na.rm = TRUE),
      avg_width = mean(width, na.rm = TRUE),
      zero_width_rate = mean(width < 1e-12, na.rm = TRUE),
      zero_fn_fallback_rate = mean(zero_fn_fallback_used, na.rm = TRUE),
      .groups = "drop"
    )

  qa <- raw %>%
    filter(metric == "sens") %>%
    group_by(design) %>%
    summarise(
      trigger_rate = mean(trigger),
      missing_variance_rate = mean(!is.finite(var_FN)),
      zero_stage2_fn_rate = mean(zero_stage2_fn, na.rm = TRUE),
      zero_width_rate = mean(width < 1e-12),
      zero_fn_fallback_rate = mean(zero_fn_fallback_used),
      minimum_n2 = suppressWarnings(min(n2_min, na.rm = TRUE)),
      mean_lambda_pool = mean(lambda_pool, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      scenario_id = scenario$scenario_id,
      allocation_topup_rate = allocation_topups / B
    )
  list(performance = performance, qa = qa)
}

make_scenarios <- function() {
  expand_grid(
    prev = c(0.05, 0.10, 0.20),
    algo_case = c("HighSe_ModSp", "ModSe_HighSp"),
    assoc_case = c("None", "Weak", "Strong"),
    mar_violation = c(FALSE, TRUE),
    m = c(150L, 275L, 400L),
    m1 = c(50L, 100L, 150L)
  ) %>%
    mutate(
      scenario_id = row_number(),
      Se_target = ifelse(algo_case == "HighSe_ModSp", 0.90, 0.80),
      Sp_target = ifelse(algo_case == "HighSe_ModSp", 0.85, 0.95),
      FPR_target = 1 - Sp_target,
      gamma3 = case_when(
        assoc_case == "None" ~ 0,
        assoc_case == "Weak" ~ -0.6,
        TRUE ~ -1.2
      )
    )
}

summarize_across_facets <- function(performance) {
  performance %>%
    group_by(prev, mar_violation, m, m1, metric, design) %>%
    summarise(
      bias = mean(bias),
      rmse = sqrt(mean(rmse^2)),
      coverage = mean(coverage),
      avg_width = mean(avg_width),
      min_coverage = min(coverage),
      max_abs_bias = max(abs(bias)),
      .groups = "drop"
    )
}

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1L) as.integer(args[1]) else 2000L
workers <- if (length(args) >= 2L) as.integer(args[2]) else {
  min(8L, max(1L, parallel::detectCores() - 1L))
}
output_dir <- if (length(args) >= 3L) args[3] else file.path(
  dirname(normalizePath(sub("^--file=", "",
                            grep("^--file=", commandArgs(), value = TRUE)[1]))),
  "outputs"
)
dir.create(file.path(output_dir, "raw"), recursive = TRUE, showWarnings = FALSE)

constants <- list(
  N = 3000L,
  H = 5L,
  alpha1 = 1,
  risk_mean = 0.30,
  beta1 = 0.8,
  R0 = 1.10,
  mnar_or = 1.20
)
scenarios <- make_scenarios()
write.csv(scenarios, file.path(output_dir, "scenario_grid.csv"), row.names = FALSE)

cat("Running", nrow(scenarios), "scenarios x", B, "replicates with",
    workers, "workers\n")
runner <- function(index) {
  run_scenario(scenarios[index, ], B, constants, output_dir)
}
results <- if (workers > 1L && .Platform$OS.type == "unix") {
  parallel::mclapply(seq_len(nrow(scenarios)), runner, mc.cores = workers,
                     mc.preschedule = FALSE, mc.set.seed = FALSE)
} else {
  lapply(seq_len(nrow(scenarios)), runner)
}

performance <- bind_rows(lapply(results, `[[`, "performance"))
qa <- bind_rows(lapply(results, `[[`, "qa"))
collapsed <- summarize_across_facets(performance)
saveRDS(performance, file.path(output_dir, "performance_by_scenario.rds"))
write.csv(performance, file.path(output_dir, "performance_by_scenario.csv"),
          row.names = FALSE)
write.csv(collapsed, file.path(output_dir, "performance_collapsed.csv"),
          row.names = FALSE)
write.csv(qa, file.path(output_dir, "implementation_qa.csv"), row.names = FALSE)

cat("Completed. Outputs:", output_dir, "\n")
