# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

###############################################################################
## Common-n subsampling sensitivity analysis for zeta diversity
## Related to Fig. S2B-D
##
## Main logic:
## - All dilution treatments are covered at common n = 8, 10, 12, 14.
## - n = 14 is the primary narrative setting because it is the largest
##   replicate number shared by all treatments.
## - n = 8, 10 and 12 are used as trend-consistency checks.
###############################################################################

library(zetadiv)
library(tidyverse)
library(ggplot2)
library(ggsci)

## Parallel backend for common-n subsampling.
## This script is tuned for a desktop with 8 cores / 16 threads and 32 GB RAM.
## On Windows, use multisession rather than multicore.
has_future = requireNamespace("future", quietly = TRUE) &&
  requireNamespace("future.apply", quietly = TRUE)
if (has_future) {
  library(future)
  library(future.apply)
}

has_openxlsx = requireNamespace("openxlsx", quietly = TRUE)

set.seed(123)

treat_levels = c("D0", "D1", "D3", "D5", "D6", "D7")
dilution_rank = c(D0 = 0, D1 = 1, D3 = 3, D5 = 5, D6 = 6, D7 = 7)

n_grid = c(8, 10, 12, 14)

## Outer subsampling iterations.
## Increase to 1000 if computation time permits.
n_iter = 100

## Monte Carlo samples inside Zeta.decline.mc.
## This follows the original Fig. S2B-D script.
zeta_sam = 1000
confint.level = 0.95

## Parallel settings.
## i7-11700T: 8 physical cores / 16 threads; 32 GB RAM.
## workers = 6 leaves resources for Windows/RStudio and reduces memory pressure.
use_parallel = TRUE
n_workers = 6

## If your computer becomes unresponsive, change n_workers to 4.
## If testing/debugging, use n_iter = 3-10 and zeta_sam = 100 first.

outdir = "zeta_commonN_outputs_parallel"
if (!dir.exists(outdir)) dir.create(outdir)

###############################################################################
## 1. Read and binarize OTU table

otu = read.table("OTUtable.txt",
                 header = TRUE, sep = "\t",
                 row.names = 1,
                 check.names = FALSE,
                 stringsAsFactors = FALSE)

data.spec = as.matrix(t(otu))
data.spec[data.spec != 0] = 1

sample_groups = lapply(treat_levels, function(g) {
  grep(paste0("^", g), rownames(data.spec), value = TRUE)
})
names(sample_groups) = treat_levels

sample_n = tibble(
  Group = names(sample_groups),
  Original_replicates = sapply(sample_groups, length)
)

print(sample_n)

if (min(sample_n$Original_replicates) < max(n_grid)) {
  stop("At least one treatment has fewer replicates than max(n_grid).")
}

###############################################################################
## 2. Helper functions

qfun = function(x, p) quantile(x, p, na.rm = TRUE)

extract_aic = function(zobj) {
  out = tibble(Model = c("Exponential", "PowerLaw"), AIC = NA_real_)

  if (!is.null(zobj$aic)) {
    aic_val = NULL
    if (is.data.frame(zobj$aic) && "AIC" %in% colnames(zobj$aic)) {
      aic_val = zobj$aic$AIC
    } else if (is.matrix(zobj$aic) && "AIC" %in% colnames(zobj$aic)) {
      aic_val = zobj$aic[, "AIC"]
    } else if (is.vector(zobj$aic)) {
      aic_val = as.numeric(zobj$aic)
    }
    if (!is.null(aic_val) && length(aic_val) >= 2) {
      out$AIC = as.numeric(aic_val[1:2])
    }
  }
  out
}

calc_curve_metrics = function(zeta_values, orders) {
  df = tibble(Order = orders, Zeta = zeta_values) %>%
    filter(!is.na(Zeta), Zeta > 0)

  if (nrow(df) < 3) {
    return(tibble(
      Exp_slope = NA_real_,
      Power_slope = NA_real_,
      Zeta_order1 = NA_real_,
      Zeta_order2 = NA_real_,
      Zeta_order3 = NA_real_,
      Zeta_last = NA_real_,
      Zeta_half_order = NA_real_,
      Last_over_first = NA_real_
    ))
  }

  exp_fit = try(lm(log(Zeta) ~ Order, data = df), silent = TRUE)
  power_fit = try(lm(log(Zeta) ~ log(Order), data = df), silent = TRUE)

  exp_slope = if (inherits(exp_fit, "try-error")) NA_real_ else coef(exp_fit)[["Order"]]
  power_slope = if (inherits(power_fit, "try-error")) NA_real_ else coef(power_fit)[["log(Order)"]]

  half_order = floor(max(orders) / 2)

  tibble(
    Exp_slope = exp_slope,
    Power_slope = power_slope,
    Zeta_order1 = zeta_values[orders == 1][1],
    Zeta_order2 = zeta_values[orders == 2][1],
    Zeta_order3 = zeta_values[orders == 3][1],
    Zeta_last = zeta_values[orders == max(orders)][1],
    Zeta_half_order = zeta_values[orders == half_order][1],
    Last_over_first = zeta_values[orders == max(orders)][1] / zeta_values[orders == 1][1]
  )
}

run_zeta_one = function(data_bin, orders, zeta_sam = 1000) {
  ## Remove all-zero OTUs within sampled subset
  data_bin = data_bin[, colSums(data_bin) > 0, drop = FALSE]
  data_df = as.data.frame(data_bin)

  Zeta.decline.mc(data_df,
                  xy = NULL,
                  orders = orders,
                  sam = zeta_sam,
                  confint.level = confint.level,
                  rescale = TRUE,
                  normalize = TRUE,
                  empty.row = "empty",
                  plot = FALSE)
}

summarise_one_zeta = function(zobj, group, common_n, iteration) {
  orders = seq_along(zobj$zeta.val)

  ## In zetadiv::Zeta.decline.mc(), zeta.val has length equal to the
  ## number of orders, whereas ratio usually has length (number of orders - 1),
  ## because it represents decline between consecutive zeta orders.
  ## To store zeta and ratio in one long table, we pad ratio to the same length.
  ## Here Ratio at Order = k is interpreted as the decline ratio from order
  ## k - 1 to order k; therefore Order = 1 has no ratio and is set to NA.
  ratio_raw = as.numeric(zobj$ratio)
  if (length(ratio_raw) == length(orders) - 1) {
    ratio_pad = c(NA_real_, ratio_raw)
  } else if (length(ratio_raw) == length(orders)) {
    ratio_pad = ratio_raw
  } else {
    warning("Unexpected ratio length: zeta length = ", length(orders),
            ", ratio length = ", length(ratio_raw),
            ". Ratio will be padded/truncated to match zeta orders.")
    ratio_pad = rep(NA_real_, length(orders))
    keep_n = min(length(ratio_raw), length(orders))
    ratio_pad[seq_len(keep_n)] = ratio_raw[seq_len(keep_n)]
  }

  curve_df = tibble(
    Common_n = common_n,
    Iteration = iteration,
    Group = group,
    Order = orders,
    Zeta = as.numeric(zobj$zeta.val),
    Ratio = ratio_pad
  )

  aic_df = extract_aic(zobj) %>%
    mutate(Common_n = common_n,
           Iteration = iteration,
           Group = group)

  metrics_df = calc_curve_metrics(as.numeric(zobj$zeta.val), orders) %>%
    mutate(Common_n = common_n,
           Iteration = iteration,
           Group = group)

  delta = NA_real_
  if (all(!is.na(aic_df$AIC))) {
    delta = aic_df$AIC[aic_df$Model == "Exponential"] -
      aic_df$AIC[aic_df$Model == "PowerLaw"]
  }

  metrics_df = metrics_df %>%
    mutate(DeltaAIC_ExpMinusPower = delta,
           PowerLaw_better = ifelse(is.na(delta), NA, delta > 0))

  list(curve = curve_df, aic = aic_df, metrics = metrics_df)
}

###############################################################################
## 3. Run common-n subsampling zeta analysis
## Parallelized version for i7-11700T / 32 GB RAM.
## Each task = one common-n x iteration x treatment zeta calculation.

run_one_task = function(task_row) {
  common_n = as.integer(task_row$Common_n)
  iter = as.integer(task_row$Iteration)
  g = as.character(task_row$Group)

  ## Deterministic seed for reproducible subsampling independent of worker order.
  seed_i = 1000000 + common_n * 10000 + iter * 100 + match(g, treat_levels)
  set.seed(seed_i)

  sample_ids = sample(sample_groups[[g]], size = common_n, replace = FALSE)
  data_g = data.spec[sample_ids, , drop = FALSE]

  zobj = try(
    run_zeta_one(data_g, orders = seq_len(common_n), zeta_sam = zeta_sam),
    silent = TRUE
  )

  if (inherits(zobj, "try-error")) {
    err = tibble(
      Common_n = common_n,
      Iteration = iter,
      Group = g,
      Error = as.character(zobj)
    )
    return(list(curve = NULL, aic = NULL, metrics = NULL, error = err))
  }

  one = summarise_one_zeta(zobj, g, common_n, iter)
  one$error = NULL
  return(one)
}

if (use_parallel && has_future) {
  workers_use = min(n_workers, max(1, future::availableCores() - 2))
  message("Parallel mode enabled: multisession workers = ", workers_use)
  future::plan(future::multisession, workers = workers_use)
  options(future.globals.maxSize = 8 * 1024^3)  # 8 GB max exported globals
} else {
  message("Parallel mode disabled or future/future.apply not installed; running sequentially.")
}

all_curve = list()
all_aic = list()
all_metrics = list()
all_errors = list()

for (common_n in n_grid) {
  message("===== Common n = ", common_n, " =====")

  task_grid_n = tibble(
    Common_n = common_n,
    Iteration = rep(seq_len(n_iter), each = length(treat_levels)),
    Group = rep(treat_levels, times = n_iter)
  )

  ## Split into row-wise task objects. This avoids nested loops and makes checkpointing easy.
  tasks = split(task_grid_n, seq_len(nrow(task_grid_n)))

  if (use_parallel && has_future) {
    res_n = future.apply::future_lapply(tasks, run_one_task, future.seed = TRUE)
  } else {
    res_n = lapply(seq_along(tasks), function(i) {
      if (i %% 50 == 0 || i == 1) {
        message("Running task ", i, " / ", length(tasks), " for common n = ", common_n)
      }
      run_one_task(tasks[[i]])
    })
  }

  curve_n = bind_rows(lapply(res_n, `[[`, "curve"))
  aic_n = bind_rows(lapply(res_n, `[[`, "aic"))
  metrics_n = bind_rows(lapply(res_n, `[[`, "metrics"))
  error_n = bind_rows(lapply(res_n, `[[`, "error"))

  ## Save checkpoint after each common-n. This is important for long runs.
  saveRDS(
    list(curve = curve_n, aic = aic_n, metrics = metrics_n, errors = error_n),
    file.path(outdir, paste0("checkpoint_common_n", common_n, ".rds"))
  )

  write.csv(curve_n, file.path(outdir, paste0("curve_common_n", common_n, ".csv")), row.names = FALSE)
  write.csv(aic_n, file.path(outdir, paste0("aic_common_n", common_n, ".csv")), row.names = FALSE)
  write.csv(metrics_n, file.path(outdir, paste0("metrics_common_n", common_n, ".csv")), row.names = FALSE)
  if (nrow(error_n) > 0) {
    write.csv(error_n, file.path(outdir, paste0("errors_common_n", common_n, ".csv")), row.names = FALSE)
  }

  all_curve[[as.character(common_n)]] = curve_n
  all_aic[[as.character(common_n)]] = aic_n
  all_metrics[[as.character(common_n)]] = metrics_n
  all_errors[[as.character(common_n)]] = error_n
}

curve_long = bind_rows(all_curve) %>%
  mutate(Group = factor(Group, levels = treat_levels))

aic_long = bind_rows(all_aic) %>%
  mutate(Group = factor(Group, levels = treat_levels))

metrics_long = bind_rows(all_metrics) %>%
  mutate(Group = factor(Group, levels = treat_levels),
         Dilution_rank = dilution_rank[as.character(Group)])

error_long = bind_rows(all_errors)
if (nrow(error_long) > 0) {
  write.csv(error_long, file.path(outdir, "Zeta_commonN_failed_tasks.csv"), row.names = FALSE)
}

if (use_parallel && has_future) {
  future::plan(future::sequential)
}

###############################################################################
## 4. Summaries

curve_summary = curve_long %>%
  group_by(Common_n, Group, Order) %>%
  summarise(
    Mean_Zeta = mean(Zeta, na.rm = TRUE),
    Median_Zeta = median(Zeta, na.rm = TRUE),
    Q2.5_Zeta = qfun(Zeta, 0.025),
    Q97.5_Zeta = qfun(Zeta, 0.975),
    Mean_Ratio = mean(Ratio, na.rm = TRUE),
    Median_Ratio = median(Ratio, na.rm = TRUE),
    Q2.5_Ratio = qfun(Ratio, 0.025),
    Q97.5_Ratio = qfun(Ratio, 0.975),
    .groups = "drop"
  )

aic_summary = aic_long %>%
  group_by(Common_n, Group, Model) %>%
  summarise(
    Mean_AIC = mean(AIC, na.rm = TRUE),
    Median_AIC = median(AIC, na.rm = TRUE),
    Q2.5_AIC = qfun(AIC, 0.025),
    Q97.5_AIC = qfun(AIC, 0.975),
    .groups = "drop"
  )

metrics_summary = metrics_long %>%
  pivot_longer(cols = c(Exp_slope, Power_slope, Zeta_order1, Zeta_order2,
                        Zeta_order3, Zeta_last, Zeta_half_order,
                        Last_over_first, DeltaAIC_ExpMinusPower),
               names_to = "Metric",
               values_to = "Value") %>%
  group_by(Common_n, Group, Metric) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Q2.5 = qfun(Value, 0.025),
    Q97.5 = qfun(Value, 0.975),
    .groups = "drop"
  )

###############################################################################
## 5. Framework-level trends by iteration

framework_by_iter = metrics_long %>%
  select(Common_n, Iteration, Group, Dilution_rank,
         Exp_slope, Power_slope, Zeta_last, Zeta_half_order,
         Last_over_first, DeltaAIC_ExpMinusPower, PowerLaw_better) %>%
  group_by(Common_n, Iteration) %>%
  summarise(
    Mean_Zeta_last_D3D7 = mean(Zeta_last[Group %in% c("D3", "D5", "D6", "D7")], na.rm = TRUE),
    Mean_Zeta_last_D5D7 = mean(Zeta_last[Group %in% c("D5", "D6", "D7")], na.rm = TRUE),
    D0_Zeta_last = Zeta_last[Group == "D0"][1],
    D1_Zeta_last = Zeta_last[Group == "D1"][1],
    D7_Zeta_last = Zeta_last[Group == "D7"][1],

    Mean_Exp_slope_D3D7 = mean(Exp_slope[Group %in% c("D3", "D5", "D6", "D7")], na.rm = TRUE),
    D0_Exp_slope = Exp_slope[Group == "D0"][1],
    D7_Exp_slope = Exp_slope[Group == "D7"][1],

    Mean_Power_slope_D3D7 = mean(Power_slope[Group %in% c("D3", "D5", "D6", "D7")], na.rm = TRUE),
    D0_Power_slope = Power_slope[Group == "D0"][1],
    D7_Power_slope = Power_slope[Group == "D7"][1],

    Spearman_ZetaLast_vs_DilutionRank =
      suppressWarnings(cor(Dilution_rank, Zeta_last, method = "spearman", use = "complete.obs")),
    Spearman_ExpSlope_vs_DilutionRank =
      suppressWarnings(cor(Dilution_rank, Exp_slope, method = "spearman", use = "complete.obs")),
    Spearman_PowerSlope_vs_DilutionRank =
      suppressWarnings(cor(Dilution_rank, Power_slope, method = "spearman", use = "complete.obs")),

    Prop_Groups_PowerLawBetter = mean(PowerLaw_better, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ZetaLast_D3D7_minus_D0 = Mean_Zeta_last_D3D7 - D0_Zeta_last,
    ZetaLast_D5D7_minus_D0 = Mean_Zeta_last_D5D7 - D0_Zeta_last,
    ZetaLast_D7_minus_D0 = D7_Zeta_last - D0_Zeta_last,

    ExpSlope_D3D7_minus_D0 = Mean_Exp_slope_D3D7 - D0_Exp_slope,
    ExpSlope_D7_minus_D0 = D7_Exp_slope - D0_Exp_slope,

    PowerSlope_D3D7_minus_D0 = Mean_Power_slope_D3D7 - D0_Power_slope,
    PowerSlope_D7_minus_D0 = D7_Power_slope - D0_Power_slope,

    Stronger_lower_high_order_zeta = ZetaLast_D3D7_minus_D0 < 0,
    D7_lower_high_order_zeta = ZetaLast_D7_minus_D0 < 0,
    Stronger_more_negative_exp_slope = ExpSlope_D3D7_minus_D0 < 0,
    D7_more_negative_exp_slope = ExpSlope_D7_minus_D0 < 0,
    ZetaLast_decreases_with_dilution = Spearman_ZetaLast_vs_DilutionRank < 0
  )

framework_summary = framework_by_iter %>%
  pivot_longer(cols = c(ZetaLast_D3D7_minus_D0, ZetaLast_D5D7_minus_D0,
                        ZetaLast_D7_minus_D0, ExpSlope_D3D7_minus_D0,
                        ExpSlope_D7_minus_D0, PowerSlope_D3D7_minus_D0,
                        PowerSlope_D7_minus_D0, Spearman_ZetaLast_vs_DilutionRank,
                        Spearman_ExpSlope_vs_DilutionRank,
                        Spearman_PowerSlope_vs_DilutionRank,
                        Prop_Groups_PowerLawBetter),
               names_to = "Framework_metric",
               values_to = "Value") %>%
  group_by(Common_n, Framework_metric) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Q2.5 = qfun(Value, 0.025),
    Q97.5 = qfun(Value, 0.975),
    .groups = "drop"
  )

framework_prob = framework_by_iter %>%
  group_by(Common_n) %>%
  summarise(
    Prob_Stronger_lower_high_order_zeta = mean(Stronger_lower_high_order_zeta, na.rm = TRUE),
    Prob_D7_lower_high_order_zeta = mean(D7_lower_high_order_zeta, na.rm = TRUE),
    Prob_Stronger_more_negative_exp_slope = mean(Stronger_more_negative_exp_slope, na.rm = TRUE),
    Prob_D7_more_negative_exp_slope = mean(D7_more_negative_exp_slope, na.rm = TRUE),
    Prob_ZetaLast_decreases_with_dilution = mean(ZetaLast_decreases_with_dilution, na.rm = TRUE),
    Mean_Prop_Groups_PowerLawBetter = mean(Prop_Groups_PowerLawBetter, na.rm = TRUE),
    .groups = "drop"
  )

###############################################################################
## 6. Trend consistency relative to n = 14

ref_n = 14

reference_metrics = metrics_summary %>%
  filter(Common_n == ref_n) %>%
  select(Group, Metric, Ref_Median = Median)

trend_consistency_vs_n14 = metrics_summary %>%
  left_join(reference_metrics, by = c("Group", "Metric")) %>%
  mutate(
    Difference_from_n14 = Median - Ref_Median,
    Ratio_to_n14 = ifelse(is.na(Ref_Median) | Ref_Median == 0, NA_real_, Median / Ref_Median),
    Same_sign_as_n14 = sign(Median) == sign(Ref_Median)
  )

###############################################################################
## 7. Save outputs

write.csv(curve_long, file.path(outdir, "Zeta_commonN_curve_long_all_iterations.csv"), row.names = FALSE)
write.csv(curve_summary, file.path(outdir, "Zeta_commonN_curve_summary.csv"), row.names = FALSE)
write.csv(aic_long, file.path(outdir, "Zeta_commonN_AIC_long_all_iterations.csv"), row.names = FALSE)
write.csv(aic_summary, file.path(outdir, "Zeta_commonN_AIC_summary.csv"), row.names = FALSE)
write.csv(metrics_long, file.path(outdir, "Zeta_commonN_metrics_long_all_iterations.csv"), row.names = FALSE)
write.csv(metrics_summary, file.path(outdir, "Zeta_commonN_metrics_summary.csv"), row.names = FALSE)
write.csv(framework_by_iter, file.path(outdir, "Zeta_commonN_framework_by_iteration.csv"), row.names = FALSE)
write.csv(framework_summary, file.path(outdir, "Zeta_commonN_framework_summary.csv"), row.names = FALSE)
write.csv(framework_prob, file.path(outdir, "Zeta_commonN_framework_probability_summary.csv"), row.names = FALSE)
write.csv(trend_consistency_vs_n14, file.path(outdir, "Zeta_commonN_trend_consistency_vs_n14.csv"), row.names = FALSE)
write.csv(sample_n, file.path(outdir, "Zeta_commonN_sample_numbers.csv"), row.names = FALSE)

parameters = tibble(
  Item = c("Purpose", "Input", "Binary transformation", "Common-n settings",
           "Primary common-n", "Subsampling iterations",
           "Monte Carlo samples in Zeta.decline.mc", "Zeta settings",
           "Main framework comparison", "Model comparison"),
  Description = c(
    "Common-n subsampling sensitivity analysis for zeta diversity underlying Fig. S2B-D.",
    "OTUtable.txt; samples as rows after transposition and OTUs as columns.",
    "All nonzero OTU abundances were converted to presence/absence.",
    paste(n_grid, collapse = ", "),
    "n = 14, the largest replicate number shared by all dilution treatments.",
    as.character(n_iter),
    as.character(zeta_sam),
    paste0("Zeta.decline.mc with sam = ", zeta_sam,
           ", confint.level = ", confint.level,
           ", rescale = TRUE, normalize = TRUE, empty.row = 'empty'."),
    "High-order zeta diversity, zeta decay slopes and zeta-ratio patterns were compared across dilution treatments.",
    "AIC values of exponential and power-law regressions were extracted from Zeta.decline.mc outputs; DeltaAIC_ExpMinusPower > 0 indicates stronger support for the power-law model."
  )
)

write.csv(parameters, file.path(outdir, "Zeta_commonN_analysis_parameters.csv"), row.names = FALSE)

if (has_openxlsx) {
  wb = openxlsx::createWorkbook()
  sheet_list = list(
    README = parameters,
    Sample_numbers = sample_n,
    Curve_summary = curve_summary,
    AIC_summary = aic_summary,
    Metrics_summary = metrics_summary,
    Framework_probability = framework_prob,
    Framework_summary = framework_summary,
    Trend_vs_n14 = trend_consistency_vs_n14
  )
  for (s in names(sheet_list)) {
    openxlsx::addWorksheet(wb, s)
    openxlsx::writeData(wb, s, sheet_list[[s]])
  }
  openxlsx::saveWorkbook(wb, file.path(outdir, "Table_Sx_Zeta_commonN_sensitivity.xlsx"), overwrite = TRUE)
}

###############################################################################
## 8. Plotting

mycols = pal_nejm()(6)

theme_zeta = theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold", size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 9, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8))

plot_n14 = curve_summary %>% filter(Common_n == 14)

p_zeta_n14 = ggplot(plot_n14, aes(x = Order, y = Median_Zeta, color = Group, fill = Group, group = Group)) +
  geom_ribbon(aes(ymin = Q2.5_Zeta, ymax = Q97.5_Zeta), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = mycols, breaks = treat_levels) +
  scale_fill_manual(values = mycols, breaks = treat_levels) +
  labs(x = "Order", y = "Zeta diversity", title = "i. Primary common-n zeta diversity decay (n = 14)") +
  theme_zeta

p_ratio_n14 = ggplot(plot_n14, aes(x = Order, y = Median_Ratio, color = Group, fill = Group, group = Group)) +
  geom_ribbon(aes(ymin = Q2.5_Ratio, ymax = Q97.5_Ratio), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = mycols, breaks = treat_levels) +
  scale_fill_manual(values = mycols, breaks = treat_levels) +
  labs(x = "Order", y = "Zeta diversity ratio", title = "ii. Primary common-n zeta ratio (n = 14)") +
  theme_zeta

aic_plot_n14 = aic_summary %>%
  filter(Common_n == 14) %>%
  mutate(Model = factor(Model, levels = c("Exponential", "PowerLaw")))

p_aic_n14 = ggplot(aic_plot_n14, aes(x = Group, y = Median_AIC, fill = Model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_errorbar(aes(ymin = Q2.5_AIC, ymax = Q97.5_AIC),
                position = position_dodge(width = 0.75), width = 0.15, linewidth = 0.3) +
  scale_fill_manual(values = c("Exponential" = "#9DBA55", "PowerLaw" = "#44A9BD")) +
  labs(x = "Dilution treatment", y = "AIC", title = "iii. Regression fitting to zeta decay (n = 14)") +
  theme_zeta

if (requireNamespace("patchwork", quietly = TRUE)) {
  fig_primary = p_zeta_n14 + p_ratio_n14 + p_aic_n14 + patchwork::plot_layout(nrow = 1, guides = "collect")
  ggsave(file.path(outdir, "FigS_Zeta_commonN_primary_n14.pdf"), fig_primary, width = 9, height = 3)
  ggsave(file.path(outdir, "FigS_Zeta_commonN_primary_n14.png"), fig_primary, width = 9, height = 3, dpi = 300)
}

framework_plot_df = framework_summary %>%
  filter(Framework_metric %in% c("ZetaLast_D3D7_minus_D0",
                                 "ZetaLast_D7_minus_D0",
                                 "ExpSlope_D3D7_minus_D0",
                                 "Spearman_ZetaLast_vs_DilutionRank",
                                 "Prop_Groups_PowerLawBetter")) %>%
  mutate(Framework_metric = recode(
    Framework_metric,
    "ZetaLast_D3D7_minus_D0" = "High-order zeta: D3-D7 minus D0",
    "ZetaLast_D7_minus_D0" = "High-order zeta: D7 minus D0",
    "ExpSlope_D3D7_minus_D0" = "Exponential decay slope: D3-D7 minus D0",
    "Spearman_ZetaLast_vs_DilutionRank" = "Spearman rho: high-order zeta vs dilution rank",
    "Prop_Groups_PowerLawBetter" = "Proportion of groups: power-law AIC < exponential AIC"
  ))

p_framework = ggplot(framework_plot_df, aes(x = Common_n, y = Median)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.25, linewidth = 0.35, color = "grey40") +
  geom_point(size = 2.4, color = "#2C7FB8") +
  facet_wrap(~ Framework_metric, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = n_grid) +
  labs(x = "Common replicate number per treatment",
       y = "Median across subsampling iterations",
       title = "Common-n sensitivity of zeta framework metrics") +
  theme_zeta +
  theme(legend.position = "none")

ggsave(file.path(outdir, "FigS_Zeta_commonN_framework_sensitivity.pdf"), p_framework, width = 5, height = 7)
ggsave(file.path(outdir, "FigS_Zeta_commonN_framework_sensitivity.png"), p_framework, width = 5, height = 7, dpi = 300)

sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Done. Outputs saved in: ", outdir)
