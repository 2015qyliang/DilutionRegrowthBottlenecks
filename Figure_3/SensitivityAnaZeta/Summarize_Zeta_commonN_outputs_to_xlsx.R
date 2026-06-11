# coding: utf-8
# Summarize existing common-n zeta outputs into Table_Sx_Zeta_commonN_sensitivity.xlsx
# This script does NOT rerun Zeta.decline.mc(). It only reads existing CSV outputs.

library(tidyverse)

has_openxlsx <- requireNamespace("openxlsx", quietly = TRUE)
if (!has_openxlsx) {
  stop("Package 'openxlsx' is required. Please install it with install.packages('openxlsx').")
}

###############################################################################
## 1. Parameters

outdir <- "zeta_commonN_outputs_parallel"
n_grid <- c(8, 10, 12, 14)
treat_levels <- c("D0", "D1", "D3", "D5", "D6", "D7")
dilution_rank <- c(D0 = 0, D1 = 1, D3 = 3, D5 = 5, D6 = 6, D7 = 7)
ref_n <- 14

xlsx_out <- file.path(outdir, "Table_Sx_Zeta_commonN_sensitivity.xlsx")

qfun <- function(x, p) {
  if (all(is.na(x))) return(NA_real_)
  as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE))
}

###############################################################################
## 2. Check and read existing CSV files

need_files <- unlist(lapply(n_grid, function(n) {
  file.path(outdir, c(
    paste0("curve_common_n", n, ".csv"),
    paste0("aic_common_n", n, ".csv"),
    paste0("metrics_common_n", n, ".csv")
  ))
}))

missing_files <- need_files[!file.exists(need_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

curve_long <- bind_rows(lapply(n_grid, function(n) {
  read.csv(file.path(outdir, paste0("curve_common_n", n, ".csv")),
           stringsAsFactors = FALSE)
})) %>%
  mutate(Group = factor(Group, levels = treat_levels),
         Common_n = as.integer(Common_n),
         Iteration = as.integer(Iteration),
         Order = as.integer(Order))

aic_long <- bind_rows(lapply(n_grid, function(n) {
  read.csv(file.path(outdir, paste0("aic_common_n", n, ".csv")),
           stringsAsFactors = FALSE)
})) %>%
  mutate(Group = factor(Group, levels = treat_levels),
         Common_n = as.integer(Common_n),
         Iteration = as.integer(Iteration),
         Model = as.character(Model))

metrics_long <- bind_rows(lapply(n_grid, function(n) {
  read.csv(file.path(outdir, paste0("metrics_common_n", n, ".csv")),
           stringsAsFactors = FALSE)
})) %>%
  mutate(Group = factor(Group, levels = treat_levels),
         Common_n = as.integer(Common_n),
         Iteration = as.integer(Iteration),
         Dilution_rank = dilution_rank[as.character(Group)])

###############################################################################
## 3. Summary tables

curve_summary <- curve_long %>%
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
    N_iterations = dplyr::n_distinct(Iteration),
    .groups = "drop"
  )

aic_summary <- aic_long %>%
  group_by(Common_n, Group, Model) %>%
  summarise(
    Mean_AIC = mean(AIC, na.rm = TRUE),
    Median_AIC = median(AIC, na.rm = TRUE),
    Q2.5_AIC = qfun(AIC, 0.025),
    Q97.5_AIC = qfun(AIC, 0.975),
    N_iterations = dplyr::n_distinct(Iteration),
    .groups = "drop"
  )

metric_cols <- c("Exp_slope", "Power_slope", "Zeta_order1", "Zeta_order2",
                 "Zeta_order3", "Zeta_last", "Zeta_half_order",
                 "Last_over_first", "DeltaAIC_ExpMinusPower")
metric_cols <- intersect(metric_cols, colnames(metrics_long))

metrics_summary <- metrics_long %>%
  pivot_longer(cols = all_of(metric_cols),
               names_to = "Metric", values_to = "Value") %>%
  group_by(Common_n, Group, Metric) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Q2.5 = qfun(Value, 0.025),
    Q97.5 = qfun(Value, 0.975),
    N_iterations = dplyr::n_distinct(Iteration),
    .groups = "drop"
  )

###############################################################################
## 4. Framework-level summaries by iteration

framework_by_iter <- metrics_long %>%
  select(Common_n, Iteration, Group, Dilution_rank,
         any_of(c("Exp_slope", "Power_slope", "Zeta_last", "Zeta_half_order",
                  "Last_over_first", "DeltaAIC_ExpMinusPower", "PowerLaw_better"))) %>%
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

    Prop_Groups_PowerLawBetter = if ("PowerLaw_better" %in% colnames(cur_data())) {
      mean(PowerLaw_better, na.rm = TRUE)
    } else {
      NA_real_
    },
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

framework_metric_cols <- c("ZetaLast_D3D7_minus_D0", "ZetaLast_D5D7_minus_D0",
                           "ZetaLast_D7_minus_D0", "ExpSlope_D3D7_minus_D0",
                           "ExpSlope_D7_minus_D0", "PowerSlope_D3D7_minus_D0",
                           "PowerSlope_D7_minus_D0", "Spearman_ZetaLast_vs_DilutionRank",
                           "Spearman_ExpSlope_vs_DilutionRank", "Spearman_PowerSlope_vs_DilutionRank",
                           "Prop_Groups_PowerLawBetter")
framework_metric_cols <- intersect(framework_metric_cols, colnames(framework_by_iter))

framework_summary <- framework_by_iter %>%
  pivot_longer(cols = all_of(framework_metric_cols),
               names_to = "Framework_metric", values_to = "Value") %>%
  group_by(Common_n, Framework_metric) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Q2.5 = qfun(Value, 0.025),
    Q97.5 = qfun(Value, 0.975),
    N_iterations = dplyr::n_distinct(Iteration),
    .groups = "drop"
  )

framework_prob <- framework_by_iter %>%
  group_by(Common_n) %>%
  summarise(
    Prob_Stronger_lower_high_order_zeta = mean(Stronger_lower_high_order_zeta, na.rm = TRUE),
    Prob_D7_lower_high_order_zeta = mean(D7_lower_high_order_zeta, na.rm = TRUE),
    Prob_Stronger_more_negative_exp_slope = mean(Stronger_more_negative_exp_slope, na.rm = TRUE),
    Prob_D7_more_negative_exp_slope = mean(D7_more_negative_exp_slope, na.rm = TRUE),
    Prob_ZetaLast_decreases_with_dilution = mean(ZetaLast_decreases_with_dilution, na.rm = TRUE),
    Mean_Prop_Groups_PowerLawBetter = mean(Prop_Groups_PowerLawBetter, na.rm = TRUE),
    N_iterations = dplyr::n_distinct(Iteration),
    .groups = "drop"
  )

###############################################################################
## 5. Trend consistency relative to n = 14

reference_metrics <- metrics_summary %>%
  filter(Common_n == ref_n) %>%
  select(Group, Metric, Ref_Median = Median)

trend_consistency_vs_n14 <- metrics_summary %>%
  left_join(reference_metrics, by = c("Group", "Metric")) %>%
  mutate(
    Difference_from_n14 = Median - Ref_Median,
    Ratio_to_n14 = ifelse(is.na(Ref_Median) | Ref_Median == 0, NA_real_, Median / Ref_Median),
    Same_sign_as_n14 = sign(Median) == sign(Ref_Median)
  )

###############################################################################
## 6. Analysis parameters / README

sample_n_path <- file.path(outdir, "Zeta_commonN_sample_numbers.csv")
if (file.exists(sample_n_path)) {
  sample_n <- read.csv(sample_n_path, stringsAsFactors = FALSE)
} else if (file.exists("OTUtable.txt")) {
  otu <- read.table("OTUtable.txt", header = TRUE, sep = "\t", row.names = 1,
                    check.names = FALSE, stringsAsFactors = FALSE)
  data.spec <- as.matrix(t(otu))
  sample_n <- tibble(
    Group = treat_levels,
    Original_replicates = sapply(treat_levels, function(g) sum(grepl(paste0("^", g), rownames(data.spec))))
  )
} else {
  sample_n <- tibble(Group = treat_levels,
                     Original_replicates = c(14, 14, 20, 20, 40, 40))
}

parameters <- tibble(
  Item = c("Purpose", "Input", "Source files", "Binary transformation", "Common-n settings",
           "Primary common-n", "Subsampling iterations", "Monte Carlo samples in Zeta.decline.mc",
           "Zeta settings", "Main framework comparison", "Model comparison"),
  Description = c(
    "Summary workbook for common-n subsampling sensitivity analysis of zeta diversity underlying Fig. S2B-D. This script does not rerun Zeta.decline.mc().",
    "Existing CSV outputs from the completed common-n zeta analysis.",
    "curve_common_n*.csv, aic_common_n*.csv, metrics_common_n*.csv.",
    "All nonzero OTU abundances were converted to presence/absence in the original zeta run.",
    paste(n_grid, collapse = ", "),
    paste0("n = ", ref_n, ", the largest replicate number shared by all dilution treatments."),
    paste(unique(metrics_long$Iteration) %>% length(), "iterations detected from CSV files."),
    "1000, according to the zeta run script.",
    "Zeta.decline.mc with confint.level = 0.95, rescale = TRUE, normalize = TRUE, empty.row = 'empty'.",
    "High-order zeta diversity, zeta decay slopes and zeta-ratio patterns were compared across dilution treatments.",
    "AIC values of exponential and power-law regressions were extracted from Zeta.decline.mc outputs; DeltaAIC_ExpMinusPower > 0 indicates stronger support for the power-law model."
  )
)

###############################################################################
## 7. Write CSV summaries and Excel workbook

write.csv(curve_summary, file.path(outdir, "Zeta_commonN_curve_summary.csv"), row.names = FALSE)
write.csv(aic_summary, file.path(outdir, "Zeta_commonN_AIC_summary.csv"), row.names = FALSE)
write.csv(metrics_summary, file.path(outdir, "Zeta_commonN_metrics_summary.csv"), row.names = FALSE)
write.csv(framework_by_iter, file.path(outdir, "Zeta_commonN_framework_by_iteration.csv"), row.names = FALSE)
write.csv(framework_summary, file.path(outdir, "Zeta_commonN_framework_summary.csv"), row.names = FALSE)
write.csv(framework_prob, file.path(outdir, "Zeta_commonN_framework_probability_summary.csv"), row.names = FALSE)
write.csv(trend_consistency_vs_n14, file.path(outdir, "Zeta_commonN_trend_consistency_vs_n14.csv"), row.names = FALSE)
write.csv(sample_n, file.path(outdir, "Zeta_commonN_sample_numbers.csv"), row.names = FALSE)
write.csv(parameters, file.path(outdir, "Zeta_commonN_analysis_parameters.csv"), row.names = FALSE)

wb <- openxlsx::createWorkbook()
sheet_list <- list(
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
  openxlsx::freezePane(wb, s, firstRow = TRUE)
}

openxlsx::saveWorkbook(wb, xlsx_out, overwrite = TRUE)

message("Done. Workbook written to: ", xlsx_out)
