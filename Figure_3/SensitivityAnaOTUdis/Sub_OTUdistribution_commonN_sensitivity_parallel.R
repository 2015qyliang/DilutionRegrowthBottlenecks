# coding: utf-8
# Fig. 3A OTU distribution-state common-n sensitivity analysis
# Purpose: lightweight replicate-equalized sensitivity for OTU distribution-state diagnostics
# Author style follows the original SubSenOTUdistribution.R workflow

library(tidyverse) 
library(ggplot2)
library(ggsci)
library(future)
library(future.apply)
library(openxlsx)

##########################################################################
## 1. Parameters

set.seed(123)

treat_levels = c("D0", "D1", "D3", "D5", "D6", "D7")
dilution_rank = c(D0 = 0, D1 = 1, D3 = 3, D5 = 5, D6 = 6, D7 = 7)

## Common-n settings. n = 14 is the primary common-n analysis.
n_grid = c(8, 10, 12, 14)
primary_n = 14

## Outer subsampling iterations.
## Use n_iter = 100 for quick run; increase to 1000 for final Table S8 if time permits.
n_iter = 1000

## Significance cutoff used in the original Fig. 3A script.
p_cut = 0.05
q_cut = 0.05

## KS tests are relatively slow and sensitive to ties/zero-inflation.
## Keep TRUE to reproduce the original Fig. 3A logic.
run_ks = TRUE

## Parallel settings for a Windows workstation such as i7-11700T / 32 GB RAM.
## Recommended final setting: use_parallel = TRUE, n_workers = 6.
## If the computer becomes unresponsive, reduce n_workers to 4.
use_parallel = TRUE
n_workers = 16
future_max_size_GB = 26

out_dir = "commonN_outputs"
# if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

##########################################################################
## 2. Read OTU table

otu_raw = read.table("OTUtable.txt", header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE)

if (!"OTUID" %in% colnames(otu_raw)) {
  stop("The first column should be named OTUID.")
}

otu_id = otu_raw$OTUID
otu_mat = as.matrix(otu_raw[, setdiff(colnames(otu_raw), "OTUID"), drop = FALSE])
mode(otu_mat) = "numeric"
rownames(otu_mat) = otu_id

sample_groups = list(
  D0 = grep("^D0", colnames(otu_mat), value = TRUE),
  D1 = grep("^D1", colnames(otu_mat), value = TRUE),
  D3 = grep("^D3", colnames(otu_mat), value = TRUE),
  D5 = grep("^D5", colnames(otu_mat), value = TRUE),
  D6 = grep("^D6", colnames(otu_mat), value = TRUE),
  D7 = grep("^D7", colnames(otu_mat), value = TRUE)
)

sample_number = tibble(
  Group = names(sample_groups),
  Original_n = as.integer(sapply(sample_groups, length))
)

print(sample_number)

if (any(sample_number$Original_n < max(n_grid))) {
  stop("At least one treatment has fewer replicates than max(n_grid).")
}

##########################################################################
## 3. Helper functions

safe_var = function(x) {
  if (length(x) <= 1) return(NA_real_)
  var(x, na.rm = TRUE)
}

## Fast chi-square test equivalent to chisq.test(x) under equal expected counts.
fast_chisq_p = function(x) {
  x = as.numeric(x)
  if (all(is.na(x)) || sum(x, na.rm = TRUE) <= 0) return(c(stat = NA_real_, p = NA_real_))
  n = length(x)
  expected = sum(x, na.rm = TRUE) / n
  if (expected <= 0) return(c(stat = NA_real_, p = NA_real_))
  stat = sum((x - expected)^2 / expected, na.rm = TRUE)
  p = stats::pchisq(stat, df = n - 1, lower.tail = FALSE)
  c(stat = stat, p = p)
}

safe_ks_p = function(x) {
  x = as.numeric(x)
  if (all(is.na(x)) || sum(x != 0, na.rm = TRUE) == 0) return(c(stat = NA_real_, p = NA_real_))
  sdx = stats::sd(x, na.rm = TRUE)
  mux = mean(x, na.rm = TRUE)
  if (is.na(sdx) || sdx == 0) return(c(stat = NA_real_, p = NA_real_))
  kt = tryCatch(
    suppressWarnings(stats::ks.test(x, "pnorm", mux, sdx)),
    error = function(e) NULL
  )
  if (is.null(kt)) return(c(stat = NA_real_, p = NA_real_))
  c(stat = unname(kt$statistic), p = kt$p.value)
}

calc_distribution_for_group = function(mat_group, group_name, common_n, iteration) {
  ## mat_group: OTUs × sampled replicates
  detected = rowSums(mat_group > 0, na.rm = TRUE) > 0
  mat_det = mat_group[detected, , drop = FALSE]

  if (nrow(mat_det) == 0) {
    return(tibble())
  }

  fM = rowMeans(mat_det, na.rm = TRUE)
  fV = apply(mat_det, 1, safe_var)
  fN = rowSums(mat_det != 0, na.rm = TRUE)
  occupancy = fN / ncol(mat_det)

  ## Index of dispersion used in the original script: log10[(variance/mean) × occurrence count]
  iod_raw = (fV / fM) * fN
  iodValue = ifelse(is.finite(iod_raw) & iod_raw > 0, log10(iod_raw), NA_real_)

  chisq_mat = t(apply(mat_det, 1, fast_chisq_p))
  colnames(chisq_mat) = c("rChisqX", "rChisqP")

  if (run_ks) {
    ks_mat = t(apply(mat_det, 1, safe_ks_p))
    colnames(ks_mat) = c("rKsD", "rKsP")
  } else {
    ks_mat = cbind(rKsD = rep(NA_real_, nrow(mat_det)),
                   rKsP = rep(NA_real_, nrow(mat_det)))
  }

  out = tibble(
    Common_n = common_n,
    Iteration = iteration,
    Group = group_name,
    OTU = rownames(mat_det),
    fV = as.numeric(fV),
    fM = as.numeric(fM),
    fN = as.numeric(fN),
    Occupancy = as.numeric(occupancy),
    iodValue = as.numeric(iodValue),
    rChisqX = as.numeric(chisq_mat[, "rChisqX"]),
    rChisqP = as.numeric(chisq_mat[, "rChisqP"]),
    rKsD = as.numeric(ks_mat[, "rKsD"]),
    rKsP = as.numeric(ks_mat[, "rKsP"])
  ) %>%
    group_by(Common_n, Iteration, Group) %>%
    mutate(
      rChisqQ = p.adjust(rChisqP, method = "BH"),
      rKsQ = p.adjust(rKsP, method = "BH"),
      ChisqRes_P = if_else(!is.na(rChisqP) & rChisqP <= p_cut, "NonUniform_P", "Uniform_P"),
      ChisqRes_Q = if_else(!is.na(rChisqQ) & rChisqQ <= q_cut, "NonUniform_Q", "Uniform_Q"),
      KsRes_P = if_else(!is.na(rKsP) & rKsP <= p_cut, "NotNorm_P", "Norm_P"),
      KsRes_Q = if_else(!is.na(rKsQ) & rKsQ <= q_cut, "NotNorm_Q", "Norm_Q"),
      OccClass = case_when(
        Occupancy == 0 ~ "Lost",
        Occupancy > 0 & Occupancy <= 0.25 ~ "Transient",
        Occupancy > 0.25 & Occupancy <= 0.50 ~ "Intermediate",
        Occupancy > 0.50 & Occupancy < 1 ~ "High",
        Occupancy == 1 ~ "Persistent",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup()

  out
}

summarise_iteration = function(otu_level_df) {
  otu_level_df %>%
    group_by(Common_n, Iteration, Group) %>%
    summarise(
      N_OTUs_tested = n(),
      Fraction_nonuniform_P = mean(rChisqP <= p_cut, na.rm = TRUE),
      Fraction_nonuniform_Q = mean(rChisqQ <= q_cut, na.rm = TRUE),
      Fraction_nonnormal_P = mean(rKsP <= p_cut, na.rm = TRUE),
      Fraction_nonnormal_Q = mean(rKsQ <= q_cut, na.rm = TRUE),
      Median_IOD = median(iodValue, na.rm = TRUE),
      Mean_IOD = mean(iodValue, na.rm = TRUE),
      Median_occupancy = median(Occupancy, na.rm = TRUE),
      Fraction_transient = mean(OccClass == "Transient", na.rm = TRUE),
      Fraction_low_occupancy = mean(Occupancy <= 0.25, na.rm = TRUE),
      Fraction_high_or_persistent = mean(Occupancy > 0.5, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_quantiles = function(iter_summary) {
  metric_cols = c(
    "N_OTUs_tested",
    "Fraction_nonuniform_P", "Fraction_nonuniform_Q",
    "Fraction_nonnormal_P", "Fraction_nonnormal_Q",
    "Median_IOD", "Mean_IOD", "Median_occupancy",
    "Fraction_transient", "Fraction_low_occupancy", "Fraction_high_or_persistent"
  )

  iter_summary %>%
    pivot_longer(cols = all_of(metric_cols), names_to = "Metric", values_to = "Value") %>%
    group_by(Common_n, Group, Metric) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      Median = median(Value, na.rm = TRUE),
      Q2.5 = quantile(Value, 0.025, na.rm = TRUE),
      Q97.5 = quantile(Value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

calc_framework_by_iteration = function(iter_summary) {
  metric_long = iter_summary %>%
    pivot_longer(cols = -c(Common_n, Iteration, Group), names_to = "Metric", values_to = "Value")

  metric_wide = metric_long %>%
    select(Common_n, Iteration, Group, Metric, Value) %>%
    pivot_wider(names_from = Group, values_from = Value)

  metric_wide %>%
    rowwise() %>%
    mutate(
      Spearman_vs_dilution_rank = suppressWarnings(cor(
        c(D0, D1, D3, D5, D6, D7),
        c(0, 1, 3, 5, 6, 7),
        method = "spearman",
        use = "complete.obs"
      )),
      D7_minus_D0 = D7 - D0,
      D6D7_mean_minus_D0 = mean(c(D6, D7), na.rm = TRUE) - D0,
      D3D7_mean_minus_D0 = mean(c(D3, D5, D6, D7), na.rm = TRUE) - D0,
      Stronger_mean_minus_D0D1 = mean(c(D3, D5, D6, D7), na.rm = TRUE) - mean(c(D0, D1), na.rm = TRUE)
    ) %>%
    ungroup()
}

##########################################################################
## 4. Main common-n subsampling loop, with optional parallel computing

## On Windows, multisession is safer than multicore. For i7-11700T / 32 GB RAM,
## n_workers = 6 is usually a good balance between speed and stability.
if (use_parallel && has_future) {
  workers_use = min(n_workers, max(1, parallel::detectCores(logical = TRUE) - 2))
  message("Using parallel computing with ", workers_use, " workers.")
  options(future.globals.maxSize = future_max_size_GB * 1024^3)
  future::plan(future::multisession, workers = workers_use)
} else {
  if (use_parallel && !has_future) {
    message("Packages future/future.apply are not available. Falling back to serial lapply().")
  } else {
    message("Using serial lapply().")
  }
}

## One iteration = randomly subsampling common_n replicates from each treatment,
## then calculating OTU distribution diagnostics for all six treatments.
run_one_iteration = function(common_n, it) {

  ## Deterministic seed for reproducibility across serial/parallel runs.
  set.seed(1000000 + common_n * 10000 + it)

  sampled_groups = lapply(sample_groups, function(cols) {
    sample(cols, size = common_n, replace = FALSE)
  })

  otu_iter_all = map_dfr(names(sampled_groups), function(g) {
    mat_g = otu_mat[, sampled_groups[[g]], drop = FALSE]
    calc_distribution_for_group(mat_g, group_name = g,
                                common_n = common_n, iteration = it)
  })

  iter_summary = summarise_iteration(otu_iter_all)

  ## Save OTU-level details only for the first 20 iterations of primary_n.
  ## This avoids very large output while retaining diagnostic/reproducibility records.
  primary_otu = NULL
  if (common_n == primary_n && it <= 20) {
    primary_otu = otu_iter_all
  }

  list(summary = iter_summary, primary_otu = primary_otu)
}

all_iteration_summary = list()
all_primary_otu_level = list()

for (common_n in n_grid) {
  message("===== Common n = ", common_n, " =====")

  iter_index = seq_len(n_iter)

  if (use_parallel && has_future) {
    res_list = future.apply::future_lapply(
      iter_index,
      function(it) run_one_iteration(common_n = common_n, it = it),
      future.seed = TRUE
    )
  } else {
    res_list = lapply(iter_index, function(it) {
      if (it %% 50 == 0 || it == 1 || it == n_iter) {
        message("Running common n = ", common_n, ", iteration ", it, " / ", n_iter)
      }
      run_one_iteration(common_n = common_n, it = it)
    })
  }

  iter_summary_n = bind_rows(lapply(res_list, `[[`, "summary"))
  all_iteration_summary[[as.character(common_n)]] = iter_summary_n

  write.csv(iter_summary_n,
            file.path(out_dir, paste0("Fig3A_distribution_iteration_summary_n", common_n, ".csv")),
            row.names = FALSE)

  ## Save checkpoint after each common-n level. If the script stops later,
  ## completed common-n outputs are still preserved.
  saveRDS(iter_summary_n,
          file.path(out_dir, paste0("checkpoint_Fig3A_iteration_summary_n", common_n, ".rds")))

  if (common_n == primary_n) {
    primary_otu_list = lapply(res_list, `[[`, "primary_otu")
    primary_otu_list = primary_otu_list[!vapply(primary_otu_list, is.null, logical(1))]
    primary_otu_df = bind_rows(primary_otu_list)

    write.csv(primary_otu_df,
              file.path(out_dir, paste0("Fig3A_distribution_primary_n", common_n,
                                        "_otu_level_first20iterations.csv")),
              row.names = FALSE)
  }
}

## Reset future plan to avoid leaving background workers active in RStudio.
if (use_parallel && has_future) {
  future::plan(future::sequential)
}

iteration_summary = bind_rows(all_iteration_summary)
write.csv(iteration_summary,
          file.path(out_dir, "Fig3A_distribution_iteration_summary_all_n.csv"),
          row.names = FALSE)

summary_quantiles = summarise_quantiles(iteration_summary)
write.csv(summary_quantiles,
          file.path(out_dir, "Fig3A_distribution_commonN_summary_quantiles.csv"),
          row.names = FALSE)

framework_by_iteration = calc_framework_by_iteration(iteration_summary)
write.csv(framework_by_iteration,
          file.path(out_dir, "Fig3A_distribution_framework_by_iteration.csv"),
          row.names = FALSE)

framework_summary = framework_by_iteration %>%
  pivot_longer(cols = c(Spearman_vs_dilution_rank, D7_minus_D0, D6D7_mean_minus_D0,
                        D3D7_mean_minus_D0, Stronger_mean_minus_D0D1),
               names_to = "Contrast", values_to = "Value") %>%
  group_by(Common_n, Metric, Contrast) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Q2.5 = quantile(Value, 0.025, na.rm = TRUE),
    Q97.5 = quantile(Value, 0.975, na.rm = TRUE),
    Prob_positive = mean(Value > 0, na.rm = TRUE),
    Prob_negative = mean(Value < 0, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(framework_summary,
          file.path(out_dir, "Fig3A_distribution_framework_summary.csv"),
          row.names = FALSE)

##########################################################################
## 5. Trend consistency relative to primary n = 14

reference = framework_summary %>%
  filter(Common_n == primary_n) %>%
  select(Metric, Contrast, Ref_Median = Median)

trend_vs_ref = framework_summary %>%
  left_join(reference, by = c("Metric", "Contrast")) %>%
  mutate(
    Difference_from_ref = Median - Ref_Median,
    Ratio_to_ref = if_else(!is.na(Ref_Median) & Ref_Median != 0, Median / Ref_Median, NA_real_),
    Same_sign_as_ref = sign(Median) == sign(Ref_Median)
  )

write.csv(trend_vs_ref,
          file.path(out_dir, "Fig3A_distribution_trend_consistency_vs_n14.csv"),
          row.names = FALSE)

##########################################################################
## 6. Optional Excel workbook

readme = tibble(
  Item = c(
    "Purpose",
    "Original Fig. 3A metrics",
    "Common-n settings",
    "Primary common-n",
    "Iterations",
    "FDR control",
    "Interpretation boundary"
  ),
  Description = c(
    "Common-n subsampling sensitivity analysis for Fig. 3A OTU distribution-state diagnostics.",
    "For each OTU within each dilution treatment, variance, mean, occurrence count, index of dispersion, Pearson chi-square P and KS P were calculated following the original script.",
    paste(n_grid, collapse = ", "),
    paste0("n = ", primary_n, ", the largest replicate number shared by all dilution treatments."),
    paste0(n_iter, " subsampling iterations for each common-n setting."),
    "Benjamini-Hochberg q-values were calculated within each dilution treatment and iteration for both chi-square and KS tests.",
    "These tests are treated as distributional diagnostics of sparsity and among-replicate heterogeneity, not as direct assembly-mechanism tests."
  )
)

if (has_openxlsx) {
  wb = openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "README")
  openxlsx::writeData(wb, "README", readme)

  openxlsx::addWorksheet(wb, "Sample_numbers")
  openxlsx::writeData(wb, "Sample_numbers", sample_number)

  openxlsx::addWorksheet(wb, "Summary_quantiles")
  openxlsx::writeData(wb, "Summary_quantiles", summary_quantiles)

  openxlsx::addWorksheet(wb, "Framework_summary")
  openxlsx::writeData(wb, "Framework_summary", framework_summary)

  openxlsx::addWorksheet(wb, "Trend_vs_n14")
  openxlsx::writeData(wb, "Trend_vs_n14", trend_vs_ref)

  openxlsx::addWorksheet(wb, "Iteration_summary")
  openxlsx::writeData(wb, "Iteration_summary", iteration_summary)

  openxlsx::saveWorkbook(wb,
                         file.path(out_dir, "Table_S8_Fig3A_distribution_commonN_sensitivity.xlsx"),
                         overwrite = TRUE)
}

##########################################################################
## 7. Lightweight plots

# plot_metrics = c(
#   "Fraction_nonuniform_P",
#   "Fraction_nonuniform_Q",
#   "Fraction_nonnormal_P",
#   "Fraction_nonnormal_Q",
#   "Median_IOD",
#   "Fraction_low_occupancy"
# )
# 
# plot_df = summary_quantiles %>%
#   filter(Metric %in% plot_metrics) %>%
#   mutate(
#     Metric = factor(Metric, levels = plot_metrics),
#     Group = factor(Group, levels = treat_levels)
#   )
# 
# p_commonN = ggplot(plot_df, aes(x = Group, y = Median, color = factor(Common_n))) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
#   geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5),
#                 width = 0.15, alpha = 0.6, linewidth = 0.35,
#                 position = position_dodge(width = 0.45)) +
#   geom_point(size = 1.8, alpha = 0.9,
#              position = position_dodge(width = 0.45)) +
#   facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
#   scale_color_brewer(palette = "Dark2", name = "Common n") +
#   labs(x = "Dilution treatment", y = "Median across subsampling iterations",
#        title = "Common-n sensitivity of Fig. 3A distribution-state diagnostics") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         plot.title = element_text(face = "bold", size = 10),
#         axis.text.x = element_text(size = 8, color = "black"),
#         axis.text.y = element_text(size = 8, color = "black"),
#         legend.position = "bottom")
# 
# ggsave(file.path(out_dir, "FigS_Fig3A_distribution_commonN_sensitivity.pdf"),
#        p_commonN, width = 7.5, height = 6)
# 


