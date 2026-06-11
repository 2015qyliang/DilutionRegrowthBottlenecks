############################################################
## Sensitivity analysis of treatment-level OTU intersections
## Common-n subsampling grid: n = 8, 10, 12, 14
## Based on the original n = 14 Fig. 2E sensitivity script
############################################################

library(tidyverse)
library(readr)
library(openxlsx)

############################################################
## 0. Parameters

set.seed(123)

otu_file = "../OTUtable.txt"
out_xlsx = "Table S7. Sensitivity analysis of OTU intersections_commonN.xlsx"
out_dir  = "OTU_intersection_commonN_outputs"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

treatment_order = c("D0", "D1", "D3", "D5", "D6", "D7")
n_grid = c(8, 10, 12, 14)
n_iter = 1000

## Treatment-level presence definition used in Fig. 2E:
## an OTU is present in a dilution treatment if detected in at least one replicate.
main_presence_mode = "min_replicates"
main_min_replicates = 1

############################################################
## 1. Read OTU table

otu = read_tsv(otu_file, show_col_types = FALSE)

if (!"OTUID" %in% colnames(otu)) {
  stop("The first column should be named OTUID.")
}

sample_groups = list(
  D0 = grep("^D0", colnames(otu), value = TRUE),
  D1 = grep("^D1", colnames(otu), value = TRUE),
  D3 = grep("^D3", colnames(otu), value = TRUE),
  D5 = grep("^D5", colnames(otu), value = TRUE),
  D6 = grep("^D6", colnames(otu), value = TRUE),
  D7 = grep("^D7", colnames(otu), value = TRUE)
)

sample_n = tibble(
  Treatment = names(sample_groups),
  No_replicates = sapply(sample_groups, length)
) %>%
  mutate(Treatment = factor(Treatment, levels = treatment_order)) %>%
  arrange(Treatment)

print(sample_n)

min_n = min(sample_n$No_replicates)
if (any(n_grid > min_n)) {
  stop("All n_grid values must be <= the minimum replicate number across treatments: ", min_n)
}

otu_num = otu %>%
  mutate(across(-OTUID, as.numeric))

############################################################
## 2. Functions

make_treatment_pa = function(otu_table, sample_groups,
                             mode = "min_replicates",
                             min_replicates = 1,
                             min_prop = NULL) {
  out = tibble(OTUID = otu_table$OTUID)
  
  for (g in names(sample_groups)) {
    cols = sample_groups[[g]]
    detected_n = rowSums(otu_table[, cols, drop = FALSE] > 0, na.rm = TRUE)
    
    if (mode == "min_replicates") {
      present = detected_n >= min_replicates
    } else if (mode == "min_prop") {
      present = detected_n / length(cols) >= min_prop
    } else {
      stop("Unknown mode: ", mode)
    }
    
    out[[g]] = present
  }
  
  out = out %>%
    filter(if_any(all_of(names(sample_groups)), ~ .x))
  
  return(out)
}

summarize_intersections = function(pa_table, treatment_cols = treatment_order) {
  ## Total OTUs per treatment
  total_metrics = map_dfr(treatment_cols, function(g) {
    tibble(
      Metric = paste0("Total OTUs in ", g),
      Value = sum(pa_table[[g]], na.rm = TRUE)
    )
  })
  
  ## Intersection membership for each OTU
  membership = pa_table %>%
    rowwise() %>%
    mutate(
      Degree = sum(c_across(all_of(treatment_cols)), na.rm = TRUE),
      Intersection_label = {
        present_groups = treatment_cols[c_across(all_of(treatment_cols))]
        if (length(present_groups) == 1) {
          paste0(present_groups, " only")
        } else {
          paste(present_groups, collapse = "&")
        }
      }
    ) %>%
    ungroup()
  
  intersection_summary = membership %>%
    count(Intersection_label, name = "Value")
  
  ## Key intersections shown or discussed in Fig. 2E and manuscript text
  key_intersections = c(
    "D0 only", "D1 only", "D3 only", "D5 only", "D6 only", "D7 only",
    "D0&D1", "D0&D1&D3", "D0&D1&D3&D5&D6&D7",
    "D1&D3&D5&D6&D7", "D3&D5&D6&D7", "D6&D7"
  )
  
  key_metrics = map_dfr(key_intersections, function(k) {
    val = intersection_summary %>%
      filter(Intersection_label == k) %>%
      pull(Value)
    if (length(val) == 0) val = 0
    tibble(
      Metric = paste0("Intersection: ", k),
      Value = val
    )
  })
  
  bind_rows(total_metrics, key_metrics)
}

## Convert one iteration's metric table into framework support indicators.
## These indicators are designed to support the manuscript framework:
## 1) D1 richness / D1-specific fraction remains detectable.
## 2) D0-D7 all-treatment shared core persists.
## 3) Multi-dilution shared fractions not including D0 persist.
## 4) Strong-dilution shared fraction persists.
framework_support_from_metrics = function(metric_df) {
  x = metric_df %>%
    select(Metric, Value) %>%
    deframe()
  
  total_metrics = paste0("Total OTUs in ", treatment_order)
  total_values = x[total_metrics]
  names(total_values) = treatment_order
  
  tibble(
    D1_has_highest_total_OTUs = total_values["D1"] == max(total_values, na.rm = TRUE),
    D1_total_rank = rank(-total_values, ties.method = "min")["D1"],
    D1_specific_positive = x["Intersection: D1 only"] > 0,
    D0_D1_shared_positive = x["Intersection: D0&D1"] > 0,
    All_treatment_shared_positive = x["Intersection: D0&D1&D3&D5&D6&D7"] > 0,
    D1_D3_D5_D6_D7_shared_positive = x["Intersection: D1&D3&D5&D6&D7"] > 0,
    D3_D5_D6_D7_shared_positive = x["Intersection: D3&D5&D6&D7"] > 0,
    D6_D7_shared_positive = x["Intersection: D6&D7"] > 0,
    D1_specific_count = x["Intersection: D1 only"],
    All_treatment_shared_count = x["Intersection: D0&D1&D3&D5&D6&D7"],
    D1_D3_D5_D6_D7_shared_count = x["Intersection: D1&D3&D5&D6&D7"],
    D3_D5_D6_D7_shared_count = x["Intersection: D3&D5&D6&D7"]
  )
}

############################################################
## 3. Original observed values from full data

pa_original = make_treatment_pa(
  otu_num,
  sample_groups,
  mode = main_presence_mode,
  min_replicates = main_min_replicates
)

observed_metrics = summarize_intersections(pa_original) %>%
  rename(Observed_value = Value)

observed_framework = framework_support_from_metrics(
  observed_metrics %>% rename(Value = Observed_value)
) %>%
  mutate(Analysis = "Original full data")

############################################################
## 4. Common-n equalized replicate-number analysis

equalized_iterations = map_dfr(n_grid, function(target_n) {
  
  map_dfr(seq_len(n_iter), function(i) {
    
    if (i %% 50 == 0 || i == 1 || i == n_iter) {
      message("Running n = ", target_n, ", iteration ", i, " / ", n_iter)
    }
    
    sampled_groups = lapply(sample_groups, function(cols) {
      sample(cols, size = target_n, replace = FALSE)
    })
    
    pa_i = make_treatment_pa(
      otu_num,
      sampled_groups,
      mode = main_presence_mode,
      min_replicates = main_min_replicates
    )
    
    summarize_intersections(pa_i) %>%
      mutate(
        Analysis = "Equalized replicate number",
        Iteration = i,
        `No. replicates per treatment` = target_n
      )
  })
})

equalized_summary = equalized_iterations %>%
  group_by(Analysis, `No. replicates per treatment`, Metric) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    `2.5% quantile` = quantile(Value, 0.025, na.rm = TRUE),
    `97.5% quantile` = quantile(Value, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(observed_metrics, by = "Metric") %>%
  mutate(
    `Observed within 95% interval` =
      Observed_value >= `2.5% quantile` &
      Observed_value <= `97.5% quantile`,
    `Mean relative to observed` = if_else(Observed_value == 0, NA_real_, Mean / Observed_value),
    `Median relative to observed` = if_else(Observed_value == 0, NA_real_, Median / Observed_value)
  )

## n = 14 as the primary reference trend
equalized_n14_reference = equalized_summary %>%
  filter(`No. replicates per treatment` == 14) %>%
  select(Metric,
         Median_n14 = Median,
         Mean_n14 = Mean,
         `2.5% quantile_n14` = `2.5% quantile`,
         `97.5% quantile_n14` = `97.5% quantile`)

trend_vs_n14 = equalized_summary %>%
  left_join(equalized_n14_reference, by = "Metric") %>%
  mutate(
    `Median relative to n14` = if_else(Median_n14 == 0, NA_real_, Median / Median_n14),
    `Mean relative to n14` = if_else(Mean_n14 == 0, NA_real_, Mean / Mean_n14),
    `Median difference from n14` = Median - Median_n14,
    `95% interval overlaps n14 interval` =
      (`2.5% quantile` <= `97.5% quantile_n14`) &
      (`97.5% quantile` >= `2.5% quantile_n14`)
  )

############################################################
## 5. Framework-support summaries by iteration

framework_by_iteration = equalized_iterations %>%
  group_by(`No. replicates per treatment`, Iteration) %>%
  group_modify(~ framework_support_from_metrics(.x)) %>%
  ungroup() %>%
  mutate(Analysis = "Equalized replicate number")

framework_summary = framework_by_iteration %>%
  group_by(Analysis, `No. replicates per treatment`) %>%
  summarise(
    Prob_D1_has_highest_total_OTUs = mean(D1_has_highest_total_OTUs, na.rm = TRUE),
    Median_D1_total_rank = median(D1_total_rank, na.rm = TRUE),
    Prob_D1_specific_positive = mean(D1_specific_positive, na.rm = TRUE),
    Prob_D0_D1_shared_positive = mean(D0_D1_shared_positive, na.rm = TRUE),
    Prob_all_treatment_shared_positive = mean(All_treatment_shared_positive, na.rm = TRUE),
    Prob_D1_D3_D5_D6_D7_shared_positive = mean(D1_D3_D5_D6_D7_shared_positive, na.rm = TRUE),
    Prob_D3_D5_D6_D7_shared_positive = mean(D3_D5_D6_D7_shared_positive, na.rm = TRUE),
    Prob_D6_D7_shared_positive = mean(D6_D7_shared_positive, na.rm = TRUE),
    Median_D1_specific_count = median(D1_specific_count, na.rm = TRUE),
    D1_specific_count_2.5 = quantile(D1_specific_count, 0.025, na.rm = TRUE),
    D1_specific_count_97.5 = quantile(D1_specific_count, 0.975, na.rm = TRUE),
    Median_all_treatment_shared_count = median(All_treatment_shared_count, na.rm = TRUE),
    All_treatment_shared_count_2.5 = quantile(All_treatment_shared_count, 0.025, na.rm = TRUE),
    All_treatment_shared_count_97.5 = quantile(All_treatment_shared_count, 0.975, na.rm = TRUE),
    Median_D1_D3_D5_D6_D7_shared_count = median(D1_D3_D5_D6_D7_shared_count, na.rm = TRUE),
    D1_D3_D5_D6_D7_shared_count_2.5 = quantile(D1_D3_D5_D6_D7_shared_count, 0.025, na.rm = TRUE),
    D1_D3_D5_D6_D7_shared_count_97.5 = quantile(D1_D3_D5_D6_D7_shared_count, 0.975, na.rm = TRUE),
    Median_D3_D5_D6_D7_shared_count = median(D3_D5_D6_D7_shared_count, na.rm = TRUE),
    D3_D5_D6_D7_shared_count_2.5 = quantile(D3_D5_D6_D7_shared_count, 0.025, na.rm = TRUE),
    D3_D5_D6_D7_shared_count_97.5 = quantile(D3_D5_D6_D7_shared_count, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
## 6. Occupancy threshold sensitivity using full data
## This retains the original threshold analysis and can be shown next to common-n results.

threshold_settings = tribble(
  ~Threshold_type, ~mode, ~min_replicates, ~min_prop,
  ">=1 replicate",   "min_replicates", 1, NA_real_,
  ">=2 replicates",  "min_replicates", 2, NA_real_,
  ">=5% replicates", "min_prop", NA_real_, 0.05,
  ">=10% replicates", "min_prop", NA_real_, 0.10
)

occupancy_summary = pmap_dfr(
  threshold_settings,
  function(Threshold_type, mode, min_replicates, min_prop) {
    pa_t = make_treatment_pa(
      otu_num,
      sample_groups,
      mode = mode,
      min_replicates = min_replicates,
      min_prop = min_prop
    )
    
    summarize_intersections(pa_t) %>%
      mutate(
        Analysis = "Occupancy threshold",
        Threshold_type = Threshold_type
      )
  }
) %>%
  left_join(observed_metrics, by = "Metric") %>%
  mutate(`Relative to original` = if_else(Observed_value == 0, NA_real_, Value / Observed_value))

############################################################
## 7. README and export

readme = tibble(
  Item = c(
    "Purpose",
    "Original treatment-level presence definition",
    "Common-n equalized replicate analysis",
    "Main narrative common-n",
    "Additional common-n settings",
    "Subsampling iterations",
    "Occupancy threshold analysis",
    "Framework-support indicators",
    "Input data"
  ),
  Description = c(
    "Sensitivity analysis of treatment-level OTU intersections underlying Fig. 2E.",
    "An OTU was considered present in a treatment if detected in at least one replicate within that treatment.",
    paste0("Replicate numbers were equalized by random subsampling without replacement within every dilution treatment. Common-n values were ", paste(n_grid, collapse = ", "), "."),
    "n = 14 was used as the primary common-n sensitivity because it is the largest replicate number shared by all dilution treatments.",
    "n = 8, 10 and 12 were used to evaluate whether the n = 14 trend depended on a single subsampling depth.",
    paste0(n_iter, " iterations for each common-n setting."),
    "Treatment-level presence was recalculated using alternative thresholds: >=1 replicate, >=2 replicates, >=5% replicates and >=10% replicates.",
    "Framework-support indicators summarize whether D1-specific OTUs, all-treatment shared OTUs and major multi-dilution shared fractions remain detectable after common-n subsampling.",
    otu_file
  )
)

analysis_parameters = tibble(
  Parameter = c("treatment_order", "n_grid", "n_iter", "presence_mode", "min_replicates", "random_seed"),
  Value = c(paste(treatment_order, collapse = ", "), paste(n_grid, collapse = ", "), n_iter, main_presence_mode, main_min_replicates, 123)
)

wb = createWorkbook()

addWorksheet(wb, "README")
writeData(wb, "README", readme)

addWorksheet(wb, "Sample_numbers")
writeData(wb, "Sample_numbers", sample_n)

addWorksheet(wb, "Analysis_parameters")
writeData(wb, "Analysis_parameters", analysis_parameters)

addWorksheet(wb, "Observed_full_data")
writeData(wb, "Observed_full_data", observed_metrics)

addWorksheet(wb, "Equalized_summary_all_n")
writeData(wb, "Equalized_summary_all_n", equalized_summary)

addWorksheet(wb, "Trend_vs_n14")
writeData(wb, "Trend_vs_n14", trend_vs_n14)

addWorksheet(wb, "Framework_by_iteration")
writeData(wb, "Framework_by_iteration", framework_by_iteration)

addWorksheet(wb, "Framework_summary")
writeData(wb, "Framework_summary", framework_summary)

addWorksheet(wb, "Occupancy_threshold_summary")
writeData(wb, "Occupancy_threshold_summary", occupancy_summary)

addWorksheet(wb, "Full_equalized_iterations")
writeData(wb, "Full_equalized_iterations", equalized_iterations)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

## Also export CSV files for easier version control
write_csv(observed_metrics, file.path(out_dir, "Observed_full_data.csv"))
write_csv(equalized_summary, file.path(out_dir, "Equalized_summary_all_n.csv"))
write_csv(trend_vs_n14, file.path(out_dir, "Trend_vs_n14.csv"))
write_csv(framework_by_iteration, file.path(out_dir, "Framework_by_iteration.csv"))
write_csv(framework_summary, file.path(out_dir, "Framework_summary.csv"))
write_csv(occupancy_summary, file.path(out_dir, "Occupancy_threshold_summary.csv"))
write_csv(equalized_iterations, file.path(out_dir, "Full_equalized_iterations.csv"))

message("Done. Results written to: ", out_xlsx)
message("CSV outputs written to: ", out_dir)

############################################################
