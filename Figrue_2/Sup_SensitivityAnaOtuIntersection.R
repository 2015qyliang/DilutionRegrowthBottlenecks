############################################################
## Sensitivity analysis of treatment-level OTU intersections

library(tidyverse)
library(readr)
library(openxlsx)

## 1. Read OTU table

otu = read_tsv("OTUtable.txt", show_col_types = FALSE)

## Expected first column: OTUID
if (!"OTUID" %in% colnames(otu)) {
  stop("The first column should be named OTUID.")
}

treatment_order = c("D0", "D1", "D3", "D5", "D6", "D7")

sample_groups = list(
  D0 = grep("^D0", colnames(otu), value = T),
  D1 = grep("^D1", colnames(otu), value = T),
  D3 = grep("^D3", colnames(otu), value = T),
  D5 = grep("^D5", colnames(otu), value = T),
  D6 = grep("^D6", colnames(otu), value = T),
  D7 = grep("^D7", colnames(otu), value = T)  )

## Check sample numbers
print(sapply(sample_groups, length))

otu_num = otu %>%
  mutate(across(-OTUID, as.numeric))

## 2. Function: make treatment-level presence/absence

make_treatment_pa = function(otu_table, sample_groups, mode = "min_replicates",
                             min_replicates = 1, min_prop = NULL) {
  
  out = tibble(OTUID = otu_table$OTUID)
  
  for (g in names(sample_groups)) {
    cols = sample_groups[[g]]
    detected_n = rowSums(otu_table[, cols] > 0, na.rm = T)
    
    if (mode == "min_replicates") {
      present = detected_n >= min_replicates
    }
    
    if (mode == "min_prop") {
      present = detected_n / length(cols) >= min_prop
    }
    
    out[[g]] = present
  }
  
  out = out %>%
    filter(if_any(all_of(names(sample_groups)), ~ .x))
  
  return(out)
}

## 3. Function: summarize key intersections

summarize_intersections = function(pa_table, treatment_cols = treatment_order) {
  
  ## total OTUs per treatment
  total_metrics = map_dfr(treatment_cols, function(g) {
    tibble(
      Metric = paste0("Total OTUs in ", g),
      Value = sum(pa_table[[g]], na.rm = T)
    )
  })
  
  ## membership summary
  membership = pa_table %>%
    rowwise() %>%
    mutate(
      Degree = sum(c_across(all_of(treatment_cols)), na.rm = T),
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
  
  ## key intersections shown or discussed in Fig. 2E
  key_intersections = c("D0 only", "D1 only", "D3 only", "D5 only", "D6 only", "D7 only",
                        "D0&D1", "D0&D1&D3", "D0&D1&D3&D5&D6&D7", "D1&D3&D5&D6&D7",
                        "D3&D5&D6&D7", "D6&D7" )
    
  
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

## 4. Original observed values

pa_original = make_treatment_pa(otu_num,  sample_groups,  
                                mode = "min_replicates",  min_replicates = 1 )

observed_metrics = summarize_intersections(pa_original) %>%
  rename(Observed_value = Value)

## 5. Equalized replicate number analysis

set.seed(123)

target_n = min(sapply(sample_groups, length))
n_iter = 1000

equalized_iterations = map_dfr(seq_len(n_iter), function(i) {
  
  ## 5 iter -> show the process
  if (i %% 5 == 0 || i == n_iter) { 
    message("Running iteration ", i, " / ", n_iter)
  }
  
  sampled_groups = lapply(sample_groups, function(cols) { 
    sample(cols, size = target_n, replace = FALSE)
  })
  
  pa_i = make_treatment_pa(
    otu_num,
    sampled_groups,
    mode = "min_replicates",
    min_replicates = 1
  )
  
  summarize_intersections(pa_i) %>%
    mutate(
      Analysis = "Equalized replicate number",
      Iteration = i,
      `No. replicates per treatment` = target_n
    )
})

equalized_summary = equalized_iterations %>%
  group_by(Analysis, `No. replicates per treatment`, Metric) %>%
  summarise(Mean = mean(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE),
            Median = median(Value, na.rm = TRUE),
            `2.5% quantile` = quantile(Value, 0.025, na.rm = TRUE),
            `97.5% quantile` = quantile(Value, 0.975, na.rm = TRUE),
            .groups = "drop"  ) %>%
  left_join(observed_metrics, by = "Metric") %>%
  mutate(
    `Observed within 95% interval` =
      Observed_value >= `2.5% quantile` &
      Observed_value <= `97.5% quantile`   )


## ---------------------------------------------------------
## 6. Occupancy threshold analysis

threshold_settings = tribble(~Threshold_type, ~mode, ~min_replicates, ~min_prop,
                             ">=1 replicate", "min_replicates", 1, NA_real_,
                             ">=2 replicates", "min_replicates", 2, NA_real_,
                             ">=5% replicates", "min_prop", NA_real_, 0.05,
                             ">=10% replicates", "min_prop", NA_real_, 0.10 )

occupancy_summary = pmap_dfr(
  threshold_settings,
  function(Threshold_type, mode, min_replicates, min_prop) {
    
    pa_t = make_treatment_pa(otu_num, sample_groups, mode = mode,
                             min_replicates = min_replicates, min_prop = min_prop ) 
    
    summarize_intersections(pa_t) %>%
      mutate(
        Analysis = "Occupancy threshold",
        Threshold_type = Threshold_type
      )
  }
) %>%
  left_join(observed_metrics, by = "Metric") %>%
  mutate(`Relative to original` = Value / Observed_value  )


## ---------------------------------------------------------
## 7. README sheet

readme = tibble(
  Item = c("Purpose", "Original treatment-level presence definition", "Equalized replicate analysis",
           "Occupancy threshold analysis",  "Metrics", "Input data"),
  Description = c(
    "Sensitivity analysis of treatment-level OTU intersections underlying Fig. 2E.",
    "An OTU was considered present in a treatment if detected in at least one replicate within that treatment.",
    paste0("Replicate numbers were equalized by randomly subsampling ", target_n,
           " replicates per treatment for ", n_iter, " iterations."), 
    "Treatment-level presence was recalculated using alternative thresholds: >=1 replicate, >=2 replicates, >=5% replicates and >=10% replicates.",
    "Metrics include total OTUs per treatment, all-treatment shared OTUs, treatment-specific OTUs and major intersections shown in Fig. 2E.",
    "Replicate-level OTU abundance table."
  )
)


## ---------------------------------------------------------
## 8. Export to Excel

wb = createWorkbook()

addWorksheet(wb, "README")
writeData(wb, "README", readme)

addWorksheet(wb, "Equalized_replicate_summary")
writeData(wb, "Equalized_replicate_summary", equalized_summary)

addWorksheet(wb, "Occupancy_threshold_summary")
writeData(wb, "Occupancy_threshold_summary", occupancy_summary)

addWorksheet(wb, "Full_equalized_iterations")
writeData(wb, "Full_equalized_iterations", equalized_iterations)

saveWorkbook(
  wb,
  "Table S7. Sensitivity analysis of OTU intersections.xlsx",
  overwrite = TRUE
)
