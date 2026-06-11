# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

############################################################################
## Fig. 2D sensitivity analysis
## Combined threshold sensitivity and replicate-number equalization
##
## Purpose:
##   1) Test whether rare/abundant OTU fate patterns are sensitive to the
##      D0 abundance threshold used to define abundant vs rare OTUs.
##   2) Test whether occupancy/retention patterns are driven by unequal
##      replicate numbers among dilution treatments.
##
## Main narrative setting:
##   dominance cutoff = 0.1% relative abundance, common n = 14 replicates
##   per dilution treatment.
##
## Additional sensitivity settings:
##   dominance cutoffs = 0.01%, 0.1%, 1%
##   common n = 8, 10, 12, 14 for all treatments: D0, D1, D3, D5, D6, D7
############################################################################

library(tidyverse)
library(scales)
library(patchwork)

############################################################################
## 1. Parameters

set.seed(123)

input_file = "RelativeAbundance.txt"
outdir = "./"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Original Fig. 2D setting: 0.1% = 0.001 if relative abundance is expressed as a proportion.
## Additional thresholds follow the review suggestion: 0.01%, 0.1%, 1%.
dominance_cutoffs = c("0.01%" = 0.0001,
                      "0.1%"  = 0.001,
                      "1%"    = 0.01)
main_cutoff_label = "0.1%"

## All dilution treatments must be covered at every common n.
n_grid = c(8, 10, 12, 14)
main_n = 14
n_iter = 1000

## Detection is defined as relative abundance greater than detect_cutoff.
detect_cutoff = 0

## Rare/abundant classification basis:
##   "subsampled_D0": recalculate D0-defined abundant/rare OTUs from the sampled D0 replicates in every iteration.
##                    This is the most conservative choice because D0 classification uncertainty is included.
##   "full_D0":       keep D0-defined abundant/rare OTUs from all D0 replicates in every iteration.
##                    This isolates treatment-level replicate equalization while holding the D0 definition fixed.
classification_basis = "subsampled_D0"

## Treatment order
treat_levels = c("D0", "D1", "D3", "D5", "D6", "D7")
dilution_numeric = c("D0" = 0, "D1" = 1, "D3" = 3, "D5" = 5, "D6" = 6, "D7" = 7)
diluted_treatments = c("D1", "D3", "D5", "D6", "D7")
strong_treatments = c("D3", "D5", "D6", "D7")

occ_levels = c("Lost (0)",
               "Transient (0-0.25]",
               "Intermediate (0.25-0.50]",
               "High (0.50-<1.00)",
               "Persistent (1.00)")

group_levels = c("Abundant OTUs", "Rare OTUs")
group_cols = c("Abundant OTUs" = "#0072B2",
               "Rare OTUs" = "#D55E00")

heat_cols = c("#f7fcf0", "#ccebc5", "#7bccc4", "#2b8cbe", "#084081")

theme_pub = theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0),
        axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(color = "black", size = 8),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold", size = 8),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))

############################################################################
## 2. Read and prepare data

ra_raw = read.table(input_file,
                    header = TRUE,
                    sep = "\t",
                    stringsAsFactors = FALSE,
                    check.names = FALSE)

if (!"OTUID" %in% colnames(ra_raw)) {
  stop("The first column should be named 'OTUID'.")
}

sample_meta = tibble(SampleID = colnames(ra_raw)[colnames(ra_raw) != "OTUID"]) %>%
  mutate(Treatment = stringr::str_extract(SampleID, "^D[0-9]")) %>%
  filter(!is.na(Treatment), Treatment %in% treat_levels) %>%
  mutate(Treatment = factor(Treatment, levels = treat_levels)) %>%
  arrange(Treatment, SampleID)

sample_groups = split(sample_meta$SampleID, sample_meta$Treatment)
sample_groups = sample_groups[treat_levels]

sample_number_table = tibble(Treatment = names(sample_groups),
                             Original_n = sapply(sample_groups, length))

print(sample_number_table)

if (any(sample_number_table$Original_n < max(n_grid))) {
  stop("At least one treatment has fewer replicates than max(n_grid). Please reduce n_grid.")
}

ra_mat = ra_raw %>%
  select(all_of(sample_meta$SampleID)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()
rownames(ra_mat) = ra_raw$OTUID

## Global pseudo-count used for log10 fold-change calculations.
## This keeps the pseudo-count identical across thresholds and subsampling iterations.
global_pseudo = min(ra_mat[ra_mat > detect_cutoff], na.rm = TRUE) / 2

############################################################################
## 3. Helper functions

mean_or_na = function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

median_or_na = function(x) {
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

q_or_na = function(x, prob) {
  if (all(is.na(x))) return(NA_real_)
  as.numeric(quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
}

all_true_or_na = function(x) {
  if (any(is.na(x))) return(NA_real_)
  as.numeric(all(x))
}

calc_median_present = function(x, detect_cutoff = 0) {
  ## x is a matrix with OTUs as rows and samples as columns
  apply(x, 1, function(v) {
    v = v[v > detect_cutoff]
    if (length(v) == 0) {
      0
    } else {
      median(v, na.rm = TRUE)
    }
  })
}

assign_occ_class = function(occupancy) {
  case_when(
    occupancy == 0 ~ "Lost (0)",
    occupancy > 0 & occupancy <= 0.25 ~ "Transient (0-0.25]",
    occupancy > 0.25 & occupancy <= 0.50 ~ "Intermediate (0.25-0.50]",
    occupancy > 0.50 & occupancy < 1 ~ "High (0.50-<1.00)",
    occupancy == 1 ~ "Persistent (1.00)",
    TRUE ~ NA_character_
  )
}

make_otu_class = function(class_samples, dominance_cutoff, cutoff_label) {
  d0_vals = ra_mat[, class_samples, drop = FALSE]
  mean_RA_D0 = rowMeans(d0_vals, na.rm = TRUE)
  occ_D0 = rowMeans(d0_vals > detect_cutoff, na.rm = TRUE)

  tibble(OTUID = rownames(ra_mat),
         mean_RA_D0 = mean_RA_D0,
         occ_D0 = occ_D0,
         Cutoff_label = cutoff_label,
         dominance_cutoff = dominance_cutoff,
         OTU_group = case_when(
           mean_RA_D0 >= dominance_cutoff ~ "Abundant OTUs",
           mean_RA_D0 > detect_cutoff & mean_RA_D0 < dominance_cutoff ~ "Rare OTUs",
           TRUE ~ "Absent in D0")) %>%
    mutate(OTU_group = factor(OTU_group,
                              levels = c("Abundant OTUs", "Rare OTUs", "Absent in D0")))
}

analyze_otu_fate = function(sampled_groups,
                            dominance_cutoff,
                            cutoff_label,
                            analysis,
                            n_subsample = NA_integer_,
                            iteration = NA_integer_) {

  if (classification_basis == "full_D0") {
    class_samples = sample_groups[["D0"]]
  } else if (classification_basis == "subsampled_D0") {
    class_samples = sampled_groups[["D0"]]
  } else {
    stop("classification_basis must be either 'subsampled_D0' or 'full_D0'.")
  }

  otu_class = make_otu_class(class_samples = class_samples,
                             dominance_cutoff = dominance_cutoff,
                             cutoff_label = cutoff_label)

  otu_class_count = otu_class %>%
    count(OTU_group, name = "n_OTUs") %>%
    complete(OTU_group = factor(c("Abundant OTUs", "Rare OTUs", "Absent in D0"),
                                levels = c("Abundant OTUs", "Rare OTUs", "Absent in D0")),
             fill = list(n_OTUs = 0)) %>%
    mutate(Analysis = analysis,
           n_subsample = n_subsample,
           Iteration = iteration,
           Cutoff_label = cutoff_label,
           dominance_cutoff = dominance_cutoff,
           .before = 1)

  focus_otus = otu_class %>%
    filter(OTU_group %in% group_levels) %>%
    pull(OTUID)

  if (length(focus_otus) == 0) {
    empty = tibble()
    return(list(otu_class_count = otu_class_count,
                otu_treat = empty,
                retention = empty,
                occupancy_class = empty,
                shift = empty,
                framework_treatment = empty,
                framework_contrast = empty))
  }

  ## OTU-level treatment summaries
  otu_treat = map_dfr(treat_levels, function(g) {
    vals = ra_mat[focus_otus, sampled_groups[[g]], drop = FALSE]
    detected = vals > detect_cutoff

    tibble(OTUID = focus_otus,
           Treatment = g,
           n_replicates = ncol(vals),
           occupancy = rowMeans(detected, na.rm = TRUE),
           retained = rowSums(detected, na.rm = TRUE) > 0,
           mean_RA = rowMeans(vals, na.rm = TRUE),
           median_RA_present = calc_median_present(vals, detect_cutoff = detect_cutoff))
  }) %>%
    left_join(otu_class %>%
                select(OTUID, OTU_group, mean_RA_D0, occ_D0),
              by = "OTUID") %>%
    mutate(OTU_group = factor(as.character(OTU_group), levels = group_levels),
           Treatment = factor(Treatment, levels = treat_levels))

  d0_base = otu_treat %>%
    filter(Treatment == "D0") %>%
    select(OTUID,
           median_RA_D0 = median_RA_present,
           occupancy_D0 = occupancy)

  otu_treat = otu_treat %>%
    left_join(d0_base, by = "OTUID") %>%
    mutate(log10_shift_vs_D0 = log10(median_RA_present + global_pseudo) -
             log10(median_RA_D0 + global_pseudo),
           Analysis = analysis,
           n_subsample = n_subsample,
           Iteration = iteration,
           Cutoff_label = cutoff_label,
           dominance_cutoff = dominance_cutoff,
           .before = 1)

  ## Panel D-i style: retained OTU proportion
  retention_df = otu_treat %>%
    group_by(Analysis, n_subsample, Iteration, Cutoff_label, dominance_cutoff,
             OTU_group, Treatment) %>%
    summarise(retained_OTUs = sum(retained, na.rm = TRUE),
              total_OTUs = n_distinct(OTUID),
              retained_prop = retained_OTUs / total_OTUs,
              mean_occupancy = mean(occupancy, na.rm = TRUE),
              .groups = "drop") %>%
    complete(Analysis, n_subsample, Iteration, Cutoff_label, dominance_cutoff,
             OTU_group = factor(group_levels, levels = group_levels),
             Treatment = factor(treat_levels, levels = treat_levels),
             fill = list(retained_OTUs = 0,
                         total_OTUs = 0,
                         retained_prop = NA_real_,
                         mean_occupancy = NA_real_))

  ## Panel D-ii style: occupancy-class proportions
  occ_count = otu_treat %>%
    mutate(Occupancy_class = assign_occ_class(occupancy),
           Occupancy_class = factor(Occupancy_class, levels = occ_levels)) %>%
    count(OTU_group, Treatment, Occupancy_class, name = "n")

  occ_total = otu_treat %>%
    count(OTU_group, Treatment, name = "total")

  occupancy_df = expand_grid(OTU_group = factor(group_levels, levels = group_levels),
                             Treatment = factor(treat_levels, levels = treat_levels),
                             Occupancy_class = factor(occ_levels, levels = occ_levels)) %>%
    left_join(occ_count, by = c("OTU_group", "Treatment", "Occupancy_class")) %>%
    left_join(occ_total, by = c("OTU_group", "Treatment")) %>%
    mutate(n = replace_na(n, 0),
           total = replace_na(total, 0),
           prop = if_else(total > 0, n / total, NA_real_),
           Analysis = analysis,
           n_subsample = n_subsample,
           Iteration = iteration,
           Cutoff_label = cutoff_label,
           dominance_cutoff = dominance_cutoff,
           .before = 1)

  ## Panel D-iii style: shift among retained OTUs
  shift_df = otu_treat %>%
    filter(Treatment != "D0", retained == TRUE, median_RA_D0 > 0) %>%
    group_by(Analysis, n_subsample, Iteration, Cutoff_label, dominance_cutoff,
             OTU_group, Treatment) %>%
    summarise(n_retained = n(),
              median_log10_shift = median(log10_shift_vs_D0, na.rm = TRUE),
              mean_log10_shift = mean(log10_shift_vs_D0, na.rm = TRUE),
              q25_log10_shift = quantile(log10_shift_vs_D0, 0.25, na.rm = TRUE),
              q75_log10_shift = quantile(log10_shift_vs_D0, 0.75, na.rm = TRUE),
              prop_increase = mean(log10_shift_vs_D0 > 0, na.rm = TRUE),
              prop_decrease = mean(log10_shift_vs_D0 < 0, na.rm = TRUE),
              .groups = "drop") %>%
    complete(Analysis, n_subsample, Iteration, Cutoff_label, dominance_cutoff,
             OTU_group = factor(group_levels, levels = group_levels),
             Treatment = factor(treat_levels[treat_levels != "D0"],
                                levels = treat_levels),
             fill = list(n_retained = 0,
                         median_log10_shift = NA_real_,
                         mean_log10_shift = NA_real_,
                         q25_log10_shift = NA_real_,
                         q75_log10_shift = NA_real_,
                         prop_increase = NA_real_,
                         prop_decrease = NA_real_))

  ## Framework-level metrics by treatment
  ret_wide = retention_df %>%
    mutate(Group_short = if_else(OTU_group == "Abundant OTUs", "Abundant", "Rare")) %>%
    select(Analysis, n_subsample, Iteration, Cutoff_label, dominance_cutoff,
           Treatment, Group_short, retained_prop, mean_occupancy, total_OTUs) %>%
    pivot_wider(names_from = Group_short,
                values_from = c(retained_prop, mean_occupancy, total_OTUs),
                names_sep = "__")

  occ_bin = occupancy_df %>%
    mutate(Group_short = if_else(OTU_group == "Abundant OTUs", "Abundant", "Rare"),
           Occ_bin = case_when(
             Occupancy_class %in% c("Lost (0)", "Transient (0-0.25]") ~ "TransientLost",
             Occupancy_class %in% c("High (0.50-<1.00)", "Persistent (1.00)") ~ "HighPersistent",
             Occupancy_class == "Intermediate (0.25-0.50]" ~ "Intermediate",
             TRUE ~ NA_character_)) %>%
    group_by(Analysis, n_subsample, Iteration, Cutoff_label, dominance_cutoff,
             Treatment, Group_short, Occ_bin) %>%
    summarise(prop = sum(prop, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = c(Group_short, Occ_bin),
                values_from = prop,
                names_prefix = "occprop__")

  shift_wide = shift_df %>%
    mutate(Group_short = if_else(OTU_group == "Abundant OTUs", "Abundant", "Rare")) %>%
    select(Analysis, n_subsample, Iteration, Cutoff_label, dominance_cutoff,
           Treatment, Group_short, median_log10_shift, prop_increase, prop_decrease) %>%
    pivot_wider(names_from = Group_short,
                values_from = c(median_log10_shift, prop_increase, prop_decrease),
                names_sep = "__")

  framework_treatment = ret_wide %>%
    left_join(occ_bin,
              by = c("Analysis", "n_subsample", "Iteration", "Cutoff_label",
                     "dominance_cutoff", "Treatment")) %>%
    left_join(shift_wide,
              by = c("Analysis", "n_subsample", "Iteration", "Cutoff_label",
                     "dominance_cutoff", "Treatment")) %>%
    mutate(Treatment = factor(Treatment, levels = treat_levels),
           dilution_rank = unname(dilution_numeric[as.character(Treatment)]),
           retention_gap_abundant_minus_rare = retained_prop__Abundant - retained_prop__Rare,
           transient_lost_gap_rare_minus_abundant = occprop__Rare_TransientLost - occprop__Abundant_TransientLost,
           high_persistent_gap_abundant_minus_rare = occprop__Abundant_HighPersistent - occprop__Rare_HighPersistent)

  ## Framework-level contrasts across treatments
  ft = framework_treatment

  rare_ret = ft %>% filter(Treatment %in% treat_levels) %>% pull(retained_prop__Rare)
  abundant_ret = ft %>% filter(Treatment %in% treat_levels) %>% pull(retained_prop__Abundant)
  dilution_vals = ft %>% filter(Treatment %in% treat_levels) %>% pull(dilution_rank)

  rare_ret_spearman = suppressWarnings(cor(dilution_vals, rare_ret, method = "spearman", use = "complete.obs"))
  abundant_ret_spearman = suppressWarnings(cor(dilution_vals, abundant_ret, method = "spearman", use = "complete.obs"))

  framework_contrast = tibble(
    Analysis = analysis,
    n_subsample = n_subsample,
    Iteration = iteration,
    Cutoff_label = cutoff_label,
    dominance_cutoff = dominance_cutoff,
    Metric = c(
      "Mean retention gap: abundant minus rare, D1-D7",
      "Mean retention gap: abundant minus rare, D3-D7",
      "Retention gap > 0 in all diluted treatments",
      "Retention gap > 0 in all stronger dilution treatments",
      "Mean rare transient/lost proportion, D1-D7",
      "Mean rare transient/lost proportion, D3-D7",
      "Mean rare-minus-abundant transient/lost gap, D1-D7",
      "Mean rare-minus-abundant transient/lost gap, D3-D7",
      "Rare transient/lost > abundant transient/lost in all diluted treatments",
      "Rare transient/lost > abundant transient/lost in all stronger dilution treatments",
      "Spearman rho: rare retained proportion vs dilution rank",
      "Spearman rho: abundant retained proportion vs dilution rank",
      "Mean abundant median log10 shift vs D0, D1-D7",
      "Mean rare median log10 shift vs D0, D1-D7",
      "Mean abundant decrease proportion, D1-D7",
      "Mean rare decrease proportion, D1-D7"),
    Value = c(
      mean_or_na(ft %>% filter(Treatment %in% diluted_treatments) %>% pull(retention_gap_abundant_minus_rare)),
      mean_or_na(ft %>% filter(Treatment %in% strong_treatments) %>% pull(retention_gap_abundant_minus_rare)),
      all_true_or_na((ft %>% filter(Treatment %in% diluted_treatments) %>% pull(retention_gap_abundant_minus_rare)) > 0),
      all_true_or_na((ft %>% filter(Treatment %in% strong_treatments) %>% pull(retention_gap_abundant_minus_rare)) > 0),
      mean_or_na(ft %>% filter(Treatment %in% diluted_treatments) %>% pull(occprop__Rare_TransientLost)),
      mean_or_na(ft %>% filter(Treatment %in% strong_treatments) %>% pull(occprop__Rare_TransientLost)),
      mean_or_na(ft %>% filter(Treatment %in% diluted_treatments) %>% pull(transient_lost_gap_rare_minus_abundant)),
      mean_or_na(ft %>% filter(Treatment %in% strong_treatments) %>% pull(transient_lost_gap_rare_minus_abundant)),
      all_true_or_na((ft %>% filter(Treatment %in% diluted_treatments) %>% pull(transient_lost_gap_rare_minus_abundant)) > 0),
      all_true_or_na((ft %>% filter(Treatment %in% strong_treatments) %>% pull(transient_lost_gap_rare_minus_abundant)) > 0),
      rare_ret_spearman,
      abundant_ret_spearman,
      mean_or_na(ft %>% filter(Treatment %in% diluted_treatments) %>% pull(median_log10_shift__Abundant)),
      mean_or_na(ft %>% filter(Treatment %in% diluted_treatments) %>% pull(median_log10_shift__Rare)),
      mean_or_na(ft %>% filter(Treatment %in% diluted_treatments) %>% pull(prop_decrease__Abundant)),
      mean_or_na(ft %>% filter(Treatment %in% diluted_treatments) %>% pull(prop_decrease__Rare))
    )
  )

  list(otu_class_count = otu_class_count,
       otu_treat = otu_treat,
       retention = retention_df,
       occupancy_class = occupancy_df,
       shift = shift_df,
       framework_treatment = framework_treatment,
       framework_contrast = framework_contrast)
}

summarise_value = function(df, group_cols, value_col = "Value") {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(Mean = mean(.data[[value_col]], na.rm = TRUE),
              SD = sd(.data[[value_col]], na.rm = TRUE),
              Median = median(.data[[value_col]], na.rm = TRUE),
              Q2.5 = q_or_na(.data[[value_col]], 0.025),
              Q97.5 = q_or_na(.data[[value_col]], 0.975),
              N_iterations = n(),
              .groups = "drop")
}

############################################################################
## 4. Full-data reference analyses for all abundance thresholds

message("Running full-data reference analyses...")

full_results = imap(dominance_cutoffs, function(cutoff, cutoff_label) {
  analyze_otu_fate(sampled_groups = sample_groups,
                   dominance_cutoff = cutoff,
                   cutoff_label = cutoff_label,
                   analysis = "Full data",
                   n_subsample = NA_integer_,
                   iteration = NA_integer_)
})

full_otu_class_count = bind_rows(lapply(full_results, `[[`, "otu_class_count"))
full_retention = bind_rows(lapply(full_results, `[[`, "retention"))
full_occupancy_class = bind_rows(lapply(full_results, `[[`, "occupancy_class"))
full_shift = bind_rows(lapply(full_results, `[[`, "shift"))
full_framework_treatment = bind_rows(lapply(full_results, `[[`, "framework_treatment"))
full_framework_contrast = bind_rows(lapply(full_results, `[[`, "framework_contrast"))

############################################################################
## 5. Combined threshold × common-n subsampling sensitivity

message("Running threshold × common-n subsampling sensitivity analyses...")

all_retention = list()
all_occupancy_class = list()
all_shift = list()
all_framework_treatment = list()
all_framework_contrast = list()
all_otu_class_count = list()

counter = 1
for (cutoff_label in names(dominance_cutoffs)) {
  cutoff = dominance_cutoffs[[cutoff_label]]

  for (n_sub in n_grid) {
    message("Cutoff = ", cutoff_label, "; common n = ", n_sub)

    for (i in seq_len(n_iter)) {
      if (i %% 50 == 0 || i == 1 || i == n_iter) {
        message("  iteration ", i, " / ", n_iter)
      }

      sampled_groups = lapply(sample_groups, function(cols) {
        sample(cols, size = n_sub, replace = FALSE)
      })
      sampled_groups = sampled_groups[treat_levels]

      res_i = analyze_otu_fate(sampled_groups = sampled_groups,
                               dominance_cutoff = cutoff,
                               cutoff_label = cutoff_label,
                               analysis = "Equalized common-n subsampling",
                               n_subsample = n_sub,
                               iteration = i)

      all_otu_class_count[[counter]] = res_i$otu_class_count
      all_retention[[counter]] = res_i$retention
      all_occupancy_class[[counter]] = res_i$occupancy_class
      all_shift[[counter]] = res_i$shift
      all_framework_treatment[[counter]] = res_i$framework_treatment
      all_framework_contrast[[counter]] = res_i$framework_contrast

      counter = counter + 1
    }
  }
}

iter_otu_class_count = bind_rows(all_otu_class_count)
iter_retention = bind_rows(all_retention)
iter_occupancy_class = bind_rows(all_occupancy_class)
iter_shift = bind_rows(all_shift)
iter_framework_treatment = bind_rows(all_framework_treatment)
iter_framework_contrast = bind_rows(all_framework_contrast)

############################################################################
## 6. Summary tables

message("Summarising sensitivity results...")

otu_class_count_summary = iter_otu_class_count %>%
  group_by(Cutoff_label, dominance_cutoff, n_subsample, OTU_group) %>%
  summarise(Mean = mean(n_OTUs, na.rm = TRUE),
            SD = sd(n_OTUs, na.rm = TRUE),
            Median = median(n_OTUs, na.rm = TRUE),
            Q2.5 = q_or_na(n_OTUs, 0.025),
            Q97.5 = q_or_na(n_OTUs, 0.975),
            N_iterations = n(),
            .groups = "drop")

retention_summary = iter_retention %>%
  group_by(Cutoff_label, dominance_cutoff, n_subsample, OTU_group, Treatment) %>%
  summarise(retained_OTUs_Mean = mean(retained_OTUs, na.rm = TRUE),
            retained_OTUs_Median = median(retained_OTUs, na.rm = TRUE),
            retained_prop_Mean = mean(retained_prop, na.rm = TRUE),
            retained_prop_Median = median(retained_prop, na.rm = TRUE),
            retained_prop_Q2.5 = q_or_na(retained_prop, 0.025),
            retained_prop_Q97.5 = q_or_na(retained_prop, 0.975),
            mean_occupancy_Median = median(mean_occupancy, na.rm = TRUE),
            N_iterations = n(),
            .groups = "drop")

occupancy_class_summary = iter_occupancy_class %>%
  group_by(Cutoff_label, dominance_cutoff, n_subsample, OTU_group, Treatment, Occupancy_class) %>%
  summarise(prop_Mean = mean(prop, na.rm = TRUE),
            prop_Median = median(prop, na.rm = TRUE),
            prop_Q2.5 = q_or_na(prop, 0.025),
            prop_Q97.5 = q_or_na(prop, 0.975),
            n_Median = median(n, na.rm = TRUE),
            N_iterations = n(),
            .groups = "drop")

shift_summary = iter_shift %>%
  group_by(Cutoff_label, dominance_cutoff, n_subsample, OTU_group, Treatment) %>%
  summarise(n_retained_Median = median(n_retained, na.rm = TRUE),
            median_log10_shift_Mean = mean(median_log10_shift, na.rm = TRUE),
            median_log10_shift_Median = median(median_log10_shift, na.rm = TRUE),
            median_log10_shift_Q2.5 = q_or_na(median_log10_shift, 0.025),
            median_log10_shift_Q97.5 = q_or_na(median_log10_shift, 0.975),
            prop_decrease_Median = median(prop_decrease, na.rm = TRUE),
            prop_increase_Median = median(prop_increase, na.rm = TRUE),
            N_iterations = n(),
            .groups = "drop")

framework_treatment_summary = iter_framework_treatment %>%
  group_by(Cutoff_label, dominance_cutoff, n_subsample, Treatment) %>%
  summarise(retention_gap_Median = median(retention_gap_abundant_minus_rare, na.rm = TRUE),
            retention_gap_Q2.5 = q_or_na(retention_gap_abundant_minus_rare, 0.025),
            retention_gap_Q97.5 = q_or_na(retention_gap_abundant_minus_rare, 0.975),
            rare_transient_lost_Median = median(occprop__Rare_TransientLost, na.rm = TRUE),
            rare_transient_lost_Q2.5 = q_or_na(occprop__Rare_TransientLost, 0.025),
            rare_transient_lost_Q97.5 = q_or_na(occprop__Rare_TransientLost, 0.975),
            TL_gap_rare_minus_abundant_Median = median(transient_lost_gap_rare_minus_abundant, na.rm = TRUE),
            TL_gap_rare_minus_abundant_Q2.5 = q_or_na(transient_lost_gap_rare_minus_abundant, 0.025),
            TL_gap_rare_minus_abundant_Q97.5 = q_or_na(transient_lost_gap_rare_minus_abundant, 0.975),
            abundant_shift_Median = median(median_log10_shift__Abundant, na.rm = TRUE),
            abundant_shift_Q2.5 = q_or_na(median_log10_shift__Abundant, 0.025),
            abundant_shift_Q97.5 = q_or_na(median_log10_shift__Abundant, 0.975),
            rare_shift_Median = median(median_log10_shift__Rare, na.rm = TRUE),
            rare_shift_Q2.5 = q_or_na(median_log10_shift__Rare, 0.025),
            rare_shift_Q97.5 = q_or_na(median_log10_shift__Rare, 0.975),
            N_iterations = n(),
            .groups = "drop")

framework_contrast_summary = iter_framework_contrast %>%
  group_by(Cutoff_label, dominance_cutoff, n_subsample, Metric) %>%
  summarise(Mean = mean(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE),
            Median = median(Value, na.rm = TRUE),
            Q2.5 = q_or_na(Value, 0.025),
            Q97.5 = q_or_na(Value, 0.975),
            N_iterations = n(),
            .groups = "drop")

## Compare all settings with the main narrative reference:
## cutoff = 0.1%, common n = 14.
reference_framework = framework_contrast_summary %>%
  filter(Cutoff_label == main_cutoff_label, n_subsample == main_n) %>%
  select(Metric, Reference_Median = Median)

framework_trend_vs_reference = framework_contrast_summary %>%
  left_join(reference_framework, by = "Metric") %>%
  mutate(Difference_from_reference = Median - Reference_Median,
         Ratio_to_reference = if_else(!is.na(Reference_Median) & Reference_Median != 0,
                                      Median / Reference_Median,
                                      NA_real_),
         Same_sign_as_reference = if_else(!is.na(Reference_Median) & !is.na(Median),
                                          sign(Median) == sign(Reference_Median),
                                          NA))

############################################################################
## 7. Write outputs

message("Writing outputs...")

write_csv(sample_number_table, file.path(outdir, "Sample_numbers.csv"))
write_csv(full_otu_class_count, file.path(outdir, "FullData_OTU_class_count.csv"))
write_csv(full_retention, file.path(outdir, "FullData_retention.csv"))
write_csv(full_occupancy_class, file.path(outdir, "FullData_occupancy_class.csv"))
write_csv(full_shift, file.path(outdir, "FullData_shift.csv"))
write_csv(full_framework_treatment, file.path(outdir, "FullData_framework_by_treatment.csv"))
write_csv(full_framework_contrast, file.path(outdir, "FullData_framework_contrast.csv"))

write_csv(iter_otu_class_count, file.path(outdir, "Equalized_OTU_class_count_iterations.csv"))
write_csv(iter_retention, file.path(outdir, "Equalized_retention_iterations.csv"))
write_csv(iter_occupancy_class, file.path(outdir, "Equalized_occupancy_class_iterations.csv"))
write_csv(iter_shift, file.path(outdir, "Equalized_shift_iterations.csv"))
write_csv(iter_framework_treatment, file.path(outdir, "Equalized_framework_by_treatment_iterations.csv"))
write_csv(iter_framework_contrast, file.path(outdir, "Equalized_framework_contrast_iterations.csv"))

write_csv(otu_class_count_summary, file.path(outdir, "Summary_OTU_class_count.csv"))
write_csv(retention_summary, file.path(outdir, "Summary_retention.csv"))
write_csv(occupancy_class_summary, file.path(outdir, "Summary_occupancy_class.csv"))
write_csv(shift_summary, file.path(outdir, "Summary_shift.csv"))
write_csv(framework_treatment_summary, file.path(outdir, "Summary_framework_by_treatment.csv"))
write_csv(framework_contrast_summary, file.path(outdir, "Summary_framework_contrast.csv"))
write_csv(framework_trend_vs_reference, file.path(outdir, "Summary_trend_vs_reference_cutoff0.1_n14.csv"))

readme = tibble(
  Item = c("Purpose",
           "Input data",
           "Detection rule",
           "Threshold sensitivity",
           "Common-n subsampling",
           "Main narrative reference",
           "Rare/abundant classification basis",
           "Iterations",
           "Key framework comparisons",
           "Interpretation boundary"),
  Description = c(
    "Combined threshold-sensitivity and replicate-equalized sensitivity analysis for Fig. 2D OTU fate patterns.",
    input_file,
    paste0("Detected if relative abundance > ", detect_cutoff, "."),
    paste0("D0-defined abundant/rare cutoffs: ", paste(names(dominance_cutoffs), dominance_cutoffs, sep = "=", collapse = "; "), "."),
    paste0("All dilution treatments were covered at common n = ", paste(n_grid, collapse = ", "), "."),
    paste0("The primary narrative setting is cutoff = ", main_cutoff_label, " and common n = ", main_n, "."),
    classification_basis,
    as.character(n_iter),
    "Retention gap between abundant and rare OTUs; rare transient/lost fraction; rare-minus-abundant transient/lost gap; median log10 abundance shift among retained OTUs.",
    "This analysis tests robustness to threshold choice and unequal replicate number. It does not model sequencing detection probability or contamination explicitly."
  )
)

analysis_parameters = tibble(
  Parameter = c("dominance_cutoffs", "main_cutoff_label", "n_grid", "main_n", "n_iter",
                "detect_cutoff", "classification_basis", "global_pseudo", "treat_levels"),
  Value = c(paste(paste(names(dominance_cutoffs), dominance_cutoffs, sep = "="), collapse = "; "),
            main_cutoff_label,
            paste(n_grid, collapse = ", "),
            as.character(main_n),
            as.character(n_iter),
            as.character(detect_cutoff),
            classification_basis,
            as.character(global_pseudo),
            paste(treat_levels, collapse = ", "))
)

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb = openxlsx::createWorkbook()

  add_sheet = function(wb, sheet, dat) {
    sheet = substr(sheet, 1, 31)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, dat)
  }

  add_sheet(wb, "README", readme)
  add_sheet(wb, "Analysis_parameters", analysis_parameters)
  add_sheet(wb, "Sample_numbers", sample_number_table)
  add_sheet(wb, "Full_framework_contrast", full_framework_contrast)
  add_sheet(wb, "Full_framework_treatment", full_framework_treatment)
  add_sheet(wb, "Summary_class_count", otu_class_count_summary)
  add_sheet(wb, "Summary_retention", retention_summary)
  add_sheet(wb, "Summary_occupancy", occupancy_class_summary)
  add_sheet(wb, "Summary_shift", shift_summary)
  add_sheet(wb, "Summary_framework_treat", framework_treatment_summary)
  add_sheet(wb, "Summary_framework_contrast", framework_contrast_summary)
  add_sheet(wb, "Trend_vs_ref_0.1_n14", framework_trend_vs_reference)

  openxlsx::saveWorkbook(wb,
                         file.path(outdir, "Table_S5_Fig2D_OTUfate_threshold_equalized_sensitivity.xlsx"),
                         overwrite = TRUE)
}

