# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(tidyverse)
library(reshape2)
library(microeco)
library(GUniFrac)
library(ape)
library(ggplot2)
library(ggsci)
library(gridExtra)

library(openxlsx)

##############################################################################

sample_info = read.table('OTUsamplesGroup.txt', header = T, 
                         sep = '\t', stringsAsFactors = F)
rownames(sample_info) = sample_info$SampleID
otu_table = read.table('OTUtable.txt', header = T, sep = '\t',
                       row.names = 1, stringsAsFactors = F)
taxonomy_table = read.table('OTUTax.txt',  header = T, sep = '\t', 
                            row.names = 1, stringsAsFactors = F)
phylo_tree = read.tree('OTUFastTree.txt')
dataset = microtable$new(sample_table = sample_info, 
                         otu_table = otu_table, 
                         tax_table = taxonomy_table, 
                         phylo_tree = phylo_tree)
dataset$cal_alphadiv(PD = T)

##############################################################################

DFalpha = cbind(dataset$alpha_diversity, sample_info)
DFalpha.melt = melt(DFalpha, measure.vars = colnames(DFalpha)[1:(ncol(DFalpha)-2)],
                    variable.name = 'Measure', value.name = 'Value')
DFalpha.melt$Value = round(DFalpha.melt$Value, digits = 3)
write.table(DFalpha.melt, 'AlphaDF.txt', 
            append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)

##############################################################################

unique(DFalpha.melt$Measure)

alphaDF = DFalpha.melt[DFalpha.melt$Measure %in% c('Observed', 'Simpson', 'Shannon', 'PD'), ]
alphaDF$Measure = factor(alphaDF$Measure, levels = c('Observed', 'Simpson', 'Shannon', 'PD'))

# boxcols = colorRampPalette(pal_npg(alpha = 0.75)(8))(length(unique(alphaDF$Group)))

p = ggplot(alphaDF, aes(x = Group, y = Value, fill = Group)) + 
  geom_boxplot(alpha = 0.9, size = 0.1, 
               outlier.color = 'grey50', outlier.size = 0.3) + 
  geom_jitter(size = 1, alpha = 0.4, width = 0.24, shape = 16, color = 'grey50') + 
  facet_wrap(facets = .~Measure, nrow = 1, scales = 'free_y') +
  # scale_fill_manual(values = rev(boxcols)) +
  scale_fill_manual(values = pal_nejm()(6),
                    breaks = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')) +
  labs(x = "", y = "") + theme_bw()+ 
  theme(panel.grid = element_blank(), 
        plot.background = element_blank(),
        strip.text = element_text(size = 10), 
        strip.background = element_blank(), 
        axis.text.y = element_text(size = 10, hjust = 0.5, vjust = 0.5, 
                                   color = 'black', angle = 90),
        axis.text.x = element_text(size = 10, hjust = 1, vjust = 0.5, 
                                   color = 'black', angle = 90),
        legend.position = "none")

ggsave("AlphaPlot.pdf", p, units = "in", width = 6, height = 2)

##############################################################################

#############
# anova TEST
for (index.alpha in c('Observed', 'Simpson', 'Shannon', 'PD')) {
  df1 = alphaDF[alphaDF$Measure == index.alpha, ]
  posthoc = TukeyHSD(aov(Value ~ Group, data = df1), 'Group', conf.level=0.95)
  write.table(data.frame(Compares = rownames(posthoc$Group), posthoc$Group), 
              paste0('Alpha_ANOVA_TurkeyHSDtest_', index.alpha, '.txt'), 
              append = F, quote = F, sep = '\t', row.names = F, col.names = T)
}

######################################################################################################
######################################################################################################

## 1. Basic settings
## ---------------------------------------------------------

treatment_order = c("D0", "D1", "D3", "D5", "D6", "D7")
alpha_metrics = c("Observed", "Simpson", "Shannon", "PD" )
metric_labels = c(Observed = "Observed OTUs", Simpson  = "Simpson diversity",
                  Shannon  = "Shannon diversity", PD= "Phylogenetic diversity" )

## 2. Check required objects
if (!exists("dataset")) {
  stop("Object 'dataset' was not found. Please run Figure3C_AlphaDiversity.R first.")
}

if (is.null(dataset$alpha_diversity)) {
  stop("dataset$alpha_diversity was not found. Please run dataset$cal_alphadiv(PD = TRUE) first.")
}

if (!exists("sample_info")) {
  stop("Object 'sample_info' was not found. Please run Figure3C_AlphaDiversity.R first.")
}

## 3. Prepare sample information

sample_info_tbl = sample_info %>%
  as.data.frame()

if (!"SampleID" %in% colnames(sample_info_tbl)) {
  sample_info_tbl = sample_info_tbl %>%
    rownames_to_column("SampleID")
}

sample_info_tbl = sample_info_tbl %>%
  as_tibble()

if (!"Group" %in% colnames(sample_info_tbl)) {
  stop("Column 'Group' was not found in sample_info. Please check sample metadata.")
}

## 4. Prepare alpha diversity values

alpha_raw = dataset$alpha_diversity %>%
  as.data.frame()

if (!"SampleID" %in% colnames(alpha_raw)) {
  alpha_raw = alpha_raw %>%
    rownames_to_column("SampleID")
}

alpha_raw = alpha_raw %>%
  as_tibble()

missing_metrics = setdiff(alpha_metrics, colnames(alpha_raw))

if (length(missing_metrics) > 0) {
  stop("The following alpha diversity metrics are missing from dataset$alpha_diversity: ",
       paste(missing_metrics, collapse = ", ")  )
}

alpha_values = alpha_raw %>%
  left_join(sample_info_tbl, by = "SampleID") %>%
  mutate( Group = factor(Group, levels = treatment_order)  ) %>% 
  select(SampleID, Group, all_of(alpha_metrics), everything() )

## Use unrounded values for all statistical tests
alpha_long = alpha_values %>%
  select(SampleID, Group, all_of(alpha_metrics)) %>%
  pivot_longer(cols = all_of(alpha_metrics),
               names_to = "Alpha diversity metric",
               values_to = "Value"  ) %>%
  mutate(`Alpha diversity metric` = factor(`Alpha diversity metric`, levels = alpha_metrics ),
         `Metric label` = recode(as.character(`Alpha diversity metric`), !!!metric_labels  ),
         Group = factor(Group, levels = treatment_order)  )

## 5. Sheet 1: sample-level alpha values

alpha_values_sheet = alpha_values %>%
  mutate(Group = as.character(Group)) %>%
  rename(`Dilution treatment` = Group,
         `Observed OTUs` = Observed,
         `Simpson diversity` = Simpson,
         `Shannon diversity` = Shannon,
         `Phylogenetic diversity` = PD  )

## 6. Sheet 2: summary by treatment

summary_by_treatment = alpha_long %>%
  group_by(`Metric label`, Group) %>%
  summarise(n = n(),
            Mean = mean(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE),
            SE = SD / sqrt(n),
            Median = median(Value, na.rm = TRUE),
            Q1 = quantile(Value, 0.25, na.rm = TRUE),
            Q3 = quantile(Value, 0.75, na.rm = TRUE),
            IQR = IQR(Value, na.rm = TRUE),
            Min = min(Value, na.rm = TRUE),
            Max = max(Value, na.rm = TRUE),
            .groups = "drop"  ) %>%
  mutate(Group = factor(Group, levels = treatment_order),
         `Metric label` = factor(`Metric label`, levels = unname(metric_labels) )  ) %>%
  arrange(`Metric label`, Group) %>%
  rename(`Alpha diversity metric` = `Metric label`, `Dilution treatment` = Group  )


## 7. Sheet 3: global ANOVA results

global_anova = map_dfr(alpha_metrics, function(metric) {
  df_metric = alpha_long %>%
    filter(`Alpha diversity metric` == metric) %>%
    mutate(Group = factor(Group, levels = treatment_order))
  
  fit = aov(Value ~ Group, data = df_metric)
  aov_tab = summary(fit)[[1]]
  
  tibble(`Alpha diversity metric` = metric_labels[[metric]],
         Test = "One-way ANOVA",
         Factor = "Dilution treatment",
         `Df factor` = aov_tab["Group", "Df"],
         `Df residual` = aov_tab["Residuals", "Df"],
         `Sum Sq factor` = aov_tab["Group", "Sum Sq"],
         `Mean Sq factor` = aov_tab["Group", "Mean Sq"],
         `F statistic` = aov_tab["Group", "F value"],
         `P value` = aov_tab["Group", "Pr(>F)"]  ) }  )


## 8. Sheet 4: pairwise Tukey HSD comparisons
## All four alpha indices are summarized in one sheet

pairwise_tukey_all = map_dfr(alpha_metrics, function(metric) {
  
  df_metric = alpha_long %>%
    filter(`Alpha diversity metric` == metric) %>%
    mutate(Group = factor(Group, levels = treatment_order))
  
  fit = aov(Value ~ Group, data = df_metric)
  
  posthoc = TukeyHSD(fit, "Group", conf.level = 0.95 )$Group
  
  posthoc_df = as.data.frame(posthoc) %>%
    rownames_to_column("Raw comparison") %>%
    as_tibble()
  
  posthoc_df %>%
    separate(`Raw comparison`, into = c("Group 1", "Group 2"),
             sep = "-", remove = FALSE ) %>%
    mutate(`Alpha diversity metric` = metric_labels[[metric]],
           Comparison = paste(`Group 1`, "vs", `Group 2`),
           Significance = case_when(`p adj` < 0.001 ~ "***",
                                    `p adj` < 0.01  ~ "**",
                                    `p adj` < 0.05  ~ "*",
                                    TRUE ~ "ns"  )  ) %>%
    transmute(`Alpha diversity metric`,
              `Group 1`,
              `Group 2`,
              Comparison,
              `Mean difference` = diff,
              `Lower 95% CI` = lwr,
              `Upper 95% CI` = upr,
              `Adjusted P value` = `p adj`,
              Significance ) }) %>%
  mutate(`Alpha diversity metric` = factor( `Alpha diversity metric`, levels = unname(metric_labels)  ),
         `Group 1` = factor(`Group 1`, levels = treatment_order),
         `Group 2` = factor(`Group 2`, levels = treatment_order)  ) %>%
  arrange(`Alpha diversity metric`,   `Group 2`,  `Group 1`  )


## 9. Sheet 5: nonparametric global sensitivity
## Kruskal-Wallis tests

kruskal_global = map_dfr(alpha_metrics, function(metric) {
  
  df_metric = alpha_long %>%
    filter(`Alpha diversity metric` == metric) %>%
    mutate(Group = factor(Group, levels = treatment_order))
  
  kt = kruskal.test(Value ~ Group, data = df_metric)
  
  tibble(`Alpha diversity metric` = metric_labels[[metric]],
         Test = "Kruskal-Wallis",
         Statistic = as.numeric(kt$statistic),
         df = as.numeric(kt$parameter),
         `P value` = kt$p.value  ) })

## 10. Sheet 6: pairwise Wilcoxon comparisons
## All four alpha indices summarized in one sheet

pairwise_wilcox_all = map_dfr(alpha_metrics, function(metric) {
  
  df_metric = alpha_long %>%
    filter(`Alpha diversity metric` == metric) %>%
    mutate(Group = factor(Group, levels = treatment_order))
  
  pw = pairwise.wilcox.test(x = df_metric$Value, g = df_metric$Group, p.adjust.method = "BH"   )
  
  as.data.frame(as.table(pw$p.value)) %>%
    filter(!is.na(Freq)) %>%
    as_tibble() %>%
    transmute(`Alpha diversity metric` = metric_labels[[metric]],
              `Group 1` = as.character(Var1),
              `Group 2` = as.character(Var2),
              Comparison = paste(Var1, "vs", Var2),
              `BH-adjusted P value` = Freq,
              Significance = case_when(Freq < 0.001 ~ "***",
                                       Freq < 0.01  ~ "**",
                                       Freq < 0.05  ~ "*",
                                       TRUE ~ "ns" )  ) }) %>%
  mutate(`Alpha diversity metric` = factor(`Alpha diversity metric`, levels = unname(metric_labels)  ),
         `Group 1` = factor(`Group 1`, levels = treatment_order),
         `Group 2` = factor(`Group 2`, levels = treatment_order)  ) %>%
  arrange(`Alpha diversity metric`, `Group 2`, `Group 1`  )

##  Export Excel workbook

wb = createWorkbook()

addWorksheet(wb, "1.Alpha_values")
writeData(wb, "1.Alpha_values", alpha_values_sheet)

addWorksheet(wb, "2.Summary_by_treatment")
writeData(wb, "2.Summary_by_treatment", summary_by_treatment)

addWorksheet(wb, "3.Global_ANOVA")
writeData(wb, "3.Global_ANOVA", global_anova)

addWorksheet(wb, "4.Pairwise_TukeyHSD_all")
writeData(wb, "4.Pairwise_TukeyHSD_all", pairwise_tukey_all)

addWorksheet(wb, "5.Kruskal_Wallis")
writeData(wb, "5.Kruskal_Wallis", kruskal_global)

addWorksheet(wb, "6.Pairwise_Wilcoxon_all")
writeData(wb, "6.Pairwise_Wilcoxon_all", pairwise_wilcox_all)

saveWorkbook(
  wb,
  "Table_Support_Fig3C.xlsx",
  overwrite = TRUE
)









##############################################################################
##############################################################################
##############################################################################

