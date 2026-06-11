# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(tidyverse)
library(ggplot2)
library(scales)
library(patchwork)

############################################################################
# Sensitivity plots

framework_treatment_summary = read.csv("Summary_framework_by_treatment.csv")

n_grid = c(8, 10, 12, 14)
treat_levels = c("D0", "D1", "D3", "D5", "D6", "D7")

# dilution_numeric = c("D0" = 0, "D1" = 1, "D3" = 3, "D5" = 5, "D6" = 6, "D7" = 7)
# diluted_treatments = c("D1", "D3", "D5", "D6", "D7")
# strong_treatments = c("D3", "D5", "D6", "D7")

dominance_cutoffs = c("0.01%" = 0.0001,
                      "0.1%"  = 0.001,
                      "1%"    = 0.01)
main_cutoff_label = "0.1%"

theme_pub = theme_bw() +
  theme(panel.background = element_blank(),
        plot.background = element_blank(), 
        plot.title = element_text(face = "bold", size = 10, hjust = 0),
        axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(color = "black", size = 8),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))

isme_cols = c(
  "8"  = "#3A7D6B",   # muted teal-green
  "10" = "#D4943A",   # warm orange
  "12" = "#6E84B7",   # muted blue
  "14" = "#B85C74"    # muted rose  "#C05A78"
)

############################################################################
############################################################################
## treatment-level metrics under the original 0.1% cutoff.

plot_treat_main = framework_treatment_summary %>%
  filter(Cutoff_label == main_cutoff_label) %>%
  mutate(n_subsample = factor(n_subsample, levels = n_grid), 
         Treatment = factor(Treatment, levels = treat_levels))  

p_retgap = ggplot(plot_treat_main,
                  aes(x = Treatment,
                      y = retention_gap_Median,
                      ymin = retention_gap_Q2.5, 
                      ymax = retention_gap_Q97.5,
                      color = n_subsample,
                      group = n_subsample)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.45), linewidth = 0.35, size = 0.3, shape = 16, alpha = 0.8) + 
  scale_color_manual(values = isme_cols, breaks = names(isme_cols) ) + 
  labs(title = "i. Retention gap: abundant - rare",
       x = NULL,
       y = "Difference in\nretained proportion",
       color = "Common n") +
  theme_pub

p_rareTL = ggplot(plot_treat_main,
                  aes(x = Treatment,
                      y = rare_transient_lost_Median,
                      ymin = rare_transient_lost_Q2.5,
                      ymax = rare_transient_lost_Q97.5,
                      color = n_subsample,
                      group = n_subsample)) +
  geom_pointrange(position = position_dodge(width = 0.45), linewidth = 0.35, size = 0.3, shape = 16, alpha = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_color_manual(values = isme_cols, breaks = names(isme_cols) ) + 
  labs(title = "ii. Rare OTUs in transient/lost states",
       x = NULL,
       y = "Proportion of\nrare OTUs",
       color = "Common n") +
  theme_pub

p_shift_abund = ggplot(plot_treat_main %>% filter(Treatment != "D0"),
                       aes(x = Treatment,
                           y = abundant_shift_Median,
                           ymin = abundant_shift_Q2.5,
                           ymax = abundant_shift_Q97.5,
                           color = n_subsample,
                           group = n_subsample)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.45), linewidth = 0.35, size = 0.3, shape = 16, alpha = 0.8) +
  scale_color_manual(values = isme_cols, breaks = names(isme_cols) ) + 
  labs(title = "iii. Abundant OTU shift among retained OTUs",
       x = "Dilution treatment",
       y = expression(Median~Delta~log[10]~RA~vs~D0),
       color = "Common n") +
  theme_pub

fig_commonN_main = (p_retgap / p_rareTL / p_shift_abund) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("Sup_Fig2D_commonN_cutoff0.1.pdf", fig_commonN_main, width = 4, height = 5)

############################################################################
############################################################################
## threshold × n framework summaries.

framework_contrast_summary = read.csv("Summary_framework_contrast.csv")

selected_framework_metrics = c(
  "Mean retention gap: abundant minus rare, D1-D7",
  "Mean retention gap: abundant minus rare, D3-D7",
  "Mean rare-minus-abundant transient/lost gap, D3-D7",
  "Spearman rho: rare retained proportion vs dilution rank",
  "Mean abundant median log10 shift vs D0, D1-D7"
)

plot_framework = framework_contrast_summary %>%
  filter(Metric %in% selected_framework_metrics) %>%
  mutate(Cutoff_label = factor(Cutoff_label, levels = names(dominance_cutoffs)),
         n_subsample = factor(n_subsample, levels = n_grid),
         Metric = factor(Metric, levels = selected_framework_metrics))

cutoff_cols = c(
  "0.01%" = "#5B5B5B",   # charcoal grey
  "0.1%"  = "#7A5C99",   # muted purple; main cutoff
  "1%"    = "#8C7A3A"    # muted olive-brown
)

p_framework = ggplot(plot_framework,
                     aes(x = n_subsample,
                         y = Median,
                         ymin = Q2.5, 
                         ymax = Q97.5,
                         color = Cutoff_label,
                         group = Cutoff_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey65") +
  geom_pointrange(position = position_dodge(width = 0.45), linewidth = 0.35, size = 0.5, shape = 17, alpha = 0.8) +
  scale_color_manual(values = cutoff_cols, breaks = names(cutoff_cols) ) + 
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  labs( #title = "Threshold × common-n sensitivity of Fig. 2D framework metrics",
       x = "Common replicate number per treatment",
       y = "Median across subsampling iterations",
       color = "D0 cutoff") + 
  theme_pub +
  theme(axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 8, face = "bold"))

ggsave("Sup_Fig2D_threshold_commonN.pdf", p_framework, width = 4, height = 5)

############################################################################


