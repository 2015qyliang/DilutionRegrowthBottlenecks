# coding: utf-8
# qfsfdxlqy@163.com
# 
# Fig. 3A OTU distribution-state common-n sensitivity analysis
# Purpose: lightweight replicate-equalized sensitivity for OTU distribution-state diagnostics


library(tidyverse)  
library(ggplot2)
library(ggsci)

summary_quantiles = read.csv('commonN_outputs/Fig3A_distribution_commonN_summary_quantiles.csv')

plot_metrics = c("Fraction_nonuniform_P", "Fraction_nonuniform_Q", "Fraction_nonnormal_P",
                 "Fraction_nonnormal_Q", "Median_IOD", "Fraction_low_occupancy" )

metric_labs = c(
  Fraction_nonuniform_P   = "Non-uniform OTUs (P)",
  Fraction_nonuniform_Q   = "Non-uniform OTUs (FDR q)",
  Fraction_nonnormal_P    = "Non-normal OTUs (P)",
  Fraction_nonnormal_Q    = "Non-normal OTUs (FDR q)",
  Median_IOD              = "Median index of dispersion",
  Fraction_low_occupancy  = "Low-occupancy OTUs"
)

treat_levels = c("D0", "D1", "D3", "D5", "D6", "D7")

isme_cols = c(
  "8"  = "#3A7D6B",   # muted teal-green
  "10" = "#D4943A",   # warm orange
  "12" = "#6E84B7",   # muted blue
  "14" = "#B85C74"    # muted rose  "#C05A78"
)

plot_df = summary_quantiles %>%
  filter(Metric %in% plot_metrics) %>%
  mutate(Metric = factor(Metric, levels = plot_metrics),
         Group = factor(Group, levels = treat_levels)   )


p_commonN = ggplot(plot_df, aes(x = Group, y = Median, color = factor(Common_n))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5),
                width = 0.15, alpha = 0.6, linewidth = 0.35,
                position = position_dodge(width = 0.45)) +
  geom_point(size = 1, shape = 16, alpha = 0.8, position = position_dodge(width = 0.45)) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2, labeller = as_labeller(metric_labs)) +
  scale_color_manual(values = isme_cols, breaks = names(isme_cols) ) + 
  labs(x = "Dilution treatments", y = "Median across subsampling iterations", color = "Common n") +
  theme_bw() +
  theme(panel.background = element_blank(), 
        plot.background = element_blank(),
        panel.grid = element_blank(), 
        strip.background = element_blank(), 
        axis.title = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"),
        legend.position = "bottom") 


ggsave("Sup_FigS2_distribution_commonN_sensitivity.pdf", p_commonN, width = 5.5, height = 5)





