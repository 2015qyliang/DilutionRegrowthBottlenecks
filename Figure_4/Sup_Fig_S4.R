
library(ggplot2)
library(ggsci)
library(reshape2)
library(ggrepel)
library(doBy)
library(tidyverse)
library(patchwork)

########################################################################################
## Treatment order
treatment_order = c("D0", "D1", "D3", "D5", "D6", "D7")

## Panel C:  Deterministic - Stochastic
dat = read.csv("Figure4C_state_space_summary.csv", check.names = F )
state_dat = dat %>%
  mutate(
    Group = factor(Group, levels = treatment_order),
    ## Deterministic contribution
    Deterministic = Mean_percent_HoS + Mean_percent_HeS,
    ## Stochastic / dispersal-related contribution
    Stochastic = Mean_percent_DL + Mean_percent_HD + Mean_percent_DR,
    ## Difference: positive values mean deterministic > stochastic
    Deterministic_minus_Stochastic = Deterministic - Stochastic,
    ## Approximate propagated SE
    ## This assumes independence among process estimates.
    ## If raw pairwise iCAMP results are available, direct aggregation is preferable.
    SE_Deterministic = sqrt(SE_percent_HoS^2 + SE_percent_HeS^2),
    SE_Stochastic = sqrt(
      SE_percent_DL^2 +
        SE_percent_HD^2 +
        SE_percent_DR^2
    ),
    Regime = ifelse(
      Deterministic > Stochastic,
      "Deterministic > stochastic",
      "Stochastic > deterministic"
    ) ) %>%
  arrange(Dilution_strength)
print(
  state_dat %>%
    select(
      Group,
      Deterministic,
      Stochastic,
      Deterministic_minus_Stochastic,
      Regime ) )

########################################################################################
## Alternative deterministic-stochastic state-space definitions

## i. HoS vs DL + DR
state_dat_s10a = state_dat %>%
  mutate(
    X = Mean_percent_DL + Mean_percent_DR,
    Y = Mean_percent_HoS,
    SE_X = sqrt(SE_percent_DL^2 + SE_percent_DR^2),
    SE_Y = SE_percent_HoS,
    Definition = "i. HoS vs DL + DR"
  )

## ii. HoS + HeS vs DL + HD + DR
state_dat_s10b = state_dat %>%
  mutate(
    X = Stochastic,
    Y = Deterministic,
    SE_X = SE_Stochastic,
    SE_Y = SE_Deterministic,
    Definition = "ii. HoS + HeS vs DL + HD + DR"
  )

## iii. HoS vs DR
state_dat_s10c = state_dat %>%
  mutate(
    X = Mean_percent_DR,
    Y = Mean_percent_HoS,
    SE_X = SE_percent_DR,
    SE_Y = SE_percent_HoS,
    Definition = "iii. HoS vs DR"
  )

state_dat_s10_all = bind_rows(
  state_dat_s10a,
  state_dat_s10b,
  state_dat_s10c
)

axis_max_s10 = max(c(state_dat_s10_all$X + state_dat_s10_all$SE_X,
                     state_dat_s10_all$Y + state_dat_s10_all$SE_Y), na.rm = T)

axis_max_s10 = ceiling(axis_max_s10 / 10) * 10
axis_max_s10 = max(axis_max_s10, 80)


########################################################################################
## A small plotting function using the same visual logic as Figure 4C

plot_state_space = function(plotdf, xlab, ylab, titlelab) {
  
  ggplot(plotdf, aes(x = X, y = Y)) +
    
    annotate("polygon", x = c(0, axis_max_s10, axis_max_s10),
             y = c(0, axis_max_s10, axis_max_s10),
             fill = "#DEEBF7", alpha = 0.45) +
    
    annotate("polygon", x = c(0, axis_max_s10, axis_max_s10),
             y = c(0, 0, axis_max_s10),
             fill = "#FEE8C8", alpha = 0.35) +
    
    geom_abline(intercept = 0, slope = 1,
                linetype = "dashed", linewidth = 0.55,
                color = "grey30") +
    
    geom_errorbar(aes(ymin = Y - SE_Y,
                      ymax = Y + SE_Y),
                  width = 0, linewidth = 0.4, color = "grey35") +
    
    geom_errorbarh(aes(xmin = X - SE_X,
                       xmax = X + SE_X),
                   height = 0, linewidth = 0.4, color = "grey35") +
    
    geom_point(aes(color = Group),
               shape = 16, size = 3.5, alpha = 0.75) +
    
    geom_text_repel(aes(label = Group),
                    size = 3.3,
                    color = "black",
                    min.segment.length = 0,
                    box.padding = 0.25,
                    point.padding = 0.2,
                    max.overlaps = 20) +
    
    annotate("text",
             x = axis_max_s10 * 0.52,
             y = axis_max_s10 * 0.55,
             label = "Equal contribution",
             angle = 45,
             color = "grey25",
             size = 2.8) +
    
    scale_color_manual(values = pal_nejm()(6),
                       breaks = c("D0", "D1", "D3", "D5", "D6", "D7")) +
    
    scale_x_continuous(limits = c(0, axis_max_s10),
                       breaks = seq(0, axis_max_s10, by = 20),
                       expand = c(0, 0)) +
    
    scale_y_continuous(limits = c(0, axis_max_s10),
                       breaks = seq(0, axis_max_s10, by = 20),
                       expand = c(0, 0)) +
    
    coord_fixed() +
    
    labs(x = xlab,
         y = ylab,
         title = titlelab) +
    
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 10, face = "bold", color = "black"),
          axis.title = element_text(color = "black", size = 9),
          axis.text = element_text(color = "black", size = 8),
          legend.position = "none")
}


########################################################################################
## Generate three panels

pS10a = plot_state_space(
  state_dat_s10a,
  xlab = "DL + DR (%)",
  ylab = "HoS (%)",
  titlelab = "i. HoS vs DL + DR"
)

pS10b = plot_state_space(
  state_dat_s10b,
  xlab = "DL + HD + DR (%)",
  ylab = "HoS + HeS (%)",
  titlelab = "ii. HoS + HeS vs DL + HD + DR"
)

pS10c = plot_state_space(
  state_dat_s10c,
  xlab = "DR (%)",
  ylab = "HoS (%)",
  titlelab = "iii. HoS vs DR"
)

pS10 = pS10a + pS10b + pS10c +
  plot_layout(nrow = 1)

pS10

ggsave("Sup_Fig_S4_alternative_state_space.pdf",
       pS10, width = 8.2, height = 3.8)



########################################################################################

state_dat_s10_out = state_dat_s10_all %>%
  select(
    Group,
    Definition,
    X,
    Y,
    SE_X,
    SE_Y,
    Regime,
    Deterministic,
    Stochastic,
    Deterministic_minus_Stochastic,
    Mean_percent_HoS,
    Mean_percent_HeS,
    Mean_percent_DL,
    Mean_percent_HD,
    Mean_percent_DR
  )

write.csv(state_dat_s10_out,
          "Sup_Fig_S4_alternative_state_space_source_data.csv",
          row.names = F)



