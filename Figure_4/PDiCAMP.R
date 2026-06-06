
library(ggplot2)
library(ggsci)
library(reshape2)
library(ggrepel)
library(doBy)
library(tidyverse)

########################################################################################

eachTurnover = read.csv('PD.ProcessImportance_EachTurnover.csv')

eachTurnover$samp1 = substr(eachTurnover$samp1, 1, 2)
eachTurnover$samp2 = substr(eachTurnover$samp2, 1, 2)

eachTurnoverSub = eachTurnover[eachTurnover$samp1 == eachTurnover$samp2, ]

eachTurnoverSub = eachTurnoverSub[, c("samp2", "HeS", "HoS", "DL", "HD", "DR")]
colnames(eachTurnoverSub)[1] = "Group"

write.csv(eachTurnoverSub, "PDiCAMP_processGroup.csv", row.names = F)

#######################################


# procdf = melt(eachTurnoverSub, variable.name = "Process", value.name = "Contribution",
#               measure.vars = c("HeS", "HoS", "DL", "HD", "DR") )

## Treatment order
treatment_order = c("D0", "D1", "D3", "D5", "D6", "D7")

dat = eachTurnoverSub %>%
  mutate(
    Group = factor(Group, levels = treatment_order),
    Dilution_strength = case_when(
      Group == "D0" ~ 0,
      Group == "D1" ~ 1,
      Group == "D3" ~ 3,
      Group == "D5" ~ 5,
      Group == "D6" ~ 6,
      Group == "D7" ~ 7
    )
  )


## Optional check: process contributions should sum to ~1
dat = dat %>%
  mutate(Process_sum = HeS + HoS + DL + HD + DR)

summary(dat$Process_sum)


## ---------------------------------------------------------
## 2. Convert to long format

process_long = dat %>%
  select(Group, Dilution_strength, HeS, HoS, DL, HD, DR) %>%
  pivot_longer(
    cols = c(HeS, HoS, DL, HD, DR),
    names_to = "Process_code",
    values_to = "Contribution"
  ) %>%
  mutate(
    Process = recode(
      Process_code,
      HoS = "Homogeneous selection",
      HeS = "Heterogeneous selection",
      DL  = "Dispersal limitation",
      HD  = "Homogenizing dispersal",
      DR  = "Drift"
    ),
    Process = factor(
      Process,
      levels = c(
        "Homogeneous selection",
        "Heterogeneous selection",
        "Dispersal limitation",
        "Homogenizing dispersal",
        "Drift"
      )
    )
  )

## ---------------------------------------------------------
## 3. Summary statistics

process_summary = process_long %>%
  group_by(Group, Dilution_strength, Process, Process_code) %>%
  summarise(
    Mean = mean(Contribution, na.rm = TRUE),
    SD = sd(Contribution, na.rm = TRUE),
    N = n(),
    SE = SD / sqrt(N),
    .groups = "drop"
  ) %>%
  mutate(
    Mean_percent = Mean * 100,
    SE_percent = SE * 100
  )

write.csv(process_summary,"Figure4ABC_iCAMP_process_summary.csv", row.names = F)

## Wide summary for panel C
state_summary = process_summary %>%
  select(Group, Dilution_strength, Process_code, Mean_percent, SE_percent) %>%
  pivot_wider(
    names_from = Process_code,
    values_from = c(Mean_percent, SE_percent)
  ) %>%
  arrange(Dilution_strength)

write.csv(state_summary,"Figure4C_state_space_summary.csv", row.names = F)


## ---------------------------------------------------------
## 4. Color palette

process_cols = c(
  "Homogeneous selection"   = "#E66101",  # orange
  "Heterogeneous selection" = "#FDB863",  # light orange
  "Dispersal limitation"    = "#1B9E77",  # teal
  "Homogenizing dispersal"  = "#7570B3",  # purple
  "Drift"                  = "#2166AC"   # blue
)

## ---------------------------------------------------------
## 5. Panel A: stacked barplot

pA = ggplot(process_summary, aes(x = Group, y = Mean_percent, fill = Process)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.35, alpha = 0.6  ) +
  scale_fill_manual(values = process_cols, name = "Assembly process"  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)   ) +
  labs(x = NULL, y = "Relative contribution (%)", title = NULL) + 
       # title = "A. Assembly-process composition"  ) +
  theme_bw() +
  theme(plot.background = element_blank(),
        # plot.title = element_text(face = "bold", size = 12),
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(color = "black", size = 8),
        legend.key.height = unit(2, 'mm'),
        legend.key.width = unit(2, 'mm'),
        legend.position = "right",
        legend.title = element_text(),
        legend.text = element_text(size = 6)  )

pA

## ---------------------------------------------------------
## 6. Panel B: HoS, DR and DL trends

trend_dat = process_summary %>%
  filter(Process_code %in% c("HoS", "DR", "DL")) %>%
  mutate(
    Process_trend = recode(
      Process_code,
      HoS = "Homogeneous selection",
      DR  = "Drift",
      DL  = "Dispersal limitation"
    ),
    Process_trend = factor(
      Process_trend,
      levels = c("Homogeneous selection", "Drift", "Dispersal limitation")
    )
  )

trend_cols = c(
  "Homogeneous selection" = "#E66101",
  "Drift" = "#2166AC",
  "Dispersal limitation" = "#1B9E77"
)

pB = ggplot(trend_dat, aes(x = Group, y = Mean_percent, group = Process_trend, color = Process_trend ) ) +
  geom_line(linewidth = 0.85, alpha = 0.6) +
  geom_point(size = 2.4, alpha = 0.6) +
  geom_errorbar(aes(ymin = Mean_percent - SE_percent, ymax = Mean_percent + SE_percent ),
                width = 0.12, linewidth = 0.45  ) +
  scale_color_manual(values = trend_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 100, 25), 
                     expand = expansion(mult = c(0.02, 0.05))  ) +
  labs(x = NULL, y = "Relative contribution (%)", title = NULL) +
  # labs(x = "Dilution treatment", y = "Relative contribution (%)",
  #      title = "B. Drift increased while homogeneous selection was retained") +
  theme_bw() +
  theme(plot.background = element_blank(),
        # plot.title = element_text(face = "bold", size = 12),
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(color = "black", size = 10),
        legend.key.height = unit(2, 'mm'),
        legend.key.width = unit(2, 'mm'),
        # legend.position = "bottom",
        legend.position = c(0.74, 0.86), 
        legend.title = element_text(),
        legend.text = element_text(size = 8)  )

pB

## ---------------------------------------------------------
## 7. Panel C:  Deterministic - Stochastic
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
    )
  ) %>%
  arrange(Dilution_strength)
print(
  state_dat %>%
    select(
      Group,
      Deterministic,
      Stochastic,
      Deterministic_minus_Stochastic,
      Regime
    )
)

write.csv(state_dat, "Figure4C_deterministic_stochastic_state_summary.csv", row.names = F)

## Axis upper limit
axis_max = max(c(state_dat$Deterministic + state_dat$SE_Deterministic,
                 state_dat$Stochastic + state_dat$SE_Stochastic),  na.rm = T )

axis_max = ceiling(axis_max / 10) * 10
axis_max = max(axis_max, 80)

pC = ggplot(state_dat,  aes(x = Stochastic,y = Deterministic  ) ) +
  
  ## Region above the 1:1 balance line
  ## Use a polygon to shade deterministic-enriched region
  annotate("polygon", x = c(0, axis_max, axis_max), y = c(0, axis_max, axis_max),
           fill = "#DEEBF7", alpha = 0.45  ) +
  
  ## Region below the 1:1 line
  annotate( "polygon", x = c(0, axis_max, axis_max), y = c(0, 0, axis_max),
            fill = "#FEE8C8", alpha = 0.35  ) +
  
  ## Balance line: deterministic = stochastic
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.6, color = "grey30"  ) +
  
  ## Trajectory along dilution order
  # geom_path(aes(group = 1), arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
  #           linewidth = 0.75, color = "grey25"  ) +
  
  ## Error bars
  geom_errorbar(aes(ymin = Deterministic - SE_Deterministic,
                    ymax = Deterministic + SE_Deterministic),
                width = 0, linewidth = 0.45, color = "grey35" ) +
  
  geom_errorbarh(aes(xmin = Stochastic - SE_Stochastic, xmax = Stochastic + SE_Stochastic),
                 height = 0, linewidth = 0.45, color = "grey35"  ) +
  
  ## Points
  geom_point(aes(color = Group), shape = 16, size = 4,alpha = 0.7) +
  
  ## Labels
  geom_text_repel(aes(label = Group), size = 3.8, 
                  color = "black", min.segment.length = 0,
                  box.padding = 0.3, point.padding = 0.25, max.overlaps = 20 ) +
  
  ## Region labels
  annotate("text", x = axis_max * 0.26, y = axis_max * 0.6,
           label = "Deterministic-enriched\nregrowth regime",
           color = "#08519C", size = 3, fontface = "italic"  ) +
  
  annotate("text", x = axis_max * 0.72, y = axis_max * 0.18,
           label = "Stochastic-dominated\nbaseline region",
           color = "#A63603", size = 3, fontface = "italic" ) +
  
  annotate("text",  x = axis_max * 0.52, y = axis_max * 0.54,
           label = "Deterministic = stochastic",
           angle = 45, color = "grey25", size = 3  ) +
  
  scale_color_manual(values = pal_nejm()(6), breaks = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')) +
  scale_x_continuous(limits = c(0, axis_max), breaks = seq(0, axis_max, by = 20), expand = c(0, 0) ) +
  scale_y_continuous(limits = c(0, axis_max), breaks = seq(0, axis_max, by = 20), expand = c(0, 0) ) +
  coord_fixed() +
  labs(x = "Stochastic / dispersal-related contribution\n(DL + HD + DR, %)",
       y = "Deterministic contribution\n(HoS + HeS, %)",
       title = NULL  ) +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.title = element_text(color = "black", size = 10 ),
        axis.text = element_text(color = "black", size = 8   ),
        legend.position = "none"   )

pC

# ggsave( "Figure4C_deterministic_stochastic_state_space.pdf",
#   pC,  width = 6.2, height = 5.4)

## ---------------------------------------------------------
## 8. Panel D:  The contribution of Bins

# filter Bins by diff-taxs
binTaxDF = read.table('PD.Bin_TopTaxon.csv', header = T, sep = ',', stringsAsFactors = F)
binsDF.f = binTaxDF[, c("Bin", "TopTaxon.Family", "TopTaxon.Phylum")]
colnames(binsDF.f) = c("Bins", "Family", "Phylum")

binContriDF = read.table('PD.BinContributeToProcess_EachGroup.csv', 
                         header = T, sep = ',', stringsAsFactors = F)
binContriDF = binContriDF[ !(str_detect(binContriDF$Group, "_vs_")) , 3:ncol(binContriDF)]
colnames(binContriDF)[3:ncol(binContriDF)] = paste0('Bin', 1:(ncol(binContriDF) - 2))
binContriDF = binContriDF[binContriDF$Process %in% c('DL', 'DR', 'HoS') , ]

binContriDFmelt = melt(binContriDF, variable.name = 'Bins',
                       measure.vars = colnames(binContriDF)[3:ncol(binContriDF)] )
binContriDFsum = summaryBy(.~Group+Bins+Process, data = binContriDFmelt, FUN = mean)

#######################################
# select Top Bins

binSum = summaryBy(.~Bins, data = binContriDFmelt, FUN = sum)
binSum = binSum[order(binSum$value.sum, decreasing = T), ]
topBins = binSum$Bins[1:10]

###################

binConDF = merge(binContriDFsum, binsDF.f, by = 'Bins', all.x = T)
binConDF = binConDF[order(binConDF$value.mean, decreasing = T), ]

# ------------------------------
binConDF$Rank_within_group_process = ave(
  binConDF$value.mean,
  binConDF$Group,
  binConDF$Process,
  FUN = function(x) rank(-x, ties.method = "min") 
)
# --------------------------------



binConDF.filt = binConDF[binConDF$Family != 'uncultured' & 
                           binConDF$Process != 'HD' & 
                           binConDF$Family != 'uncultured__bacterium' & 
                           binConDF$Bins %in% topBins, ]
###################

binConDF.filtCut = binConDF.filt[binConDF.filt$value.mean > 0.01 , ]
binConDF.filtCut$Family = gsub('__', '_',binConDF.filtCut$Family)
binConDF.filtCut$value.mean = log10(100*binConDF.filtCut$value.mean)
binConDF.filtCut$Family = paste0(binConDF.filtCut$Family, " ", binConDF.filtCut$Bins)

pD = ggplot(binConDF.filtCut, aes(x = Family, y = value.mean, fill = Process)) + 
  geom_col(width = 0.7, alpha = 0.6, position = 'stack') + 
  scale_fill_manual(values = c('#1B9E77', '#2166AC',  '#E66101'), 
                    breaks = c('DL', 'DR', 'HoS')) +
  facet_grid(Group~., scales = 'free', space = 'free') +
  labs(x = NULL, y = 'log10(Relative importance\nof ecological process)') + 
  theme_bw() +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        strip.background = element_blank(), 
        legend.position = 'none', 
        axis.title.x = element_text(size = 8,  color = 'black'),
        axis.text.y = element_text(size = 8,face = 'italic', color = 'black',
                                   hjust = 1, vjust = 0.5 ),
        axis.text.x = element_text(size = 8, color = 'black') )  + coord_flip()

colnames(binConDF)[4] = 'Relative importance'

write.csv(binConDF, "Figure4D_TheRelativeImportanceOfBins_summary_V2.csv", row.names = F)

##############################################

ggsave('Figure4A.pdf', pA, width = 3.5, height = 1.7 )
ggsave('Figure4B.pdf', pB, width = 2.1, height = 1.7 )
ggsave('Figure4C.pdf', pC, width = 4.5, height = 2.8 )
ggsave('Figure4D_TopBins_RelativeImportance.pdf', pD, width = 3, height = 5 )

##############################################


# ----------------------------------
shownKeys = unique(binConDF.filtCut[, c("Bins", "Group", "Process")])
shownKeys$Shown_in_Fig4D = "Yes"

binConDF = merge(
  binConDF,
  shownKeys,
  by = c("Bins", "Group", "Process"),
  all.x = TRUE )

binConDF$Shown_in_Fig4D[is.na(binConDF$Shown_in_Fig4D)] = "No"

binConDF$Family = gsub("__", "_", binConDF$Family)

binConDF.out = binConDF

colnames(binConDF.out)[colnames(binConDF.out) == "Bins"] = "Bin ID"
colnames(binConDF.out)[colnames(binConDF.out) == "Group"] = "Dilution treatment"
colnames(binConDF.out)[colnames(binConDF.out) == "Process"] = "Assembly process"
colnames(binConDF.out)[colnames(binConDF.out) == "value.mean"] = "Relative importance"
colnames(binConDF.out)[colnames(binConDF.out) == "Rank_within_group_process"] = "Rank within group-process"

binConDF.out$`Relative importance (%)` = binConDF.out$`Relative importance` * 100

binConDF.out = binConDF.out[, c(
  "Bin ID",
  "Dilution treatment",
  "Assembly process",
  "Relative importance",
  "Relative importance (%)",
  "Rank within group-process",
  "Family",
  "Phylum"  )]
write.csv(binConDF.out, "Figure4D_TheRelativeImportanceOfBins_summary_V1.csv",  row.names = F)
# -----------------------------------


