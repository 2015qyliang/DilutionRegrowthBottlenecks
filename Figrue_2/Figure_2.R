library(reshape2)
library(ggplot2)
library(ggsci)
# library(ggrepel)
library(vegan)
library(sads)
library(ape)
library(stringr)
# library(dplyr)
library(tidyverse)
library(viridis) # scale_fill_viridis()
library(ggpointdensity) 
library(scales)

library(patchwork)

############################################################################
############################################################################
############################################################################

# # Figure 2A

otudf = read.table('OTUtable.txt', header = TRUE, sep = '\t', row.names = 1)
zero.otudf = otudf[, str_detect(colnames(otudf), 'D0')]
otuSeeds = rowSums(zero.otudf) / ncol(zero.otudf)
otuSeeds = otuSeeds[otuSeeds > 0]
otuSeeds = sort(as.integer(otuSeeds), decreasing = TRUE)
otuSeeds[otuSeeds == 0] = 1 

# -----------------------------------

# SADdf = data.frame(rankA = log10(1:length(otuSeeds)), lgA = log10(otuSeeds))

# package - sads
spe.vec.rad = rad(otuSeeds)
rk.abund = spe.vec.rad[,2]
names(rk.abund) = rownames(spe.vec.rad)

# package - vegan
# The null model (Brokenstick)
rk.Broken = as.numeric(rad.null(rk.abund)$fitted.values)
# The niche preemption model (A purely nonlinear model.)
rk.preempt = as.numeric(rad.preempt(rk.abund)$fitted.values)
# The log-Normal model (Generalized linear models)
rk.lognormal = as.numeric(rad.lognormal(rk.abund)$fitted.values)
# The Zipf model (kind of The Power Distribution; Generalized linear models )
rk.zipf = as.numeric(rad.zipf(rk.abund)$fitted.values)
# The Zipf--Mandelbrot model (Nonlinear model, adding the nonlinear parameter)
rk.zipfbrot = as.numeric(rad.zipfbrot(rk.abund)$fitted.values)

rk.lens = c(length(rk.abund), length(rk.Broken),
            length(rk.preempt), length(rk.lognormal),
            length(rk.zipf), length(rk.zipfbrot))
if (length(unique(rk.lens)) == 1) {
  orders = rep(1:length(rk.abund), 6)
  values = c(rk.abund, rk.Broken, rk.preempt, rk.lognormal, rk.zipf, rk.zipfbrot)
  kinds = c(rep('Abundance-Null', length(rk.abund)),
            rep('Brokenstick', length(rk.Broken)),
            rep('Niche-preemption', length(rk.preempt)),
            rep('Log-Normal', length(rk.lognormal)),
            rep('Zipf', length(rk.zipf)),
            rep('Zipf-Mandelbrot', length(rk.zipfbrot)))
  rk.df = data.frame(orders, values, kinds)
}

#####################

otuSeeds = as.integer(rowSums(zero.otudf)/ncol(zero.otudf))
names(otuSeeds) = names(rowSums(zero.otudf))
otuSeeds[otuSeeds == 0] = 1

spe.radfit = radfit(otuSeeds)
srmodels = spe.radfit$models
paraAICs = c(as.integer(srmodels$Null$aic),
             as.integer(srmodels$Preemption$aic),
             as.integer(srmodels$Lognormal$aic),
             as.integer(srmodels$Zipf$aic),
             as.integer(srmodels$Mandelbrot$aic))
names(paraAICs) = c('Brokenstick', 'Niche-Preemption',
                    'Lognormal', 'Zipf', 'Zipf-Mandelbrot')

write.table(data.frame(Model = names(paraAICs), AICs = paraAICs), 
            'DistributionModel.txt', append = F, sep = '\t', quote = F, 
            row.names = F, col.names = T)

#####################

p = ggplot(rk.df, aes(x = log10(orders), y = log10(values), color = kinds)) + 
  geom_line(linewidth = .5, alpha = 0.8 ) + 
  scale_color_nejm() +
  scale_y_continuous(limits = c(0, 3.7)) + 
  scale_x_continuous(limits = c(0, 3.3)) + 
  labs(x = expression(paste('log'[10] ,'(rank)')), 
       y = expression(paste('log'[10] ,'(abundance)')),
       title = paste0('> The most best model: ', names(which.min(paraAICs))),
       colour = NULL) +
  theme( panel.background = element_rect(fill="white",color="black", linewidth = 0.5),
         panel.grid = element_blank(), 
         plot.title = element_text(size = 10),
         # legend.position = "right",
         legend.position = c(.35, .5), 
         legend.background = element_blank(), 
         legend.text = element_text(size = 8),
         legend.key = element_blank())

ggsave('Figure_2A_all.pdf', p , width = 3, height = 2.2 )

rk.dfsub = rk.df[rk.df$kinds %in% c("Abundance-Null", "Zipf-Mandelbrot"), ]
p = ggplot(rk.dfsub, aes(x = log10(orders), y = log10(values), color = kinds)) + 
  geom_line(linewidth = 1, alpha = 0.6 ) + 
  scale_color_nejm() +
  scale_y_continuous(limits = c(0, 3.7)) + 
  scale_x_continuous(limits = c(0, 3.3)) + 
  labs(x = expression(paste('log'[10] ,'(rank)')), 
       y = expression(paste('log'[10] ,'(abundance)')),
       colour = NULL) +
  theme( panel.background = element_rect(fill="white",color="black", linewidth = 0.5),
         panel.grid = element_blank(), 
         plot.title = element_text(size = 10),
         # legend.position = "right",
         legend.position = c(.35, .5), 
         legend.background = element_blank(), 
         legend.text = element_text(size = 8),
         legend.key = element_blank())

ggsave('Figure_2A.pdf', p , width = 3, height = 2.2 )



############################################################################
############################################################################

# Supplyment Figure 

newdf = data.frame(OTUID = names(sort(rowSums(zero.otudf), decreasing = T)), 
                   observed = spe.vec.rad$abund, 
                   Brokenstick = rk.Broken, 
                   preemption = rk.preempt, 
                   lognormal = rk.lognormal, 
                   zipf = rk.zipf,
                   Mandelbrot = rk.zipfbrot)

write.table(newdf, 'RankAbundanceFit.txt', 
            append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)

#################################

rankFitDf = read.table('RankAbundanceFit.txt', header = T, sep = '\t')
smpgp = read.table('samplesGroup.txt', header = T, sep = '\t')
otudf = read.table('OTUtable.txt', header = T, sep = '\t', row.names = 1)
otudf = otudf[rankFitDf$OTUID, ]
otudf.m = as.matrix(otudf)
otudf.m[otudf.m != 0] = 1

for (dl in c('D1', 'D3', 'D5', 'D6', 'D7')) { 
  diluDf = otudf.m[ , smpgp$SampleID[smpgp$Group == dl]]
  diluDf = data.frame(OTUID = rownames(diluDf), diluDf)
  gps = rep('Obs', ncol(diluDf) - 1)
  
  for (ci in c('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7')) {
    simuFn = paste0('/SimulationFiles/', ci, '-simuDF_', gsub('D', '', dl) , '.txt')
    simuDf = read.table(simuFn, header = T, sep = '\t')
    diluDf = merge(diluDf, simuDf, by = 'OTUID', all.x = T)
    gps = append(gps, rep(ci, ncol(simuDf) - 1))
  }
  
  diluDf.mg = diluDf[ , 2:ncol(diluDf)]
  rownames(diluDf.mg) = diluDf$OTUID
  diluDf.mgt = t(diluDf.mg)
  diluDf.mgtDF = data.frame(Group = gps, diluDf.mgt)
  write.table(diluDf.mgtDF, paste0(dl, 'merged.txt'), 
              append = F, quote = F, sep = '\t', na = '0', 
              row.names = F, col.names = T)
} 

######################################

for (dl in c('D1', 'D3', 'D5', 'D6', 'D7')) { 
  mergeDf = read.table(paste0( dl, 'merged.txt'), 
                       header = T, sep = '\t')
  mergeDf.m = mergeDf[ , 2:ncol(mergeDf)]
  pcoa = pcoa(vegdist(mergeDf.m, 'jac', binary = T))
  dt.scores = data.frame(PCoA1 = pcoa$vectors[ , 1], 
                         PCoA2 = pcoa$vectors[ , 2],
                         gps = factor(gsub('Obs', 'Aobs', mergeDf$Group)))
  
  p = ggplot(dt.scores, aes(x = PCoA1, y = PCoA2, color = gps)) + 
    geom_point(size = 0.75, alpha = 0.5, shape = 16) + 
    scale_color_d3() +
    geom_hline(yintercept=0, colour="#000000", linetype="dashed")+
    geom_vline(xintercept=0, colour="#000000", linetype="dashed")+
    theme( panel.background = element_rect(fill="white",color="black"),
           panel.grid = element_blank(), 
           axis.title = element_text(size = 10, color = 'black'),
           axis.text.y = element_text(size = 8, hjust = 0.5, vjust = 0.5, 
                                      color = 'black', angle = 90),
           axis.text.x = element_text(size = 8, hjust = 0.5, color = 'black'),
           legend.title = element_blank(), 
           legend.text = element_text(size = 8, color = 'black'),
           # legend.position = "right",
           legend.position = "none",
           legend.key = element_blank())
  ggsave(paste0('PCoA_jac_', dl, '.pdf'), p, width = 1.6, height = 1.6)
}

############################################################################
############################################################################

# Figure 2C
# 
# compare the difference among observed and simulated strategies
# 

dilution = simulation = permanova_R2 = permanova_F = permanova_p = anosim_R = anosim_p = c()
mrpp_delta = mrpp_E_delta = mrpp_A = mrpp_p = betadisper_F = betadisper_p = c()

for (dl in c('D1', 'D3', 'D5', 'D6', 'D7')) { 
  
  for (strategy in c( 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7') ){
    mergeDf = read.table(paste0( dl, 'merged.txt'),  header = T, sep = '\t')
    mergeDf = mergeDf[mergeDf$Group %in% c('Obs', strategy) , ]
    mergeDf$Group = gsub('Obs', dl, mergeDf$Group)
    rownames(mergeDf) = paste0('S', 1:nrow(mergeDf))
    
    mergeDf.m = mergeDf[ , 2:ncol(mergeDf)]
    mergeDf.m = mergeDf.m[, colSums(mergeDf.m) != 0 ]    
    
    dfGroup = data.frame(SampleID = rownames(mergeDf.m), Group = mergeDf$Group)
    
    ## Binary Jaccard distance
    dist_obj = vegdist(mergeDf.m, method = "jaccard", binary = TRUE)
    
    ## -------------------------
    ## PERMANOVA
    ## -------------------------
    ad = adonis2(dist_obj ~ Group, data = dfGroup, permutations = 999)
    permanova_R2 =  append(permanova_R2, ad[1, "R2"])
    permanova_F =   append(permanova_F, ad[1, "F"])
    permanova_p =   append(permanova_p, ad[1, "Pr(>F)"])
    
    ## -------------------------
    ## ANOSIM
    ## -------------------------
    an = anosim(dist_obj, grouping = dfGroup$Group, permutations = 999)
    anosim_R =   append(anosim_R, an$statistic)
    anosim_p =   append(anosim_p, an$signif)
    
    ## -------------------------
    ## MRPP
    ## -------------------------
    mr = mrpp(mergeDf.m, grouping = dfGroup$Group, 
              distance = "jaccard", permutations = 999,  weight.type = 1)
    
    mrpp_delta =   append(mrpp_delta, mr$delta)
    mrpp_E_delta =   append(mrpp_E_delta, mr$E.delta)
    mrpp_A =   append(mrpp_A, mr$A)
    mrpp_p =   append(mrpp_p, mr$Pvalue)
    
    ## -------------------------
    ## Betadisper
    ## check “whether PERMANOVA be influenced by the high distribution within group"
    ## -------------------------
    bd = betadisper(dist_obj, group = dfGroup$Group)
    bd_perm = permutest(bd, permutations = 999)
    
    betadisper_F =   append(betadisper_F, bd_perm$tab[1, "F"])
    betadisper_p =  append(betadisper_p, bd_perm$tab[1, "Pr(>F)"])

    dilution = append(dilution, dl)
    simulation = append(simulation, strategy)
  }
}

dfStat = data.frame(Group = dilution, code = simulation, 
                    permanova_R2, permanova_F, permanova_p,
                    anosim_R, anosim_p, 
                    mrpp_delta, mrpp_E_delta, mrpp_A, mrpp_p,
                    betadisper_F, betadisper_p)

strategy_key = data.frame(code = c('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'), 
                          Simulation_Strategy = c('Random_select', "Observed_probability", 
                                                  "Brokenstick", "Preemption", "Lognormal", "Zipf", "Mandelbrot"))
dfStat = merge(dfStat, strategy_key, by = 'code', all.x = T )
dfStat = dfStat[ , c(2, 14, 3:13)]
dfStat = dfStat[order(dfStat$Group, decreasing = F), ]

dfStat = dfStat %>%
  mutate( Group = factor(Group, levels = treatment_order),
          Simulation_Strategy = factor(Simulation_Strategy, levels = strategy_order)  ) %>%
  group_by(Group) %>%
  mutate(Rank_by_PERMANOVA_R2 = rank(permanova_R2, ties.method = "min"),
         Closest_strategy = Rank_by_PERMANOVA_R2 == 1 ) %>%
  ungroup()

write.table(dfStat, 'SupplyTableDistance_Obs_VS_Sim.txt', append = F, 
            quote = F, sep = '\t', row.names = F, col.names = T)

#############################################

# visualize permanova_R2 among observed and dif_simulation
# 

treatment_order = c("D1", "D3", "D5", "D6", "D7")
strategy_order = c("Random_select", "Observed_probability", "Brokenstick", 
                   "Preemption", "Lognormal", "Zipf", "Mandelbrot")

strategy_labels = c( "Random\nselect", "Obs.\nprob.",
                     "Broken-\nstick", "Preemption", "Lognormal", "Zipf", "Mandelbrot")


p = ggplot(dfStat, aes(x = Simulation_Strategy, y = Group, fill = permanova_R2 )) +
  geom_tile( color = "white", linewidth = 0.5 ) +
  ## black outline for closest strategy within each dilution treatment
  # geom_tile( data = dfStat %>% filter(Closest_strategy), aes( x = Simulation_Strategy, y = Group  ),
  #            fill = NA, color = "black", linewidth = 1.2 ) +
  geom_tile( data = dfStat %>% filter(Closest_strategy), fill = NA, 
             color = "black", linewidth = 0.8, width = 0.92, height = 0.85  ) +
  
  ## show R2 values
  # geom_text( aes(label = sprintf("%.3f", permanova_R2)), size = 3.4, color = "black" ) +
  geom_text(aes(  label = sprintf("%.3f", permanova_R2), colour = permanova_R2 > 0.28),
            size = 2.6,show.legend = FALSE ) +
  labs(x = NULL, y = NULL) + 
  # labs( x = "Simulated entry scenario", y = "Dilution"  ) +
  scale_x_discrete( limits = strategy_order, labels = strategy_labels ) +
  scale_y_discrete( limits = rev(treatment_order) ) +
  
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white") ) +
  scale_fill_gradientn(colours = c("#fff7ec", "#fdd49e", "#fc8d59", "#d7301f", "#7f0000"),
                       values  = scales::rescale(c(0.08, 0.15, 0.23, 0.32, 0.40)),
                       limits  = c(0.08, 0.40),
                       breaks  = c(0.10, 0.20, 0.30, 0.40),
                       name    = expression(PERMANOVA~R^2),
                       guide = guide_colorbar(title.position = "top",  title.hjust = 0.5,  
                                              barheight = unit(28, "mm"),  barwidth  = unit(3.2, "mm"),  
                                              ticks.colour = "black",  frame.colour = "black")  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.background = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black" , size = 10),
        axis.text.y = element_text(color = "black", size = 10 ),
        axis.title = element_text(color = "black", size = 12),
        legend.title = element_text(face = "bold", size = 10),
        legend.key.width = unit(2, 'mm'),
        legend.direction = 'vertical',
        legend.background = element_blank(), 
        legend.box.background = element_blank(), 
        legend.position = "bottom" )

ggsave('Figure_2C_permanovaR2.pdf', p, width = 3.2, height = 3.5)


############################################################################
############################################################################
############################################################################

# Figure 2D
 
#  Parameters
dominance_cutoff = 0.001
detect_cutoff = 0

treat_levels = c("D0", "D1", "D3", "D5", "D6", "D7")

group_cols = c( "Abundant OTUs" = "#0072B2", "Rare OTUs"     = "#D55E00" )

heat_cols = c( "#f7fcf0", "#ccebc5", "#7bccc4", "#2b8cbe", "#084081" )

theme_pub = theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0),
        axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(color = "black", size = 8),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold", size = 8),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8)  )


# =========================
# Read and reshape data
ra_raw = read.table( "RelativeAbundance.txt", header = T, sep = '\t')

ra_long = ra_raw %>%
  pivot_longer(cols = -OTUID,
               names_to = "SampleID",
               values_to = "RelAbund" ) %>%
  mutate(RelAbund = as.numeric(RelAbund),
         Treatment = stringr::str_extract(SampleID, "^D[0-9]"),
         Treatment = factor(Treatment, levels = treat_levels) ) %>%
  filter(!is.na(Treatment))

sample_info = ra_long %>%
  distinct(SampleID, Treatment) %>%
  count(Treatment)

print(sample_info)

###########################################

otu_class = ra_long %>%
  filter(Treatment == "D0") %>%
  group_by(OTUID) %>%
  summarise(mean_RA_D0 = mean(RelAbund, na.rm = TRUE),
            occ_D0 = mean(RelAbund > detect_cutoff, na.rm = TRUE),
            .groups = "drop"  ) %>%
  mutate(OTU_group = case_when(mean_RA_D0 >= dominance_cutoff ~ "Abundant OTUs",
                               mean_RA_D0 > detect_cutoff & mean_RA_D0 < dominance_cutoff ~ "Rare OTUs",
                               TRUE ~ "Absent in D0"  )  )

otu_class_count = otu_class %>%
  count(OTU_group)

print(otu_class_count)

# -----------------------------------

ra_focus = ra_long %>%
  inner_join(otu_class %>%
               filter(OTU_group %in% c("Abundant OTUs", "Rare OTUs")) %>%
               select(OTUID, OTU_group, mean_RA_D0, occ_D0),
             by = "OTUID"  ) %>%
  mutate(OTU_group = factor(OTU_group,levels = c("Abundant OTUs", "Rare OTUs")  )  )

group_n = ra_focus %>%
  distinct(OTUID, OTU_group) %>%
  count(OTU_group)

group_labs = setNames( paste0(group_n$OTU_group, " (", group_n$n, ")"),  group_n$OTU_group )

print(group_n)

# ---------------------------------------

pseudo = min(ra_focus$RelAbund[ra_focus$RelAbund > detect_cutoff], na.rm = TRUE) / 2

otu_treat = ra_focus %>%
  group_by(OTUID, OTU_group, Treatment) %>%
  summarise(n_replicates = n(),
            occupancy = mean(RelAbund > detect_cutoff, na.rm = TRUE),
            retained = any(RelAbund > detect_cutoff, na.rm = TRUE),
            mean_RA = mean(RelAbund, na.rm = TRUE),
            median_RA_present = if (sum(RelAbund > detect_cutoff, na.rm = TRUE) > 0) {
              median(RelAbund[RelAbund > detect_cutoff], na.rm = TRUE)
            } else {  0  },
            .groups = "drop" )

d0_base = otu_treat %>%
  filter(Treatment == "D0") %>%
  select( OTUID, median_RA_D0 = median_RA_present, occupancy_D0 = occupancy  )

otu_treat = otu_treat %>%
  left_join(d0_base, by = "OTUID") %>%
  mutate(log10_shift_vs_D0 = log10(median_RA_present + pseudo) - log10(median_RA_D0 + pseudo)   )


# -------------------------------------------
# Panel D-i：retained OTU proportion
retention_df = otu_treat %>%
  group_by(OTU_group, Treatment) %>%
  summarise( retained_OTUs = sum(retained), total_OTUs = n_distinct(OTUID), 
             retained_prop = retained_OTUs / total_OTUs, 
             mean_occupancy = mean(occupancy), .groups = "drop")

p_retention = ggplot(retention_df, aes(x = Treatment, y = retained_prop, group = OTU_group, color = OTU_group ) ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_color_manual(values = group_cols, name = NULL) +
  scale_y_continuous( labels = percent_format(accuracy = 1), 
                      limits = c(0, 1), 
                      expand = expansion(mult = c(0.02, 0.05))) +
  labs(title = "i. Retention of D0-defined OTUs", x = NULL, y = "Retained OTUs (%)" ) +
  theme_pub

# -------------------------------------------------
# Panel D-ii：occupancy class heatmap
occ_levels = c( "Lost\n(0)", "Transient\n(0–0.25]", "Intermediate\n(0.25–0.50]", 
                "High\n(0.50–<1.00)", "Persistent\n(1.00)")

heat_occ = otu_treat %>%
  mutate(Occupancy_class = case_when(occupancy == 0 ~ "Lost\n(0)",
                                     occupancy > 0 & occupancy <= 0.25 ~ "Transient\n(0–0.25]",
                                     occupancy > 0.25 & occupancy <= 0.50 ~ "Intermediate\n(0.25–0.50]",
                                     occupancy > 0.50 & occupancy < 1 ~ "High\n(0.50–<1.00)",
                                     occupancy == 1 ~ "Persistent\n(1.00)"),
         Occupancy_class = factor(Occupancy_class, levels = occ_levels) ) %>%
  count(OTU_group, Treatment, Occupancy_class, name = "n") %>%
  group_by(OTU_group, Treatment) %>%
  mutate( total = sum(n), prop = n / total) %>%
  ungroup() %>%
  complete( OTU_group, Treatment, Occupancy_class, fill = list(n = 0, prop = 0) )

p_occ_heat = ggplot( heat_occ, aes(x = Treatment, y = Occupancy_class, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.35) +
  facet_wrap( ~ OTU_group, nrow = 1, labeller = as_labeller(group_labs)) +
  scale_fill_gradientn( colors = heat_cols, 
                        labels = percent_format(accuracy = 1), 
                        name = "Proportion\nof OTUs") +
  scale_y_discrete(drop = FALSE) +
  labs( title = "ii. Shift from persistent to transient/lost states", x = NULL, y = "Occupancy class") +
  theme_pub +
  theme(panel.spacing = unit(0.8, "lines") )

# --------------------------------------
# Panel D-iii：retained OTUs
shift_retained = otu_treat %>%
  filter( Treatment != "D0", retained == TRUE, median_RA_D0 > 0)

p_shift = ggplot( shift_retained, aes(x = Treatment, y = log10_shift_vs_D0, fill = OTU_group) ) +
  geom_hline( yintercept = 0, linetype = "dashed", 
              linewidth = 0.4, color = "grey40" ) +
  geom_violin( position = position_dodge(width = 0.78), 
               width = 0.72, trim = TRUE, scale = "width", 
               linewidth = 0.25, alpha = 0.75) +
  geom_boxplot( position = position_dodge(width = 0.78), 
                width = 0.16, outlier.shape = NA, 
                linewidth = 0.25, alpha = 0.85) +
  scale_fill_manual(values = group_cols, name = NULL) +
  labs( title = "iii. Regrowth shift among retained OTUs", x = "Dilution treatment", 
        y = expression(Delta~log[10]~"relative abundance vs D0")) +
  theme_pub

# --------------------------------------------

fig2D_main = (p_retention /
                p_occ_heat /
                p_shift ) +
  plot_layout( heights = c(0.9, 1.35, 1.15), guides = "collect" ) &
  theme(legend.position = "bottom" )

fig2D_main

ggsave("Figure2D_compact_main_vertical.pdf", fig2D_main, width = 4.4, height = 5.7)

write_csv(retention_df, "Fig2D_retention_summary.csv")
write_csv(heat_occ, "Fig2D_occupancy_class_summary.csv")
write_csv(otu_treat, "Fig2D_OTU_treatment_summary.csv")

# ------------------------------------------
# 

p_supp_occ_abund = ggplot( otu_treat, aes(   x = median_RA_present + pseudo,   y = occupancy ) ) +
  geom_point( aes(color = OTU_group), size = 0.45, alpha = 0.35  ) +
  facet_grid(OTU_group ~ Treatment, labeller = labeller(OTU_group = as_labeller(group_labs)) ) +
  scale_x_log10(labels = scientific_format() ) +
  scale_color_manual(values = group_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1) ) +
  labs(x = "Median relative abundance when present",
       y = "Occupancy",
       title = "Supplementary Fig. Sx. Occupancy–abundance relationships"
  ) +
  theme_pub +
  theme( legend.position = "bottom", 
         axis.text = element_text(color = "black", size = 10),
         axis.title = element_text(color = "black", size = 10),
         strip.text = element_text(size = 10, face = "bold"), 
         axis.text.x = element_text(angle = 45, hjust = 1))

p_supp_occ_abund

ggsave("Figure2D_Sup_occupancy_abundance.pdf", p_supp_occ_abund, width = 8, height = 4.7)

############################################################################
############################################################################
############################################################################

# Figure 2E
# 
# the distribution of shared/unique OTUs via UpSet
# 

otu = read_tsv("OTUtable.txt", show_col_types = FALSE)

sample_groups = list(D0 = grep("^D0", colnames(otu), value = TRUE),
                     D1 = grep("^D1", colnames(otu), value = TRUE),
                     D3 = grep("^D3", colnames(otu), value = TRUE),
                     D5 = grep("^D5", colnames(otu), value = TRUE),
                     D6 = grep("^D6", colnames(otu), value = TRUE),
                     D7 = grep("^D7", colnames(otu), value = TRUE) )
set_names = names(sample_groups)

## Convert OTU abundance to treatment-level presence/absence
## min_replicates = 1 ; min_replicates = 2, test it in Supplementary sensitivity analysis。
min_replicates = 1

otu_num = otu %>%
  mutate( across( .cols = -OTUID, .fns = ~ as.numeric(.x) ) )
otu_sets = tibble(OTUID = otu_num$OTUID)

for (g in set_names) {
  cols = sample_groups[[g]]
  otu_sets[[g]] = rowSums(otu_num[, cols] > 0, na.rm = TRUE) >= min_replicates
}

## delete absent OTU in all treatment 
otu_sets = otu_sets %>%
  filter(if_any(all_of(set_names), ~ .x))

## treatment-level presence/absence table
write_csv( otu_sets, "Figure2E_OTU_treatment_presence_absence.csv" )

## Set sizes: number of OTUs detected in each treatment
set_sizes = otu_sets %>%
  summarise(across(all_of(set_names), ~ sum(.x))) %>%
  pivot_longer( cols = everything(), names_to = "Treatment", values_to = "Set_size"
  ) %>%
  mutate(Treatment = factor(Treatment, levels = rev(set_names)) )
write_csv(set_sizes, "Figure2E_set_sizes.csv")

## Compute exact intersection sizes
intersection_table = otu_sets %>%
  group_by(across(all_of(set_names))) %>%
  summarise( Intersection_size = n(), OTUs = paste(OTUID, collapse = ";"), .groups = "drop"  )

## remove Null intersection
intersection_table = intersection_table %>%
  mutate( Degree = rowSums(across(all_of(set_names)))  ) %>%
  filter(Degree > 0)

## create intersection label
intersection_table$Intersection_label = apply( intersection_table[, set_names], 1, 
                                               function(x) { paste(set_names[as.logical(x)], collapse = "&") } )

## mark shared and specific
intersection_table = intersection_table %>%
  mutate( Shared_by_all = Degree == length(set_names), Unique_to_one = Degree == 1 ) %>%
  arrange(desc(Intersection_size))

## export intersection table
write_csv( intersection_table, "Figure2E_all_intersections.csv" )

## ---------------------------------------------------------
## Select intersections for main Figure 2E

## Main show top_n intersection，
## 1)  shared in six
## 2) specific respectively
top_n = 20

top_intersections = intersection_table %>%
  arrange(desc(Intersection_size)) %>%
  slice_head(n = top_n)

must_keep = intersection_table %>%
  filter(Shared_by_all | Unique_to_one)

plot_intersections = bind_rows(top_intersections, must_keep) %>%
  distinct(Intersection_label, .keep_all = TRUE) %>%
  arrange(desc(Intersection_size), desc(Degree))

plot_intersections = plot_intersections %>%
  mutate( Intersection_id = row_number(), Intersection_id = factor(Intersection_id, levels = Intersection_id) )

## export the showed intersection in visualization
write_csv( plot_intersections, "Figure2E_main_plot_intersections.csv" )


## ---------------------------------------------------------
## Prepare matrix data for UpSet plot

matrix_dat = plot_intersections %>%
  select(Intersection_id, Intersection_label, all_of(set_names)) %>%
  pivot_longer( cols = all_of(set_names), names_to = "Treatment", values_to = "Present"  ) %>%
  mutate(Treatment = factor(Treatment, levels = rev(set_names)),
         Intersection_id = factor( Intersection_id, levels = plot_intersections$Intersection_id ) )

## link present point
line_dat = matrix_dat %>%
  filter(Present) %>%
  group_by(Intersection_id) %>%
  mutate(n_present = n()) %>%
  ungroup() %>%
  filter(n_present > 1)

## Plot: intersection size barplot
p_intersection = ggplot( plot_intersections, aes(x = Intersection_id, y = Intersection_size) ) +
  geom_col( fill = "#4C78A8", width = 0.75  ) +
  # geom_text( aes(label = Intersection_size), vjust = -0.25, size = 2.8  ) +
  labs( y = NULL, x = NULL  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))  ) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1, color = "black"),
        axis.ticks.x = element_blank(),
        # axis.title.y = element_text(size = 10),
        plot.margin = margin(2, 2, 0, 2, unit = "mm"))

            
## Plot: UpSet combination matrix
p_matrix = ggplot( matrix_dat, aes(x = Intersection_id, y = Treatment) ) +
  ## lines connecting present sets
  geom_line( data = line_dat, aes(group = Intersection_id), color = "black", linewidth = 0.45  ) +
  ## absent points
  geom_point( data = matrix_dat %>% filter(!Present), shape = 21, size = 2.6, 
              fill = "#E0E0E0", color = "#BDBDBD", stroke = 0.4) +
  ## present points
  geom_point( data = matrix_dat %>% filter(Present), shape = 21, size = 2.8, 
              fill = "#1F78B4", color = "black", stroke = 0.45) +
  labs( x = NULL, y = NULL ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        plot.margin = margin(0, 2, 2, 2, unit = "mm"))


## Plot: total OTUs per treatment as inset
## put into right-top
set_sizes_inset = set_sizes %>%
  mutate(Treatment = factor(Treatment, levels = set_names)  )

p_setsize_inset = ggplot( set_sizes_inset, aes(x = Treatment, y = Set_size )) +
  geom_col( fill = "#9ECAE1", width = 0.7   ) +
  # geom_text( aes(label = Set_size), angle = 90, vjust = -0.25, size = 2.2  ) +
  labs( x = NULL, y = "Total OTUs\nper dilution group"  ) +
  scale_y_continuous( position = "left") +
  theme_bw(base_size = 7) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text( size = 10, angle = 0, color = "black" ),
        axis.text.y = element_text(size = 8, color = "black"),
        # axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1,size = 10, color = "black"),
        axis.ticks = element_line(linewidth = 0.25),
        axis.title.y = element_text(size = 10, color = "black" ) ,
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)   )

## Add inset into intersection-size barplot
p_intersection_inset = p_intersection +
  inset_element( p_setsize_inset, left = 0.2, bottom = 0.20, right = 0.98, top = 0.95, align_to = "panel"   )

## Combine intersection bars and matrix
p_upset_inset = p_intersection_inset / p_matrix +
  plot_layout(heights = c(2.6, 1.15) ) +
  theme_bw( )

ggsave("Figure2E_UpSet_main_with_inset.pdf", p_upset_inset, width = 3.2, height = 3)



















