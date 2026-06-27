###############################################################################
# 01-FAPROTAX_supp4panel_figures_only.R
#
# Purpose:
#   microeco + FAPROTAX, compare predicted functional differences
#   among D0, D1, D3, D5, D6, D7 groups.
#
# Output:
#   Only figures are exported in the current working directory.
#   No intermediate .txt result files are written.
#
# Figures:
#   panelA_FAPROTAX_coverage.pdf
#   panelB_FAPROTAX_PCoA_coverageNorm.pdf
#   panelC_FAPROTAX_focal_guilds_coverageNorm.pdf
#   panelD_FAPROTAX_dispersion_coverageNorm.pdf
#   FAPROTAX_supplementary_figure_current_path.pdf
###############################################################################

library(microeco)
library(ape)
library(vegan)
library(reshape2)
library(doBy)
library(ggplot2)
library(ggsci)
library(patchwork)

# ------------------------------------------------------------
# 1. Read data + build microtable
# ------------------------------------------------------------

sample_info = read.table('OTUsamplesGroup.txt', header = T,
                         sep = '\t', stringsAsFactors = F)
rownames(sample_info) = sample_info$SampleID

otu_table = read.table('OTUtable.txt', header = T, sep = '\t',
                       row.names = 1, stringsAsFactors = F,
                       check.names = F)

taxonomy_table = read.table('OTUTax.txt',  header = T, sep = '\t',
                            row.names = 1, stringsAsFactors = F,
                            check.names = F)

phylo_tree = read.tree('OTUFastTree.txt')

dataset = microtable$new(sample_table = sample_info,
                         otu_table = otu_table,
                         tax_table = taxonomy_table,
                         phylo_tree = phylo_tree)

dataset$tidy_dataset()

# ------------------------------------------------------------
# 2. Group information and color setting
# ------------------------------------------------------------

# 如果你的分组列叫 sample_group，就改成：group_col = 'sample_group'
group_col = 'Group'

group_levels = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')

sample_info = dataset$sample_table
sample_info$SampleID = rownames(sample_info)
sample_info[, group_col] = factor(sample_info[, group_col], levels = group_levels)

meta_group = data.frame(SampleID = sample_info$SampleID,
                        Group = sample_info[, group_col],
                        stringsAsFactors = F)
rownames(meta_group) = meta_group$SampleID

group_pal = pal_nejm()(6)
names(group_pal) = group_levels

# ------------------------------------------------------------
# 3. Run FAPROTAX
# ------------------------------------------------------------

func = trans_func$new(dataset = dataset)
func$for_what = 'prok'

func$cal_spe_func(prok_database = 'FAPROTAX')
func$cal_spe_func_perc(abundance_weighted = TRUE, perc = TRUE, dec = 4)

# ------------------------------------------------------------
# 4. Annotation coverage
#    Coverage = read fraction from OTUs with at least one FAPROTAX function
# ------------------------------------------------------------

otu_anno = func$res_spe_func

if(!any(rownames(otu_anno) %in% rownames(dataset$otu_table)) &&
   any(colnames(otu_anno) %in% rownames(dataset$otu_table))){
  otu_anno = t(otu_anno)
}

otu_anno = as.matrix(otu_anno)

annotated_otus = rownames(otu_anno)[rowSums(otu_anno > 0) > 0]

otu_counts = as.matrix(dataset$otu_table)
otu_relabun = sweep(otu_counts, 2, colSums(otu_counts), '/')

annotated_otus = intersect(annotated_otus, rownames(otu_relabun))
coverage_sample = colSums(otu_relabun[annotated_otus, , drop = F])

coverage_df = data.frame(SampleID = names(coverage_sample),
                         Coverage = as.numeric(coverage_sample),
                         stringsAsFactors = F)

coverage_df$Group = meta_group[coverage_df$SampleID, 'Group']
coverage_df$Group = factor(coverage_df$Group, levels = group_levels)

coverage_group = summaryBy(Coverage ~ Group,
                           data = coverage_df,
                           FUN = c(mean, median, sd))

coverage_kw = kruskal.test(Coverage ~ Group, data = coverage_df)
coverage_p = coverage_kw$p.value

# ------------------------------------------------------------
# 5. Sample x function matrix
#    func_mat     : % of total reads
#    func_mat_cov : % of annotated reads, normalized by annotation coverage
# ------------------------------------------------------------

func_table = func$res_spe_func_perc
func_table$SampleID = rownames(func_table)

func_group = merge(meta_group[, c('SampleID', 'Group')],
                   func_table,
                   by = 'SampleID')

func_mat = func_group[, !(colnames(func_group) %in% c('SampleID', 'Group'))]
rownames(func_mat) = func_group$SampleID

for(i in 1:ncol(func_mat)){
  func_mat[, i] = as.numeric(func_mat[, i])
}

func_mat[is.na(func_mat)] = 0
func_mat = func_mat[, colSums(func_mat) > 0]

meta_func = meta_group[rownames(func_mat), ]
meta_func$Group = factor(meta_func$Group, levels = group_levels)

cov_frac = coverage_sample[rownames(func_mat)]

func_mat_cov = as.data.frame(sweep(func_mat, 1, cov_frac, '/'))
func_mat_cov[!is.finite(as.matrix(func_mat_cov))] = 0

# Long table for plotting and summary
func_df_cov = as.data.frame(func_mat_cov)
func_df_cov$SampleID = rownames(func_df_cov)
func_df_cov$Group = meta_func$Group

func_long_cov = melt(func_df_cov,
                     id.vars = c('SampleID', 'Group'),
                     variable.name = 'Function',
                     value.name = 'Abundance')

func_long_cov$Abundance = as.numeric(func_long_cov$Abundance)
func_long_cov$Group = factor(func_long_cov$Group, levels = group_levels)

group_mean_cov = summaryBy(Abundance ~ Group + Function,
                           data = func_long_cov,
                           FUN = c(mean, median, sd))

# ------------------------------------------------------------
# 6. Bray-Curtis + PERMANOVA for functional structure
# ------------------------------------------------------------

func_bray_cov = vegdist(func_mat_cov, method = 'bray')

permanova_cov = adonis2(func_bray_cov ~ Group,
                        data = meta_func,
                        permutations = 999)

permanova_R2 = round(permanova_cov$R2[1], 3)
permanova_P = permanova_cov$`Pr(>F)`[1]

permanova_label = paste('PERMANOVA: R2 =', permanova_R2,
                        ', P =', signif(permanova_P, 3))

# ------------------------------------------------------------
# 7. Among-replicate dispersion
# ------------------------------------------------------------

dispersion_cov = betadisper(func_bray_cov,
                            group = meta_func$Group)

dispersion_anova = anova(dispersion_cov)
dispersion_perm = permutest(dispersion_cov, permutations = 999)

dispersion_P = dispersion_perm$tab$`Pr(>F)`[1]
dispersion_label = paste('betadisper: P =', signif(dispersion_P, 3))

# ------------------------------------------------------------
# 8. Focal guilds and Kruskal-Wallis test
# ------------------------------------------------------------

focal_levels = c('chemoheterotrophy',
                 'aerobic_chemoheterotrophy',
                 'anaerobic_chemoheterotrophy',
                 'sulfate_respiration',
                 'sulfite_respiration',
                 'respiration_of_sulfur_compounds')

focal_levels = intersect(focal_levels, colnames(func_mat_cov))

focal_mean = group_mean_cov[group_mean_cov$Function %in% focal_levels, ]

focal_mean_wide = dcast(focal_mean,
                        Function ~ Group,
                        value.var = 'Abundance.mean')

focal_kw = data.frame()

for(fn in focal_levels){

  d = data.frame(Abundance = func_mat_cov[, fn],
                 Group = meta_func$Group)

  d$Group = factor(d$Group, levels = group_levels)

  kw = kruskal.test(Abundance ~ Group, data = d)

  H = as.numeric(kw$statistic)
  k = length(unique(d$Group))
  n = nrow(d)

  epsilon_sq = (H - k + 1) / (n - k)

  temp = data.frame(Function = fn,
                    EpsilonSq = epsilon_sq,
                    Pvalue = kw$p.value,
                    stringsAsFactors = F)

  focal_kw = rbind(focal_kw, temp)
}

focal_kw$FDR = p.adjust(focal_kw$Pvalue, method = 'BH')

focal_kw$Label = paste(focal_kw$Function,
                       '\nKW FDR =', signif(focal_kw$FDR, 3),
                       sep = '')

focal_box = func_long_cov[func_long_cov$Function %in% focal_levels, ]

focal_box = merge(focal_box,
                  focal_kw[, c('Function', 'Label')],
                  by = 'Function')

label_levels = focal_kw$Label[match(focal_levels, focal_kw$Function)]

focal_box$Label = factor(focal_box$Label, levels = label_levels)
focal_box$Group = factor(focal_box$Group, levels = group_levels)

# ------------------------------------------------------------
# 9. Build panel A: annotation coverage
# ------------------------------------------------------------

pA = ggplot(coverage_df,
            aes(x = Group, y = Coverage, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, linewidth = 0.1) +
  # geom_jitter(aes(color = Group), width = 0.15,
  #             size = 1.2, alpha = 0.55) +
  scale_fill_manual(values = group_pal,
                    breaks = group_levels,
                    drop = F) +
  scale_colour_manual(values = group_pal,
                      breaks = group_levels,
                      drop = F) +
  theme_bw() +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(),
        plot.subtitle = element_text(size = 8),
        axis.text = element_text(size = 8, color = 'black'), 
        axis.title = element_text(size = 10),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(x = 'Group',
       y = 'Annotation coverage\n(fraction of reads)',
       # title = 'A. FAPROTAX annotation coverage',
       subtitle = paste('Kruskal-Wallis: P =', signif(coverage_p, 3)))

# ------------------------------------------------------------
# 10. Build panel B: coverage-normalized PCoA
# ------------------------------------------------------------

pcoa_cov = cmdscale(func_bray_cov, eig = T, k = 2)

eig_cov = pcoa_cov$eig
eig_cov = eig_cov[eig_cov > 0]

pcoa_cov_df = data.frame(SampleID = rownames(pcoa_cov$points),
                         PCoA1 = pcoa_cov$points[, 1],
                         PCoA2 = pcoa_cov$points[, 2],
                         Group = meta_func$Group)

pcoa_cov_df$Group = factor(pcoa_cov_df$Group, levels = group_levels)

pB = ggplot(pcoa_cov_df,
            aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(size = 1, alpha = 0.6) +
  stat_ellipse(level = 0.95) +
  scale_colour_manual(values = group_pal,
                      breaks = group_levels,
                      drop = F) +
  theme_bw() +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(),
        plot.subtitle = element_text(size = 8),
        axis.text = element_text(size = 8, color = 'black'), 
        axis.title = element_text(size = 10),
        legend.position = 'none') +
  labs(x = paste('PCoA1 (',
                 round(pcoa_cov$eig[1] / sum(eig_cov) * 100, 2),
                 '%)', sep = ''),
       y = paste('PCoA2\n(',
                 round(pcoa_cov$eig[2] / sum(eig_cov) * 100, 2),
                 '%)', sep = ''),
       color = 'Group',
       # title = 'B. Functional structure (coverage-normalized)',
       subtitle = permanova_label)

# ------------------------------------------------------------
# 11. Build panel C: focal guilds
# ------------------------------------------------------------

pC = ggplot(focal_box,
            aes(x = Group, y = Abundance, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, linewidth = 0.1) +
  # geom_jitter(aes(color = Group), width = 0.15,
  #             size = 0.8, alpha = 0.45) +
  facet_wrap(~ Label, scales = 'free_y', ncol = 2) +
  scale_fill_manual(values = group_pal,
                    breaks = group_levels,
                    drop = F) +
  scale_colour_manual(values = group_pal,
                      breaks = group_levels,
                      drop = F) +
  theme_bw() +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(),
        plot.subtitle = element_text(size = 8),
        axis.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8, color = 'black'), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 8, color = 'black'), 
        legend.position = 'none') + 
  labs(x = 'Group', y = 'Functional abundance (% of annotated reads)'
       # y = 'Functional abundance (% of annotated reads)',
       # title = 'C. Focal FAPROTAX guilds'
       )

# ------------------------------------------------------------
# 12. Build panel D: dispersion
# ------------------------------------------------------------

disp_cov_df = data.frame(SampleID = names(dispersion_cov$distances),
                         DistToCentroid = as.numeric(dispersion_cov$distances),
                         stringsAsFactors = F)

disp_cov_df$Group = meta_func[disp_cov_df$SampleID, 'Group']
disp_cov_df$Group = factor(disp_cov_df$Group, levels = group_levels)

pD = ggplot(disp_cov_df,
            aes(x = Group, y = DistToCentroid, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, linewidth = 0.1) +
  # geom_jitter(aes(color = Group), width = 0.15,
  #             size = 1.2, alpha = 0.55) +
  scale_fill_manual(values = group_pal,
                    breaks = group_levels,
                    drop = F) +
  scale_colour_manual(values = group_pal,
                      breaks = group_levels,
                      drop = F) +
  theme_bw() +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(),
        plot.subtitle = element_text(size = 8),
        axis.text = element_text(size = 8, color = 'black'), 
        axis.title = element_text(size = 10),
        axis.title.x = element_blank(), 
        legend.position = 'none') +
  labs(x = 'Group',
       y = 'Distance to centroid\n(Bray-Curtis)',
       # title = 'D. Among-replicate functional dispersion',
       subtitle = dispersion_label)

# ------------------------------------------------------------
# 13. Save only figures in current directory
# ------------------------------------------------------------

library(grid)
# set the layout of result
pdf(width = 8, height = 5.5, file = 'FAPROTAX_supplementary_figure_current_path.pdf' )
grid.newpage()  
pushViewport(viewport(layout = grid.layout(3, 3)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(pA, vp = vplayout(1,   1))
print(pB, vp = vplayout(2,   1))
print(pC, vp = vplayout(1:3, 2:3))
print(pD, vp = vplayout(3,   1))
dev.off() 

# combined = (pA | pB) / pC / pD +
#   plot_layout(heights = c(1, 1.6, 0.9))
# 
# combined = (pA / pB / pD) | pC +
#   plot_layout(heights = c(1, 1.6, 0.9))
# 
# ggsave('FAPROTAX_supplementary_figure_current_path.pdf',
#        combined, width = 8.2, height = 6)
