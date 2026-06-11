
# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(vegan)
library(ggsci)
library(ggplot2)
library(tidyverse)

#############################################################
# Supplementary Fig. S2E-F
# E: Bray-Curtis NMDS
# F: CLR/Aitchison PCoA
#############################################################

set.seed(123)

#############################################################
# 1. Data input

dclassDF = read.table('OTUsamplesGroup.txt', header = T,
                      sep = '\t', stringsAsFactors = F)

otuSampleDF = t(read.table('OTUtable.txt', header = T,
                           sep = '\t', row.names = 1,
                           check.names = FALSE))

# Match sample order between OTU table and metadata
common.samples = rownames(otuSampleDF)[rownames(otuSampleDF) %in% dclassDF$SampleID]

otuSampleDF = otuSampleDF[common.samples, , drop = FALSE]
dclassDF = dclassDF[match(common.samples, dclassDF$SampleID), , drop = FALSE]

# Remove all-zero samples and OTUs, if any
otuSampleDF = otuSampleDF[rowSums(otuSampleDF) > 0, , drop = FALSE]
dclassDF = dclassDF[match(rownames(otuSampleDF), dclassDF$SampleID), , drop = FALSE]

otuSampleDF = otuSampleDF[, colSums(otuSampleDF) > 0, drop = FALSE]

dclassDF$Group = factor(dclassDF$Group,
                        levels = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7'))

#############################################################
# 2. Common theme, following Figure3D_NMDS_AbsT.R

mycols = pal_nejm()(6)

# mytheme = theme(panel.background = element_rect(fill = "white", color = "black"),
#                 panel.grid = element_blank(),
#                 plot.caption = element_text(size = 10, color = 'black'),
#                 axis.title = element_text(size = 10, color = 'black'),
#                 axis.text.y = element_text(size = 8, hjust = 0.5, vjust = 0.5,
#                                            angle = 90, color = 'black'),
#                 axis.text.x = element_text(size = 8, color = 'black'),
#                 legend.position = "right",
#                 legend.title = element_blank(),
#                 legend.text = element_text(size = 8),
#                 legend.key = element_blank())

#############################################################
# 3. Helper functions for statistics

get_global_stats = function(dis, meta, distance.name) {
  
  ad = adonis2(dis ~ Group, data = meta, permutations = 999)
  an = anosim(dis, grouping = meta$Group, permutations = 999)
  mr = mrpp(dis, grouping = meta$Group, permutations = 999)
  
  bd = betadisper(dis, group = meta$Group)
  bd.test = permutest(bd, permutations = 999)
  
  out = data.frame(
    Distance = distance.name,
    Test = c('PERMANOVA', 'ANOSIM', 'MRPP', 'Beta-dispersion'),
    Statistic = c(ad$F[1],
                  an$statistic,
                  mr$A,
                  bd.test$tab[1, 'F']),
    R2 = c(ad$R2[1],
           NA,
           NA,
           NA),
    P = c(ad$`Pr(>F)`[1],
          an$signif,
          mr$Pvalue,
          bd.test$tab[1, 'Pr(>F)']),
    stringsAsFactors = FALSE
  )
  
  return(out)
}

pairwise_permanova = function(dis, meta, distance.name) {
  
  dis.mat = as.matrix(dis)
  gps = levels(meta$Group)
  pair.list = combn(gps, 2, simplify = FALSE)
  
  out = map_df(pair.list, function(x) {
    
    idx = meta$Group %in% x
    sub.dis = as.dist(dis.mat[idx, idx])
    sub.meta = droplevels(meta[idx, , drop = FALSE])
    
    ad = adonis2(sub.dis ~ Group, data = sub.meta, permutations = 999)
    
    data.frame(
      Distance = distance.name,
      Group1 = x[1],
      Group2 = x[2],
      F = ad$F[1],
      R2 = ad$R2[1],
      P = ad$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    )
  })
  
  out$P.adj.BH = p.adjust(out$P, method = 'BH')
  
  return(out)
}

#############################################################
# 4. Fig. S2E: Bray-Curtis NMDS
# Abundance-sensitive ordination

otuRelDF = decostand(otuSampleDF, method = 'total')

bray.dis = vegdist(otuRelDF, method = 'bray')

nmds.bray = metaMDS(bray.dis,
                    k = 2,
                    trymax = 100,
                    wascores = TRUE)

bray.scores = as.data.frame(scores(nmds.bray))
bray.df = data.frame(bray.scores,
                     dclass = dclassDF$Group,
                     SampleID = dclassDF$SampleID)

head(bray.df)

bray.p = ggplot(bray.df, aes(x = NMDS1, y = NMDS2, color = dclass)) +
  geom_point(size = 1.5, alpha = 0.8, shape = 16) +
  labs(caption = paste0('Stress: ',
                        paste(round(nmds.bray$stress, digits = 2), collapse = ', '))) +
  scale_colour_manual(values = mycols,
                      breaks = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')) +
  geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank(),
        plot.caption = element_text(size = 8, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        axis.text.y = element_text(size = 8, hjust = 0.5, vjust = 0.5,
                                   angle = 90, color = 'black'),
        axis.text.x = element_text(size = 8, color = 'black'),
        legend.position = "none" )

ggsave('Sup_Fig_S3_B_NMDS_BrayCurtis.pdf',
       bray.p, units = 'in', width = 2.5, height = 2.5)

stressplot

#############################################################
# 5. Fig. S2F: CLR/Aitchison PCoA
# Compositional-aware ordination
#
# Aitchison distance = Euclidean distance on CLR-transformed abundance.
# Zero values are handled by adding a pseudocount before CLR transformation.
#############################################################

pseudo.count = 0.5

clr_fun = function(x) {
  x = x + pseudo.count
  x = x / sum(x)
  lx = log(x)
  return(lx - mean(lx))
}

otuCLRDF = t(apply(otuSampleDF, 1, clr_fun))

ait.dis = dist(otuCLRDF, method = 'euclidean')

pcoa.ait = cmdscale(ait.dis,
                    k = 2,
                    eig = TRUE)

ait.scores = as.data.frame(pcoa.ait$points)
colnames(ait.scores) = c('PCoA1', 'PCoA2')

eig = pcoa.ait$eig
eig.pos = eig[eig > 0]
eig.percent = round(eig.pos / sum(eig.pos) * 100, 1)

ait.df = data.frame(ait.scores,
                    dclass = dclassDF$Group,
                    SampleID = dclassDF$SampleID)

head(ait.df)

ait.p = ggplot(ait.df, aes(x = PCoA1, y = PCoA2, color = dclass)) +
  geom_point(size = 1.5, alpha = 0.8, shape = 16) +
  labs(x = paste0('PCoA1 (', eig.percent[1], '%)'),
       y = paste0('PCoA2 (', eig.percent[2], '%)'),
       caption = paste0('CLR/Aitchison PCoA; pseudocount = ', pseudo.count)) +
  scale_colour_manual(values = mycols,
                      breaks = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')) +
  geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank(),
        plot.caption = element_text(size = 8, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        axis.text.y = element_text(size = 8, hjust = 0.5, vjust = 0.5,
                                   angle = 90, color = 'black'),
        axis.text.x = element_text(size = 8, color = 'black'),
        legend.position = "none" )

ggsave('Sup_Fig_S3_C_PCoA_CLR_Aitchison.pdf',
       ait.p, units = 'in', width = 2.5, height = 2.5)

#############################################################
# 6. Global tests and pairwise PERMANOVA
# These results can be added to Table S11.

global.stats = bind_rows(
  get_global_stats(bray.dis, dclassDF, 'Bray-Curtis'),
  get_global_stats(ait.dis, dclassDF, 'CLR/Aitchison')
)

pairwise.stats = bind_rows(
  pairwise_permanova(bray.dis, dclassDF, 'Bray-Curtis'),
  pairwise_permanova(ait.dis, dclassDF, 'CLR/Aitchison')
)

ordination.info = data.frame(
  Distance = c('Bray-Curtis', 'CLR/Aitchison'),
  Ordination = c('NMDS', 'PCoA'),
  Diagnostic1 = c(paste0('Stress = ', round(nmds.bray$stress, 4)),
                  paste0('PCoA1 = ', eig.percent[1], '%')),
  Diagnostic2 = c('trymax = 100',
                  paste0('PCoA2 = ', eig.percent[2], '%')),
  ZeroHandling = c('Not applicable',
                   paste0('Pseudocount = ', pseudo.count, ' before CLR')),
  stringsAsFactors = FALSE
)

write.csv(global.stats,
          'Sup_Fig_S3_GlobalStats.csv',
          row.names = FALSE)

write.csv(pairwise.stats,
          'Sup_Fig_S3_PairwisePERMANOVA.csv',
          row.names = FALSE)

write.csv(ordination.info,
          'Sup_Fig_S3_OrdinationInfo.csv',
          row.names = FALSE)

#############################################################
# 7. Print summary

print(global.stats)
print(pairwise.stats)
print(ordination.info)

#############################################################
