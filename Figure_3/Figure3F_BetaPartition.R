# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang
# http://webspersoais.usc.es/export9/sites/persoais/persoais/andres.baselga/pdfs/Baselga-Orme2012.pdf
# https://github.com/pnnl/brislawn-2018-founders-species
# https://cran.r-project.org/web/packages/betapart/betapart.pdf

library(betapart)
library(ggplot2)
library(ggsci)
library(grid)
library(viridis)
library(tidyverse)
library(reshape2)
library(vegan)

################################################

otuSampleDF = t(read.table('OTUtable.txt', header = T,sep = '\t',row.names = 1))
otuSampleDF = as.matrix(otuSampleDF)
otuSampleDF[otuSampleDF != 0] = 1

beta = beta.pair(otuSampleDF, index.family="jaccard")

saveRDS(beta, 'Figure3F-betaPart.rds')

################################################

distmelt = function(d){
  d.melt = d %>% as.matrix %>% as.data.frame %>% tibble::rownames_to_column() %>% melt()
  d.melt$variable = as.character(d.melt$variable)
  return(d.melt[d.melt$variable > d.melt$rowname, ])
}

# Better combined graphs
plot_dm = function(df, ylab){
  return(
    ggplot(df, aes(y = variable, x = rowname, fill = value)) + geom_raster() +
      theme(panel.background = element_rect(fill="white", color="white"),
            strip.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 10),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = c(0.98, 0.4), 
            legend.background = element_blank(), 
            legend.key.width = unit(1.5, 'mm'), 
            legend.key = element_blank(),
            legend.title = element_text(size = 6), 
            legend.text = element_text(size = 6)) +
      labs(y = ylab, fill = "Binary\nJaccard\nDistance") +
      facet_grid(var_cat ~ Type + row_cat, scales = "free", switch = "y")
  )
}

###########################################################

# Add columns listing the categories that the samples are in. Match them using grep.
fix_rows = function(df){
  df$Type = factor(df$Type, levels = c("Total", "Nestedness", "Turnover"))
  df$row_cat = 'D0'
  df$row_cat[grepl("D1", df$rowname, fixed = T)] = 'D1'
  df$row_cat[grepl("D3", df$rowname, fixed = T)] = 'D3'
  df$row_cat[grepl("D5", df$rowname, fixed = T)] = 'D5'
  df$row_cat[grepl("D6", df$rowname, fixed = T)] = 'D6'
  df$row_cat[grepl("D7", df$rowname, fixed = T)] = 'D7'
  df$row_cat = factor(df$row_cat, levels = c('D0', 'D1', 'D3','D5','D6', 'D7'))
  df$var_cat = 'D0'
  df$var_cat[grepl("D1", df$variable, fixed = T)] = 'D1'
  df$var_cat[grepl("D3", df$variable, fixed = T)] = 'D3'
  df$var_cat[grepl("D5", df$variable, fixed = T)] = 'D5'
  df$var_cat[grepl("D6", df$variable, fixed = T)] = 'D6'
  df$var_cat[grepl("D7", df$variable, fixed = T)] = 'D7'
  df$var_cat = factor(df$var_cat, levels = rev(c('D0', 'D1', 'D3','D5','D6', 'D7')))
  return(df)
}

bp = readRDS('Figure3F-betaPart.rds')
jaccard.melt = bp$beta.jac %>% distmelt()
turnover.melt = bp$beta.jtu %>% distmelt()
nestedness.melt = bp$beta.jne %>% distmelt()

# Merge for combined graph
jaccard.melt$Type = "Total"
turnover.melt$Type = "Turnover"
nestedness.melt$Type = "Nestedness"
all.melt = rbind(jaccard.melt, turnover.melt, nestedness.melt)

all.melt = all.melt %>% fix_rows

p1 = all.melt %>% plot_dm(ylab = "Dilutions") + 
  scale_fill_viridis(limits=c(0,1)) 

ggsave('Figure3F-betapart.pdf', p1, width = 6.5, height = 2.6)

###########################################################
###########################################################
###########################################################
# summary Info

bp = readRDS('Figure3F-betaPart.rds')

groupDF = read.table('OTUtable.txt', header = T, sep = '\t', stringsAsFactors = F)
group = paste0('D', substr(colnames(groupDF)[2:ncol(groupDF)], 1, 2))

aovResult = adonis(bp$beta.jac ~ group)$aov.tab
jac.r = aovResult$R2[1]
jac.p = aovResult$`Pr(>F)`[1]

aovResult = adonis(bp$beta.jne ~ group)$aov.tab
nes.r = aovResult$R2[1]
nes.p = aovResult$`Pr(>F)`[1]

aovResult = adonis(bp$beta.jtu ~ group)$aov.tab
tur.r = aovResult$R2[1]
tur.p = aovResult$`Pr(>F)`[1]

c.values = c(jac.r, jac.p, nes.r, nes.p, tur.r, tur.p)

#########################

indexs = c('jac.r', 'jac.p', 'nes.r', 'nes.p', 'tur.r', 'tur.p')
newdf = data.frame(indexs, c.values)

write.table(newdf, 'Figure3F-adnois.txt', append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)

