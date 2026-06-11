# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(vegan)
library(ggsci)
library(ggplot2)
library(tidyverse)

#############################################################

dclassDF = read.table('OTUsamplesGroup.txt', header = T, 
                      sep = '\t', stringsAsFactors = F)

otuSampleDF = t(read.table('OTUtable.txt', header = T,sep = '\t',row.names = 1))
 
nmds = metaMDS(vegdist(otuSampleDF, method = 'jaccard', binary = T),  
               k = 2, trymax = 100, wascores = TRUE)
# data.scores = as.data.frame(nmds$points)
data.scores = as.data.frame(scores(nmds))
newdf = data.frame(data.scores, dclass = dclassDF$Group )

head(newdf)

nmds.p = ggplot(newdf, aes(x = NMDS1, y = NMDS2, color = dclass)) +
  geom_point(size = 1.5, alpha = 0.8, shape = 16) +
  labs(caption = paste0('Stress: ', 
                        paste(round(nmds$stress, digits = 2), collapse = ', '))) +
  scale_colour_manual(values = pal_nejm()(6),
                      breaks = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')) +
  geom_hline(yintercept=0, colour="grey60", linetype="dashed")+
  geom_vline(xintercept=0, colour="grey60", linetype="dashed")+
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        plot.caption =  element_text(size = 10, color = 'black'), 
        axis.title = element_text(size = 10, color = 'black'),
        axis.text.y = element_text(size = 8, hjust = 0.5, vjust = 0.5, 
                                   angle = 90, color = 'black'),
        axis.text.x = element_text(size = 8, color = 'black'),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key = element_blank())
ggsave('Figure3D_NMDS.pdf',  nmds.p, units = 'in', width = 3, height = 2.5)

#############################################################






