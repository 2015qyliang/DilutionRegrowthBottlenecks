

library(ggplot2)
library(ggsci)
library(tidyverse)
# library(reshape2)

##################################################################

reabun = read.table('OTUtable.txt', header = T, sep = '\t')
Group = otus = fV = fM = fN = rChisqX = rChisqP = rKsD = rKsP = c()
for (otu in sort(unique(reabun$OTUID))) {
  dldfotu = reabun[reabun$OTUID == otu, ]
  for (dl in c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')) {
    dldf = dldfotu[ , str_detect(colnames(dldfotu), dl)]
    if (sum(dldf == 0) != ncol(dldf)) {
      fV = append(fV, var(as.numeric(dldf)))
      fM = append(fM, mean(as.numeric(dldf)))
      fN = append(fN, sum(as.numeric(dldf) != 0))
      Group = append(Group, dl)
      otus = append(otus, otu)
      rChisqX = append(rChisqX, chisq.test(as.numeric(dldf))$statistic)
      rChisqP = append(rChisqP, chisq.test(as.numeric(dldf))$p.value)
      rKsD = append(rKsD, ks.test(as.numeric(dldf), "pnorm", 
                                  mean(as.numeric(dldf)), 
                                  sd(as.numeric(dldf)))$statistic)
      rKsP = append(rKsP, ks.test(as.numeric(dldf), "pnorm", 
                                  mean(as.numeric(dldf)), 
                                  sd(as.numeric(dldf)))$p.value)
    }
  }
}

IoDs = data.frame(Group, OTU = otus, fV, fM, fN, rChisqX, rChisqP, rKsD, rKsP)
write.table(IoDs[order(IoDs$Group), ], 'IndexOfDispersion.txt',
            append = F, quote = F, sep = '\t', row.names = F, col.names = T)

####################################

colnames(IoDs)
IoDs$iodValue = log10((IoDs$fV/IoDs$fM)*IoDs$fN)
IoDs$ChisqRes = 'Uniform'
IoDs$ChisqRes[IoDs$rChisqP <= 0.05] = 'NonUniform'

p = ggplot(IoDs, aes(x = fN, y = iodValue, color = ChisqRes)) + 
  geom_point(shape = 16, alpha = 0.2, size = 0.5) + 
  scale_color_manual(values = c('#f36e42', '#2e8abf'), 
                     breaks = c('Uniform', 'NonUniform')) + 
  labs(x = 'Occurrence (Times)', y = 'Log10 (Index\nof dispersion)') + 
  facet_grid(.~Group, scales = 'free', space = 'free_y') + 
  theme_bw() + 
  theme(#plot.margin = unit(c(0.3,0.2,0,0.2),"in"),
        plot.background = element_blank(), 
        panel.grid = element_blank(), 
        # panel.grid.minor = element_blank(), 
        strip.background = element_blank(), 
        strip.text.x = element_text(size = 10, color = 'black'),
        legend.position = 'bottom', 
        axis.title = element_text(size = 10, color = 'black'),
        axis.text = element_text(size = 8, color = 'black'))
  
ggsave('IndexOfDispersionChisq.pdf', p, width = 6, height = 2)

######################################

IoDs$KsRes = 'Norm'
IoDs$KsRes[IoDs$rKsP <= 0.05] = 'NotNorm'

p = ggplot(IoDs, aes(x = fN, y = iodValue, color = KsRes)) + 
  geom_point(shape = 16, alpha =  0.2, size = 0.5) + 
  scale_color_manual(values = c('#f36e42', '#2e8abf'), 
                     breaks = c('Norm', 'NotNorm')) + 
  labs(x = 'Occurrence (Times)', y = 'Log10 (Index\nof dispersion)') + 
  facet_grid(.~Group, scales = 'free', space = 'free_y') + 
  theme_bw() + 
  theme(#plot.margin = unit(c(0.3,0.2,0,0.2),"in"),
        plot.background = element_blank(), 
        panel.grid = element_blank(), 
        # panel.grid.minor = element_blank(), 
        strip.background = element_blank(), 
        strip.text.x = element_text(size = 10, color = 'black'),
        legend.position = 'bottom', 
        axis.title = element_text(size = 10, color = 'black'),
        axis.text = element_text(size = 8, color = 'black'))

ggsave('IndexOfDispersionNorm.pdf', p, width = 6, height = 2)


# Pearson's Chi-squared Test
# Kolmogorov-Smirnov Tests

######################################################################################
######################################################################################







######################################################################################
######################################################################################




