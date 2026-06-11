
library(reshape2)
library(ggplot2)
library(ggsci)
library(vegan)
library(sads)
library(dplyr)
library(tidyverse)

######################################

# cs = c('C2', 'C3', 'C4', 
#        'C5', 'C6', 'C7')
# cfit = c("observed", "Brokenstick", "preemption", 
#          "lognormal", "zipf", "Mandelbrot")

otudf = read.table('OTUtable.txt', header = TRUE, sep = '\t', row.names = 1)
zero.otudf = otudf[, str_detect(colnames(otudf), 'D0')]
otuSeeds = rowSums(zero.otudf) / ncol(zero.otudf)
otuSeeds = otuSeeds[otuSeeds > 0]
otuSeeds = sort(as.integer(otuSeeds), decreasing = TRUE)
otuSeeds[otuSeeds == 0] = 1 

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

#####################

# p = ggplot(rk.df, aes(x = log10(orders), y = log10(values), color = kinds)) + 
#   geom_line(linewidth = .5) + 
#   scale_color_j() +
#   scale_y_continuous(limits = c(0, 3.7)) + 
#   scale_x_continuous(limits = c(0, 3.3)) + 
#   labs(x = expression(paste('log'[10] ,'(rank)')), 
#        y = expression(paste('log'[10] ,'(abundance)')),
#        title = paste0('> The most best model: ', names(which.min(paraAICs))),
#        colour = NULL) +
#   theme( panel.background = element_rect(fill="white",color="black"),
#          panel.grid = element_blank(), 
#          plot.title = element_text(size = 10),
#          # legend.position = "right",
#          legend.position = c(.25, .25), 
#          legend.background = element_blank(), 
#          legend.text = element_text(size = 8),
#          legend.key = element_blank())
# 
# ggsave('21-RankAbundanceFit.pdf', p , width = 4, height = 4, units = 'in')

############################################################################
############################################################################
############################################################################

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

############################################################################
############################################################################
############################################################################
############################################################################
# Simulation of seven conditions

rdcount = 100

######################################################
# C1: Random select, equal probility

otudf = read.table('OTUtable.txt', header = T, sep = '\t', row.names = 1)
smpgp = read.table('samplesGroup.txt', header = T, sep = '\t', row.names = 1)

otus.rd1 = otus.rd3 = otus.rd5 = otus.rd6 = otus.rd7 = list()

zero.otudf = otudf[ , str_detect(colnames(otudf), 'D0') ]
otuSeeds = sort(rowSums(zero.otudf)/ncol(zero.otudf), decreasing = T)
# otuSeeds[otuSeeds == 0] = 1

for (i in 1:rdcount) {
  print(paste0('Simulation - ', i))
  # print(paste0('Start - ', Sys.time()))
  # Dilution -1
  otuSeeds1 = rep(names(otuSeeds), times = 10)
  set.seed(i)
  smpdf = sample_frac(data.frame(otuSeeds1),
                      size = 0.1, replace = F)
  otus.rd1[[i]] = sort(unique(smpdf$otuSeeds1))
  
  # Dilution -3
  otuSeeds1 = rep(sort(unique(smpdf$otuSeeds1)), times = 100)
  set.seed(i+rdcount*2)
  smpdf = sample_frac(data.frame(otuSeeds1),
                      size = 0.01, replace = F)
  otus.rd3[[i]] = sort(unique(smpdf$otuSeeds1))
  
  # Dilution -5
  otuSeeds1 = rep(sort(unique(smpdf$otuSeeds1)), times = 100)
  set.seed(i+rdcount*3)
  smpdf = sample_frac(data.frame(otuSeeds1),
                      size = 0.01, replace = F)
  otus.rd5[[i]] = sort(unique(smpdf$otuSeeds1))
  
  # Dilution -6
  otuSeeds1 = rep(sort(unique(smpdf$otuSeeds1)), times = 10)
  set.seed(i+rdcount*4)
  smpdf = sample_frac(data.frame(otuSeeds1),
                      size = 0.1, replace = F)
  otus.rd6[[i]] = sort(unique(smpdf$otuSeeds1))
  
  # Dilution -7
  otuSeeds1 = rep(sort(unique(smpdf$otuSeeds1)), times = 10)
  set.seed(i+rdcount*5)
  smpdf = sample_frac(data.frame(otuSeeds1),
                      size = 0.1, replace = F)
  otus.rd7[[i]] = sort(unique(smpdf$otuSeeds1))
  
  # print(paste0('Done - ', Sys.time()))
}

saveRDS(otus.rd1, 'SimulationFiles/C1-simuRandom_1.rds')
saveRDS(otus.rd3, 'SimulationFiles/C1-simuRandom_3.rds')
saveRDS(otus.rd5, 'SimulationFiles/C1-simuRandom_5.rds')
saveRDS(otus.rd6, 'SimulationFiles/C1-simuRandom_6.rds')
saveRDS(otus.rd7, 'SimulationFiles/C1-simuRandom_7.rds')

####################

for (dilu in c(1, 3, 5, 6, 7)) {
  otulist = readRDS(paste0('SimulationFiles/C1-simuRandom_', dilu, '.rds'))
  
  otus = c()
  for (i in 1:rdcount) { otus = append(otus, otulist[[i]]) }
  
  newdf = data.frame(OTUID = sort(unique(otus)))
  for (i in 1:rdcount) {
    rdf = data.frame(OTUID = otulist[[i]], cnt = 1)
    newdf = merge(newdf, rdf, by = 'OTUID', all.x = T)
  }
  colnames(newdf)[2:(rdcount + 1)] = paste0('R', 1:rdcount)
  
  write.table(newdf, paste0('SimulationFiles/C1-simuDF_', dilu, '.txt'), 
              append = F, quote = F, sep = '\t', na = '0', 
              row.names = F, col.names = T)
}

####################

# for (dilu in c(1, 3, 5, 6, 7)) {
#   dfSimu = read.table(paste0('SimulationFiles/C1-simuDF_', dilu, '.txt'), header = T, sep = '\t')
#   htdf = melt(dfSimu, measure.vars = colnames(dfSimu)[2:ncol(dfSimu)], 
#               variable.name = 'simuRd')
#   htdf$simuRd = factor(htdf$simuRd, levels = colnames(dfSimu)[2:ncol(dfSimu)])
#   
#   p = ggplot(htdf, aes(x = simuRd, y = OTUID, fill = value)) +
#     geom_tile() + 
#     scale_fill_gradient(low = "white", high = "black") + 
#     theme( panel.background = element_rect(fill="white",color="black", linewidth = 0.1),
#            panel.grid = element_blank(), 
#            plot.title = element_blank(), 
#            axis.title = element_blank(), 
#            axis.ticks = element_blank(), 
#            axis.text.y = element_blank(), 
#            axis.text.x = element_blank(), 
#            legend.position = "none")
#   
#   ggsave(paste0('SimulationFiles/C1-simuGrid_', dilu, '.jpeg'), p, 
#          width = 1, height = 1.5, units = 'in', dpi = 600)
# }


######################################################
# C2-7: "observed", "Brokenstick", "preemption", 
#       "lognormal", "zipf", "Mandelbrot"

cs = c('C2', 'C3', 'C4', 
       'C5', 'C6', 'C7')
cfit = c("observed", "Brokenstick", "preemption", 
         "lognormal", "zipf", "Mandelbrot")

for (ci in 1:6) {
  newdf = read.table('RankAbundanceFit.txt', header = T, sep = '\t')
  
  otus.rd1 = otus.rd3 = otus.rd5 = otus.rd6 = otus.rd7 = list()
  
  otuSeeds = as.numeric(newdf[ , cfit[ci]])
  names(otuSeeds) = newdf$OTUID
  
  otuSeeds.freq = otuSeeds/sum(otuSeeds)
  names(otuSeeds.freq) = names(otuSeeds)
  
  for (i in 1:rdcount) {
    print(paste0(cs[ci], ' - Simulation - ', i))
    # print(paste0('Start - ', Sys.time()))
    # Dilution -1
    otuSeeds1 = rep(names(otuSeeds), times = 10)
    otuSeeds.freq1 = rep(otuSeeds.freq, times = 10)
    set.seed(i)
    smpdf = sample_frac(data.frame(otuSeeds1), weight = otuSeeds.freq1,
                        size = 0.1, replace = F)
    otus.rd1[[i]] = sort(unique(smpdf$otuSeeds1))
    
    # Dilution -3
    otuSeeds1 = rep(sort(unique(smpdf$otuSeeds1)), times = 100)
    otuSeeds.freq1 = rep(otuSeeds.freq[sort(unique(smpdf$otuSeeds1))], times = 100)
    set.seed(i+rdcount*2)
    smpdf = sample_frac(data.frame(otuSeeds1), weight = otuSeeds.freq1,
                        size = 0.01, replace = F)
    otus.rd3[[i]] = sort(unique(smpdf$otuSeeds1))
    
    # Dilution -5
    otuSeeds1 = rep(sort(unique(smpdf$otuSeeds1)), times = 100)
    otuSeeds.freq1 = rep(otuSeeds.freq[sort(unique(smpdf$otuSeeds1))], times = 100)
    set.seed(i+rdcount*3)
    smpdf = sample_frac(data.frame(otuSeeds1), weight = otuSeeds.freq1,
                        size = 0.01, replace = F)
    otus.rd5[[i]] = sort(unique(smpdf$otuSeeds1))
    
    # Dilution -6
    otuSeeds1 = rep(sort(unique(smpdf$otuSeeds1)), times = 10)
    otuSeeds.freq1 = rep(otuSeeds.freq[sort(unique(smpdf$otuSeeds1))], times = 10)
    set.seed(i+rdcount*4)
    smpdf = sample_frac(data.frame(otuSeeds1), weight = otuSeeds.freq1,
                        size = 0.1, replace = F)
    otus.rd6[[i]] = sort(unique(smpdf$otuSeeds1))
    
    # Dilution -7
    otuSeeds1 = rep(sort(unique(smpdf$otuSeeds1)), times = 10)
    otuSeeds.freq1 = rep(otuSeeds.freq[sort(unique(smpdf$otuSeeds1))], times = 10)
    set.seed(i+rdcount*5)
    smpdf = sample_frac(data.frame(otuSeeds1), weight = otuSeeds.freq1,
                        size = 0.1, replace = F)
    otus.rd7[[i]] = sort(unique(smpdf$otuSeeds1))
    
    # print(paste0('Done - ', Sys.time()))
  }
  
  saveRDS(otus.rd1, paste0('SimulationFiles/', cs[ci], '-simuRandom_1.rds'))
  saveRDS(otus.rd3, paste0('SimulationFiles/', cs[ci], '-simuRandom_3.rds'))
  saveRDS(otus.rd5, paste0('SimulationFiles/', cs[ci], '-simuRandom_5.rds'))
  saveRDS(otus.rd6, paste0('SimulationFiles/', cs[ci], '-simuRandom_6.rds'))
  saveRDS(otus.rd7, paste0('SimulationFiles/', cs[ci], '-simuRandom_7.rds'))
  
  ####################
  
  for (dilu in c(1, 3, 5, 6, 7)) {
    otulist = readRDS(paste0('SimulationFiles/', cs[ci],'-simuRandom_', dilu, '.rds'))
    otus = c()
    for (i in 1:rdcount) { otus = append(otus, otulist[[i]]) }
    newdf = data.frame(OTUID = sort(unique(otus)))
    for (i in 1:rdcount) {
      rdf = data.frame(OTUID = otulist[[i]], cnt = 1)
      newdf = merge(newdf, rdf, by = 'OTUID', all.x = T)
    }
    colnames(newdf)[2:(rdcount + 1)] = paste0('R', 1:rdcount)
    
    write.table(newdf, paste0('SimulationFiles/', cs[ci],'-simuDF_', dilu, '.txt'), 
                append = F, quote = F, sep = '\t', na = '0', 
                row.names = F, col.names = T)
  }
  
  ####################
#   
#   totaldf = data.frame(OTUID = names(otuSeeds))
#   
#   for (dilu in c(1, 3, 5, 6, 7)) {
#     dfSimu = read.table(paste0('SimulationFiles/', cs[ci],'-simuDF_', dilu, '.txt'), 
#                         header = T, sep = '\t')
#     dfSimu = merge(totaldf, dfSimu, by = 'OTUID', all.x = T)
#     dfSimu$OTUID = factor(dfSimu$OTUID, names(otuSeeds))
#     
#     htdf = melt(dfSimu, measure.vars = colnames(dfSimu)[2:ncol(dfSimu)],
#                 variable.name = 'simuRd')
#     htdf$simuRd = factor(htdf$simuRd, levels = colnames(dfSimu)[2:ncol(dfSimu)])
#     htdf$value[is.na(htdf$value)] = 0
#     
#     p = ggplot(htdf, aes(x = simuRd, y = OTUID, fill = value)) +
#       geom_tile() +
#       scale_fill_gradient(low = "white", high = "black") +
#       theme( panel.background = element_rect(fill="white",color="black", linewidth = 0.1),
#              panel.grid = element_blank(),
#              plot.title = element_blank(),
#              axis.title = element_blank(),
#              axis.ticks = element_blank(),
#              axis.text.y = element_blank(),
#              axis.text.x = element_blank(),
#              legend.position = "none")
#     
#     ggsave(paste0('SimulationFiles/',  cs[ci],'-simuGrid_', dilu, '.jpeg'), p,
#            width = 1, height = 1.5, units = 'in', dpi = 600) 
#   }  
# }


############################################################################
############################################################################
############################################################################
# plot Observed 

# rankFitDf = read.table('22-RankAbundanceFit.txt', header = T, sep = '\t')
# smpgp = read.table('samplesGroup.txt', header = T, sep = '\t')
# otudf = read.table('OTUtable.txt', header = T, sep = '\t', row.names = 1)
# otudf = otudf[rankFitDf$OTUID, ]
# otudf.m = as.matrix(otudf)
# otudf.m[otudf.m != 0] = 1
# 
# for (dl in c('D1', 'D3', 'D5', 'D6', 'D7')) {
#   diluDf = otudf.m[ , smpgp$SampleID[smpgp$Group == dl]]
#   diluDf = data.frame(OTUID = rownames(diluDf), diluDf)
#   diluDf$OTUID = factor(diluDf$OTUID, levels = diluDf$OTUID)
#   
#   htdf = melt(diluDf, measure.vars = colnames(diluDf)[2:ncol(diluDf)],
#               variable.name = 'smpid')
#   htdf$smpid = factor(htdf$smpid, levels = colnames(diluDf)[2:ncol(diluDf)])
#   htdf$value[is.na(htdf$value)] = 0
#   
#   p = ggplot(htdf, aes(x = smpid, y = OTUID, fill = value)) +
#     geom_tile() +
#     scale_fill_gradient(low = "white", high = "black") +
#     theme( panel.background = element_rect(fill="white",color="black", linewidth = 0.1),
#            panel.grid = element_blank(),
#            plot.title = element_blank(),
#            axis.title = element_blank(),
#            axis.ticks = element_blank(),
#            axis.text.y = element_blank(),
#            axis.text.x = element_blank(),
#            legend.position = "none")
#   
#   ggsave(paste0('SimulationFiles/', dl,'-Grid.jpeg'), p, 
#          # width = 0.3, height = 1.5, units = 'in', dpi = 600)
#          width = 2*( (ncol(diluDf)-1)/100 ), height = 1.5, units = 'in', dpi = 600)
# }

######################################


######################################




######################################


######################################































