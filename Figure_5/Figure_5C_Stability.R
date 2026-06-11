# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(ggplot2)
library(ggsci)
library(igraph)
library(brainGraph)
library(tidyverse)
library(NetSwan)
library(randomForest)
library(doBy)
library(scales)

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

# simulation of net_knock

#consider cascade effects
rand.remov1.once<-function(netRaw, rm.num, keystonelist, sp.ra, abundance.weighted=T){
  rm.num2<-ifelse(rm.num > length(keystonelist), length(keystonelist), rm.num)
  id.rm<-sample(keystonelist, rm.num2) 
  net.Raw=netRaw #don't want change netRaw
  
  net.new=net.Raw[!names(sp.ra) %in% id.rm, !names(sp.ra) %in% id.rm]   ##remove all the links to these species
  if (nrow(net.new)<2){
    0
  } else {
    sp.ra.new=sp.ra[!names(sp.ra) %in% id.rm]
    
    if (abundance.weighted){
      net.stength= net.new*sp.ra.new
    } else {
      net.stength= net.new
    }
    
    sp.meanInteration<-colMeans(net.stength)
    
    
    while ( length(sp.meanInteration)>1 & min(sp.meanInteration) <=0){
      id.remain<- which(sp.meanInteration>0) 
      net.new=net.new[id.remain,id.remain]
      sp.ra.new=sp.ra.new[id.remain]
      
      if (abundance.weighted){
        net.stength= net.new*sp.ra.new
      } else {
        net.stength= net.new
      }
      
      if (length(net.stength)>1){
        sp.meanInteration<-colMeans(net.stength)
      } else{
        sp.meanInteration<-0
      }
      
    }
    
    remain.percent<-length(sp.ra.new)/length(sp.ra)
    
    remain.percent}
}

#not consider cascade effects
rand.remov1.once_2extinctOnly<-function(netRaw, rm.num, keystonelist, sp.ra, abundance.weighted=T){
  rm.num2<-ifelse(rm.num > length(keystonelist), length(keystonelist), rm.num)
  id.rm<-sample(keystonelist, rm.num2)
  net.Raw=netRaw #don't want change netRaw
  
  net.new=net.Raw[!names(sp.ra) %in% id.rm, !names(sp.ra) %in% id.rm]   ##remove all the links to these species
  if (nrow(net.new)<2){
    0
  } else {
    sp.ra.new=sp.ra[!names(sp.ra) %in% id.rm]
    
    if (abundance.weighted){
      net.stength= net.new*sp.ra.new
    } else {
      net.stength= net.new
    }
    
    sp.meanInteration<-colMeans(net.stength)
    
    id.remain<- which(sp.meanInteration>0) 
    sp.ra.new=sp.ra.new[id.remain]
    
    remain.percent<-length(sp.ra.new)/length(sp.ra)
    
    remain.percent}
}

#rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu1<-function(netRaw, rm.p.list, keystonelist,sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov1.once(netRaw=netRaw, rm.num=x, keystonelist=keystonelist, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}
rmsimu1_2extinctOnly<-function(netRaw, rm.p.list, keystonelist,sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov1.once_2extinctOnly(netRaw=netRaw, rm.num=x, keystonelist=keystonelist, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}



#####################random deletion######################

#consider cascade effects: removed species will further influence the remaining nodes

rand.remov2.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  
  net.new=net.Raw[-id.rm, -id.rm]   ##remove all the links to these species
  if (nrow(net.new)<2){
    0
  } else {
    sp.ra.new=sp.ra[-id.rm]
    
    if (abundance.weighted){
      net.stength= net.new*sp.ra.new
    } else {
      net.stength= net.new
    }
    
    sp.meanInteration<-colMeans(net.stength)
    
    
    while ( length(sp.meanInteration)>1 & min(sp.meanInteration) <=0){
      id.remain<- which(sp.meanInteration>0) 
      net.new=net.new[id.remain,id.remain]
      sp.ra.new=sp.ra.new[id.remain]
      
      if (abundance.weighted){
        net.stength= net.new*sp.ra.new
      } else {
        net.stength= net.new
      }
      
      if (length(net.stength)>1){
        sp.meanInteration<-colMeans(net.stength)
      } else{
        sp.meanInteration<-0
      }
      
    }
    
    remain.percent<-length(sp.ra.new)/length(sp.ra)
    
    remain.percent}
}

rand.remov2.once_2extinctOnly<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  
  net.new=net.Raw[-id.rm, -id.rm]   ##remove all the links to these species
  if (nrow(net.new)<2){
    0
  } else {
    sp.ra.new=sp.ra[-id.rm]
    
    if (abundance.weighted){
      net.stength= net.new*sp.ra.new
    } else {
      net.stength= net.new
    }
    
    sp.meanInteration<-colMeans(net.stength)
    id.remain<- which(sp.meanInteration>0) 
    sp.ra.new=sp.ra.new[id.remain]
    remain.percent<-length(sp.ra.new)/length(sp.ra)
    remain.percent}
}


#rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov2.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}
rmsimu_2extinctOnly<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov2.once_2extinctOnly(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}

##########################################################################

nts = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')

for (netTag in nts) {
  corrFile = paste0('KnockFns/', netTag, '_Netcorr_full.txt')
  otuFile = paste0('KnockFns/', netTag, '_Netotutable.txt')
  modFile = paste0('KnockFns/', netTag, '_NetModule.txt')
  
  cormatrix<-read.csv(corrFile,sep="\t",header = F)
  otutab<-read.csv(otuFile,sep="\t",header = F,row.names = 1)
  
  otutab[is.na(otutab)]<-0
  comm<-t(otutab)
  comm<-comm/rowSums(comm)   
  
  sp.ra<-colMeans(comm)  #relative abundance of each species
  
  
  row.names(cormatrix)<-colnames(cormatrix)<-colnames(comm)
  
  
  cormatrix2<-cormatrix*(abs(cormatrix)>=0.80)  #only keep links above the cutoff point
  cormatrix2[is.na(cormatrix2)]<-0
  diag(cormatrix2)<-0    #no links for self-self    
  sum(abs(cormatrix2)>0)/2  #this should be the number of links. 
  
  sum(colSums(abs(cormatrix2))>0)  #?? species have at least one linkage with others.
  
  network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
  sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
  sum(row.names(network.raw)==names(sp.ra2))  #check if matched
  
  #input network matrix, percentage of randomly removed species, and ra of all species
  #return the proportion of species remained
  
  node.attri<-read.csv(modFile,sep="\t",skip = 2,row.names = 1)
  module.hub<-as.character(node.attri$ID[node.attri$Zi > 2.5 & node.attri$Pi <= 0.62])
  length(module.hub)
  
  
  print(paste0(netTag, ' --- 1 ---'))
  
  Weighted.simu<-rmsimu1(netRaw=network.raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
  Unweighted.simu<-rmsimu1(netRaw=network.raw, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)
  
  Weighted.simu_2extinctOnly<-rmsimu1_2extinctOnly(netRaw=network.raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
  Unweighted.simu_2extinctOnly<-rmsimu1_2extinctOnly(netRaw=network.raw, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)
  
  dat1<-data.frame(Number.hub.removed=rep(1:length(module.hub),4),rbind(Weighted.simu,Unweighted.simu,Weighted.simu_2extinctOnly,Unweighted.simu_2extinctOnly),
                   weighted=rep(c("weighted","unweighted","weighted","unweighted"),each=length(module.hub)),consider_cascade=rep(c("Yes","No"),each=2*length(module.hub)),
                   network=rep(netTag,4*length(module.hub)))
  
  currentdat_target<-dat1
  currentdat_target<-rbind(dat1,currentdat_target)
  write.table(currentdat_target, paste0('KnockRes/', netTag, "_target_deletion.txt"), 
              append = F, quote = F, sep = '\t', row.names = F, col.names = T)
  
  print(paste0(netTag, ' --- 2 ---'))
  
  Weighted.simu2<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
  Unweighted.simu2<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)
  
  Weighted.simu2_2extinctOnly<-rmsimu_2extinctOnly(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
  Unweighted.simu2_2extinctOnly<-rmsimu_2extinctOnly(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)
  
  
  dat2<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),4),rbind(Weighted.simu2,Unweighted.simu2,Weighted.simu2_2extinctOnly,Unweighted.simu2_2extinctOnly),
                   weighted=rep(c("weighted","unweighted","weighted","unweighted"),each=20),
                   consider_cascade=rep(c("Yes","No"),each=40),
                   network=rep(netTag,80))
  
  #random deletion
  currentdat<-dat2
  currentdat<-rbind(dat2,currentdat)
  write.table(currentdat, paste0('KnockRes/', netTag, "_random_deletion.txt"), 
              append = F, quote = F, sep = '\t', row.names = F, col.names = T)
  
}

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

# output the data of Figure_5C i

dfFiltGet = function(fn) {
  df.r = read.table(file.path('KnockRes', fn), header = T,
                    sep = '\t', stringsAsFactors = F)
  df = df.r[df.r$consider_cascade == 'No' & df.r$weighted == 'unweighted',]
  df = df[df$Proportion.removed == 0.5, ]
  return(df[1, ])
}


fns.random = list.files('KnockRes/', pattern = '_random_deletion.txt')
fns.target = list.files('KnockRes/', pattern = '_target_deletion.txt')


randomDF = data.frame()
for (fn in fns.random) {
  randomDF = rbind(randomDF, dfFiltGet(fn))
}
write.table(randomDF, 'Figure_5C-simuRandomDelt.txt', append = F, quote = F, 
            sep = '\t', row.names = F, col.names = T)

# ————————————————————————————————————

dfFiltGet = function(fn) {
  df.r = read.table(file.path('KnockRes', fn), header = T,
                    sep = '\t', stringsAsFactors = F)
  df = df.r[df.r$consider_cascade == 'No' & df.r$weighted == 'unweighted',]
  df = df[df$Number.hub.removed == 5, ]
  return(df[1, ])
}

randomDF = data.frame()
for (fn in fns.target) {
  randomDF = rbind(randomDF, dfFiltGet(fn))
}
write.table(randomDF, 'Figure_5C-simuTargetDelt.txt', append = F, quote = F, 
            sep = '\t', row.names = F, col.names = T)

#######################################################################################################

# output the data of Figure_5C ii

getMinE = function(fn){
  nodes = read.table(paste0('./NetInfo/', fn,' node_attribute.txt'), 
                     header = T, sep = '\t', stringsAsFactors = F)
  nodes$No..module = paste0('M_', nodes$No..module)
  rownames(nodes) = nodes$id
  links = read.table(paste0('./NetInfo/', fn,' edge_attribute.txt'),
                     header = T, sep = '\t', stringsAsFactors = F)
  net.graph = graph_from_data_frame( d = links, vertices = nodes, directed = F)
  
  e.global = efficiency(net.graph, type = 'global')
  
  e.nodes = c()
  for (nd in nodes$id) {
    nodes.nd = nodes[nodes$id != nd, ]
    links.nd = links[(links$from != nd) & (links$to != nd) , ]
    net.graph.nd = graph_from_data_frame( d = links.nd, vertices = nodes.nd, directed = F)
    e.nd = efficiency(net.graph.nd, type = 'global')
    e.nodes = append(e.nodes, e.nd)
  }
  
  return(c(e.global, round(max((e.global-e.nodes)/e.global), digits = 4)))
}

# ————————————————————————————————————————

sampleGroups = read.table('RMT_value.txt', header = T, sep = '\t', stringsAsFactors = F)
vulnerability = c()
efficient = c()
for (i in 1:length(sampleGroups$groups)) {
  fn = paste0(sampleGroups$groups[i], ' ', sprintf('%0.2f', sampleGroups$value[i]))
  efficient = append(efficient, getMinE(fn)[1])
  vulnerability = append(vulnerability, getMinE(fn)[2]) 
}

vulDF = data.frame(Group = sampleGroups$groups, efficient, vulnerability)

write.table(vulDF, 'Figure_5C-EffiVulnDF.txt', append = F, sep = '\t', row.names = F, col.names = T)

#######################################################################################################

# Figure_5C i

randomDF = read.table('Figure_5C-simuRandomDelt.txt', header = T, sep = '\t', stringsAsFactors = F)

# Kolmogorov-Smirnov test
com.pair = c("D0", "D1", "D3", "D5", "D6", "D7")
cp = statisticD = Pvalue = c()
for (i in 1:5) {
  set.seed(618)
  c = randomDF[randomDF$network == com.pair[i], 'remain.mean']
  m = randomDF[randomDF$network == com.pair[i+1], 'remain.mean']
  sd.c = randomDF[randomDF$network == com.pair[i], 'remain.sd']
  sd.m = randomDF[randomDF$network == com.pair[i+1], 'remain.sd']
  kr = ks.test(rnorm(100, mean = c, sd = sd.c), rnorm(100, mean = m, sd = sd.m))
  cp = append(cp, paste0(com.pair[i], '_vs_', com.pair[i+1]))
  statisticD = append(statisticD, kr$statistic)
  Pvalue = append(Pvalue, kr$p.value)
}

write.table(data.frame(ComparePairs = cp, statisticD, Pvalue), 'Figure_5C-KStestRandom.txt', 
            append = F, quote = F, sep = '\t', row.names = F, col.names = T)

# ##########################
# # Kolmogorov-Smirnov test
# com.pair = c("D0", "D1", "D3", "D5", "D6", "D7")
# cp = statisticD = Pvalue = c()
# for (i in 1:5) {
#   set.seed(618)
#   c = randomDF[randomDF$network == com.pair[i], 'remain.mean']
#   m = randomDF[randomDF$network == com.pair[i+1], 'remain.mean']
#   sd.c = randomDF[randomDF$network == com.pair[i], 'remain.sd']
#   sd.m = randomDF[randomDF$network == com.pair[i+1], 'remain.sd']
#   kr = wilcox.test(rnorm(100, mean = c, sd = sd.c), rnorm(100, mean = m, sd = sd.m))
#   cp = append(cp, paste0(com.pair[i], '_vs_', com.pair[i+1]))
#   statisticD = append(statisticD, kr$statistic)
#   Pvalue = append(Pvalue, kr$p.value)
# }
# 
# write.table(data.frame(ComparePairs = cp, statisticD, Pvalue), 'Figure_5C-WMtestRandom.txt', 
#             append = F, quote = F, sep = '\t', row.names = F, col.names = T)
# 
# # ————————————————————————————————————————
# 
# randomDF = read.table('Figure_5C-simuTargetDelt.txt', header = T, sep = '\t', stringsAsFactors = F)
# 
# # Kolmogorov-Smirnov test
# com.pair = c("D0", "D1", "D3", "D5", "D6", "D7")
# cp = statisticD = Pvalue = c()
# for (i in 1:5) {
#   set.seed(618)
#   c = randomDF[randomDF$network == com.pair[i], 'remain.mean']
#   m = randomDF[randomDF$network == com.pair[i+1], 'remain.mean']
#   sd.c = randomDF[randomDF$network == com.pair[i], 'remain.sd']
#   sd.m = randomDF[randomDF$network == com.pair[i+1], 'remain.sd']
#   kr = ks.test(rnorm(100, mean = c, sd = sd.c), rnorm(100, mean = m, sd = sd.m))
#   cp = append(cp, paste0(com.pair[i], '_vs_', com.pair[i+1]))
#   statisticD = append(statisticD, kr$statistic)
#   Pvalue = append(Pvalue, kr$p.value)
# }
# 
# write.table(data.frame(ComparePairs = cp, statisticD, Pvalue), 'Figure_5C-KStestTarget.txt', 
#             append = F, quote = F, sep = '\t', row.names = F, col.names = T)
# 
# #########################
# # Kolmogorov-Smirnov test
# com.pair = c("D0", "D1", "D3", "D5", "D6", "D7")
# cp = statisticD = Pvalue = c()
# for (i in 1:5) {
#   set.seed(618)
#   c = randomDF[randomDF$network == com.pair[i], 'remain.mean']
#   m = randomDF[randomDF$network == com.pair[i+1], 'remain.mean']
#   sd.c = randomDF[randomDF$network == com.pair[i], 'remain.sd']
#   sd.m = randomDF[randomDF$network == com.pair[i+1], 'remain.sd']
#   kr = wilcox.test(rnorm(100, mean = c, sd = sd.c), rnorm(100, mean = m, sd = sd.m))
#   cp = append(cp, paste0(com.pair[i], '_vs_', com.pair[i+1]))
#   statisticD = append(statisticD, kr$statistic)
#   Pvalue = append(Pvalue, kr$p.value)
# }
# 
# write.table(data.frame(ComparePairs = cp, statisticD, Pvalue), 'Figure_5C-WMtestTarget.txt', 
#             append = F, quote = F, sep = '\t', row.names = F, col.names = T)

# ————————————————————————————————————————

dfrandom = read.table('Figure_5C-simuRandomDelt.txt', header = T, sep = '\t')
dfrandom$Type = 'Random'
dftarget = read.table('Figure_5C-simuTargetDelt.txt', header = T, sep = '\t')
dftarget$Type = 'Target'
colnames(dfrandom)[1] = colnames(dftarget)[1] = 'removed'
randomDF = rbind(dfrandom, dftarget)
randomDF = randomDF[!(is.na(randomDF$network)), ]

p = ggplot(randomDF, aes(x = network, y = remain.mean, fill = Type)) +
  geom_col(position = 'dodge', alpha = 0.6) +
  geom_errorbar(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd), 
                width=0.2, position=position_dodge(.9)) + 
  scale_fill_manual(values = c('#2F5D50', '#B7C9D3'), 
                    breaks = c('Random', 'Target'))+ 
  labs(x= NULL, y = 'Robustness') +
  theme( panel.background = element_rect(fill = NA, color = 'black'), 
         plot.background = element_blank(), 
         panel.grid = element_blank(), 
         axis.text.x = element_text(size = 8, color = "black", 
                                    hjust = 0.5, vjust = 0.5),
         axis.text.y = element_text(size = 8, angle = 90, color = "black",
                                    hjust = 0.5, vjust = 0.5),
         axis.title.y = element_text(size = 10), 
         strip.background = element_blank(),
         strip.text = element_blank(),
         legend.position = 'none')

ggsave('Figure_5C-Simulation.pdf', p, units = 'in', width = 1.6, height = 1.5)

#######################################################################################################

# Figure_5C ii

vulDF = read.table('Figure_5C-EffiVulnDF.txt', header = T, sep = '\t', stringsAsFactors = F)

p = ggplot(vulDF, aes(x = Group, y = vulnerability, fill = Group)) +
  geom_col(width = 0.6, alpha = 0.6) +
  scale_fill_manual(values = pal_nejm()(6), 
                    breaks = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')) + 
  labs(x= NULL, y = 'Vulnerability') +
  theme( panel.background = element_rect(fill = NA, color = 'black'), 
         plot.background = element_blank(), 
         panel.grid = element_blank(), 
         axis.text.x = element_text(size = 8, color = "black",),
         axis.text.y = element_text(size = 8, angle = 90, color = "black", 
                                    hjust = 0.5, vjust = 0.5),
         axis.title.y = element_text(size = 10), 
         legend.position = 'none')

ggsave('Figure_5C-SimulationVulnetability.pdf', p, width = 1.6, height = 1.5 )


#######################################################################################################

# Figure_5C iii

otusDF = read.table('OTUtableRela.txt', header = T, sep = '\t', stringsAsFactors = F, row.names = 1)
ndDF.0 = read.table('NetInfo/D0 0.95 node_attribute.txt', header = T, sep = '\t', row.names = 1)
nodesFns = list.files('NetInfo/', 'node')

persistence = c()
for (hfn in c('D1', 'D3', 'D5', 'D6', 'D7')) {
  fns = nodesFns[str_detect(nodesFns, hfn)]
  Ndf = read.table(paste0('NetInfo/', fns), header = T, sep = '\t', row.names = 1)
  totoalNDs = length(unique(c(rownames(ndDF.0), rownames(Ndf))))
  otudf.nd = otusDF[, str_detect(colnames(otusDF), hfn)]
  otudf.nd[otudf.nd > 0] = 1
  persistence = append(persistence, sum(rowSums(otudf.nd) == ncol(otudf.nd))/totoalNDs)
}

DFpersist = data.frame(Persistence = persistence, 
                       Group = c('D1', 'D3', 'D5', 'D6', 'D7'))
# write.table(DFpersist, 'Figure_5C-Persistence.txt', append = F, sep = '\t', 
#             row.names = F, col.names = T, quote = F)

p = ggplot(DFpersist, aes(x = Group, y = Persistence, fill = Group)) +
  geom_col(width = 0.6, alpha = 0.6) + 
  scale_fill_manual(values = pal_nejm()(6)[2:6], 
                    breaks = c('D1', 'D3', 'D5', 'D6', 'D7'))+
  labs(x= NULL, y = 'Persistence of nodes') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        plot.background = element_blank(),  
        panel.grid = element_blank(), 
        axis.text.x = element_text(size = 8, color="black", 
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 8, angle = 90, color="black", 
                                   vjust = 0, hjust = 0.5),
        axis.title.y = element_text(size = 10), 
        legend.position = 'none')

ggsave('Figure_5C-PersistenceNodes.pdf', p, width = 1.6, height = 1.5)

#######################################################################################################




#######################################################################################################

# Figure_5D

gini = data.frame(NetIndex = c('node.degree', 'node.betw', 'node.stress', 
                               'node.evcent', 'Zi', 'Pi'))

for (dl in list.files(path = 'NetInfo/', pattern = 'node_attribute.txt')) {
  nodes = read.table(paste0('NetInfo/', dl), header = T, sep = '\t')
  nodes.f = nodes[ , c(-6,-7)]
  print('-----------')
  print(substr(dl, 1, 2))
  resRF = randomForest(scale(nodes.f[ , 2:7], center = T, scale = T),  
                       mtry = 2, ntree = 5000, proximity = T)
  print(resRF$importance)
  
  gini = cbind(gini, as.data.frame(resRF$importance))
  colnames(gini)[ncol(gini)] = substr(dl, 1, 2)
}

write.table(gini, 'Figure_5D-MeanDecreaseGini.txt', append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)

# ------------------------------

nd = data.frame()

for (dl in list.files(path = 'NetInfo/', pattern = 'node_attribute.txt')) {
  nodes.f = read.table(paste0('NetInfo/', dl), header = T, sep = '\t')
  nodes.f$No..module = paste0('M', nodes.f$No..module + 1)
  mtable = data.frame(table(nodes.f$No..module))
  mtable = mtable[mtable$Freq >= 5 , ]
  nodes.f = nodes.f[nodes.f$No..module %in% mtable$Var1, ]
  
  nbetw = nodes.f$node.betw
  nstre = nodes.f$node.stress
  names(nbetw) = names(nstre) = nodes.f$id
  
  ids.nbetw = names(sort(nbetw, decreasing = T)[1:(0.2*(nrow(nodes.f)))])
  ids.nstre = names(sort(nstre, decreasing = T)[1:(0.2*(nrow(nodes.f)))])
  
  nodes.f2 = nodes.f[nodes.f$id %in% intersect(ids.nbetw, ids.nstre), ]
  
  
  nodes.f$node.betw = scale(nodes.f$node.betw, center = T, scale = T)
  nodes.f$node.stress = scale(nodes.f$node.stress, center = T, scale = T)
  
  nodes.f2 = nodes.f[nodes.f$node.betw >= 0. , ]
  nodes.f2$Group = substr(dl, 1, 2)
  
  nd = rbind(nd, nodes.f2)
}

tax = read.table('OTUTax.txt', header = T, sep = '\t')
colnames(tax)[1] = 'id'

nd = merge(nd, tax, by = 'id', all.x = T)

write.table(nd, 'Figure_5D-filterNodes.txt', append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)

# -------------------------------

gini = read.table('Figure_5D-filterNodes.txt', header = T, sep = '\t')

sort(table(gini$Family))

table(gini$Group)
# D0 D1 D3 D5 D6 D7 
# 80 85 36 30 36 16 

gini.f = gini[gini$Family %in% names(table(gini$Family) )[table(gini$Family) > 5], ]

table(gini.f$Family)   

fams = c("Anaerolineaceae", "Dethiosulfatibacteraceae", "Fusobacteriaceae", 
         "Izemoplasmataceae","Lentimicrobiaceae","Marinifilaceae", 
         "Marinilabiliaceae","Prolixibacteraceae","Spirochaetaceae","Vibrionaceae")

# fams = c('Prolixibacteraceae','Vibrionaceae','Marinifilaceae',
#          'Izemoplasmataceae','Dethiosulfatibacteraceae','Marinilabiliaceae')

gini.f = gini[gini$Family %in% fams, c('Group', 'Family')] 
gini.f$Count = 1

df = summaryBy(Count ~ Group + Family, gini.f, FUN = sum)
colnames(df)[3] = 'Count'

df$bubble_size = rescale(sqrt(df$Count), to = c(3.2, 14.5))
df$label_size  = rescale(sqrt(df$Count), to = c(1.9, 4.7))

p = ggplot(df, aes(x = Group, y = Family, size = Count, label = Count)) +
  geom_point(aes(size = bubble_size), color = '#DDE7F0') +
  geom_text(aes(size = label_size), color = '#3A0D14') +
  theme_bw() + 
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8, color = 'black', face = 'italic',
                                   hjust = 1, vjust = 0.5),
        axis.text.x = element_text(size = 8, color = 'black' ))

ggsave('Figure_5D-largeFamily.pdf', p, width = 3, height = 1.5)















