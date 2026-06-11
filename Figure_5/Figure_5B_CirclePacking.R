# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(tidyverse)
library(packcircles)
library(ggplot2)
library(viridis)
library(jcolors)

##########################################################################

otusRela.r = read.table('OTUtableRela.txt', header = T, sep = '\t', row.names = 1)

otusRela = data.frame()
for (dl in paste0('D', c(0, 1, 3, 5, 6, 7))) {
  otur = rowMeans(otusRela.r[ , str_detect(colnames(otusRela.r), dl)])
  newdf = data.frame(OTUID = names(otur), RelAbund = otur, Group = dl)
  otusRela = rbind(otusRela, newdf)
}

##########################################################################

sampleGroups = read.table('RMT_value.txt', header = T, sep = '\t', stringsAsFactors = F)
ndOrder = readLines('RMT_AllNodesLeafNodes.txt') 

# mdcolors = colorRampPalette(pal_npg()(10))(50)
mdcolors = colorRampPalette(jcolors('pal8'))(50) 

##########################################################################
##########################################################################
##########################################################################

i = 1
rankotus.filt = otusRela[otusRela$Group == 'D0', ]
fn = paste0(sampleGroups$groups[i], ' ', sprintf('%0.2f', sampleGroups$value[i]))
nodes = read.table(paste0('./NetInfo/', fn,' node_attribute.txt'),
                   header = T, sep = '\t', stringsAsFactors = F)
colnames(nodes)[1] = 'OTUID'
nodes = merge(nodes, rankotus.filt, by = 'OTUID', all.x = T) 
nodes$No..module = paste0('M', nodes$No..module+1)
lgmd = names(table(nodes$No..module)[table(nodes$No..module) >=5])
nodes$No..module[!(nodes$No..module %in% lgmd)] = NA
# nodes = nodes[nodes$No..module %in% lgmd, ]

nodes$rareAbun = NA
nodes$rareAbun[nodes$RelAbund >= 0.001] = 'Abundant'
nodes$rareAbun[nodes$RelAbund < 0.001] = 'Rare'
nodes$sizeAbun[nodes$RelAbund >= 0.001] = 'Size1'
nodes$sizeAbun[nodes$RelAbund >= 0.0001 & nodes$RelAbund < 0.001] = 'Size2'
nodes$sizeAbun[nodes$RelAbund >= 0.00001 & nodes$RelAbund < 0.0001] = 'Size3'
nodes$sizeAbun[nodes$RelAbund < 0.00001] = 'Size4'

nodes$RelAbund[nodes$sizeAbun == 'Size4'] = 0.2
nodes$RelAbund[nodes$sizeAbun == 'Size3'] = 1
nodes$RelAbund[nodes$sizeAbun == 'Size2'] = 3
nodes$RelAbund[nodes$sizeAbun == 'Size1'] = 6

packing = circleProgressiveLayout(nodes$RelAbund, sizetype='area')
data = cbind(nodes, packing)
dat.gg = circleLayoutVertices(packing, npoints=100)
nodes$id = 1:length(unique(nodes$OTUID))
nodes0 = merge(nodes, dat.gg, by = 'id', all.x = T )

mdcls = mdcolors[1:11]
names(mdcls) = sort(unique(nodes0$No..module))

set.seed(123) 
p0 = ggplot() + 
  geom_polygon(data = nodes0, aes(x, y, group = id, fill = No..module), 
               colour = NA, alpha = 0.8) +
  scale_fill_manual(values = mdcls[sort(unique(nodes0$No..module))], 
                    breaks = sort(unique(nodes0$No..module)), na.value = 'grey85') +
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

ggsave('CirclePacking-D0.pdf', p0, width = 2, height = 2) 

##########################################################################
##########################################################################
##########################################################################

i = 2
rankotus.filt = otusRela[otusRela$Group == 'D1', ]
fn = paste0(sampleGroups$groups[i], ' ', sprintf('%0.2f', sampleGroups$value[i]))
nodes = read.table(paste0('./NetInfo/', fn,' node_attribute.txt'),
                   header = T, sep = '\t', stringsAsFactors = F)
colnames(nodes)[1] = 'OTUID'
nodes = merge(nodes, rankotus.filt, by = 'OTUID', all.x = T) 
nodes$No..module = paste0('M', nodes$No..module+1)
lgmd = names(table(nodes$No..module)[table(nodes$No..module) >=5])
nodes$No..module[!(nodes$No..module %in% lgmd)] = NA
# nodes = nodes[nodes$No..module %in% lgmd, ]

nodes$rareAbun = NA
nodes$rareAbun[nodes$RelAbund >= 0.001] = 'Abundant'
nodes$rareAbun[nodes$RelAbund < 0.001] = 'Rare'
nodes$sizeAbun[nodes$RelAbund >= 0.001] = 'Size1'
nodes$sizeAbun[nodes$RelAbund >= 0.0001 & nodes$RelAbund < 0.001] = 'Size2'
nodes$sizeAbun[nodes$RelAbund >= 0.00001 & nodes$RelAbund < 0.0001] = 'Size3'
nodes$sizeAbun[nodes$RelAbund < 0.00001] = 'Size4'

nodes$RelAbund[nodes$sizeAbun == 'Size4'] = 0.2
nodes$RelAbund[nodes$sizeAbun == 'Size3'] = 1
nodes$RelAbund[nodes$sizeAbun == 'Size2'] = 3
nodes$RelAbund[nodes$sizeAbun == 'Size1'] = 6

packing = circleProgressiveLayout(nodes$RelAbund, sizetype='area')
data = cbind(nodes, packing)
dat.gg = circleLayoutVertices(packing, npoints=100)
nodes$id = 1:length(unique(nodes$OTUID))
nodes1 = merge(nodes, dat.gg, by = 'id', all.x = T )

nodes1$mtOTUs = NA
nodes1$mtOTUsWd = NA
nodes1$mtOTUs[nodes1$OTUID %in% nodes0$OTUID] = 'black'
nodes1$mtOTUsWd[nodes1$OTUID %in% nodes0$OTUID] = 0.3

mdcls = mdcolors[12:(12+12)]
names(mdcls) = sort(unique(nodes1$No..module))

set.seed(123) 
p1 = ggplot() + 
  geom_polygon(data = nodes1, aes(x, y, group = id, fill = No..module), 
               colour = nodes1$mtOTUs, linewidth = nodes1$mtOTUsWd, alpha = 0.8) +
  scale_fill_manual(values = mdcls[sort(unique(nodes1$No..module))], 
                    breaks = sort(unique(nodes1$No..module)), na.value = 'grey85') +
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

ggsave('CirclePacking-D1.pdf', p1, width = 2, height = 2) 


##########################################################################
##########################################################################
##########################################################################

i = 3
rankotus.filt = otusRela[otusRela$Group == 'D3', ]
fn = paste0(sampleGroups$groups[i], ' ', sprintf('%0.2f', sampleGroups$value[i]))
nodes = read.table(paste0('./NetInfo/', fn,' node_attribute.txt'),
                   header = T, sep = '\t', stringsAsFactors = F)
colnames(nodes)[1] = 'OTUID'
nodes = merge(nodes, rankotus.filt, by = 'OTUID', all.x = T) 
nodes$No..module = paste0('M', nodes$No..module+1)
lgmd = names(table(nodes$No..module)[table(nodes$No..module) >=5])
nodes$No..module[!(nodes$No..module %in% lgmd)] = NA
# nodes = nodes[nodes$No..module %in% lgmd, ]

nodes$rareAbun = NA
nodes$rareAbun[nodes$RelAbund >= 0.001] = 'Abundant'
nodes$rareAbun[nodes$RelAbund < 0.001] = 'Rare'
nodes$sizeAbun[nodes$RelAbund >= 0.001] = 'Size1'
nodes$sizeAbun[nodes$RelAbund >= 0.0001 & nodes$RelAbund < 0.001] = 'Size2'
nodes$sizeAbun[nodes$RelAbund >= 0.00001 & nodes$RelAbund < 0.0001] = 'Size3'
nodes$sizeAbun[nodes$RelAbund < 0.00001] = 'Size4'

nodes$RelAbund[nodes$sizeAbun == 'Size4'] = 0.2
nodes$RelAbund[nodes$sizeAbun == 'Size3'] = 1
nodes$RelAbund[nodes$sizeAbun == 'Size2'] = 3
nodes$RelAbund[nodes$sizeAbun == 'Size1'] = 6

packing = circleProgressiveLayout(nodes$RelAbund, sizetype='area')
data = cbind(nodes, packing)
dat.gg = circleLayoutVertices(packing, npoints=100)
nodes$id = 1:length(unique(nodes$OTUID))
nodes3 = merge(nodes, dat.gg, by = 'id', all.x = T )

nodes3$mtOTUs = NA
nodes3$mtOTUsWd = NA
nodes3$mtOTUs[nodes3$OTUID %in% nodes0$OTUID] = 'black'
nodes3$mtOTUsWd[nodes3$OTUID %in% nodes0$OTUID] = 0.3

mdcls = mdcolors[25:(25+8)]
names(mdcls) = sort(unique(nodes3$No..module))

set.seed(123) 
p3 = ggplot() + 
  geom_polygon(data = nodes3, aes(x, y, group = id, fill = No..module), 
               colour = nodes3$mtOTUs, linewidth = nodes3$mtOTUsWd, alpha = 0.8) +
  scale_fill_manual(values = mdcls[sort(unique(nodes3$No..module))], 
                    breaks = sort(unique(nodes3$No..module)), na.value = 'grey85') +
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

ggsave('CirclePacking-D3.pdf', p3, width = 2, height = 2) 

##########################################################################
##########################################################################
##########################################################################

i = 4
rankotus.filt = otusRela[otusRela$Group == 'D5', ]
fn = paste0(sampleGroups$groups[i], ' ', sprintf('%0.2f', sampleGroups$value[i]))
nodes = read.table(paste0('./NetInfo/', fn,' node_attribute.txt'),
                   header = T, sep = '\t', stringsAsFactors = F)
colnames(nodes)[1] = 'OTUID'
nodes = merge(nodes, rankotus.filt, by = 'OTUID', all.x = T) 
nodes$No..module = paste0('M', nodes$No..module+1)
lgmd = names(table(nodes$No..module)[table(nodes$No..module) >=5])
nodes$No..module[!(nodes$No..module %in% lgmd)] = NA
# nodes = nodes[nodes$No..module %in% lgmd, ]

nodes$rareAbun = NA
nodes$rareAbun[nodes$RelAbund >= 0.001] = 'Abundant'
nodes$rareAbun[nodes$RelAbund < 0.001] = 'Rare'
nodes$sizeAbun[nodes$RelAbund >= 0.001] = 'Size1'
nodes$sizeAbun[nodes$RelAbund >= 0.0001 & nodes$RelAbund < 0.001] = 'Size2'
nodes$sizeAbun[nodes$RelAbund >= 0.00001 & nodes$RelAbund < 0.0001] = 'Size3'
nodes$sizeAbun[nodes$RelAbund < 0.00001] = 'Size4'

nodes$RelAbund[nodes$sizeAbun == 'Size4'] = 0.2
nodes$RelAbund[nodes$sizeAbun == 'Size3'] = 1
nodes$RelAbund[nodes$sizeAbun == 'Size2'] = 3
nodes$RelAbund[nodes$sizeAbun == 'Size1'] = 6

packing = circleProgressiveLayout(nodes$RelAbund, sizetype='area')
data = cbind(nodes, packing)
dat.gg = circleLayoutVertices(packing, npoints=100)
nodes$id = 1:length(unique(nodes$OTUID))
nodes5 = merge(nodes, dat.gg, by = 'id', all.x = T )

nodes5$mtOTUs = NA
nodes5$mtOTUsWd = NA
nodes5$mtOTUs[nodes5$OTUID %in% nodes0$OTUID] = 'black'
nodes5$mtOTUsWd[nodes5$OTUID %in% nodes0$OTUID] = 0.3

mdcls = mdcolors[34:(34+6)]
names(mdcls) = sort(unique(nodes5$No..module))

set.seed(123) 
p5 = ggplot() + 
  geom_polygon(data = nodes5, aes(x, y, group = id, fill = No..module), 
               colour = nodes5$mtOTUs, linewidth = nodes5$mtOTUsWd, alpha = 0.8) +
  scale_fill_manual(values = mdcls[sort(unique(nodes5$No..module))], 
                    breaks = sort(unique(nodes5$No..module)),  na.value = 'grey85') +
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

ggsave('CirclePacking-D5.pdf', p5, width = 2, height = 2) 

##########################################################################
##########################################################################
##########################################################################

i = 5
rankotus.filt = otusRela[otusRela$Group == 'D6', ]
fn = paste0(sampleGroups$groups[i], ' ', sprintf('%0.2f', sampleGroups$value[i]))
nodes = read.table(paste0('./NetInfo/', fn,' node_attribute.txt'),
                   header = T, sep = '\t', stringsAsFactors = F)
colnames(nodes)[1] = 'OTUID'
nodes = merge(nodes, rankotus.filt, by = 'OTUID', all.x = T) 
nodes$No..module = paste0('M', nodes$No..module+1)
lgmd = names(table(nodes$No..module)[table(nodes$No..module) >=5])
nodes$No..module[!(nodes$No..module %in% lgmd)] = NA
# nodes = nodes[nodes$No..module %in% lgmd, ]

nodes$rareAbun = NA
nodes$rareAbun[nodes$RelAbund >= 0.001] = 'Abundant'
nodes$rareAbun[nodes$RelAbund < 0.001] = 'Rare'
nodes$sizeAbun[nodes$RelAbund >= 0.001] = 'Size1'
nodes$sizeAbun[nodes$RelAbund >= 0.0001 & nodes$RelAbund < 0.001] = 'Size2'
nodes$sizeAbun[nodes$RelAbund >= 0.00001 & nodes$RelAbund < 0.0001] = 'Size3'
nodes$sizeAbun[nodes$RelAbund < 0.00001] = 'Size4'

nodes$RelAbund[nodes$sizeAbun == 'Size4'] = 0.2
nodes$RelAbund[nodes$sizeAbun == 'Size3'] = 1
nodes$RelAbund[nodes$sizeAbun == 'Size2'] = 3
nodes$RelAbund[nodes$sizeAbun == 'Size1'] = 6

packing = circleProgressiveLayout(nodes$RelAbund, sizetype='area')
data = cbind(nodes, packing)
dat.gg = circleLayoutVertices(packing, npoints=100)
nodes$id = 1:length(unique(nodes$OTUID))
nodes6 = merge(nodes, dat.gg, by = 'id', all.x = T )

nodes6$mtOTUs = NA
nodes6$mtOTUsWd = NA
nodes6$mtOTUs[nodes6$OTUID %in% nodes0$OTUID] = 'black'
nodes6$mtOTUsWd[nodes6$OTUID %in% nodes0$OTUID] = 0.3

mdcls = mdcolors[41:(41+6)]
names(mdcls) = sort(unique(nodes6$No..module))

set.seed(123) 
p6 = ggplot() + 
  geom_polygon(data = nodes6, aes(x, y, group = id, fill = No..module), 
               colour = nodes6$mtOTUs, linewidth = nodes6$mtOTUsWd, alpha = 0.8) +
  scale_fill_manual(values = mdcls[sort(unique(nodes6$No..module))], 
                    breaks = sort(unique(nodes6$No..module)), na.value = 'grey85') +
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

ggsave('CirclePacking-D6.pdf', p6, width = 2, height = 2) 

##########################################################################
##########################################################################
##########################################################################

i = 6
rankotus.filt = otusRela[otusRela$Group == 'D7', ]
fn = paste0(sampleGroups$groups[i], ' ', sprintf('%0.2f', sampleGroups$value[i]))
nodes = read.table(paste0('./NetInfo/', fn,' node_attribute.txt'),
                   header = T, sep = '\t', stringsAsFactors = F)
colnames(nodes)[1] = 'OTUID'
nodes = merge(nodes, rankotus.filt, by = 'OTUID', all.x = T) 
nodes$No..module = paste0('M', nodes$No..module+1)
lgmd = names(table(nodes$No..module)[table(nodes$No..module) >=5])
nodes$No..module[!(nodes$No..module %in% lgmd)] = NA
# nodes = nodes[nodes$No..module %in% lgmd, ]

nodes$rareAbun = NA
nodes$rareAbun[nodes$RelAbund >= 0.001] = 'Abundant'
nodes$rareAbun[nodes$RelAbund < 0.001] = 'Rare'
nodes$sizeAbun[nodes$RelAbund >= 0.001] = 'Size1'
nodes$sizeAbun[nodes$RelAbund >= 0.0001 & nodes$RelAbund < 0.001] = 'Size2'
nodes$sizeAbun[nodes$RelAbund >= 0.00001 & nodes$RelAbund < 0.0001] = 'Size3'
nodes$sizeAbun[nodes$RelAbund < 0.00001] = 'Size4'

nodes$RelAbund[nodes$sizeAbun == 'Size4'] = 0.2
nodes$RelAbund[nodes$sizeAbun == 'Size3'] = 1
nodes$RelAbund[nodes$sizeAbun == 'Size2'] = 3
nodes$RelAbund[nodes$sizeAbun == 'Size1'] = 6

packing = circleProgressiveLayout(nodes$RelAbund, sizetype='area')
data = cbind(nodes, packing)
dat.gg = circleLayoutVertices(packing, npoints=100)
nodes$id = 1:length(unique(nodes$OTUID))
nodes7 = merge(nodes, dat.gg, by = 'id', all.x = T )

nodes7$mtOTUs = NA
nodes7$mtOTUsWd = NA
nodes7$mtOTUs[nodes7$OTUID %in% nodes0$OTUID] = 'black'
nodes7$mtOTUsWd[nodes7$OTUID %in% nodes0$OTUID] = 0.3

mdcls = mdcolors[48:50]
names(mdcls) = sort(unique(nodes7$No..module))

set.seed(123) 
p7 = ggplot() + 
  geom_polygon(data = nodes7, aes(x, y, group = id, fill = No..module), 
               colour = nodes7$mtOTUs, linewidth = nodes7$mtOTUsWd, alpha = 0.8) +
  scale_fill_manual(values = mdcls[sort(unique(nodes7$No..module))], 
                    breaks = sort(unique(nodes7$No..module)), na.value = 'grey85') +
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

ggsave('CirclePacking-D7.pdf', p7, width = 2, height = 2) 

##########################################################################




