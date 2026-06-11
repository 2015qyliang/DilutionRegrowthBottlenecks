# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggraph)
library(igraph)
library(colormap)
library(hrbrthemes)
library(viridis)
library(patchwork)

##########################################################################

cultDF2 = read.table('03-subAllOtusCult.txt', header = F, sep = '\t')
ndLargeMd = read.table('03-filterNodes.txt', header = T, sep = '\t')
colnames(ndLargeMd)[1] = 'OTUID'
ndLargeMd$isCult = NA
ndLargeMd$isCult[ndLargeMd$OTUID %in% cultDF2$V1] = 'Cultured' 

# rankotus.raw = read.table('03-AllOtuRelAbundRankTax.txt', header = T, sep = '\t')
# rownames(rankotus.raw) = rankotus.raw$OTUID

sampleGroups = read.table('02-RMTvalue.txt', header = T, sep = '\t', stringsAsFactors = F)

# cultDF = read.table('03-subAllOtusCult.txt', header = F, sep = '\t')
# cultDF = data.frame(OTUID = cultDF$V1, isCult = 'CuledYes')
# rankotus = merge(rankotus.raw, cultDF, by = 'OTUID', all.x = T)
# rownames(rankotus) = rankotus$OTUID
# rankotus = rankotus[rownames(rankotus.raw), c('OTUID', 'RelAbund')]
# # rankotus$OTUID = factor(rankotus$OTUID, levels = rankotus$OTUID)
# 
ndOrder = readLines('00-AllNodesLeafNodes.txt')

##########################################################################
##########################################################################
##########################################################################

i = 1

rankotus.raw = read.table('03-AllOtuRelAbundRankTax.txt', header = T, sep = '\t')
rankotus.raw = rankotus.raw[rankotus.raw$Group == 'D0', ]
rownames(rankotus.raw) = rankotus.raw$OTUID
cultDF = read.table('03-subAllOtusCult.txt', header = F, sep = '\t')
cultDF = data.frame(OTUID = cultDF$V1, isCult = 'CuledYes')
rankotus = merge(rankotus.raw, cultDF, by = 'OTUID', all.x = T)
rownames(rankotus) = rankotus$OTUID
rankotus = rankotus[rownames(rankotus.raw), c('OTUID', 'RelAbund')]


fn = paste0(sampleGroups$groups[i], ' ', sprintf('%0.2f', sampleGroups$value[i]))
edges = read.table(paste0('../NetInfo/', fn,' edge_attribute.txt'),
                   header = T, sep = '\t', stringsAsFactors = F)
nds = ndLargeMd$OTUID[ndLargeMd$Group == sampleGroups$groups[i]]
links = edges[edges$from %in% nds & edges$to %in% nds, ]

unique(links$type)

nodes = ndLargeMd[ndLargeMd$Group == sampleGroups$groups[i], c('OTUID', 'isCult', 'No..module')]
nodes = merge(rankotus, nodes,by = 'OTUID', all.x = T) 

nodes$rareAbun = NA
nodes$rareAbun[nodes$RelAbund >= 0.001] = 'Abundant'
nodes$rareAbun[nodes$RelAbund < 0.001] = 'Rare'
nodes$sizeAbun[nodes$RelAbund >= 0.001] = 'Size1'
nodes$sizeAbun[nodes$RelAbund >= 0.0001 & nodes$RelAbund < 0.001] = 'Size2'
nodes$sizeAbun[nodes$RelAbund >= 0.00001 & nodes$RelAbund < 0.0001] = 'Size3'
nodes$sizeAbun[nodes$RelAbund < 0.00001] = 'Size4'

rownames(nodes) = nodes$OTUID
nodes = nodes[ndOrder, ]

net.graph = graph_from_data_frame( d = links, vertices = nodes, directed = F)

p0 = ggraph(net.graph, layout="linear") + 
  geom_edge_arc(aes(color = type),alpha = 0.4, width = 0.6, fold=F) +
  scale_edge_color_gradient(low = '#2e8abf',  high = '#f36e42') + # low - blue; high - orange
  geom_node_point(aes(color = rareAbun,  size = sizeAbun), alpha = 0.85, shape = 16) + 
  # geom_node_text(aes(label = OTUIDlabel), size = 2, angle = 90) +
  scale_color_manual(values = c('#4daf4a', '#984ea3'), breaks = c('Abundant', 'Rare')) +
  scale_size_manual(breaks = paste0('Size', 4:1), values = c(1.0, 1.4, 1.8, 2.2)) + 
  theme_void() +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null") )  

ggsave(paste0('06-Arc-', sampleGroups$groups[i], '.pdf'), 
       p0, width = 9.5, height = 0.83, units = 'in')

#################################

p0 = ggraph(net.graph, layout="linear") + 
  geom_node_point(aes(color = No..module), size = 0.78 , alpha = 0.85, shape = 15) + 
  scale_color_manual(na.value = 'grey70',
                     values = colorRampPalette(pal_npg()(8))(length(unique(nodes$No..module)))) +
  theme_void() +
  theme(legend.position = "none",
        plot.margin=unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null") )  

ggsave(paste0('06-Arc-', sampleGroups$groups[i], '-Module.pdf'), 
       p0, width = 9.5, height = 0.1, units = 'in')

#################################

p0 = ggraph(net.graph, layout="linear") + 
  geom_node_point(aes(color = isCult), size = 0.7 , alpha = 0.85, shape = 17) + 
  scale_color_manual(na.value = 'grey70', values = '#2e8abf') + 
  theme_void() +
  theme(legend.position = "none",
        plot.margin=unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null") )  

ggsave(paste0('06-Arc-', sampleGroups$groups[i], '-Culed.pdf'), 
       p0, width = 9.5, height = 0.1, units = 'in')

##########################################################################
##########################################################################
#################################

i = 2

rankotus.raw = read.table('03-AllOtuRelAbundRankTax.txt', header = T, sep = '\t')
rankotus.raw = rankotus.raw[rankotus.raw$Group == 'D1', ]
rownames(rankotus.raw) = rankotus.raw$OTUID
cultDF = read.table('03-subAllOtusCult.txt', header = F, sep = '\t')
cultDF = data.frame(OTUID = cultDF$V1, isCult = 'CuledYes')
rankotus = merge(rankotus.raw, cultDF, by = 'OTUID', all.x = T)
rownames(rankotus) = rankotus$OTUID
rankotus = rankotus[rownames(rankotus.raw), c('OTUID', 'RelAbund')]


fn = paste0(sampleGroups$groups[i], ' ', sprintf('%0.2f', sampleGroups$value[i]))
edges = read.table(paste0('../NetInfo/', fn,' edge_attribute.txt'),
                   header = T, sep = '\t', stringsAsFactors = F)
nds = ndLargeMd$OTUID[ndLargeMd$Group == sampleGroups$groups[i]]
links = edges[edges$from %in% nds & edges$to %in% nds, ]

unique(links$type)

nodes = ndLargeMd[ndLargeMd$Group == sampleGroups$groups[i], c('OTUID', 'isCult', 'No..module')]
nodes = merge(rankotus, nodes,by = 'OTUID', all.x = T) 

nodes$rareAbun = NA
nodes$rareAbun[nodes$RelAbund >= 0.001] = 'Abundant'
nodes$rareAbun[nodes$RelAbund < 0.001] = 'Rare'
nodes$sizeAbun[nodes$RelAbund >= 0.001] = 'Size1'
nodes$sizeAbun[nodes$RelAbund >= 0.0001 & nodes$RelAbund < 0.001] = 'Size2'
nodes$sizeAbun[nodes$RelAbund >= 0.00001 & nodes$RelAbund < 0.0001] = 'Size3'
nodes$sizeAbun[nodes$RelAbund < 0.00001] = 'Size4'

rownames(nodes) = nodes$OTUID
nodes = nodes[ndOrder, ]

net.graph = graph_from_data_frame( d = links, vertices = nodes, directed = F)

p1 = ggraph(net.graph, layout="linear") + 
  # geom_edge_arc(color = '#f36e42', alpha = 0.4, width = 0.6, fold=F) +
  geom_edge_arc(aes(color = type),alpha = 0.4, width = 0.6, fold=F) +
  scale_edge_color_gradient(low = '#2e8abf',  high = '#f36e42') + # low - blue; high - orange
  geom_node_point(aes(color = rareAbun,  size = sizeAbun), alpha = 0.85, shape = 16) + 
  scale_color_manual(values = c('#4daf4a', '#984ea3'), breaks = c('Abundant', 'Rare')) +
  scale_size_manual(breaks = paste0('Size', 4:1), values = c(1.0, 1.4, 1.8, 2.2)) + 
  theme_void() +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null") )  

ggsave(paste0('06-Arc-', sampleGroups$groups[i], '.pdf'), 
       p1, width = 9.5, height = 0.83, units = 'in')

#################################


p1 = ggraph(net.graph, layout="linear") + 
  geom_node_point(aes(color = No..module), size = 0.78 , alpha = 0.85, shape = 15) + 
  scale_color_manual(na.value = 'grey70',
                     values = colorRampPalette(pal_npg()(8))(length(unique(nodes$No..module)))) +
  theme_void() +
  theme(legend.position = "none",
        plot.margin=unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null") )  

ggsave(paste0('06-Arc-', sampleGroups$groups[i], '-Module.pdf'), 
       p1, width = 9.5, height = 0.1, units = 'in')

#################################

# p1 = ggraph(net.graph, layout="linear") + 
#   geom_node_point(aes(color = isCult), size = 1.5 , alpha = 0.85, shape = 16) + 
#   scale_color_manual(na.value = 'grey70', values = '#2e8abf') + 
#   theme_void() +
#   theme(legend.position = "none",
#         plot.margin=unit(c(0,0,0,0), "null"),
#         panel.spacing=unit(c(0,0,0,0), "null") )  
# 
# ggsave(paste0('06-Arc-', sampleGroups$groups[i], '-Culed.pdf'), 
#        p1, width = 7, height = 0.1, units = 'in')

##########################################################################
##########################################################################
#################################

i = 3

rankotus.raw = read.table('03-AllOtuRelAbundRankTax.txt', header = T, sep = '\t')
rankotus.raw = rankotus.raw[rankotus.raw$Group == 'D3', ]
rownames(rankotus.raw) = rankotus.raw$OTUID
cultDF = read.table('03-subAllOtusCult.txt', header = F, sep = '\t')
cultDF = data.frame(OTUID = cultDF$V1, isCult = 'CuledYes')
rankotus = merge(rankotus.raw, cultDF, by = 'OTUID', all.x = T)
rownames(rankotus) = rankotus$OTUID
rankotus = rankotus[rownames(rankotus.raw), c('OTUID', 'RelAbund')]


fn = paste0(sampleGroups$groups[i], ' ', sprintf('%0.2f', sampleGroups$value[i]))
edges = read.table(paste0('../NetInfo/', fn,' edge_attribute.txt'),
                   header = T, sep = '\t', stringsAsFactors = F)
nds = ndLargeMd$OTUID[ndLargeMd$Group == sampleGroups$groups[i]]
links = edges[edges$from %in% nds & edges$to %in% nds, ]

unique(links$type)

nodes = ndLargeMd[ndLargeMd$Group == sampleGroups$groups[i], c('OTUID', 'isCult', 'No..module')]
nodes = merge(rankotus, nodes,by = 'OTUID', all.x = T) 

nodes$rareAbun = NA
nodes$rareAbun[nodes$RelAbund >= 0.001] = 'Abundant'
nodes$rareAbun[nodes$RelAbund < 0.001] = 'Rare'
nodes$sizeAbun[nodes$RelAbund >= 0.001] = 'Size1'
nodes$sizeAbun[nodes$RelAbund >= 0.0001 & nodes$RelAbund < 0.001] = 'Size2'
nodes$sizeAbun[nodes$RelAbund >= 0.00001 & nodes$RelAbund < 0.0001] = 'Size3'
nodes$sizeAbun[nodes$RelAbund < 0.00001] = 'Size4'

rownames(nodes) = nodes$OTUID
nodes = nodes[ndOrder, ]

net.graph = graph_from_data_frame( d = links, vertices = nodes, directed = F)

p3 = ggraph(net.graph, layout="linear") + 
  # geom_edge_arc(color = '#2e8abf', alpha = 0.5, width = 0.7, fold=F) +
  geom_edge_arc(aes(color = type),alpha = 0.5, width = 0.7, fold=F) +
  scale_edge_color_gradient(low = '#2e8abf',  high = '#f36e42') + # low - blue; high - orange
  geom_node_point(aes(color = rareAbun,  size = sizeAbun), alpha = 0.85, shape = 16) + 
  scale_color_manual(values = c('#4daf4a', '#984ea3'), breaks = c('Abundant', 'Rare')) +
  scale_size_manual(breaks = paste0('Size', 4:1), values = c(1.0, 1.4, 1.8, 2.2)) + 
  theme_void() +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null") ) 

ggsave(paste0('06-Arc-', sampleGroups$groups[i], '.pdf'), 
       width = 9.5, height = 0.83, units = 'in')

#################################


p3 = ggraph(net.graph, layout="linear") + 
  geom_node_point(aes(color = No..module), size = 0.78 , alpha = 0.85, shape = 15) + 
  scale_color_manual(na.value = 'grey70',
                     values = colorRampPalette(pal_npg()(8))(length(unique(nodes$No..module)))) +
  theme_void() +
  theme(legend.position = "none",
        plot.margin=unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null") )  

ggsave(paste0('06-Arc-', sampleGroups$groups[i], '-Module.pdf'), 
       p3, width = 9.5, height = 0.1, units = 'in')

#################################

# p3 = ggraph(net.graph, layout="linear") + 
#   geom_node_point(aes(color = isCult), size = 1.5 , alpha = 0.85, shape = 16) + 
#   scale_color_manual(na.value = 'grey70', values = '#2e8abf') + 
#   theme_void() +
#   theme(legend.position = "none",
#         plot.margin=unit(c(0,0,0,0), "null"),
#         panel.spacing=unit(c(0,0,0,0), "null") )  
# 
# ggsave(paste0('06-Arc-', sampleGroups$groups[i], '-Culed.pdf'), 
#        p3, width = 7, height = 0.1, units = 'in')


##########################################################################
##########################################################################
#################################

i = 4

rankotus.raw = read.table('03-AllOtuRelAbundRankTax.txt', header = T, sep = '\t')
rankotus.raw = rankotus.raw[rankotus.raw$Group == 'D5', ]
rownames(rankotus.raw) = rankotus.raw$OTUID
cultDF = read.table('03-subAllOtusCult.txt', header = F, sep = '\t')
cultDF = data.frame(OTUID = cultDF$V1, isCult = 'CuledYes')
rankotus = merge(rankotus.raw, cultDF, by = 'OTUID', all.x = T)
rownames(rankotus) = rankotus$OTUID
rankotus = rankotus[rownames(rankotus.raw), c('OTUID', 'RelAbund')]


fn = paste0(sampleGroups$groups[i], ' ', sprintf('%0.2f', sampleGroups$value[i]))
edges = read.table(paste0('../NetInfo/', fn,' edge_attribute.txt'),
                   header = T, sep = '\t', stringsAsFactors = F)
nds = ndLargeMd$OTUID[ndLargeMd$Group == sampleGroups$groups[i]]
links = edges[edges$from %in% nds & edges$to %in% nds, ]

unique(links$type)

nodes = ndLargeMd[ndLargeMd$Group == sampleGroups$groups[i], c('OTUID', 'isCult', 'No..module')]
nodes = merge(rankotus, nodes,by = 'OTUID', all.x = T) 

nodes$rareAbun = NA
nodes$rareAbun[nodes$RelAbund >= 0.001] = 'Abundant'
nodes$rareAbun[nodes$RelAbund < 0.001] = 'Rare'
nodes$sizeAbun[nodes$RelAbund >= 0.001] = 'Size1'
nodes$sizeAbun[nodes$RelAbund >= 0.0001 & nodes$RelAbund < 0.001] = 'Size2'
nodes$sizeAbun[nodes$RelAbund >= 0.00001 & nodes$RelAbund < 0.0001] = 'Size3'
nodes$sizeAbun[nodes$RelAbund < 0.00001] = 'Size4'

rownames(nodes) = nodes$OTUID
nodes = nodes[ndOrder, ]

net.graph = graph_from_data_frame( d = links, vertices = nodes, directed = F)

p5 = ggraph(net.graph, layout="linear") + 
  geom_edge_arc(color = '#2e8abf', alpha = 0.5, width = 0.7, fold=F) +
  geom_node_point(aes(color = rareAbun,  size = sizeAbun), alpha = 0.85, shape = 16) + 
  scale_color_manual(values = c('#4daf4a', '#984ea3'), breaks = c('Abundant', 'Rare')) +
  scale_size_manual(breaks = paste0('Size', 4:1), values = c(1.0, 1.4, 1.8, 2.2)) + 
  theme_void() +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null") )   

ggsave(paste0('06-Arc-', sampleGroups$groups[i], '.pdf'), 
       p5, width = 9.5, height = 0.83, units = 'in')

#################################


p5 = ggraph(net.graph, layout="linear") + 
  geom_node_point(aes(color = No..module), size = 0.78 , alpha = 0.85, shape = 15) + 
  scale_color_manual(na.value = 'grey70',
                     values = colorRampPalette(pal_npg()(8))(length(unique(nodes$No..module)))) +
  theme_void() +
  theme(legend.position = "none",
        plot.margin=unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null") )  

ggsave(paste0('06-Arc-', sampleGroups$groups[i], '-Module.pdf'), 
       p5, width = 9.5, height = 0.1, units = 'in')

#################################

# p5 = ggraph(net.graph, layout="linear") + 
#   geom_node_point(aes(color = isCult), size = 1.5 , alpha = 0.85, shape = 16) + 
#   scale_color_manual(na.value = 'grey70', values = '#2e8abf') + 
#   theme_void() +
#   theme(legend.position = "none",
#         plot.margin=unit(c(0,0,0,0), "null"),
#         panel.spacing=unit(c(0,0,0,0), "null") )  
# 
# ggsave(paste0('06-Arc-', sampleGroups$groups[i], '-Culed.pdf'), 
#        p5, width = 7, height = 0.1, units = 'in')


##########################################################################
##########################################################################
#################################

i = 5

rankotus.raw = read.table('03-AllOtuRelAbundRankTax.txt', header = T, sep = '\t')
rankotus.raw = rankotus.raw[rankotus.raw$Group == 'D6', ]
rownames(rankotus.raw) = rankotus.raw$OTUID
cultDF = read.table('03-subAllOtusCult.txt', header = F, sep = '\t')
cultDF = data.frame(OTUID = cultDF$V1, isCult = 'CuledYes')
rankotus = merge(rankotus.raw, cultDF, by = 'OTUID', all.x = T)
rownames(rankotus) = rankotus$OTUID
rankotus = rankotus[rownames(rankotus.raw), c('OTUID', 'RelAbund')]


fn = paste0(sampleGroups$groups[i], ' ', sprintf('%0.2f', sampleGroups$value[i]))
edges = read.table(paste0('../NetInfo/', fn,' edge_attribute.txt'),
                   header = T, sep = '\t', stringsAsFactors = F)
nds = ndLargeMd$OTUID[ndLargeMd$Group == sampleGroups$groups[i]]
links = edges[edges$from %in% nds & edges$to %in% nds, ]

unique(links$type)

nodes = ndLargeMd[ndLargeMd$Group == sampleGroups$groups[i], c('OTUID', 'isCult', 'No..module')]
nodes = merge(rankotus, nodes,by = 'OTUID', all.x = T) 

nodes$rareAbun = NA
nodes$rareAbun[nodes$RelAbund >= 0.001] = 'Abundant'
nodes$rareAbun[nodes$RelAbund < 0.001] = 'Rare'
nodes$sizeAbun[nodes$RelAbund >= 0.001] = 'Size1'
nodes$sizeAbun[nodes$RelAbund >= 0.0001 & nodes$RelAbund < 0.001] = 'Size2'
nodes$sizeAbun[nodes$RelAbund >= 0.00001 & nodes$RelAbund < 0.0001] = 'Size3'
nodes$sizeAbun[nodes$RelAbund < 0.00001] = 'Size4'

rownames(nodes) = nodes$OTUID
nodes = nodes[ndOrder, ]

net.graph = graph_from_data_frame( d = links, vertices = nodes, directed = F)

p6 = ggraph(net.graph, layout="linear") + 
  geom_edge_arc(color = '#2e8abf', alpha = 0.5, width = 0.7, fold=F) +
  geom_node_point(aes(color = rareAbun,  size = sizeAbun), alpha = 0.85, shape = 16) + 
  scale_color_manual(values = c('#4daf4a', '#984ea3'), breaks = c('Abundant', 'Rare')) +
  scale_size_manual(breaks = paste0('Size', 4:1), values = c(1.0, 1.4, 1.8, 2.2)) + 
  theme_void() +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null") ) 

ggsave(paste0('06-Arc-', sampleGroups$groups[i], '.pdf'), 
       p6, width = 9.5, height = 0.83, units = 'in')

#################################


p6 = ggraph(net.graph, layout="linear") + 
  geom_node_point(aes(color = No..module), size = 0.78 , alpha = 0.85, shape = 15) + 
  scale_color_manual(na.value = 'grey70',
                     values = colorRampPalette(pal_npg()(8))(length(unique(nodes$No..module)))) +
  theme_void() +
  theme(legend.position = "none",
        plot.margin=unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null") )  

ggsave(paste0('06-Arc-', sampleGroups$groups[i], '-Module.pdf'), 
       p6, width = 9.5, height = 0.1, units = 'in')

#################################

# p6 = ggraph(net.graph, layout="linear") + 
#   geom_node_point(aes(color = isCult), size = 1.5 , alpha = 0.85, shape = 16) + 
#   scale_color_manual(na.value = 'grey70', values = '#2e8abf') + 
#   theme_void() +
#   theme(legend.position = "none",
#         plot.margin=unit(c(0,0,0,0), "null"),
#         panel.spacing=unit(c(0,0,0,0), "null") )  
# 
# ggsave(paste0('06-Arc-', sampleGroups$groups[i], '-Culed.pdf'), 
#        p6, width = 7, height = 0.1, units = 'in')

##########################################################################
##########################################################################
#################################

i = 6

rankotus.raw = read.table('03-AllOtuRelAbundRankTax.txt', header = T, sep = '\t')
rankotus.raw = rankotus.raw[rankotus.raw$Group == 'D7', ]
rownames(rankotus.raw) = rankotus.raw$OTUID
cultDF = read.table('03-subAllOtusCult.txt', header = F, sep = '\t')
cultDF = data.frame(OTUID = cultDF$V1, isCult = 'CuledYes')
rankotus = merge(rankotus.raw, cultDF, by = 'OTUID', all.x = T)
rownames(rankotus) = rankotus$OTUID
rankotus = rankotus[rownames(rankotus.raw), c('OTUID', 'RelAbund')]


fn = paste0(sampleGroups$groups[i], ' ', sprintf('%0.2f', sampleGroups$value[i]))
edges = read.table(paste0('../NetInfo/', fn,' edge_attribute.txt'),
                   header = T, sep = '\t', stringsAsFactors = F)
nds = ndLargeMd$OTUID[ndLargeMd$Group == sampleGroups$groups[i]]
links = edges[edges$from %in% nds & edges$to %in% nds, ]

unique(links$type)

nodes = ndLargeMd[ndLargeMd$Group == sampleGroups$groups[i], c('OTUID', 'isCult', 'No..module')]
nodes = merge(rankotus, nodes,by = 'OTUID', all.x = T) 

nodes$rareAbun = NA
nodes$rareAbun[nodes$RelAbund >= 0.001] = 'Abundant'
nodes$rareAbun[nodes$RelAbund < 0.001] = 'Rare'
nodes$sizeAbun[nodes$RelAbund >= 0.001] = 'Size1'
nodes$sizeAbun[nodes$RelAbund >= 0.0001 & nodes$RelAbund < 0.001] = 'Size2'
nodes$sizeAbun[nodes$RelAbund >= 0.00001 & nodes$RelAbund < 0.0001] = 'Size3'
nodes$sizeAbun[nodes$RelAbund < 0.00001] = 'Size4'

rownames(nodes) = nodes$OTUID
nodes = nodes[ndOrder, ]

net.graph = graph_from_data_frame( d = links, vertices = nodes, directed = F)

p7 = ggraph(net.graph, layout="linear") + 
  geom_edge_arc(color = '#2e8abf', alpha = 0.5, width = 0.7, fold=F) +
  geom_node_point(aes(color = rareAbun,  size = sizeAbun), alpha = 0.85, shape = 16) + 
  scale_color_manual(values = c('#4daf4a', '#984ea3'), breaks = c('Abundant', 'Rare')) +
  scale_size_manual(breaks = paste0('Size', 4:1), values = c(1.0, 1.4, 1.8, 2.2)) + 
  theme_void() +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null") ) 

ggsave(paste0('06-Arc-', sampleGroups$groups[i], '.pdf'), 
       width = 9.5, height = 0.83, units = 'in')

#################################


p7 = ggraph(net.graph, layout="linear") + 
  geom_node_point(aes(color = No..module), size = 0.78 , alpha = 0.85, shape = 15) + 
  scale_color_manual(na.value = 'grey70',
                     values = colorRampPalette(pal_npg()(8))(length(unique(nodes$No..module)))) +
  theme_void() +
  theme(legend.position = "none",
        plot.margin=unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null") )  

ggsave(paste0('06-Arc-', sampleGroups$groups[i], '-Module.pdf'), 
       p7, width = 9.5, height = 0.1, units = 'in')

#################################

# p7 = ggraph(net.graph, layout="linear") + 
#   geom_node_point(aes(color = isCult), size = 1.5 , alpha = 0.85, shape = 16) + 
#   scale_color_manual(na.value = 'grey70', values = '#2e8abf') + 
#   theme_void() +
#   theme(legend.position = "none",
#         plot.margin=unit(c(0,0,0,0), "null"),
#         panel.spacing=unit(c(0,0,0,0), "null") )  
# 
# ggsave(paste0('06-Arc-', sampleGroups$groups[i], '-Culed.pdf'), 
#        p7, width = 7, height = 0.1, units = 'in')

##########################################################################
##########################################################################










