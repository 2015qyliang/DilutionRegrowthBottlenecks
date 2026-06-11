
library(grid)
library(ggplot2)
library(ggsci)

#####################################################################

ztOrder = ztDiversity = Group = c()
for (dl in c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')) {
  dtzt = readRDS(paste0('Sup_Fig_S3-', dl, 'zeta.decline.mc.rds'))
  ztDiversity = append(ztDiversity, dtzt$zeta.val)
  ztOrder = append(ztOrder, 1:length(dtzt$zeta.val))
  Group = append(Group, rep(dl, length(dtzt$zeta.val)))
}
ztDiversityDF = data.frame(ztOrder, ztDiversity, Group)

write.table(ztDiversityDF, '33-zetaDiv.txt', append = F,quote = F,
            sep = '\t', row.names = F,col.names = T)

###################################

ztOrder = ztRatio = Group = c()
for (dl in c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')) {
  dtzt = readRDS(paste0('Sup_Fig_S3-', dl, 'zeta.decline.mc.rds'))
  ztRatio = append(ztRatio, dtzt$ratio)
  ztOrder = append(ztOrder, 1:length(dtzt$ratio))
  Group = append(Group, rep(dl, length(dtzt$ratio)))
}
ztRatioDF = data.frame(ztOrder, ztRatio, Group)

###################################

ztAIC = ztRegr = Group = c()
for (dl in c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')) {
  dtzt = readRDS(paste0('Sup_Fig_S3-', dl, 'zeta.decline.mc.rds'))
  ztAIC = append(ztAIC, dtzt$aic$AIC)
  ztRegr = append(ztRegr, c('Exponential', 'PowerLaw'))
  Group = append(Group, rep(dl, 2))
}
ztRegrDF = data.frame(ztAIC, ztRegr, Group)

#####################################################################

pDiversity = ggplot(ztDiversityDF, aes(x = ztOrder, 
                                       y = ztDiversity, 
                                       group = Group, 
                                       color = Group)) + 
  geom_point(shape = 16, ) + geom_line() + 
  # scale_color_npg() +
  scale_colour_manual(values = pal_nejm()(6),
                      breaks = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')) +
  labs(x = 'Order', y = 'Zeta diversity') + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = NULL, color = 'black'),
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.position = 'none', # right  none
        axis.text.x = element_text(size = 10, color = 'black'), 
        axis.text.y = element_text(size = 10, color = 'black', 
                                   angle = 90, hjust = 0.5) )

pRatio = ggplot(ztRatioDF, aes(x = ztOrder, y = ztRatio,
                               group = Group, color = Group)) + 
  geom_point(shape = 16, ) + geom_line() + 
  # scale_color_npg() +
  scale_colour_manual(values = pal_nejm()(6),
                      breaks = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')) +
  labs(x = 'Order', y = 'Zeta diversity ratio') + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = NULL, color = 'black'),
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.position = 'none', # right  none
        axis.text.x = element_text(size = 10, color = 'black'), 
        axis.text.y = element_text(size = 10, color = 'black', 
                                   angle = 90, hjust = 0.5) )

ztRegr = ggplot(ztRegrDF, aes(x = Group, y = ztAIC, fill = ztRegr)) + 
  geom_col(position = 'dodge', alpha = 1) + 
  scale_fill_manual(values = c('#f36e42', '#2e8abf'), 
                    breaks = c('Exponential', 'PowerLaw')) +
  labs(x = 'Dilutions', y = 'Akaike information criterion (AIC)') + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = NULL, color = 'black'),
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.position = 'none', # right  none
        axis.text.x = element_text(size = 10, color = 'black'), 
        axis.text.y = element_text(size = 10, color = 'black', 
                                   angle = 90, hjust = 0.5) )

pdf(width = 8, height = 3, file = 'Sup_Fig_S3-zeta.pdf' )
grid.newpage()  
pushViewport(viewport(layout = grid.layout(1, 3)))
vplayout = function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(pDiversity, vp = vplayout(1, 1))
print(pRatio, vp = vplayout(1, 2))
print(ztRegr, vp = vplayout(1, 3))
dev.off()



