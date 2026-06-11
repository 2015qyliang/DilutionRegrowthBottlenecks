
library(tidyverse)
library(ggfortify)
library(lamW)
library(GGally)
library(ggcorrplot)
library(pracma)
library(ggsci)
library(grid)

source("./MacLaw_theme.R")

##################################

load( "./MacLawDataestimate.Rdata" )

##### LOGNORMAL
nbin = 20
statp = gamma_pars %>%  mutate( lf = log(f ) , cutoff = -100  ) %>% 
  group_by( idall  ) %>% mutate( df = (max(lf)-min(lf))/nbin ) %>% 
  mutate( b = as.integer( (lf - min(lf) )/df )  ) %>% 
  ungroup() %>%  group_by( idall, b, df ) %>% 
  summarise(  lf = mean(lf), cutoff = mean(cutoff), n = n() ) %>% 
  ungroup() %>% group_by(idall) %>%  
  mutate( p = n / sum(n) / df ) %>%  ungroup()

lognlog  = function(x) { 10^( -0.5 * x^2 - log10(2 * pi)/2 ) }

pldf = statp %>%  mutate( xx = runif(dim(statp)[1])  )  %>%  
  left_join(mean_pars) %>% filter(lf > c, n > 5) %>% 
  arrange(xx, sname) 
plnind = ggplot(pldf, aes( x =  (lf-mu)/sqrt(2*sigma^2)  , y = p, color = as.factor(sname) , shape = as.factor(sname) )  ) + 
  geom_point( size = 1, alpha = 0.5 ) +
  stat_function( fun = lognlog, color = "black", size = 0.5 ) +
  scalecols + scaleshapes +
  scale_x_continuous( "Rescaled log average\nrelative abundance", limits = c(-2, 2.5) ) +
  scale_y_log10( "Probability\ndensity", limits = c(0.009, 0.8), labels = fancy_scientific  )  +
  theme_bw() + 
  theme( panel.background = element_blank(),
         plot.background = element_blank(),
         panel.grid = element_blank(),
         legend.position = "none", 
         axis.text = element_text(size = 8, color = 'black'),
         axis.title = element_text(size = 10, color = 'black')  ) 
# ggsave( filename =  "26-logn_ind.pdf", plnind, width = 2.5, height = 2.4 )

##################################

####### TAYLOR
nbins = 20
taylor_binnes = gamma_pars %>% group_by( idall, sname ) %>% 
  mutate( lf = log(f), dlf = (max(lf)-min(lf))/nbins  ) %>% 
  mutate(b = as.integer( (lf-min(lf))/dlf )  ) %>% 
  ungroup() %>% group_by(idall, sname,b) %>% 
  summarise( vf = mean(vf), f = mean(f)  ) %>% 
  ungroup() %>% mutate(xx = runif(n() ) ) %>% arrange(xx) 

ptaylor_all = taylor_binnes  %>% 
  ggplot() + 
  aes( x = f,  y = vf, color = as.factor(sname), shape = as.factor(sname) ) +
  scalecols + scaleshapes +
  geom_point( size = 1, alpha = 0.5 ) +   
  geom_abline( slope = 2, size = 0.5, intercept = 1.  ) +
  scale_x_log10( "Average\nrelative abundance", 
                 limits = c(-2.5,2.5), 
                 labels = fancy_scientificb ) +
  scale_y_log10( "Variance of\nrelative abundance" , 
                 labels = fancy_scientific  )  +
  theme_bw() + 
  theme( panel.background = element_blank(),
         plot.background = element_blank(),
         panel.grid = element_blank(),
         legend.position = "none", 
         axis.text = element_text(size = 8, color = 'black'),
         axis.title = element_text(size = 10, color = 'black'))  
# ggsave( filename =  "26-taylor_ind.pdf", ptaylor_all, width = 2.5, height = 2.4  )
# ptaylor_all

##################################

######### GAMMA
nbin = 20
dh = proj %>% filter( o > 0.9 ) %>% 
  mutate( l = log(count/nreads), f = count/nreads ) %>% 
  group_by( otu_id, idall, sname ) %>% 
  mutate( ml = mean(l), sl = sd(l), k = mean(f)^2/var(f) ) %>% ungroup() %>% 
  mutate( lf = (l-ml)/sl ) %>% 
  group_by( idall , sname ) %>% mutate( df = (max(lf)-min(lf))/nbin ) %>% 
  mutate( b = as.integer( (lf - min(lf) )/df )  ) %>% 
  ungroup() %>%  group_by( idall, b, df, sname ) %>% 
  summarise(  lf = mean(lf),  n = n() , k = mean(k)) %>% 
  ungroup() %>% group_by(idall) %>%  
  mutate( p = n / sum(n) / df ) %>%  ungroup()

gammalog  = function(x, k = 1.7) { 10^( ( k*trigamma(k)*x - exp( sqrt(trigamma(k))*x+ digamma(k)) ) - log(gamma(k)) + k*digamma(k) + log10(exp(1)) ) }
p2 = dh %>%  mutate( xx = runif(dim(dh)[1])  ) %>%   
  arrange(xx, sname) %>% ggplot() + 
  aes(x =  lf/ sqrt(2) , y = 10^( log10(p)  ),
      color = as.factor(sname) , shape = as.factor(sname) ) + 
  geom_point( size = 1, alpha = 0.5 ) +   
  geom_function( fun = gammalog, color = "black", size = 0.5 ) +
  scalecols + scaleshapes +
  scale_x_continuous( "Rescaled log\nrelative abundance") +
  scale_y_log10( "Probability\ndensity", 
                 limits = c(0.001,0.8), 
                 labels = fancy_scientific )  +
  theme_bw() + 
  theme( panel.background = element_blank(),
         plot.background = element_blank(),
         panel.grid = element_blank(),
         legend.position = "none", 
         axis.text = element_text(size = 8, color = 'black'),
         axis.title = element_text(size = 10, color = 'black')) 
# ggsave( filename =  "26-gamma_fluct.pdf", p2, width = 2.5, height = 2.4  )

#########################################################################

# plegend = dh %>%  mutate( xx = runif(dim(dh)[1])  ) %>%   
#   arrange(xx, sname) %>% ggplot() + 
#   aes(x =  lf/ sqrt(2) , y = 10^( log10(p)  ),
#       color = as.factor(sname) , shape = as.factor(sname) ) + 
#   geom_point( size = 1, alpha = 0.5 ) +   
#   geom_function( fun = gammalog, color = "black", size = 0.5 ) +
#   scalecols + scaleshapes +
#   scale_x_continuous( "Rescaled log\nrelative abundance") +
#   scale_y_log10( "Probability density", 
#                  limits = c(0.001,0.8), 
#                  labels = fancy_scientific )  +
#   theme_bw() + 
#   theme( panel.background = element_blank(),
#          plot.background = element_blank(),
#          panel.grid = element_blank(),
#          legend.position = "right", 
#          axis.text = element_text(size = 8, color = 'black'),
#          axis.title = element_text(size = 10, color = 'black')) 
# ggsave( filename =  "26-Legend.pdf", plegend, width = 2.5, height = 2.4  )

#########################################################################

pdf(width = 2, height = 4, file = 'MacLaw-gamma_talor_logn.pdf' )
grid.newpage()  
pushViewport(viewport(layout = grid.layout(3, 1))) 
vplayout = function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p2 , vp = vplayout(1, 1))
print(ptaylor_all , vp = vplayout(2, 1))
print(plnind , vp = vplayout(3, 1))
dev.off()





