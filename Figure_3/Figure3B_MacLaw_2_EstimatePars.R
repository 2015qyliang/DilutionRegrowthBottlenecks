
library(tidyverse)
library(ggfortify)
library(lamW)
library(GGally)
library(ggcorrplot)
library(pracma)
library(ggsci)
source("./MacLaw_theme.R")

##################################

datatax = read.table('MacLaw.txt', header = T, sep = '\t')

summarydata = datatax %>% 
  group_by(project_id, classification) %>%
  summarise( n_of_runs = n_distinct(run_id), 
             mean_nreads = mean(nreads)) %>%  
  ungroup() %>%  as.data.frame() %>%  
  filter( n_of_runs > 0  ) %>% 
  mutate(idall = paste(project_id, classification))

summarydata

##################################

names = scat = c("D0", "D1", "D3", "D5", "D6", "D7")
shortnames = as.data.frame(list( idall = summarydata$idall, 
                                 sname = names, scat = scat  ))
proj = datatax %>% 
  group_by( project_id, classification, otu_id ) %>% 
  mutate( tf = mean(count/nreads), o = n(), 
          tvpf = mean( (count^2 - count)/nreads^2 ) ) %>% 
  ungroup() %>%
  group_by( project_id, classification ) %>%  
  mutate(o = o / n_distinct(run_id) ) %>% 
  mutate( f = o*tf, vf = o*tvpf ) %>% 
  mutate(vf = vf - f^2 ) %>%  ungroup() %>% 
  mutate(idall = paste(project_id, classification)) %>% 
  select( -tf, -tvpf )

proj = proj %>% left_join( shortnames)  

##################################

############ PARAMETERS OF THE GAMMA DISTRIBUTION
gamma_pars = proj %>% select( idall, sname, otu_id, o, f, vf ) %>% 
  mutate( cv = sqrt(vf/f^2) ) %>% distinct() %>% 
  mutate( beta = 1./cv^2, theta = f/beta )

gamma_pars %>%  ggplot() + mytheme +
  aes( x = f, y = vf ) + 
  geom_point( alpha = 0.5 ) + 
  facet_wrap(  ~ sname ) +
  scale_x_log10() + scale_y_log10()


gamma_pars %>% filter(f > 2*10^-5) %>% 
  group_by(sname) %>% 
  summarise( cor(beta, log(f)) ) %>%  ungroup()
### think about this correlation (sampling) and variability of cv/beta 
# (could be able to say something from the model)

##################################################################

############# LOGNORMAL PARAMETERS
nbin = 20
cutoffs = list( sname = c("D0", "D1", "D3", "D5", "D6", "D7"),
                 c = rep(-18, 6) ) %>% 
  as.data.frame()

fun_erf = function( mu, c, m1, m2  ){
  sigma <- sqrt(-c*m1 + m2 + c*mu - m1*mu )
  x <- (c-mu)/sigma/sqrt(2.)
  f <- (mu-m1)*erfc(x) + exp(- x^2) * sqrt(2/pi) * sigma 
  return(f)
}

estimatemean_f = function( c, m1, m2 ){
  mumin <- (c*m1 - m2)/(c - m1)
  muest <- uniroot(fun_erf, c = c, m1 = m1, m2 = m2, 
                   interval = c(-50,  mumin-0.001), tol = 0.0001 )$root
  return( muest )
}

estimate_mean = gamma_pars %>% 
  mutate(lf = log(f)) %>%
  left_join(cutoffs) %>%  
  ungroup() %>% filter( lf > c ) %>% 
  group_by(idall, sname) %>% 
  summarise( c = mean(c),  m1 = mean(lf) ,m2 = mean(lf^2), 
             ns_obs = n_distinct(otu_id), nf = sum(f)  )  %>%   
  ungroup() %>% rowwise() %>% 
  mutate( mu = estimatemean_f(c,m1,m2) ) %>% 
  mutate( sigma = sqrt(-c*m1 + m2 + c*mu - m1*mu ) ) %>% 
  ungroup() %>%  as.data.frame() %>% 
  mutate( stot = 2*ns_obs / erfc( (c-mu)/sigma/sqrt(2)  ) ) # , np_tot_r = 1./exp(mu+sigma^2/2.)/nf )

statp =  gamma_pars %>% 
  mutate( lf = log(f ) , cutoff = -100  ) %>% 
  group_by( idall  ) %>% 
  mutate( df = (max(lf)-min(lf))/nbin ) %>% 
  mutate( b = as.integer( (lf - min(lf) )/df )  ) %>% 
  ungroup() %>%  group_by( idall, b, df ) %>% 
  summarise(  lf = mean(lf), cutoff = mean(cutoff), n = n() ) %>% 
  ungroup() %>% group_by(idall) %>%  
  mutate( p = n / sum(n) / df ) %>%  ungroup()

simple_par = function(x) { - x^2   }  

p1 = statp %>%  mutate( xx = runif(dim(statp)[1])  )  %>% 
  filter(lf > cutoff) %>% 
  left_join(estimate_mean) %>% 
  arrange(xx, sname) %>% 
  ggplot() + mytheme +
  stat_function( fun = simple_par, color = "black", size = 1 ) +
  aes(x = ( (lf-mu)/sqrt(2*sigma^2) ) ,
    y = 10^( log(p) - log( 0.5 * erfc( (mu-c)/sqrt(2)/sigma ) ) + 0.5*log(2.*pi) ),
    color = as.factor(sname) , shape = as.factor(sname) ) + 
  geom_point( size = 3 , stroke = 1) +
  geom_vline( aes( xintercept = (c-mu)/sqrt(2*sigma^2)  )  ) +
  geom_vline( aes( xintercept = 0  )  ) +
  scalecols + scaleshapes +
  scale_x_continuous( "Rescaled log average abundance" ) +
  scale_y_log10( "Fraction of Species"   )  +
  theme( legend.position = "none" ) +    facet_wrap( ~ sname )
p1

mean_pars  = estimate_mean %>% 
  select( idall, sname, nf, c, mu, sigma, stot  ) %>% 
  left_join( gamma_pars %>% 
               filter(f > 10^-5, vf > 0 ) %>% 
               group_by(sname) %>% 
               summarise( mbeta = mean(beta) ) %>%  
               ungroup() )

save( proj, gamma_pars, mean_pars , file = "MacLawDataestimate.Rdata" )

########################################################################
