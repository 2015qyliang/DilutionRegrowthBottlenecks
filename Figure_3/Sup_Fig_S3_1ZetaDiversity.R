

library(zetadiv)
library(tidyverse)

####################################################################

data.spec = as.matrix(t(read.table('OTUtable.txt', 
                                   header = T, sep = '\t', row.names = 1)))
data.spec[data.spec != 0] = 1

D0 = as.data.frame(data.spec[str_detect(rownames(data.spec), 'D0'), 
                             colSums(data.spec) != 0])
D0zeta = Zeta.decline.mc(D0,  
                         xy = NULL,  
                         orders = 1:14, 
                         sam = 1000,  
                         confint.level = 0.95, 
                         rescale = T, 
                         normalize = T, 
                         empty.row = "empty", 
                         plot = TRUE 
                         )
saveRDS(D0zeta, 'Sup_Fig_S3-D0zeta.decline.mc.rds')

D1 = as.data.frame(data.spec[str_detect(rownames(data.spec), 'D1'), 
                             colSums(data.spec) != 0])
D1zeta = Zeta.decline.mc(D1,  
                         xy = NULL,  
                         orders = 1:14, 
                         sam = 1000,  
                         confint.level = 0.95, 
                         rescale = T, 
                         normalize = T, 
                         empty.row = "empty", 
                         plot = TRUE 
)
saveRDS(D1zeta, 'Sup_Fig_S3-D1zeta.decline.mc.rds')

D3 = as.data.frame(data.spec[str_detect(rownames(data.spec), 'D3'), 
                             colSums(data.spec) != 0])
D3zeta = Zeta.decline.mc(D3,  
                         xy = NULL,  
                         orders = 1:20, 
                         sam = 1000,  
                         confint.level = 0.95, 
                         rescale = T, 
                         normalize = T, 
                         empty.row = "empty", 
                         plot = TRUE 
)
saveRDS(D3zeta, 'Sup_Fig_S3-D3zeta.decline.mc.rds')

D5 = as.data.frame(data.spec[str_detect(rownames(data.spec), 'D5'), 
                             colSums(data.spec) != 0])
D5zeta = Zeta.decline.mc(D5,  
                         xy = NULL,  
                         orders = 1:20, 
                         sam = 1000,  
                         confint.level = 0.95, 
                         rescale = T, 
                         normalize = T, 
                         empty.row = "empty", 
                         plot = TRUE 
)
saveRDS(D5zeta, 'Sup_Fig_S3-D5zeta.decline.mc.rds')

D6 = as.data.frame(data.spec[str_detect(rownames(data.spec), 'D6'), 
                             colSums(data.spec) != 0])
D6zeta = Zeta.decline.mc(D6,  
                         xy = NULL,  
                         orders = 1:40, 
                         sam = 1000,  
                         confint.level = 0.95, 
                         rescale = T, 
                         normalize = T, 
                         empty.row = "empty", 
                         plot = TRUE 
)
saveRDS(D6zeta, 'Sup_Fig_S3-D6zeta.decline.mc.rds')

D7 = as.data.frame(data.spec[str_detect(rownames(data.spec), 'D7'), 
                             colSums(data.spec) != 0])
D7zeta = Zeta.decline.mc(D7,  
                         xy = NULL,  
                         orders = 1:40, 
                         sam = 1000,  
                         confint.level = 0.95, 
                         rescale = T, 
                         normalize = T, 
                         empty.row = "empty", 
                         plot = TRUE 
)
saveRDS(D7zeta, 'Sup_Fig_S3-D7zeta.decline.mc.rds')

####################################################################

summary(D0zeta$zeta.pl)

summary(D0zeta$zeta.exp)
