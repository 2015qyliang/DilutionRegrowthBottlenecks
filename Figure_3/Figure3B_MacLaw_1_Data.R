
library(tidyverse)
library(reshape2)

otudf = read.table('OTUtable.txt', header = T, sep = '\t', row.names = 1)
nreads = colSums(otudf)

otudf = read.table('OTUtable.txt', header = T, sep = '\t')

lawdf = melt(otudf, measure.vars = colnames(otudf)[2:ncol(otudf)])
colnames(lawdf) = c('otu_id', 'sample_id', 'count')
lawdf$run_id = lawdf$sample_id
lawdf$project_id = substr(lawdf$sample_id, 1, 2)
lawdf$classification = lawdf$project_id

lawdf = lawdf[lawdf$count != 0, ]

dfcount = data.frame(sample_id = names(nreads), nreads)

lawdf = merge(lawdf, dfcount, by = 'sample_id', all.x = T)

colnames(lawdf)
lawdf = lawdf[, c("otu_id", "count", "project_id", "sample_id", 
                  "run_id", "nreads", "classification" )]

write.table(lawdf, 'MacLaw.txt', append = F, quote = F, 
            sep = '\t', row.names = F, col.names = T)







