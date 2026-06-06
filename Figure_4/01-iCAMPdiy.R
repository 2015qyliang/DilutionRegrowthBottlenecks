# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

##########################################################################
# 
library(NST)
library(iCAMP)
library(ape)

rm(list = ls())
# 
rand.time = 100
nworker = 8 #
memory.G = 26 # 

# 
comm = t(read.table("OTUtable.txt", header = T, sep = "\t",
                    row.names = 1, as.is = T, stringsAsFactors = F ))
treat = read.table("OTUsamplesGroup.txt", header = T, sep = "\t",
                   row.names = 1, stringsAsFactors = F )
tree = read.tree('OTUFastTree.txt')
taxon = read.table('OTUTax.txt', header = TRUE, sep = "\t", 
                   row.names = 1, stringsAsFactors = FALSE)

# 
prefix = save.wd = "PD"
# calculate pairwise phylogenetic distance matrix
if(!dir.exists(save.wd)){dir.create(save.wd)}
if(!file.exists(paste0(save.wd, "/pd.desc"))) {
  pd.big = pdist.big(tree = tree, wd = save.wd, nworker = nworker, memory.G = memory.G)
} else {
  pd.big = list()
  pd.big$tip.label = read.csv(paste0(save.wd, "/pd.taxon.name.csv"),
                              row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd = save.wd
  pd.big$pd.file = 'pd.desc'
  pd.big$pd.name.file = 'pd.taxon.name.csv'
}

##########################################################################
# without omitting small bins.
bin.size.limit = 12  
sig.index="Confidence" 
icres = icamp.big(comm = comm,
                  pd.desc = pd.big$pd.file,
                  pd.spname = pd.big$tip.label,
                  pd.wd = pd.big$pd.wd,
                  rand = rand.time,
                  tree = tree,
                  prefix = prefix,
                  output.wd = save.wd,
                  nworker = nworker,
                  memory.G = memory.G,
                  unit.sum = rowSums(comm),
                  ds = 0.2, pd.cut = NA, 
                  sp.check = TRUE,
                  phylo.rand.scale = "within.bin",
                  taxa.rand.scale = "across.all",
                  phylo.metric = "bMPD",
                  sig.index = 'Confidence',
                  bin.size.limit = 12,
                  rtree.save = FALSE,
                  detail.save = TRUE,
                  qp.save = FALSE,
                  detail.null = FALSE,
                  ignore.zero = TRUE,
                  correct.special = TRUE,
                  special.method = "depend",
                  ses.cut = 1.96, rc.cut = 0.95, conf.cut = 0.975,
                  omit.option = "no", meta.ab = NULL)

saveRDS(icres, 'icres.rds')

##########################################################################
# booting
icamp.result = icres$CbMPDiCBraya
icboot = icamp.boot(icamp.result = icamp.result,
                    treat = treat,
                    rand.time = rand.time,
                    compare = TRUE,
                    silent = FALSE,
                    between.group = TRUE,
                    ST.estimation = TRUE)
saveRDS(icboot, 'PD.Boot.rds')
save(icboot, file = paste0(prefix,".Boot.",colnames(treat),".rda"))
write.csv(icboot$summary, 
          file = paste0(prefix,".BootSummary.",colnames(treat),".csv"),
          row.names = FALSE)
write.csv(icboot$compare, 
          file = paste0(prefix, ".Compare.", colnames(treat),".csv"),
          row.names = FALSE)

##########################################################################
# iCAMP bin level statistics
icbin = icamp.bins(icamp.detail = icres$detail,
                   treat = treat,
                   clas = taxon,
                   silent = FALSE,
                   boot = TRUE,
                   rand.time = rand.time,
                   between.group = TRUE)
save(icbin,file = paste0(prefix,".Summary.rda")) 
write.csv(icbin$Pt,file = paste0(prefix,".ProcessImportance_EachGroup.csv"),
          row.names = FALSE)
write.csv(icbin$Ptk,file = paste0(prefix,".ProcessImportance_EachBin_EachGroup.csv"),
          row.names = FALSE)
write.csv(icbin$Ptuv,file = paste0(prefix,".ProcessImportance_EachTurnover.csv"),
          row.names = FALSE)
write.csv(icbin$BPtk,file = paste0(prefix,".BinContributeToProcess_EachGroup.csv"),
          row.names = FALSE)
write.csv(data.frame(ID = rownames(icbin$Class.Bin),
                     icbin$Class.Bin,
                     stringsAsFactors = FALSE),
          file = paste0(prefix,".Taxon_Bin.csv"), row.names = FALSE)
write.csv(icbin$Bin.TopClass,file = paste0(prefix,".Bin_TopTaxon.csv"),
          row.names = FALSE)

##########################################################################
