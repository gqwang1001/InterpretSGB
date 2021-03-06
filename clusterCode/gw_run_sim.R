setwd("/SFS/user/ctc/waguoqi2/InterpretSGB/")
dirResults = "/SFS/scratch/waguoqi2/InterPretSGBResults/"

source("simplified.Tree.R")
source("SubgroupBoost_local.R")

temp = commandArgs(trailingOnly = TRUE)
# temp=c(1, 1,1)
seeds = as.integer(temp[1])
isim = as.integer(temp[2])
nrep = as.integer(temp[3])

simList = paste0("sim", c(105, 106, 107, 101, 102, 103, 104, 32, 42, 43, 62, 7, 8))
simIdx = simList[isim]
print(c(seeds, simIdx))

saveDir = paste0(dirResults, simIdx)
if (!dir.exists(saveDir))
  dir.create(saveDir, recursive = T)
#saveRDS(object = seeds, file = paste0(saveDir, "/test.rds"))

library(xgboost)
library(caret)

for (iseeds in (seeds - 1) * nrep + (1:nrep)) {
  
  # if (!file.exists(paste0(saveDir, "/", iseeds, ".rds"))) {
    set.seed(iseeds)
    source(paste0("EricCode/MySimData_", simIdx, ".R"))
   # model <- SubgroupBoost.RMST_local(data.simulation[[1]])
    model = readRDS(paste0(saveDir, "/NoTop_",iseeds,".rds"))$xgbModel

    spltree = simplified.Tree(
      model = model,
      datalist = data.simulation,
      cutoff = 0.95,
      # plot.name = paste0("Results/", simIdx,"_simplified_tree.png"),
      seed = iseeds,
      top3=F
    )
    saveRDS(object = spltree, file = paste0(saveDir, "/NoTop_", iseeds, ".rds"))
  # }
}

rm(list = ls())

iseeds=1
spltree = readRDS(paste0(saveDir, "/NoTop_",iseeds,".rds"))
spltree$evalLoss.lm.train
spltree$evalLoss.spltree.train

