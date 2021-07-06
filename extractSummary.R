setwd("/SFS/scratch/waguoqi2/InterPretSGBResults")
# in the Results folder
library(parallel)
simList = paste0("sim", c(101,102,105,106)); scen = paste0("scenario",c(7:10)); cases="MoreFeatures"
simList = paste0("sim", c(32, 62, 43, 42, 7, 8)); scen =paste0("scenario",c(1:6)); cases="OrginalPaper"
varNames = paste0("z", 1:5)
NREP =500
source("/SFS/user/ctc/waguoqi2/InterpretSGB/var.rank.R")
source("/SFS/user/ctc/waguoqi2/InterpretSGB/var.rank.xgb.all.R")

savedir = "/SFS/scratch/waguoqi2/InterPretSGBResults/NoTopResults/"
# savedir = "/SFS/scratch/waguoqi2/InterPretSGBResults/Top3Results/"
if(!dir.exists(savedir)) {
  dir.create(paste0(savedir, "results/"), recursive = T)
  }

DataReshape = function(isim, savedir){
  
  acc.xgboost = rep(0, NREP)
  acc.splTree = rep(0, NREP)
  sen.xgboost = rep(0, NREP)
  sen.splTree = rep(0, NREP)  
  spe.xgboost = rep(0, NREP)
  spe.splTree = rep(0, NREP)
  acc.splTree.gain = rep(0, NREP)
  sen.splTree.gain = rep(0, NREP)  
  spe.splTree.gain = rep(0, NREP)
  
  RMSE.lmVsfull = RMSE.spltreeVsFull = rep(0, NREP)
  dice.lmVSfull = dice.spltreeVsFULL = rep(0, NREP)
  
  ORT.lm = ORT.spltree = rep(0,NREP)
  
  ranking.s1 = matrix(rep(0, NREP*4), ncol=4)
  ranking.s2 = matrix(rep(0, NREP*4), ncol=4)
  
  ranking.all = matrix(NaN, nrow = NREP, ncol = length(varNames))
  ranking.all.spltree = matrix(NaN, nrow = NREP, ncol = length(varNames))
  ranking.all.spltree.gain = matrix(NaN, nrow = NREP, ncol = length(varNames))
  
  for (ir in 1:NREP) {
    spltree = readRDS(paste0("/SFS/scratch/waguoqi2/InterPretSGBResults/", isim, "/NoTop_", ir, ".rds"))
    acc.xgboost[ir] = spltree$summaryPred.xgb$accuracy
    sen.xgboost[ir] = spltree$summaryPred.xgb$sensitivity 
    spe.xgboost[ir] = spltree$summaryPred.xgb$specificity
    
    acc.splTree[ir] = spltree$summaryPred.splTree.shap$accuracy
    sen.splTree[ir] = spltree$summaryPred.splTree.shap$sensitivity
    spe.splTree[ir] = spltree$summaryPred.splTree.shap$specificity
    
    acc.splTree.gain[ir] = spltree$summaryPred.splTree.gain$accuracy
    sen.splTree.gain[ir] = spltree$summaryPred.splTree.gain$sensitivity
    spe.splTree.gain[ir] = spltree$summaryPred.splTree.gain$specificity
    
    
    RMSE.lmVsfull[ir] = spltree$RMSE.lm
    RMSE.spltreeVsFull[ir] = spltree$RMSE.spltree
    
    dice.lmVSfull[ir] = spltree$acc.comp.lm
    dice.spltreeVsFULL[ir] = spltree$acc.comp.splTreeVsxgb
    
    ORT.lm[ir] = spltree$evalLoss.lm
    ORT.spltree[ir] = spltree$evalLoss.spltree
    
    ranking.s1[ir,] = do.call("c", var.rank(spltree,"s1"))
    ranking.s2[ir,] = do.call("c", var.rank(spltree,"s2"))
    
    ranking.all[ir,] = var.rank.xgb.all(spltree, varNames)$xgb.shap
    ranking.all.spltree[ir, ] = var.rank.xgb.all(spltree, varNames)$sptrees.xgb
    ranking.all.spltree.gain[ir, ] = var.rank.xgb.all(spltree, varNames)$sptrees.gain
    
  }
  
  acc.df = data.frame(acc.xgboost=acc.xgboost,
                      acc.splTree=acc.splTree,
                      acc.splTree.gain = acc.splTree.gain)
  sen.df = data.frame(sen.xgboost=sen.xgboost,
                      sen.splTree=sen.splTree,
                      sen.splTree.gain=sen.splTree.gain)  
  spe.df = data.frame(spe.xgboost=spe.xgboost,
                      spe.splTree=spe.splTree,
                      spe.splTree.gain=spe.splTree.gain) 
  rmse.df = data.frame(RMSE.lm = RMSE.lmVsfull,
                       RMSE.spltree = RMSE.spltreeVsFull)
  dice.df = data.frame(dice.lm = dice.lmVSfull,
                       dice.spltree = dice.spltreeVsFULL)
  
  ORT.df = data.frame(ORT.lm = ORT.lm, ORT.spltree = ORT.spltree)
  
  colnames(ranking.s1) = c("splTrees_SHAP","SplTrees_GAIN", "XGB_SHAP", "XGB_GAIN")
  colnames(ranking.s2) = c("splTrees_SHAP","SplTrees_GAIN", "XGB_SHAP", "XGB_GAIN")
  colnames(ranking.all) = c(varNames)
  colnames(ranking.all.spltree) = c(varNames)
  colnames(ranking.all.spltree.gain) = c(varNames)
  
  write.csv(acc.df, file = paste0(savedir,"results/acc_", isim,".csv"))
  write.csv(sen.df, file = paste0(savedir,"results/sen_", isim,".csv"))
  write.csv(spe.df, file = paste0(savedir,"results/spe_", isim,".csv"))
  
  write.csv(as.data.frame(ranking.s1), paste0(savedir,"results/Ranking_S1_", isim,".csv"))
  write.csv(as.data.frame(ranking.s2), paste0(savedir,"results/Ranking_S2_", isim,".csv"))
  
  write.csv(rmse.df, paste0(savedir,"results/rmse_comparison_", isim,".csv"))
  write.csv(dice.df, paste0(savedir,"results/dice_comparison_", isim,".csv"))
  write.csv(ORT.df, paste0(savedir,"results/ORT_comparison_", isim,".csv"))
  
  
  write.csv(as.data.frame(ranking.all), paste0(savedir,"results/rankingAll_", isim,".csv"))
  write.csv(as.data.frame(ranking.all.spltree), paste0(savedir,"results/rankingAll_spltree_", isim,".csv"))
  write.csv(as.data.frame(ranking.all.spltree.gain), paste0(savedir,"results/rankingAll_spltree_gain_", isim,".csv"))
  
}


mclapply(simList, DataReshape, savedir = savedir, mc.cores = length(simList))


for (isim in simList){
  print(isim)
  acc.xgboost = rep(0, NREP)
  acc.splTree = rep(0, NREP)
  sen.xgboost = rep(0, NREP)
  sen.splTree = rep(0, NREP)  
  spe.xgboost = rep(0, NREP)
  spe.splTree = rep(0, NREP)
  acc.splTree.gain = rep(0, NREP)
  sen.splTree.gain = rep(0, NREP)  
  spe.splTree.gain = rep(0, NREP)
  
  RMSE.lmVsfull = RMSE.spltreeVsFull = rep(0, NREP)
  dice.lmVSfull = dice.spltreeVsFULL = rep(0, NREP)
  
  ORT.lm = ORT.spltree = rep(0,NREP)

  ranking.s1 = matrix(rep(0, NREP*4), ncol=4)
  ranking.s2 = matrix(rep(0, NREP*4), ncol=4)
  
  ranking.all = matrix(NaN, nrow = NREP, ncol = length(varNames))
  ranking.all.spltree = matrix(NaN, nrow = NREP, ncol = length(varNames))
  ranking.all.spltree.gain = matrix(NaN, nrow = NREP, ncol = length(varNames))
  
  for (ir in 1:NREP) {
    spltree = readRDS(paste0("/SFS/scratch/waguoqi2/InterPretSGBResults/", isim, "/NoTop_", ir, ".rds"))
    acc.xgboost[ir] = spltree$summaryPred.xgb$accuracy
    sen.xgboost[ir] = spltree$summaryPred.xgb$sensitivity 
    spe.xgboost[ir] = spltree$summaryPred.xgb$specificity
    
    acc.splTree[ir] = spltree$summaryPred.splTree.shap$accuracy
    sen.splTree[ir] = spltree$summaryPred.splTree.shap$sensitivity
    spe.splTree[ir] = spltree$summaryPred.splTree.shap$specificity
    
    acc.splTree.gain[ir] = spltree$summaryPred.splTree.gain$accuracy
    sen.splTree.gain[ir] = spltree$summaryPred.splTree.gain$sensitivity
    spe.splTree.gain[ir] = spltree$summaryPred.splTree.gain$specificity
    
    
    RMSE.lmVsfull[ir] = spltree$RMSE.lm
    RMSE.spltreeVsFull[ir] = spltree$RMSE.spltree
    
    dice.lmVSfull[ir] = spltree$acc.comp.lm
    dice.spltreeVsFULL[ir] = spltree$acc.comp.splTreeVsxgb

    ORT.lm[ir] = spltree$evalLoss.lm
    ORT.spltree[ir] = spltree$evalLoss.spltree
    
    ranking.s1[ir,] = do.call("c", var.rank(spltree,"s1"))
    ranking.s2[ir,] = do.call("c", var.rank(spltree,"s2"))
    
    ranking.all[ir,] = var.rank.xgb.all(spltree, varNames)$xgb.shap
    ranking.all.spltree[ir, ] = var.rank.xgb.all(spltree, varNames)$sptrees.xgb
    ranking.all.spltree.gain[ir, ] = var.rank.xgb.all(spltree, varNames)$sptrees.gain
    
  }
  
  acc.df = data.frame(acc.xgboost=acc.xgboost,
                      acc.splTree=acc.splTree,
                      acc.splTree.gain = acc.splTree.gain)
  sen.df = data.frame(sen.xgboost=sen.xgboost,
                      sen.splTree=sen.splTree,
                      sen.splTree.gain=sen.splTree.gain)  
  spe.df = data.frame(spe.xgboost=spe.xgboost,
                      spe.splTree=spe.splTree,
                      spe.splTree.gain=spe.splTree.gain) 
  rmse.df = data.frame(RMSE.lm = RMSE.lmVsfull,
                       RMSE.spltree = RMSE.spltreeVsFull)
  dice.df = data.frame(dice.lm = dice.lmVSfull,
                       dice.spltree = dice.spltreeVsFULL)
  
  ORT.df = data.frame(ORT.lm = ORT.lm, ORT.spltree = ORT.spltree)
  
  colnames(ranking.s1) = c("splTrees_SHAP","SplTrees_GAIN", "XGB_SHAP", "XGB_GAIN")
  colnames(ranking.s2) = c("splTrees_SHAP","SplTrees_GAIN", "XGB_SHAP", "XGB_GAIN")
  colnames(ranking.all) = c(varNames)
  colnames(ranking.all.spltree) = c(varNames)
  colnames(ranking.all.spltree.gain) = c(varNames)
  
  write.csv(acc.df, file = paste0(savedir,"results/acc_", isim,".csv"))
  write.csv(sen.df, file = paste0(savedir,"results/sen_", isim,".csv"))
  write.csv(spe.df, file = paste0(savedir,"results/spe_", isim,".csv"))
  
  write.csv(as.data.frame(ranking.s1), paste0(savedir,"results/Ranking_S1_", isim,".csv"))
  write.csv(as.data.frame(ranking.s2), paste0(savedir,"results/Ranking_S2_", isim,".csv"))
  
  write.csv(rmse.df, paste0(savedir,"results/rmse_comparison_", isim,".csv"))
  write.csv(dice.df, paste0(savedir,"results/dice_comparison_", isim,".csv"))
  write.csv(ORT.df, paste0(savedir,"results/ORT_comparison_", isim,".csv"))
  
  
  write.csv(as.data.frame(ranking.all), paste0(savedir,"results/rankingAll_", isim,".csv"))
  write.csv(as.data.frame(ranking.all.spltree), paste0(savedir,"results/rankingAll_spltree_", isim,".csv"))
  write.csv(as.data.frame(ranking.all.spltree.gain), paste0(savedir,"results/rankingAll_spltree_gain_", isim,".csv"))
  
}



mean.acc.xgboost = rep(0, length(simList))
mean.acc.splTree = rep(0, length(simList))
sd.acc.xgboost = rep(0, length(simList))
sd.acc.splTree = rep(0, length(simList))
acc.df.xgboost = matrix(nrow = NREP, ncol = length(simList))
acc.df.splTree = matrix(nrow = NREP, ncol = length(simList))
acc.df.splTree.gain = matrix(nrow = NREP, ncol = length(simList))

sen.df.xgboost = matrix(nrow = NREP, ncol = length(simList))
sen.df.splTree = matrix(nrow = NREP, ncol = length(simList))
sen.df.splTree.gain = matrix(nrow = NREP, ncol = length(simList))

spe.df.xgboost = matrix(nrow = NREP, ncol = length(simList))
spe.df.splTree = matrix(nrow = NREP, ncol = length(simList))
spe.df.splTree.gain = matrix(nrow = NREP, ncol = length(simList))


rmse.df.lm = rmse.df.splTree = matrix(nrow = NREP, ncol = length(simList))
dice.df.lm = dice.df.splTree = matrix(nrow = NREP, ncol = length(simList))
ORT.df.lm = ORT.df.splTree = matrix(nrow = NREP, ncol = length(simList))

for (i in 1:length(simList)) {
  isim = simList[i]
  acc.df = read.csv(paste0(savedir,"results/acc_", isim, ".csv"))
  acc.df.xgboost[,i] = acc.df$acc.xgboost
  acc.df.splTree[,i] = acc.df$acc.splTree
  acc.df.splTree.gain[,i] = acc.df$acc.splTree.gain
  
  mean.acc.xgboost[i] = mean(acc.df$acc.xgboost)
  mean.acc.splTree[i] = mean(acc.df$acc.splTree)
  sd.acc.xgboost[i] = sd(acc.df$acc.xgboost)
  sd.acc.splTree[i] = sd(acc.df$acc.splTree)
  
  sen.df = read.csv(paste0(savedir,"results/sen_", isim,".csv"))
  spe.df = read.csv(paste0(savedir,"results/spe_", isim,".csv"))
  sen.df.xgboost[,i] = sen.df$sen.xgboost
  sen.df.splTree[,i] = sen.df$sen.splTree
  sen.df.splTree.gain[,i] = sen.df$sen.splTree.gain
  
  spe.df.xgboost[,i] = spe.df$spe.xgboost
  spe.df.splTree[,i] = spe.df$spe.splTree
  spe.df.splTree.gain[,i] = spe.df$spe.splTree.gain
  
  rmse.df = read.csv(paste0(savedir,"results/rmse_comparison_", isim,".csv"))
  rmse.df.lm[,i] = rmse.df$RMSE.lm
  rmse.df.splTree[,i] = rmse.df$RMSE.spltree
  
  dice.df = read.csv(paste0(savedir,"results/dice_comparison_", isim,".csv"))
  dice.df.lm[,i] = dice.df$dice.lm
  dice.df.splTree[,i] = dice.df$dice.spltree
  
  ORT.df = read.csv(paste0(savedir,"results/ORT_comparison_", isim,".csv"))
  ORT.df.lm[,i] = ORT.df$ORT.lm
  ORT.df.splTree[,i] = ORT.df$ORT.spltree
  
}

comparison.df = rbind(mean.acc.xgboost, sd.acc.xgboost, mean.acc.splTree,sd.acc.splTree)
colnames(comparison.df)=simList
write.csv(comparison.df, paste0(savedir,"results/comparison_splTreeVSXgboost.csv"))

colnames(acc.df.xgboost)=scen
colnames(acc.df.splTree)=scen
colnames(acc.df.splTree.gain) = scen

acc.df.xgboost.long = reshape2::melt(acc.df.xgboost,value.name = "Accuracy", variable.name="Simulation")
acc.df.xgboost.long$Method = "XGBoost"
head(acc.df.xgboost.long)
nrow(acc.df.xgboost.long)
acc.df.splTree.long = reshape2::melt(acc.df.splTree,value.name = "Accuracy", variable.name="Simulation")
acc.df.splTree.long$Method = "Simplified Tree By SHAP"

acc.df.splTree.gain.long = reshape2::melt(acc.df.splTree.gain,value.name = "Accuracy", variable.name="Simulation")
acc.df.splTree.gain.long$Method = "Simplified Tree By GAIN"

acc.df.long=rbind(acc.df.xgboost.long, acc.df.splTree.long, acc.df.splTree.gain.long)
tail(acc.df.long)

colnames(acc.df.long) = c("id", "Simulation", "Accuracy", "Method")
write.csv(acc.df.long, paste0(savedir,"results/Accuracy_comparision_long_", cases, ".csv"))

library(ggplot2)
acc.df.long = read.csv(paste0(savedir,"results/Accuracy_comparision_long_", cases, ".csv"))
acc.df.long$Simulation = factor(acc.df.long$Simulation, levels = scen)

p1 = ggplot(as.data.frame(acc.df.long),aes(x=Simulation, y=Accuracy))+
  geom_boxplot(aes(fill=Method))
ggsave(paste0(savedir,"results/Accuracy_comparison_XgboostVsSplTree_", cases, ".png"), p1, width = 14)

colnames(sen.df.xgboost)=scen
colnames(sen.df.splTree)=scen
colnames(sen.df.splTree.gain) = scen

sen.df.xgboost.long = reshape2::melt(sen.df.xgboost,value.name = "Sensitivity", variable.name="Simulation")
sen.df.xgboost.long$Method = "XGBoost"
head(sen.df.xgboost.long)
nrow(sen.df.xgboost.long)
sen.df.splTree.long = reshape2::melt(sen.df.splTree,value.name = "Sensitivity", variable.name="Simulation")
sen.df.splTree.long$Method = "Simplified Tree By SHAP"
sen.df.splTree.gain.long = reshape2::melt(sen.df.splTree.gain,value.name = "Sensitivity", variable.name="Simulation")
sen.df.splTree.gain.long$Method = "Simplified Tree By GAIN"

sen.df.long=rbind(sen.df.xgboost.long, sen.df.splTree.long, sen.df.splTree.gain.long)
head(sen.df.long)
colnames(sen.df.long) = c("id", "Simulation", "Sensitivity", "Method")
write.csv(sen.df.long, paste0(savedir,"results/Sensitivity_comparision_long_", cases, ".csv"))

library(ggplot2)
sen.df.long = read.csv(paste0(savedir,"results/Sensitivity_comparision_long_", cases, ".csv"))
sen.df.long$Simulation = factor(sen.df.long$Simulation, levels = scen)

p1 = ggplot(as.data.frame(sen.df.long),aes(x=Simulation, y=Sensitivity))+
  geom_boxplot(aes(fill=Method))
ggsave(paste0(savedir,"results/Sensitivity_comparison_XgboostVsSplTree_", cases, ".png"),p1, width = 14)


colnames(spe.df.xgboost)=scen
colnames(spe.df.splTree)=scen
colnames(spe.df.splTree.gain) = scen

spe.df.xgboost.long = reshape2::melt(spe.df.xgboost,value.name = "Specificity", variable.name="Simulation")
spe.df.xgboost.long$Method = "XGBoost"
head(spe.df.xgboost.long)
nrow(spe.df.xgboost.long)
spe.df.splTree.long = reshape2::melt(spe.df.splTree,value.name = "Specificity", variable.name="Simulation")
spe.df.splTree.long$Method = "Simplified Tree By SHAP"
spe.df.splTree.gain.long = reshape2::melt(spe.df.splTree.gain,value.name = "Specificity", variable.name="Simulation")
spe.df.splTree.gain.long$Method = "Simplified Tree By GAIN"

spe.df.long=rbind(spe.df.xgboost.long, spe.df.splTree.long,spe.df.splTree.gain.long)
head(spe.df.long)
colnames(spe.df.long) = c("id", "Simulation", "Specificity", "Method")
write.csv(spe.df.long, paste0(savedir,"results/Specificity_comparision_long_", cases, ".csv"))

library(ggplot2)
spe.df.long = read.csv(paste0(savedir,"results/Specificity_comparision_long_", cases, ".csv"))
spe.df.long$Simulation = factor(spe.df.long$Simulation, levels = scen)

p1 = ggplot(as.data.frame(spe.df.long),aes(x=Simulation, y=Specificity))+
  geom_boxplot(aes(fill=Method))
ggsave(paste0(savedir,"results/Specificity_comparison_XgboostVsSplTree_", cases, ".png"),p1, width = 14)

# comparing linear model and decision trees -------------------------------

colnames(rmse.df.lm)=scen
colnames(rmse.df.splTree)=scen

rmse.df.lm.long = reshape2::melt(rmse.df.lm,value.name = "RMSE", variable.name="Simulation")
rmse.df.lm.long$Method = "Linear Model"

rmse.df.splTree.long = reshape2::melt(rmse.df.splTree,value.name = "RMSE", variable.name="Simulation")
rmse.df.splTree.long$Method = "Decision Tree"

rmse.df.long=rbind(rmse.df.lm.long, rmse.df.splTree.long)
colnames(rmse.df.long) = c("id", "Simulation", "RMSE", "Method")
write.csv(rmse.df.long, paste0(savedir,"results/RMSE_comparision_long_", cases, ".csv"))

library(ggplot2)
rmse.df.long = read.csv(paste0(savedir,"results/RMSE_comparision_long_", cases, ".csv"))
rmse.df.long$Simulation = factor(rmse.df.long$Simulation, levels = scen)

p1 = ggplot(as.data.frame(rmse.df.long),aes(x=Simulation, y=RMSE))+
  geom_boxplot(aes(fill=Method))
ggsave(paste0(savedir,"results/RMSE_comparison_", cases, ".png"),p1, width = 14)


colnames(dice.df.lm)=scen
colnames(dice.df.splTree)=scen

dice.df.lm.long = reshape2::melt(dice.df.lm,value.name = "DICE", variable.name="Simulation")
dice.df.lm.long$Method = "Linear Model"

dice.df.splTree.long = reshape2::melt(dice.df.splTree,value.name = "DICE", variable.name="Simulation")
dice.df.splTree.long$Method = "Decision Tree"

dice.df.long=rbind(dice.df.lm.long, dice.df.splTree.long)
colnames(dice.df.long) = c("id", "Simulation", "DICE", "Method")
write.csv(dice.df.long, paste0(savedir,"results/dice_comparision_long_", cases, ".csv"))

library(ggplot2)
dice.df.long = read.csv(paste0(savedir,"results/dice_comparision_long_", cases, ".csv"))
dice.df.long$Simulation = factor(dice.df.long$Simulation, levels = scen)

p1 = ggplot(as.data.frame(dice.df.long),aes(x=Simulation, y=DICE))+
  geom_boxplot(aes(fill=Method))
ggsave(paste0(savedir,"results/dice_comparison_", cases, ".png"),p1, width = 14)



colnames(ORT.df.lm)=scen
colnames(ORT.df.splTree)=scen

ORT.df.lm.long = reshape2::melt(ORT.df.lm,value.name = "ORT", variable.name="Simulation")
ORT.df.lm.long$Method = "Linear Model"

ORT.df.splTree.long = reshape2::melt(ORT.df.splTree,value.name = "ORT", variable.name="Simulation")
ORT.df.splTree.long$Method = "Decision Tree"

ORT.df.long=rbind(ORT.df.lm.long, ORT.df.splTree.long)
colnames(ORT.df.long) = c("id", "Simulation", "ORT", "Method")
write.csv(ORT.df.long, paste0(savedir,"results/ORT_comparision_long_", cases, ".csv"))

library(ggplot2)
ORT.df.long = read.csv(paste0(savedir,"results/ORT_comparision_long_", cases, ".csv"))
ORT.df.long$Simulation = factor(ORT.df.long$Simulation, levels = scen)

p1 = ggplot(as.data.frame(ORT.df.long),aes(x=Simulation, y=ORT))+
  geom_boxplot(aes(fill=Method))
ggsave(paste0(savedir,"results/ORT_comparison_", cases, ".png"),p1, width = 14)

# rankings ----------------------------------------------------------------

prop.misclass.s1 = matrix(
  nrow = 4,
  ncol = length(simList),
  dimnames = list(c(
    "SplTrees_SHAP","SplTrees_Gain", "XGBoost_SHAP", "XGBoost_Gain"
  ),
  scen)
)
prop.misclass.s2 = prop.misclass.s1

for (i in 1:length(simList)) {
  isim = simList[i]
  ranking.s1 = read.csv(paste0(savedir,"results/Ranking_S1_", isim, ".csv"))
  prop.misclass.s1[, i] = colMeans(ranking.s1 == 0 | ranking.s1 > 2)[-1]
  
  ranking.s2 = read.csv(paste0(savedir,"results/Ranking_S2_", isim, ".csv"))
  prop.misclass.s2[, i] = colMeans(ranking.s2 == 0 | ranking.s2 > 2)[-1]
  prop.misclass.s2[prop.misclass.s2==1]=NaN
}

write.csv(as.data.frame(prop.misclass.s1[,-1]), paste0(savedir,"results/Ranking_summary_S1.csv"))
write.csv(as.data.frame(prop.misclass.s2[,-1]), paste0(savedir,"results/Ranking_summary_S2.csv"))


# rankings xgb all --------------------------------------------------------
prop.rank = matrix(
  nrow = length(varNames),
  ncol = length(simList),
  dimnames = list(varNames,
                  scen)
)
prop.rank.spltree = prop.rank.spltree.gain = prop.rank
trueranking = 1:5

for (i in 1:length(simList)) {
  isim = simList[i]
  ranking = read.csv(paste0(savedir,"results/rankingAll_", isim, ".csv"))[,-1]
  prop.rank[, i] = sapply(trueranking,function(i) mean(ranking[,i]==i))
  
  ranking2 = read.csv(paste0(savedir,"results/rankingAll_spltree_", isim, ".csv"))[,-1]
  prop.rank.spltree[, i] = sapply(trueranking,function(i) mean(ranking2[,i]==i))
  
  ranking3 = read.csv(paste0(savedir,"results/rankingAll_spltree_gain_", isim, ".csv"))[,-1]
  prop.rank.spltree.gain[, i] = sapply(trueranking,function(i) mean(ranking3[,i]==i))
}

write.csv(as.data.frame(prop.rank), paste0(savedir,"results/Ranking_summary_all.csv"))
write.csv(as.data.frame(prop.rank.spltree), paste0(savedir,"results/Ranking_summary_spltree.csv"))
write.csv(as.data.frame(prop.rank.spltree.gain), paste0(savedir,"results/Ranking_summary_spltree_gain.csv"))

# rankings xgb average --------------------------------------------------------
prop.rank = matrix(
  nrow = length(varNames),
  ncol = length(simList),
  dimnames = list(varNames,
                  scen)
)
prop.rank.spltree = prop.rank.spltree.gain = prop.rank
trueranking = 1:5

for (i in 1:length(simList)) {
  isim = simList[i]
  ranking = read.csv(paste0(savedir,"results/rankingAll_", isim, ".csv"))[,-1]
  prop.rank[, i] = colMeans(ranking)
  
  ranking2 = read.csv(paste0(savedir,"results/rankingAll_spltree_", isim, ".csv"))[,-1]
  prop.rank.spltree[, i] = colMeans(ranking2)
  
  ranking3 = read.csv(paste0(savedir,"results/rankingAll_spltree_gain_", isim, ".csv"))[,-1]
  prop.rank.spltree.gain[, i] = colMeans(ranking3)
}

write.csv(as.data.frame(prop.rank), paste0(savedir,"results/Ranking_average_summary_all.csv"))
write.csv(as.data.frame(prop.rank.spltree), paste0(savedir,"results/Ranking_average_summary_spltree.csv"))
write.csv(as.data.frame(prop.rank.spltree.gain), paste0(savedir,"results/Ranking_average_summary_spltree_gain.csv"))



