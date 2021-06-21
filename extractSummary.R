# in the Results folder

simList = paste0("sim", c(101,102,103,104)); scen = paste0("scenario",c(7:10)); cases="MoreFeatures"
# simList = paste0("sim", c(101, 32, 62, 43, 32, 7, 8)); scen =paste0("scenario",c(7, 1:6)); cases="OrginalPaper"

NREP =500
source("/SFS/user/ctc/waguoqi2/InterpretSGB/var.rank.R")

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
  
  ranking.s1 = matrix(rep(0, NREP*4), ncol=4)
  ranking.s2 = matrix(rep(0, NREP*4), ncol=4)

  for (ir in 1:NREP) {
    spltree = readRDS(paste0(isim, "/", ir, ".rds"))
    acc.xgboost[ir] = spltree$accuracy_xgb
    acc.splTree[ir] = spltree$accuracy_simplifiedTree
    sen.xgboost[ir] = spltree$sensitivity_xgb
    sen.splTree[ir] = spltree$sensitivity_simplifiedTree
    spe.xgboost[ir] = spltree$specificity_xgb
    spe.splTree[ir] = spltree$specificity_simplifiedTree
    
    acc.splTree.gain[ir] = spltree$accuracy_simplifiedTree_gain
    sen.splTree.gain[ir] = spltree$sensitivity_simplifiedTree_gain
    spe.splTree.gain[ir] = spltree$specificity_simplifiedTree_gain
    
    ranking.s1[ir,] = do.call("c", var.rank(spltree,"s1"))
    ranking.s2[ir,] = do.call("c", var.rank(spltree,"s2"))
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
  
  colnames(ranking.s1) = c("splTrees_SHAP","SplTrees_GAIN", "XGB_SHAP", "XGB_GAIN")
  colnames(ranking.s2) = c("splTrees_SHAP","SplTrees_GAIN", "XGB_SHAP", "XGB_GAIN")
  
  write.csv(acc.df, file = paste0("results/acc_", isim,".csv"))
  write.csv(sen.df, file = paste0("results/sen_", isim,".csv"))
  write.csv(spe.df, file = paste0("results/spe_", isim,".csv"))
  
  write.csv(as.data.frame(ranking.s1), paste0("results/Ranking_S1_", isim,".csv"))
  write.csv(as.data.frame(ranking.s2), paste0("results/Ranking_S2_", isim,".csv"))
  
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

for (i in 1:length(simList)) {
  isim = simList[i]
  acc.df = read.csv(paste0("results/acc_", isim, ".csv"))
  acc.df.xgboost[,i] = acc.df$acc.xgboost
  acc.df.splTree[,i] = acc.df$acc.splTree
  acc.df.splTree.gain[,i] = acc.df$acc.splTree.gain
  
  mean.acc.xgboost[i] = mean(acc.df$acc.xgboost)
  mean.acc.splTree[i] = mean(acc.df$acc.splTree)
  sd.acc.xgboost[i] = sd(acc.df$acc.xgboost)
  sd.acc.splTree[i] = sd(acc.df$acc.splTree)
  
  sen.df = read.csv(paste0("results/sen_", isim,".csv"))
  spe.df = read.csv(paste0("results/spe_", isim,".csv"))
  sen.df.xgboost[,i] = sen.df$sen.xgboost
  sen.df.splTree[,i] = sen.df$sen.splTree
  sen.df.splTree.gain[,i] = sen.df$sen.splTree.gain
  
  spe.df.xgboost[,i] = spe.df$spe.xgboost
  spe.df.splTree[,i] = spe.df$spe.splTree
  spe.df.splTree.gain[,i] = spe.df$spe.splTree.gain
}

comparison.df = rbind(mean.acc.xgboost, sd.acc.xgboost, mean.acc.splTree,sd.acc.splTree)
colnames(comparison.df)=simList
write.csv(comparison.df, "results/comparison_splTreeVSXgboost.csv")

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
head(acc.df.long)
colnames(acc.df.long) = c("id", "Simulation", "Accuracy", "Method")
write.csv(acc.df.long, paste0("results/Accuracy_comparision_long_", cases, ".csv"))

library(ggplot2)
acc.df.long = read.csv(paste0("results/Accuracy_comparision_long_", cases, ".csv"))
p1 = ggplot(as.data.frame(acc.df.long),aes(x=Simulation, y=Accuracy))+
  geom_boxplot(aes(fill=Method))
ggsave(paste0("results/Accuracy_comparison_XgboostVsSplTree_", cases, ".png"),p1, width = 14)

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
write.csv(sen.df.long, paste0("results/Sensitivity_comparision_long_", cases, ".csv"))

library(ggplot2)
sen.df.long = read.csv(paste0("results/Sensitivity_comparision_long_", cases, ".csv"))
p1 = ggplot(as.data.frame(sen.df.long),aes(x=Simulation, y=Sensitivity))+
  geom_boxplot(aes(fill=Method))
ggsave(paste0("results/Sensitivity_comparison_XgboostVsSplTree_", cases, ".png"),p1, width = 14)


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
write.csv(spe.df.long, paste0("results/Specificity_comparision_long_", cases, ".csv"))

library(ggplot2)
spe.df.long = read.csv(paste0("results/Specificity_comparision_long_", cases, ".csv"))
p1 = ggplot(as.data.frame(spe.df.long),aes(x=Simulation, y=Specificity))+
  geom_boxplot(aes(fill=Method))
ggsave(paste0("results/Specificity_comparison_XgboostVsSplTree_", cases, ".png"),p1, width = 14)


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
  ranking.s1 = read.csv(paste0("results/Ranking_S1_", isim, ".csv"))
  prop.misclass.s1[, i] = colMeans(ranking.s1 == 0 | ranking.s1 > 2)[-1]
  
  ranking.s2 = read.csv(paste0("results/Ranking_S2_", isim, ".csv"))
  prop.misclass.s2[, i] = colMeans(ranking.s2 == 0 | ranking.s2 > 2)[-1]
  prop.misclass.s2[prop.misclass.s2==1]=NaN
}

write.csv(as.data.frame(prop.misclass.s1[,-1]), paste0("results/Ranking_summary_S1.csv"))
write.csv(as.data.frame(prop.misclass.s2[,-1]), paste0("results/Ranking_summary_S2.csv"))



