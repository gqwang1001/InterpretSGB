# devtools::install_github("liupeng2117/SubgroupBoost")


# real data reproducing ---------------------------------------------------
library(SubgroupBoost)
library(dplyr)
library(DiagrammeR)
dat = speff2trial::ACTG175
datReal = dat %>% filter(arms %in% c(1, 3))
info.dat = datReal %>% select(arms, days, cens)
dat.coviates = datReal %>% select(age,
                                  wtkg,
                                  karnof,
                                  cd40,
                                  cd80,
                                  hemo,
                                  homo,
                                  drugs,
                                  race,
                                  gender,
                                  str2,
                                  symptom)
dat.real = data.frame(
  dat.coviates,
  trt01p = (info.dat$arms - 1) / 2,
  #0,1
  evnt = info.dat$cens,
  #0,1
  aval = info.dat$days
)

idx_train = 1:800
dtrain = dat.real[idx_train,]
dtest  = dat.real[-idx_train,]

set.seed(123)
model_RMST <- SubgroupBoost.RMST(dtrain)
importanceReal = xgb.importance(model_RMST$feature_names, model_RMST)
importanceReal
importance_plot = xgb.plot.importance(importanceReal)
export_graph(
  importance_plot,
  file_name = "train_xgboost_importance_plot.png",
  file_type = "png",
  width = 1000,
  height = 1000
)


tree_plot = xgb.plot.tree(model = model_RMST,
                          trees = 0:3,
                          render = F)
export_graph(tree_plot,
             "./Results/train_xgboost_tree_plot.pdf",
             width = 2000,
             height = 3000)


# xgboostExplainer --------------------------------------------------------
# devtools::install_github("AppliedDataSciencePartners/xgboostExplainer")
library(xgboostExplainer)
datTrain = prepareData(dtrain)
datTest = prepareData(dtest)

explainer = buildExplainer(
  model_RMST,
  datTrain$DatMat,
  type = "binary",
  base_score = 0.5,
  trees_idx = NULL
)
pred.breakdown = explainPredictions(model_RMST, explainer, datTest$DatMat)

idx_to_get = 100
showWaterfall(model_RMST,
              explainer,
              datTest$DatMat,
              data.matrix(dtest[,1:12]),
              idx_to_get,
              type = "binary")

cr <- colorRamp(c("blue", "red"))

plot(
  dtest$cd80,
  pred.breakdown$cd80,
  cex = 0.4,
  pch = 16,
  xlab = "CD8",
  ylab = "CD8 impact on log-odds"
)
plot(
  dtest$cd40,
  pred.breakdown$cd40,
  cex = 0.4,
  pch = 16,
  col = rgb(cr(dtest$str2), max=255),
  xlab = "CD4",
  ylab = "CD4 impact on log-odds"
)
legend(
  "topright",
  title = "Antiretroviral History",
  legend = c("naive", "experienced"),
  col = c("blue", "red"),
  cex = 1,
  pch = 16
)

plot(
  dtest$str2,
  pred.breakdown$str2,
  cex = 0.4,
  pch = 16,
  xlab = "ST2",
  ylab = "ST2 impact on log-odds"
)

plot(
  dtest$age,
  pred.breakdown$age,
  cex = 0.4,
  pch = 16,
  xlab = "AGE",
  ylab = "AGE impact on log-odds"
)
plot(
  dtest$wtkg,
  pred.breakdown$wtkg,  
  col = rgb(cr(dtest$str2), max=255),
  cex = 0.4,
  pch = 16,
  xlab = "Weight",
  ylab = "Weight impact on log-odds"
)

legend(
  "topright",
  title = "Antiretroviral History",
  legend = c("naive", "experienced"),
  col = c("blue", "red"),
  cex = 1,
  pch = 16
)

plot(
  dtest$age,
  pred.breakdown$age,
  col = rgb(cr(dtest$str2), max=255),
  cex = 0.4,
  pch = 16,
  xlab = "AGE",
  ylab = "AGE impact on log-odds"
)
legend(
  "topright",
  title = "Antiretroviral History",
  legend = c("naive", "experienced"),
  col = c("blue", "red"),
  cex = 1,
  pch = 16
)

dtest$wtstrata = ifelse(dtest$wtkg<60,1,ifelse(dtest$wtkg>100,16,2))
par(mar = c(5, 5, 4, 9))
plot(
  dtest$cd40,
  pred.breakdown$cd40,
  cex = 1,
  pch = dtest$wtstrata,
  col = rgb(cr(dtest$str2), max=255),
  xlab = "CD4",
  ylab = "CD4 impact on log-odds"
)

legend(
  "topright",inset = c(-0.4, 0),
  title = "Antiretroviral History",
  legend = c("naive", "experienced"),
  col = c("blue", "red"),
  cex = 1,
  xpd = T,
  pch = 16
)
legend(
  "right",inset = c(-0.4, 0),
  title = "Weight(Kg)",
  legend = c("weight<60", "60<=weight<100", "weight>=100"),
  col = 1,
  cex = 1,
  xpd = T,
  pch = c(1,2,16)
)

# ggplot results ----------------------------------------------------------
library(ggplot2)
library(xgboostExplainer)
logOdds = rowSums(pred.breakdown)
dat.plot = data.frame(CD4 = dtest$cd40,
                      Antiretroviral_history = factor(as.character(dtest$str2), labels = c("Naive", "Experienced")),
                      weight = dtest$wtkg,
                      logOdds = logOdds,
                      subgroup = as.factor(ifelse(logOdds>0, "postive treatment effect", "no effect")))
g1 <- ggplot(dat.plot, aes(x=CD4, y=weight, color =subgroup, shape = Antiretroviral_history))+
  geom_point(size=2)+theme(legend.position = 'bottom')
g1

ggsave("./Results/logOdds_test.png",g1)

g7 <- ggplot(dat.plot, aes(x=CD4, y = logOdds, color=Antiretroviral_history))+
  geom_point(size=2, alpha=0.5)+theme(legend.position = 'bottom')
g7
ggsave("./Results/logOdds_test_first2.png",g7)

explainer = buildExplainer(
  model_RMST,
  datTrain$DatMat,
  type = "binary",
  base_score = 0.5,
  trees_idx = NULL
)
pred.breakdown.train = explainPredictions(model_RMST, explainer, datTrain$DatMat)

logOdds.train = rowSums(pred.breakdown.train)
dat.plot.train = data.frame(CD4 = dtrain$cd40,
                            Antiretroviral_history = factor(as.character(dtrain$str2), labels = c("Naive", "Experienced")),
                      weight = dtrain$wtkg,
                      logOdds = logOdds.train,
                      subgroup = as.factor(ifelse(logOdds.train>0, "postive treatment effect", "no effect")))
g2 <- ggplot(dat.plot.train, aes(x=CD4, y=weight, color =subgroup, shape = Antiretroviral_history))+
  geom_point(size=2)+theme(legend.position = 'bottom')
g2
ggsave("./Results/logOdds_train.png",g2)

g8 <- ggplot(dat.plot.train, aes(x=CD4, y = logOdds, color=Antiretroviral_history))+
  geom_point(size=2, alpha=0.5)+theme(legend.position = 'bottom')
g8
ggsave("./Results/logOdds_train_first2.png",g8)

# remove high leverage data -----------------------------------------------
dat.robust = dat.real[dat.real$cd40<1000,]

idx_train = 1:800
dtrain = dat.robust[idx_train,]
dtest  = dat.robust[-idx_train,]

set.seed(123)
model_RMST_robust <- SubgroupBoost.RMST(dtrain)
importanceReal_robust = xgb.importance(model_RMST_robust$feature_names, model_RMST_robust)
importanceReal_robust
importance_plot = xgb.plot.importance(importanceReal_robust)

library(xgboostExplainer)
datTrain = prepareData(dtrain)
datTest = prepareData(dtest)

explainer = buildExplainer(
  model_RMST_robust,
  datTrain$DatMat,
  type = "binary",
  base_score = 0.5,
  trees_idx = NULL
)
pred.breakdown.robust = explainPredictions(model_RMST_robust, explainer, datTest$DatMat)
logOdds = rowSums(pred.breakdown.robust)
dat.plot.train.robust = data.frame(CD4 = dtest$cd40,
                            Antiretroviral_history = factor(as.character(dtest$str2), labels = c("Naive", "Experienced")),
                            weight = dtest$wtkg,
                            logOdds = logOdds,
                            subgroup = as.factor(ifelse(logOdds>0, "postive treatment effect", "no effect")))
g3 <- ggplot(dat.plot.train.robust, aes(x=CD4, y=weight, color =subgroup, shape = Antiretroviral_history))+
  geom_point(size=2)+theme(legend.position = 'bottom')
g3
ggsave("./Results/logOdds_test_robust.png",g3)
g6 <- ggplot(dat.plot.train.robust, aes(x=CD4, y = logOdds, color=Antiretroviral_history))+
  geom_point(size=2, alpha=0.5)+theme(legend.position = 'bottom')
g6
ggsave("./Results/logOdds_test_robust_first2.png",g6)

pred.breakdown.robust = explainPredictions(model_RMST_robust, explainer, datTrain$DatMat)
logOdds = rowSums(pred.breakdown.robust)
dat.plot.train.robust = data.frame(CD4 = dtrain$cd40,
                                   Antiretroviral_history = factor(as.character(dtrain$str2), labels = c("Naive", "Experienced")),
                                   weight = dtrain$wtkg,
                                   logOdds = logOdds,
                                   subgroup = as.factor(ifelse(logOdds>0, "postive treatment effect", "no effect")))
g4 <- ggplot(dat.plot.train.robust, aes(x=CD4, y=weight, color =subgroup, shape = Antiretroviral_history))+
  geom_point(size=2)+theme(legend.position = 'bottom')
g4
ggsave("./Results/logOdds_train_robust.png",g4)

g5 <- ggplot(dat.plot.train.robust, aes(x=CD4, y = logOdds, color=Antiretroviral_history))+
  geom_point(size=2, alpha=0.5)+theme(legend.position = 'bottom')
g5
ggsave("./Results/logOdds_train_robust_first2.png",g5)



# traditional approach ----------------------------------------------------

preds = predict(model_RMST_robust, datTrain$DatMat)
preds
