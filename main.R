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
saveRDS(model_RMST, "realDataResults.rds")
importanceReal = xgb.importance(model_RMST$feature_names, model_RMST)
importanceReal
importance_plot = xgb.plot.importance(importanceReal)
export_graph(
  importance_plot,
  file_name = "Results/train_xgboost_importance_plot.png",
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

cd4range = range(dat.coviates$cd40)
wtrange = range(dat.coviates$wtkg)

RealNew = expand.grid(cd40 = seq(cd4range[1], cd4range[2], by=10),
                      wtkg = seq(wtrange[1], wtrange[2], by=1),
                      str2 = c(0,1))
ndata = nrow(RealNew)
RealNew$age = mean(dat.coviates$age)
RealNew$karnof = mean(dat.coviates$karnof)
RealNew$gender = 1
RealNew$cd80 = mean(dat.coviates$cd80)
RealNew$drugs = 1
RealNew$race = 1
RealNew$symptom=1
RealNew$hemo = 1
RealNew$homo = 1

RealNew = RealNew[, model_RMST_robust$feature_names]
preds.realnew = predict(model_RMST_robust, as.matrix(RealNew))

dat.plot = data.frame(RealNew,
                      logOdds=preds.realnew,
                      subGroup = ifelse(preds.realnew>0, "+","-"),
                      subGroups = ifelse(preds.realnew>0, 1,0),
                      Antiretroviral_history = factor(as.character(RealNew$str2), 
                                                      labels = c("Naive", "Experienced")))

p1.realdata <- ggplot(dat.plot, aes(cd40, wtkg)) +
  geom_raster(aes(fill=logOdds)) +
  scale_fill_distiller(palette = "Spectral", direction = -1)+
  geom_contour(aes(z=subGroups),breaks = 0.1, size=2, color = 'black')+
  geom_point(data = dat.plot.train.robust, aes(CD4, weight, color = subgroup),alpha=0.5)+
  facet_wrap(~Antiretroviral_history)+
  xlab("CD4")+
  ylab("Weight(kg)")
p1.realdata
ggsave("Results/realdata_interpret_plots.png",p1.realdata)

p2.realdata <-
  ggplot(dat.plot, aes(cd40, wtkg))+
  geom_contour(aes(z=subGroups),breaks = .5, size=2, color = 'black')+
  geom_point(data = dat.plot.train.robust, aes(CD4, weight, color=subgroup),alpha=0.5)+
  facet_wrap(~Antiretroviral_history)+
  xlab("CD4")+
  ylab("Weight(kg)")
p2.realdata


# SHAP approach -----------------------------------------------------------

source("InterpretSGB/shap.R")
library(SHAPforxgboost)
library(dplyr)
# data.shap = as.matrix(sim62[[2]][idx, m62$feature_names])
data.shap = as.matrix(dtrain[,model_RMST$feature_names])
shap_values = shap.values(model_RMST, data.shap)
shap_int <- shap.prep.interaction(xgb_mod = model_RMST, X_train = data.shap)
shap_results = shap.score.rank(
  xgb_model = model_RMST,
  X_train = data.shap,
  shap_approx = F
)
var_importance(shap_results, top_n = 10)+ggtitle("Importance by SHAP")
ggsave("Results/aids_SHAP_importancePlot.png")

shap_long = shap.prep(shap = shap_results,
                      X_train = data.shap,
                      top_n = 5)
shap.summaryPlot = plot.shap.summary(data_long = shap_long)+ggtitle("SHAP summary")
ggsave("Results/aids_SHAP_summaryplot.png", shap.summaryPlot)

g2 <- shap.plot.dependence(data_long = shap_long,
                           x= "str2", y = "cd40", 
                           color_feature = "cd40")
g2
g3 <- shap.plot.dependence(data_long = shap_long,
                           x= "cd40", y = "str2", 
                           color_feature = "str2")+
  ggtitle("SHAP values of str2 vs CD4")
g3
ggsave("Results/aids_shap_plots.png", g3)

g4 <- shap.plot.dependence.local(data_long = shap_long,
                                 data_int = shap_int,
                                 x= "s1", y = "s2", 
                                 color_feature = "s2")
g4

g5 <- shap.plot.dependence.local(data_long = shap_long,
                                 data_int = shap_int,
                                 x= "cd40", y = "str2", 
                                 color_feature = "str2") + 
  ggtitle("SHAP Interaction")
ggsave("Results/aids_shapInteraction.png", g5)
g6 <- shap.plot.dependence.local(data_long = shap_long,
                                 data_int = shap_int,
                                 x= "cd80", y = "str2", 
                                 color_feature = "str2") + 
  ggtitle("SHAP Interaction")
ggsave("Results/aids_shapInteraction_cd80_str2.png", g6)

plot_data <- shap.prep.stack.data(shap_contrib = shap_values$shap_score, top_n = 10, n_groups = 4, cluster_method = "ward.D")

p62.1.shap.fp=shap.plot.force_plot(plot_data, zoom_in = F)
# ggplot(test62, aes(x = 1:4000,y=1, color=trueLabel))+geom_point()
ggsave("Results/aids_test_results_shap_forceplot.png", p62.1.shap.fp)

p62.1.shap.fp_by_group=shap.plot.force_plot_bygroup(plot_data)+ggtitle("forceplot by group")
ggsave("Results/aids_test_results_shap_forceplot_by_group.png", p62.1.shap.fp_by_group)

sorted_plot_data <- plot_data[order(plot_data$ID),]
# test62$subGroup_SHAP = ifelse(sorted_plot_data$group==1,"+","-")
# test62$subGroup_SHAP = as.factor(sorted_plot_data$group)
# 
# acc = with(test62, mean(1*(subGroup_SHAP==trueLabel)))
# acc
# acc.sd = with(test62, sd(1*(subGroup_SHAP==trueLabel)))
# acc.sd/sqrt(500)
# 
# p62.1.shap <- ggplot(test62, aes(s1, s2))+
#   geom_point(aes(color=subGroup_SHAP), alpha=0.5)+
#   geom_abline(slope = 1)+
#   ggtitle("SHAP clustering")
# p62.1.shap
# 
# ggsave("Results/sim62_test_results_shap_clustering.png", p62.1.shap)
# 
# # use xgboost package
# png("Results/ShapplotPKG.png")
# shapPlotPKG = xgb.plot.shap(data = as.matrix(sim62[[2]][, m62$feature_names]), model = m62)
# dev.off()
# 
# 
# 






