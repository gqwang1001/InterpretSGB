source("InterpretSGB/EricCode/MySimData_sim62.R")
library(SubgroupBoost)

expit = function(logodds) 1/(1+exp(-logodds))
set.seed(123)
# m62 <- SubgroupBoost.RMST_local(sim62[[2]][1:1000,])
m62 <- SubgroupBoost.RMST_local(sim62[[1]])

m62$params

imp62 <- xgb.importance(m62$feature_names, m62)
imp62 # finds s1 s2
xgb.plot.importance(imp62, main = "Importance Plot by Gain")
xgb.plot.tree(model = m62, trees = 0:5, show_node_id = TRUE)

# predicts on test set ----------------------------------------------------
idx = -c(1:1e3)
test62 <- sim62[[2]][idx,m62$feature_names]
pred62.test <- predict(m62, as.matrix(test62))
test62$subGroup = ifelse(pred62.test>0, "+","-")
test62$subGroupNumber = ifelse(pred62.test>0, 1,0)
test62$logOdds = pred62.test
test62$trueLabel = ifelse(sim62[[4]][idx],"+","-")

acc = with(test62, mean(1*(subGroup==trueLabel)))
acc
acc.sd = with(test62, sd(1*(subGroup==trueLabel)))
acc.sd/sqrt(500)

p62.1.truth <- ggplot(test62, aes(s1, s2))+
  geom_point(aes(color=trueLabel), alpha=0.5)+
  ggtitle("True Label")
p62.1.truth
ggsave("Results/sim62_test_truth.png", p62.1.truth)

p62.1 <- ggplot(test62, aes(s1, s2))+
  geom_point(aes(color=subGroup), alpha=0.5)+
  geom_abline(slope = 1)+
  ggtitle("Predicted Label")
p62.1
ggsave("Results/sim62_test_results.png", p62.1)

#generate data grid to show the decision boundary
xnew62 <- with(test62, expand.grid(s1 = seq(range(s1)[1], range(s1)[2], by=0.1),
                                   s2 = seq(range(s2)[1], range(s2)[2], by=0.1)))
xmeans <- matrix(rep(t(colMeans(sim62[[2]][,m62$feature_names])), nrow(xnew62)), nrow = nrow(xnew62), byrow = T)
colnames(xmeans) = m62$feature_names
dat62.new = as.data.frame(xmeans)
dat62.new[,c("s1", "s2")] = xnew62

pred62.grid <- predict(m62, as.matrix(dat62.new))
dat62.grid = data.frame(dat62.new,
                        logOdds = pred62.grid,
                        prob = expit(pred62.grid),
                        subGroup = ifelse(pred62.grid>0, "+","-"),
                        subGroups = ifelse(pred62.grid>0, 1,0))

p62.4 <- ggplot(dat62.grid, aes(s1, s2)) +
  geom_raster(aes(fill=prob)) +
  scale_fill_distiller(palette = "Spectral", direction = 1)+
  geom_point(data = test62, aes(s1, s2, color = subGroup),alpha=0.5)+
  geom_contour(aes(z=subGroups),color='black',breaks = .5, size=2)+
  geom_abline(slope = 1)
p62.4

ggsave("Results/sim62_plot.png", p62.4)

# SHAP value --------------------------------------------------------------
source("InterpretSGB/shap.R")
library(SHAPforxgboost)
library(dplyr)
# data.shap = as.matrix(sim62[[2]][idx, m62$feature_names])
data.shap = as.matrix(sim62[[1]][, m62$feature_names])
shap_values = shap.values(m62, data.shap)
shap_int <- shap.prep.interaction(xgb_mod = m62, X_train = data.shap)
shap_results = shap.score.rank(
  xgb_model = m62,
  X_train = data.shap,
  shap_approx = F
)
var_importance(shap_results, top_n = 5)+ggtitle("Importance by SHAP")
ggsave("Results/sim62_SHAP_importancePlot.png")

shap_long = shap.prep(shap = shap_results,
                      X_train = data.shap,
                      top_n = 5)
shap.summaryPlot = plot.shap.summary(data_long = shap_long)+ggtitle("SHAP summary")
ggsave("Results/sim62_SHAP_summaryplot.png", shap.summaryPlot)

g3 <- shap.plot.dependence(data_long = shap_long,
                           data_int = shap_int,
                           x= "s2", y = "s2", 
                           color_feature = "s1")
g3

g4 <- shap.plot.dependence.local(data_long = shap_long,
                           data_int = shap_int,
                           x= "s1", y = "s2", 
                           color_feature = "s2")
g4

g5 <- shap.plot.dependence.local(data_long = shap_long,
                                 data_int = shap_int,
                                 x= "s2", y = "s1", 
                                 color_feature = "s1") + 
  ggtitle("SHAP Interaction")
ggsave("Results/shapInteraction.png", g5)

plot_data <- shap.prep.stack.data(shap_contrib = shap_values$shap_score, top_n = 10, n_groups = 4, cluster_method = "ward.D")

shap.plot.force_plot(plot_data, zoom_in = F)
# ggplot(test62, aes(x = 1:4000,y=1, color=trueLabel))+geom_point()

p62.1.shap.fp_by_group=shap.plot.force_plot_bygroup(plot_data)+ggtitle("forceplot by group")
ggsave("Results/sim62_test_results_shap_forceplot_by_group.png", p62.1.shap.fp_by_group)

sorted_plot_data <- plot_data[order(plot_data$ID),]
# test62$subGroup_SHAP = ifelse(sorted_plot_data$group==1,"+","-")
test62$subGroup_SHAP = as.factor(sorted_plot_data$group)

acc = with(test62, mean(1*(subGroup_SHAP==trueLabel)))
acc
acc.sd = with(test62, sd(1*(subGroup_SHAP==trueLabel)))
acc.sd/sqrt(500)

p62.1.shap <- ggplot(test62, aes(s1, s2))+
  geom_point(aes(color=subGroup_SHAP), alpha=0.5)+
  geom_abline(slope = 1)+
  ggtitle("SHAP clustering")
p62.1.shap

ggsave("Results/sim62_test_results_shap_clustering.png", p62.1.shap)

# use xgboost package
png("Results/ShapplotPKG.png")
shapPlotPKG = xgb.plot.shap(data = as.matrix(sim62[[2]][, m62$feature_names]), model = m62, top_n = 20, plot = F)
dev.off()

shapPlotPKG = xgb.plot.shap(data = as.matrix(sim62[[2]][, m62$feature_names]), model = m62, top_n = 25, plot = F)
shapresults = apply(shapPlotPKG$shap_contrib, 2, function(c) mean(abs(c)))
sort(shapresults, decreasing = T)
cumsum(sort(shapresults, decreasing = T)/sum(shapresults))


