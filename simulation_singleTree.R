
simIdx = "sim7" # sim42 sim43 sim62 sim7
source(paste0("InterpretSGB/EricCode/MySimData_", simIdx,".R"))
library(SubgroupBoost)

expit = function(logodds) 1/(1+exp(-logodds))
set.seed(2021)
# idx = 1:1e3
# model <- SubgroupBoost.RMST(data.simulation[[2]][idx,])
model <- SubgroupBoost.RMST(data.simulation[[1]])

imp <- xgb.importance(model$feature_names, model)
imp # finds s1 s2
xgb.plot.importance(imp, main = "Importance Plot by Gain")
xgb.plot.tree(model = model, trees = 0:5, show_node_id = TRUE)

# predicts on test set ----------------------------------------------------
idx = 1:5e3
test.data <- data.simulation[[2]][idx,model$feature_names]
pred.test <- predict(model, as.matrix(test.data))
test.data$subGroup = ifelse(pred.test>0, "+","-")
test.data$subGroupNumber = ifelse(pred.test>0, 1,0)
test.data$logOdds = pred.test
test.data$trueLabel = ifelse(data.simulation[[4]][idx],"+","-")

acc = with(test.data, mean(1*(subGroup==trueLabel)))
acc

acc.sd = with(test.data, sd(1*(subGroup==trueLabel)))
acc.sd

p.truth <- ggplot(test.data, aes(s1, z1))+
  geom_point(aes(color=trueLabel), alpha=0.5)+
  ggtitle("True Label")
p.truth
ggsave(paste0("Results/", simIdx,"_test_truth.png"), p.truth)

p.fitted <- ggplot(test.data, aes(s1, z1))+
  geom_point(aes(color=subGroup), alpha=0.5)+
  ggtitle("Predicted Label")
p.fitted
ggsave(paste0("Results/", simIdx,"_test_results.png"), p.fitted)

# simplified trees --------------------------------------------------------
source("InterpretSGB/simplified.Tree.R")
spltree = simplified.Tree(model=model, datalist = data.simulation, cutoff = 0.8, plot = T,seed = 1) 
test.data$spltree_pred = ifelse(spltree$stree_pred>0,"+","-")
p.spltree <- ggplot(test.data, aes(s1, z1))+
  geom_point(aes(color=spltree_pred), alpha=0.5)+
  ggtitle("Predicted Label")
p.spltree
ggsave(paste0("Results/", simIdx,"_test_splfTree.png"), p.spltree)
