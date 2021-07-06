library(SubgroupBoost)
expit = function(logodds)
  1 / (1 + exp(-logodds))
source("InterpretSGB/simplified.Tree.R")
simList = paste0("sim", c(108,105, 106,107,101,102,103,104, 32, 42, 43, 62, 7, 8))
seeds = 2
# Train model -------------------------------------------------------------

for (i in 1:6) {
  simIdx = simList[i]
  set.seed(seeds)
  source(paste0("InterpretSGB/EricCode/MySimData_", simIdx, ".R"))
  # idx = 1:1e3
  # model <- SubgroupBoost.RMST(data.simulation[[2]][idx,])
  model <- SubgroupBoost.RMST_local(data.simulation[[1]])
  saveRDS(model, file = paste0("Results/modelFitting_", simIdx, ".rds"))
}

# View results ------------------------------------------------------------
for (i in 1:6) {
  simIdx = simList[i]
  print(simIdx)
  set.seed(seeds)
  source(paste0("InterpretSGB/EricCode/MySimData_", simIdx, ".R"))
  model = readRDS(paste0("Results/modelFitting_", simIdx, ".rds"))
  imp <- xgb.importance(model$feature_names, model)
  imp
  # imp # finds s1 s2
  # xgb.plot.importance(imp, main = "Importance Plot by Gain")
  # xgb.plot.tree(model = model,
  #               trees = 0:5,
  #               show_node_id = TRUE)
  spltree = simplified.Tree(
    model = model,
    datalist = data.simulation,
    cutoff = 0.95,
    # plot.name = paste0("Results/", simIdx,"_simplified_tree.png"),
    seed = seeds,
    top3 = F
  )
  
  # saveRDS(spltree, file = paste0("Results/", simIdx,"/", seeds,".rds"))
  # saveRDS(spltree, file = paste0("Results/", seeds,".rds"))
  spltree$evalLoss.lm
  spltree$evalLoss.spltree
  spltree$evalLoss.lm.train
  spltree$evalLoss.spltree.train
  
  
  spltree$imporant_variables
  spltree$imporant_variables_sptree
  spltree$imporant_variables_sptree_gain
  
  spltree$RMSE.lm
  spltree$RMSE.spltree
  
  spltree$RMSE.lm.train
  spltree$RMSE.spltree.train

  spltree$summaryPred.lm$accuracy
  spltree$summaryPred.splTree.shap$accuracy
  spltree$summaryPred.splTree.gain$accuracy
  spltree$summaryPred.xgb$accuracy
  spltree$acc.comp.lm
  spltree$acc.comp.splTreeVsxgb

  spltree$imporant_variables_xgb_all
  
    finaltree = spltree$simplified_tree$finalModel
  split.var = unique(finaltree$frame$var[finaltree$frame$var!="<leaf>"])
  
  rattle::fancyRpartPlot(finaltree)
  rattle::fancyRpartPlot(spltree$simplified_tree_gain$finalModel)
  
  sum.finaltree = rpart::printcp(finaltree)

    summary(finaltree)
  finaltree$splits[1,]
  varImp(spltree$simplified_tree$finalModel)
  
  vi.shap = vip::vi_shap(finaltree, pred_wrapper = predict) %>% filter(Importance!=0) %>% arrange(desc(Importance))
  
 ranks = var.rank(spltree,"s2")
  do.call("c",var.rank(spltree,"s2"))
  
  impvars = spltree$imporant_variables
  impvars
  
  spltree$imporant_variables_sptree
  spltree$imporant_variables_xgb_all
  
  # predicts on test set ----------------------------------------------------
  idx = 1:5e3
  test.data <- data.simulation[[2]][idx, model$feature_names]
  pred.test <- predict(model, as.matrix(test.data))
  test.data$subGroup = ifelse(pred.test > 0, "+", "-")
  test.data$subGroupNumber = ifelse(pred.test > 0, 1, 0)
  test.data$logOdds = pred.test
  test.data$trueLabel = ifelse(data.simulation[[4]][idx], "+", "-")
  
  acc = with(test.data, mean(1 * (subGroup == trueLabel)))
  acc
  
  acc.sd = with(test.data, sd(1 * (subGroup == trueLabel)))
  acc.sd/sqrt(5000)
  
  p.truth <- ggplot(test.data, aes_string(impvars[1], impvars[2])) +
    geom_point(aes(color = trueLabel), alpha = 0.5) +
    ggtitle("True Label")
  # p.truthc
  # ggsave(paste0("Results/", simIdx,"_test_truth.png"), p.truth)
  
  p.fitted <- ggplot(test.data, aes_string(impvars[1], impvars[2])) +
    geom_point(aes(color = subGroup), alpha = 0.5) +
    ggtitle("Predicted Label by XGboost")+
    geom_label(aes(x = mean(test.data[,impvars[1]]), y = min(test.data[,impvars[2]])), label = paste("Accuracy =", acc))
  
  # p.fitted.comparison <-
  #   gridExtra::grid.arrange(p.truth, p.fitted, nrow = 1)
  # ggsave(
  #   paste0("Results/", simIdx, "_test_results.png"),
  #   p.fitted.comparison,
  #   width = 15,
  #   height = 6.5
  # )
  
  # simplified trees --------------------------------------------------------

  # ggsave(paste0("Results/", simIdx, "_test_splfTree_plot.png"),
  #        spltree$tree_figure)
  test.data$spltree_pred = ifelse(spltree$stree_pred > 0, "+", "-")
  
  p.spltree <- ggplot(test.data, aes_string(impvars[1], impvars[2])) +
    geom_point(aes(color = spltree_pred), alpha = 0.5) +
    ggtitle("Predicted Label by Simplifed Tree with Top 3 Features")+
    geom_label(aes(x = mean(test.data[,impvars[1]]), y = min(test.data[,impvars[2]])), label = paste("Accuracy =", spltree$accuracy_simplifiedTree))
  # p.spltree
  
  p.fitted.comparison.spltree <-
    gridExtra::grid.arrange(p.truth, p.fitted, p.spltree, nrow = 1)
  ggsave(
    paste0("Results/", simIdx, "_test_splfTree_top3.png"),
    p.fitted.comparison.spltree,
    width = 21,
    height = 6.5
  )
}
