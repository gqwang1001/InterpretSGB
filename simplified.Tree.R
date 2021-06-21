simplified.Tree = function(model,
                           datalist,
                           cutoff = 0.8,
                           plot.name = NaN,
                           seed = 1,
                           top3 = T) {
  require(caret)
  # require(rpart.plot)
  # require(rattle)
  require(xgboost)
  library(xgboost)
  library(caret)
  
  options(warn = -1)
  shapPlotPKG <- xgb.plot.shap(
    data = as.matrix(datalist[[1]][, model$feature_names]),
    model = model,
    top_n = 25,
    plot = F
  )
  shapresults <- colMeans(abs(shapPlotPKG$shap_contrib))
  shap_cumsum <- cumsum(sort(shapresults, decreasing = T) / sum(shapresults))
  valNames <- names(shap_cumsum)
  imp.val <- valNames[which(shap_cumsum < cutoff)]
  
  imp.gain <- xgb.importance(model$feature_names, model)
  imp.gain.top = imp.gain$Feature[1:3]
  
  if (top3 & length(imp.val)>=3) imp.val=imp.val[1:3]
  if (top3 &!is.na(plot.name)) plot.name = paste0("Top3_",plot.name)
  if (length(imp.val)==1) imp.val = valNames[1:2]
  
  # fit decision tree -------------------------------------------------------
  dataInPred <- data.frame(shapPlotPKG$data[, imp.val],
                          logOdds = rowSums(shapPlotPKG$shap_contrib))
  dataInPred.gain <- data.frame(shapPlotPKG$data[, imp.gain.top],
                           logOdds = rowSums(shapPlotPKG$shap_contrib))
  trctrl <-
    trainControl(method = "repeatedcv",
                 number = 10,
                 repeats = 3)
  set.seed(seed)
  dtree_fit <- train(
    logOdds ~ .,
    data = dataInPred,
    method = "rpart",
    trControl = trctrl,
    tuneLength = 10
  )
  
  dtree_fit_gain <- train(
    logOdds ~ .,
    data = dataInPred.gain,
    method = "rpart",
    trControl = trctrl,
    tuneLength = 10
  )
  
  
  # summary(dtree_fit)
  if (!is.na(plot.name)) {
    library(rattle)
    png(filename = plot.name)
    rattle::fancyRpartPlot(dtree_fit$finalModel)
    dev.off()
  }
  
  pred.simple = predict(dtree_fit, datalist[[2]])
  pred.simple.gain = predict(dtree_fit_gain, datalist[[2]])
  
  pred.xgb <- predict(model, as.matrix(datalist[[2]][,model$feature_names]))
  
  rmse =  RMSE(pred.simple, pred.xgb)
  acc.simplifiedTree = mean(1 * ((pred.simple > 0) == datalist[[4]]))
  acc.simplifiedTree.gain = mean(1 * ((pred.simple.gain > 0) == datalist[[4]]))
  
  acc.xgb = mean(1 * ((pred.xgb > 0) == datalist[[4]]))
  acc.comparison=mean(1*(((pred.simple > 0) == (pred.xgb > 0))))
  
  
  sen.spl <- sum(pred.simple>=0 & datalist[[4]]==1)/sum(datalist[[4]]==1)
  spe.spl <- sum(pred.simple<0 & datalist[[4]]==0)/sum(datalist[[4]]==0)
  sen.xgb <- sum(pred.xgb>=0 & datalist[[4]]==1)/sum(datalist[[4]]==1)
  spe.xgb <- sum(pred.xgb<0 & datalist[[4]]==0)/sum(datalist[[4]]==0)
  sen.spl.gain <- sum(pred.simple.gain>=0 & datalist[[4]]==1)/sum(datalist[[4]]==1)
  spe.spl.gain <- sum(pred.simple.gain<0 & datalist[[4]]==0)/sum(datalist[[4]]==0)
  
  return(
    list(
      xgbModel = model,
      imporant_variables = imp.val,
      imporant_variables_xgb_all = valNames,
      imporant_variables_sptree = dtree_fit$finalModel$variable.importance,
      imporant_variables_sptree_gain = dtree_fit_gain$finalModel$variable.importance,
      
      simplified_tree = dtree_fit,
      simplified_tree_gain = dtree_fit_gain,
      
      stree_pred = pred.simple,
      RMSE = rmse,    
      accuracy_xgb = acc.xgb,
      accuracy_simplifiedTree = acc.simplifiedTree,
      accuracy_simplifiedTree_gain = acc.simplifiedTree.gain,
      
      sensitivity_xgb = sen.xgb,
      sensitivity_simplifiedTree = sen.spl,   
      sensitivity_simplifiedTree_gain = sen.spl.gain,     
      
      specificity_xgb = spe.xgb,
      specificity_simplifiedTree = spe.spl,
      specificity_simplifiedTree_gain = spe.spl.gain,
      
      acc_comparison = acc.comparison
    )
  )
}
