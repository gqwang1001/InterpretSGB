simplified.Tree = function(model,
                           datalist,
                           cutoff = 0.8,
                           plot.name = NaN,
                           seed = 1) {
  require(caret)
  require(rpart.plot)
  require(rattle)
  require(xgboost)
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
  
  # fit decision tree -------------------------------------------------------
  dataInPred <- data.frame(shapPlotPKG$data[, imp.val],
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
  # summary(dtree_fit)
  if (!is.na(plot.name)) {
    png(filename = plot.name)
    fancyRpartPlot(dtree_fit$finalModel)
    dev.off()
  }
  
  pred.simple = predict(dtree_fit, datalist[[2]])
  pred.xgb <- predict(model, as.matrix(datalist[[2]][,model$feature_names]))
  
  rmse =  RMSE(pred.simple, pred.xgb)
  acc.simplifiedTree = mean(1 * ((pred.simple > 0) == datalist[[4]]))
  acc.xgb = mean(1 * ((pred.xgb > 0) == datalist[[4]]))
  acc.comparison=mean(1*(((pred.simple > 0) == (pred.xgb > 0))))
  
  
  return(
    list(
      imporant_variables = imp.val,
      simplified_tree = dtree_fit,
      stree_pred = pred.simple,
      RMSE = rmse,    
      accuracy_xgb = acc.xgb,
      accuracy_simplifiedTree = acc.simplifiedTree,
      acc_comparison = acc.pred
    )
  )
}
