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
  require(vip)
  library(xgboost)
  library(caret)
  library(dplyr)
  
  options(warn = -1)
  shapPlotPKG <- xgb.plot.shap(
    data = as.matrix(datalist[[1]][, model$feature_names]),
    model = model,
    top_n = 25,
    plot = F
  )
  shapresults <- colMeans(abs(shapPlotPKG$shap_contrib))
  shap_cumsum <-
    cumsum(sort(shapresults, decreasing = T) / sum(shapresults))
  valNames <- names(shap_cumsum)
  imp.val <- valNames[which(shap_cumsum < cutoff)]
  
  imp.gain <- xgb.importance(model$feature_names, model)
  imp.gain.top = imp.gain$Feature[1:min(10,length(imp.gain$Feature))]#[1:5]
  
  if (top3 & length(imp.val) >= 3)
    imp.val = imp.val[1:3]
  if (top3 & !is.na(plot.name))
    plot.name = paste0("Top3_", plot.name)
  if (length(imp.val) == 1)
    imp.val = valNames[1:2]
  
  dataInPred <- data.frame(shapPlotPKG$data[, imp.val],
                           logOdds = rowSums(shapPlotPKG$shap_contrib))
  # fit decision tree -------------------------------------------------------
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
    method = "rpart2",
    trControl = trctrl,
    tuneLength = 3
  )
  
  dtree_fit_gain <- train(
    logOdds ~ .,
    data = dataInPred.gain,
    method = "rpart2",
    trControl = trctrl,
    tuneLength = 3
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
  
  pred.xgb <-
    predict(model, as.matrix(datalist[[2]][, model$feature_names]))
  
  rmse.splTreeVsxgb =  RMSE(pred.simple, pred.xgb)
  summaryPred.splTree.shap = summaryPred(pred.simple, datalist)
  summaryPred.splTree.gain = summaryPred(pred.simple.gain, datalist)
  summaryPred.xgb = summaryPred(pred.xgb, datalist)
  
  acc.comp.splTreeVsxgb = 2 * sum(1 * (((pred.simple > 0) == (pred.xgb > 0)))) /
    (length(pred.simple) + length(pred.xgb))
  
  # fit linear regression model ---------------------------------------------
  fit.lm.full = lm(data = dataInPred, logOdds ~ .)
  fit.lm = step(fit.lm.full, trace = 0)
  
  pred.lm = predict(fit.lm, datalist[[2]])
  acc.comp.lm = 2 * sum(1 * (((pred.lm > 0) == (pred.xgb > 0)))) /
    (length(pred.lm) + length(pred.xgb))
  rmse.lmVsxgb =  RMSE(pred.lm, pred.xgb)
  summaryPred.lm = summaryPred(pred.lm, datalist)
  
  
  return(
    list(
      xgbModel = model,
      imporant_variables = imp.val,
      imporant_variables_xgb_all = valNames,
      imporant_variables_sptree = vip::vi_shap(dtree_fit$finalModel, pred_wrapper = predict) %>% filter(Importance!=0) %>% arrange(desc(Importance)),
      imporant_variables_sptree_gain = vip::vi_shap(dtree_fit_gain$finalModel, pred_wrapper = predict) %>% filter(Importance!=0) %>% arrange(desc(Importance)),
      
      simplified_tree = dtree_fit,
      simplified_tree_gain = dtree_fit_gain,
      stree_pred = pred.simple,
      
      summaryPred.xgb = summaryPred.xgb,
      summaryPred.splTree.shap = summaryPred.splTree.shap,
      summaryPred.splTree.gain = summaryPred.splTree.gain,
      
      acc.comp.splTreeVsxgb = acc.comp.splTreeVsxgb,
      RMSE.spltree = rmse.splTreeVsxgb,
      
      lm_pred = pred.lm,
      acc.comp.lm = acc.comp.lm,
      RMSE.lm = rmse.lmVsxgb,
      summaryPred.lm = summaryPred.lm,
      lmModel = fit.lm
    )
  )
}


summaryPred = function(pred, datalist) {
  return(list(
    accuracy = mean(1 * ((pred > 0) == datalist[[4]])),
    sensitivity = sum(pred >= 0 & datalist[[4]] == 1) / sum(datalist[[4]] == 1),
    specificity = sum(pred < 0 & datalist[[4]] == 0) / sum(datalist[[4]] == 0)
  ))
}
