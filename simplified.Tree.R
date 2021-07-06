
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
  shap_cumsum <- cumsum(sort(shapresults, decreasing = T) / sum(shapresults))
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
  
  pred.xgb.train =predict(model, as.matrix(datalist[[1]][, model$feature_names]), outputmargin = TRUE)
  dataInPred <- data.frame(datalist[[1]][, imp.val], logOdds = pred.xgb.train)
  # fit decision tree -------------------------------------------------------
  dataInPred.gain <- data.frame(datalist[[1]][, imp.gain.top], logOdds = pred.xgb.train)
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
  set.seed(seed)
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
  pred.xgb =predict(model, as.matrix(datalist[[2]][, model$feature_names]), outputmargin = T)

  # prediction on training data
  pred.simple.train = predict(dtree_fit, datalist[[1]])
  pred.simple.gain.train  = predict(dtree_fit_gain, datalist[[1]])
  
    
  rmse.splTreeVsxgb =  Rsquare(pred.simple, pred.xgb)
  rmse.splTreeVsxgb.train = Rsquare(pred.simple.train, pred.xgb.train)
  
  summaryPred.splTree.shap = summaryPred(pred.simple, datalist)
  summaryPred.splTree.gain = summaryPred(pred.simple.gain, datalist)
  summaryPred.xgb = summaryPred(pred.xgb, datalist)
  
  acc.comp.splTreeVsxgb = 2 * sum(1 * (((pred.simple > 0) == (pred.xgb > 0)))) / (length(pred.simple) + length(pred.xgb))
  acc.comp.splTreeVsxgb.train = 2 * sum(1 * (((pred.simple.train> 0) == (pred.xgb.train > 0)))) / (length(pred.simple.train) + length(pred.xgb.train))
  
  # fit linear regression model ---------------------------------------------
  fit.lm.full = lm(data = dataInPred, logOdds ~ .)
  fit.lm = step(fit.lm.full, trace = 0)
  pred.lm.train = predict(fit.lm, datalist[[1]])
  pred.lm = predict(fit.lm, datalist[[2]])
  acc.comp.lm = 2 * sum(1 * (((pred.lm > 0) == (pred.xgb > 0)))) / (length(pred.lm) + length(pred.xgb))
  acc.comp.lm.train = 2 * sum(1 * (((pred.lm.train > 0) == (pred.xgb.train > 0)))) / (length(pred.lm.train) + length(pred.xgb.train))
  
  rmse.lmVsxgb = Rsquare(pred.lm, pred.xgb)
  rmse.lmVsxgb.train = Rsquare(pred.lm.train, pred.xgb.train)
  summaryPred.lm = summaryPred(pred.lm, datalist)
  
  evalLoss.lm.train = lossRatio(pred.lm.train, pred.xgb.train, datalist[[1]])
  evalLoss.spltree.train = lossRatio(pred.simple.train, pred.xgb.train, datalist[[1]])
  
  evalLoss.lm = lossRatio(pred.lm,pred.xgb, datalist[[2]])
  evalLoss.spltree = lossRatio(pred.simple,pred.xgb, datalist[[2]])
  
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
      acc.comp.splTreeVsxgb.train = acc.comp.splTreeVsxgb.train,
      
      RMSE.spltree = rmse.splTreeVsxgb,
      RMSE.spltree.train = rmse.splTreeVsxgb.train,
      
      lm_pred = pred.lm,
      acc.comp.lm = acc.comp.lm,
      acc.comp.lm.train = acc.comp.lm.train,
      
      RMSE.lm = rmse.lmVsxgb,
      RMSE.lm.train = rmse.lmVsxgb.train,
      summaryPred.lm = summaryPred.lm,
      lmModel = fit.lm,
      
      evalLoss.lm = evalLoss.lm,
      evalLoss.spltree=evalLoss.spltree,
      evalLoss.lm.train = evalLoss.lm.train,
      evalLoss.spltree.train = evalLoss.spltree.train
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

Rsquare = function(m.s, m.full){
  return( sum((m.s-m.full)^2)/sum((m.full-mean(m.full))^2) )  
}
  
lossRatio = function(preds.s,preds.full, dat){
  return(evalueLoss(as.vector(preds.s), dat)/evalueLoss(as.vector(preds.full), dat))
}

evalueLoss <- function(preds, dat) {
  N <- nrow(dat)
  labels <- rep(NA, N)
  labels[dat$trt01p == 1 &
           dat$evnt == 1] <- -1000 - dat$aval[dat$trt01p == 1 & dat$evnt == 1]
  labels[dat$trt01p == 1 &
           dat$evnt == 0] <- -1 - dat$aval[dat$trt01p == 1 & dat$evnt == 0]
  labels[dat$trt01p == 0 &
           dat$evnt == 1] <- 1000 + dat$aval[dat$trt01p == 0 & dat$evnt == 1]
  labels[dat$trt01p == 0 &
           dat$evnt == 0] <- dat$aval[dat$trt01p == 0 & dat$evnt == 0]
  dat$evnt <- NULL
  dat$aval <- NULL
  dat$trt01p <- NULL
  dtrain <- xgb.DMatrix(as.matrix(dat), label = labels)
  
  
  labels <- getinfo(dtrain, "label")
  trt01p <- rep(NA, length(labels))
  evnt <- rep(NA, length(labels))
  aval <- rep(NA, length(labels))
  trt01p[labels < 0] <- 1
  trt01p[labels >= 0] <- 0
  evnt[abs(labels) >= (1000)] <- 1
  evnt[abs(labels) < (1000)] <- 0
  aval[abs(labels) >= (1000)] <-
    abs(labels[abs(labels) >= (1000)]) - 1000
  aval[labels < 0 &
         labels > -1000] <- -labels[labels < 0 & labels > -1000] - 1
  aval[labels >= 0 &
         labels < 1000] <- labels[labels >=  0 & labels < 1000]
  arm.val <- c(1, 0)
  km.dat <- data.frame(trt01p, evnt, aval)
  km.dat$pred <- 1 / (1 + exp(-preds))
  km.dat <- km.dat[order(km.dat$aval), ]
  n <- dim(km.dat)[1]
  utime <- unique(km.dat$aval[km.dat$evnt == 1])
  dt <- utime - c(0, utime[1:length(utime) - 1])
  rmst.diff.r1 <- 0
  rmst.diff.r2 <- 0
  for (i in 0:(length(utime) - 1)) {
    if (i == 0) {
      H1.r1 <- 0
      H0.r1 <- 0
      H1.r2 <- 0
      H0.r2 <- 0
    }
    else {
      denom <- subset(km.dat, aval >= utime[i])
      nume <- subset(km.dat, aval == utime[i] & evnt == 1)
      gH1.r1.denom <-
        sum((denom$trt01p == arm.val[1]) * denom$pred)
      gH0.r1.denom <-
        sum((denom$trt01p == arm.val[2]) * denom$pred)
      gH1.r2.denom <-
        sum((denom$trt01p == arm.val[1]) * (1 - denom$pred))
      gH0.r2.denom <-
        sum((denom$trt01p == arm.val[2]) * (1 - denom$pred))
      if (gH1.r1.denom > 0) {
        H1.r1 <-
          H1.r1 + sum((nume$trt01p == arm.val[1]) * (nume$evnt == 1) * nume$pred) / gH1.r1.denom
      }
      if (gH0.r1.denom > 0) {
        H0.r1 <-
          H0.r1 + sum((nume$trt01p == arm.val[2]) * (nume$evnt == 1) * nume$pred) / gH0.r1.denom
      }
      if (gH1.r2.denom > 0) {
        H1.r2 <-
          H1.r2 + sum((nume$trt01p == arm.val[1]) * (nume$evnt == 1) * (1 - nume$pred)) / gH1.r2.denom
      }
      if (gH0.r2.denom > 0) {
        H0.r2 <-
          H0.r2 + sum((nume$trt01p == arm.val[2]) * (nume$evnt == 1) * (1 - nume$pred)) / gH0.r2.denom
      }
    }
    rmst.diff.r1 <-
      rmst.diff.r1 + (exp(-H1.r1) - exp(-H0.r1)) * dt[i + 1]
    rmst.diff.r2 <-
      rmst.diff.r2 + (exp(-H1.r2) - exp(-H0.r2)) * dt[i + 1]
  }
  err <-
    (-1) * (sum(km.dat$pred) * rmst.diff.r1 - sum(1 -  km.dat$pred) * rmst.diff.r2)
  return(err)
}
