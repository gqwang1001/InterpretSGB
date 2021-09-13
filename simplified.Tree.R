simplified.Tree = function(model,
                           datalist,
                           cutoff = 0.8,
                           plot.name = NaN,
                           seed = 1,
                           top3 = T) {
  
  # model: the model output from the XGboost
  # datalist: the data used for fitting the XGboost model.
  # plot.name: if NaN, no figure is plotted, else show figure with the given name.
  # top3: logic: whether the top 3 features are restricted. 
  
  require(caret)
  # require(rpart.plot)
  # require(rattle)
  require(xgboost)
  require(vip)
  library(xgboost)
  library(caret)
  library(dplyr)
  
  options(warn = -1)
  print(datalist[[3]] %>% head())
  shapvi <- vip::vi_shap(model, train=as.matrix(datalist[[1]][, model$feature_names]), pred_wrapper = predict) %>% 
    filter(Importance!=0) %>% 
    arrange(desc(Importance)) %>% 
    mutate(CumSum.Imp = cumsum(Importance)/sum(Importance))
  valNames <- shapvi$Variable
  imp.val <- shapvi$Variable[shapvi$CumSum.Imp<cutoff]
  imp.gain <- xgb.importance(model$feature_names, model)
  imp.gain.top = imp.gain$Feature[1:min(10,length(imp.gain$Feature))]#[1:5]
  
  if (top3 & length(imp.val) >= 3)
    imp.val = imp.val[1:3]
  if (top3 & !is.na(plot.name))
    plot.name = paste0("Top3_", plot.name)
  if (length(imp.val) == 1)
    imp.val = valNames[1:2]
  
  pred.xgb.train <- predict(model, as.matrix(datalist[[1]][, model$feature_names]), outputmargin = TRUE)
  print(head(pred.xgb.train))
  pred.xgb = predict(model, as.matrix(datalist[[2]][, model$feature_names]), outputmargin = T)
  
  dataInPred <- data.frame(datalist[[1]][, imp.val], logOdds = pred.xgb.train)
  dataInPred.gain <- data.frame(datalist[[1]][, imp.gain.top], logOdds = pred.xgb.train)
  summaryPred.xgb = summaryPred(pred.xgb, datalist)

  # fit decision tree -------------------------------------------------------
  spl.shap <- simplified.models(dat.train = dataInPred, dat.test = datalist[[2]], modelname = "tree", seed,pred.xgb.train = pred.xgb.train, pred.xgb = pred.xgb)
  spl.gain <- simplified.models(dat.train = dataInPred.gain, dat.test = datalist[[2]], modelname = "tree", seed, pred.xgb.train = pred.xgb.train, pred.xgb = pred.xgb)
  
  if (!is.na(plot.name)) {
    library(rattle)
    png(filename = plot.name)
    rattle::fancyRpartPlot(spl.shap$model$finalModel)
    dev.off()
  }

  summaryPred.splTree.shap = summaryPred(spl.shap$pred, datalist)
  summaryPred.splTree.gain = summaryPred(spl.gain$pred, datalist)
  
  # fit linear regression model ---------------------------------------------
  fit.lm = simplified.models(dat.train = dataInPred, dat.test = datalist[[2]], modelname = "lm", seed=seed, pred.xgb.train = pred.xgb.train, pred.xgb = pred.xgb)
  fit.lm.gain = simplified.models(dat.train = dataInPred.gain, dat.test = datalist[[2]], modelname = "lm", seed=seed, pred.xgb.train = pred.xgb.train, pred.xgb = pred.xgb)
  
  summaryPred.lm = summaryPred(fit.lm$pred, datalist)
  summaryPred.lm.gain = summaryPred(fit.lm.gain$pred, datalist)
  
  # print(fit.lm$vi$data)
  
  # print(summaryPred.lm)
  # metric 2: use the loss function -----------------------------------------------
  # pred.cv.lm = pred.cv.simplied.models(modelname = "lm", dat.train = dataInPred, seed=seed)
  # pred.cv.spltree = pred.cv.simplied.models(modelname = "tree", dat.train = dataInPred, seed=seed)
  # pred.cv.lm.gain = pred.cv.simplied.models(modelname = "lm", dat.train = dataInPred.gain, seed=seed)
  # pred.cv.spltree.gain = pred.cv.simplied.models(modelname = "tree", dat.train = dataInPred.gain, seed=seed)
  # 
  # print(cbind(pred.cv.lm[1:30], pred.xgb.train[1:30]))
  
  # evalLoss.lm.train = pred.cv.simplied.models(modelname = "lm",pred.xgb.train = pred.xgb.train, dat.train = dataInPred,datalist = datalist, seed=seed)
  # evalLoss.spltree.train = pred.cv.simplied.models(modelname = "tree",pred.xgb.train = pred.xgb.train, dat.train = dataInPred, datalist = datalist,seed=seed)
  # evalLoss.lm.train.gain = pred.cv.simplied.models(modelname = "lm",pred.xgb.train = pred.xgb.train, dat.train = dataInPred.gain, datalist = datalist,seed=seed)
  # evalLoss.spltree.train.gain = pred.cv.simplied.models(modelname = "tree",pred.xgb.train = pred.xgb.train, dat.train = dataInPred.gain,datalist = datalist, seed=seed)
  
  # evalLoss.lm.train = lossRatio(fit.lm$pred.train, pred.xgb.train, datalist[[1]])
  # evalLoss.spltree.train = lossRatio(spl.shap$pred.train, pred.xgb.train, datalist[[1]])
  # evalLoss.lm.train.gain = lossRatio(pred.cv.lm.gain, pred.xgb.train, datalist[[1]])
  # evalLoss.spltree.train.gain = lossRatio(pred.cv.spltree.gain, pred.xgb.train, datalist[[1]])
  # 
  testIdxs = 1:100
  evalLoss.lm.test = lossRatio(fit.lm$pred[testIdxs], pred.xgb[testIdxs], datalist[[2]][testIdxs,])
  evalLoss.spltree.test = lossRatio(spl.shap$pred[testIdxs], pred.xgb[testIdxs], datalist[[2]][testIdxs,])
  evalLoss.lm.test.gain = lossRatio(fit.lm.gain$pred[testIdxs], pred.xgb[testIdxs], datalist[[2]][testIdxs,])
  evalLoss.spltree.test.gain = lossRatio(spl.gain$pred[testIdxs], pred.xgb[testIdxs], datalist[[2]][testIdxs,])
  
  # print(evalLoss.lm.train$loss)
  # print(evalLoss.spltree.train$loss)
  # print(evalLoss.lm.train$lossvec)
  # print(evalLoss.spltree.train$lossvec)
  print(evalLoss.lm.test)
  print(evalLoss.spltree.test)
  # print(lossRatio(fit.lm$pred.train, pred.xgb.train, datalist[[1]]))
  # print(lossRatio(spl.shap$pred.train, pred.xgb.train, datalist[[1]]))
  # 
  return(
    list(
      xgbModel = model,
      imporant_variables = imp.val,
      imporant_variables_xgb_all_shap = shapvi,
      imporant_variables_gain = imp.gain,
      
      simplified_tree = spl.shap,
      simplified_tree_gain = spl.gain,
      
      summaryPred.xgb = summaryPred.xgb,
      summaryPred.splTree.shap = summaryPred.splTree.shap,
      summaryPred.splTree.gain = summaryPred.splTree.gain,

      summaryPred.lm = summaryPred.lm,
      summaryPred.lm.gain = summaryPred.lm.gain,
      
      lmModel = fit.lm,
      lmModel_gain = fit.lm.gain,
    
      # evalLoss.lm.train = evalLoss.lm.train,
      # evalLoss.spltree.train = evalLoss.spltree.train,   
      # evalLoss.lm.train.gain = evalLoss.lm.train.gain,
      # evalLoss.spltree.train.gain = evalLoss.spltree.train.gain,
      
      evalLoss.lm.test = evalLoss.lm.test,
      evalLoss.spltree.test = evalLoss.spltree.test,
      evalLoss.lm.test.gain = evalLoss.lm.test.gain,
      evalLoss.spltree.test.gain = evalLoss.spltree.test.gain
      )
  )
}

pred.cv.simplied.models = function(modelname = "lm", dat.train, datalist, pred.xgb.train, seed, nfold = 10){
  
  preds = rep(NA, length(dat.train[,1]))
  set.seed(seed)
  testSet = caret::createFolds(dat.train[,1], k=nfold, list = T, returnTrain = F)
  lossvec = rep(NA, nfold)
  
  # print(testSet[[1]])
  for (k in 1:nfold){
    Traindata = dat.train[-testSet[[k]],]
    Testdata = dat.train[testSet[[k]],]
    # Testdata = datalist[[1]][testSet[[k]],]
    splModel = simplified.models(dat.train = Traindata, dat.test = Testdata, modelname = modelname, seed = seed,
                                 pred.xgb.train = Traindata$logOdds, pred.xgb = Testdata$logOdds)
    preds[testSet[[k]]] = splModel$pred
    
    lossvec[k]=lossRatio(preds.s = splModel$pred, preds.full = pred.xgb.train[testSet[[k]]], dat = datalist[[1]][testSet[[k]],])
    # if (k==1) print(splModel$pred)
  }
  return(list(preds = preds,lossvec=lossvec, loss = mean(lossvec)))
}

simplified.models <- function(dat.train, dat.test, modelname = "lm", seed, pred.xgb.train=NaN, pred.xgb=NaN) {
  vi = list()
  vi$data$Variable = NA
  if (modelname == "tree") { 
    trctrl <-
      trainControl(method = "repeatedcv",
                   number = 10,
                   repeats = 3)
    
    set.seed(seed)
    model.fit <- train(
      logOdds ~ .,
      data = dat.train,
      method = "rpart2",
      trControl = trctrl,
      tuneLength = 3
    )
    pred = predict(model.fit, dat.test)
    pred.train = predict(model.fit, dat.train)
    vi = vip::vi_shap(model.fit$finalModel, pred_wrapper = predict) %>% filter(Importance!=0) %>% arrange(desc(Importance))
    
  } else if (modelname == "lm") {
    fit.lm.full = lm(data = dat.train, logOdds ~ .)
    model.fit = step(fit.lm.full, trace = 0)
    pred = predict(model.fit, dat.test)
    pred.train = predict(model.fit, dat.train)
    if (length(coef(model.fit))>1) vi = vip::vip(model.fit)
  }
  
  if (is.nan(pred.xgb) && is.nan(pred.xgb.train)){
    
    return(list(model = model.fit,
                pred = pred,
                pred.train = pred.train,
                vi = vi))
    
  }else{
    acc.comp.train = 2 * sum(1 * (((pred.train > 0) == (pred.xgb.train > 0)))) / (length(pred.train) + length(pred.xgb.train))
    rmse.splVsxgb.train = Rsquare(pred.train, pred.xgb.train)
    
    acc.comp = 2 * sum(1 * (((pred > 0) == (pred.xgb > 0)))) / (length(pred) + length(pred.xgb))
    rmse.splVsxgb = Rsquare(pred, pred.xgb)
    
    return(list(model = model.fit,
                pred = pred,
                pred.train = pred.train,
                vi=vi,
                acc.comp.train = acc.comp.train,
                acc.comp = acc.comp,
                rmse.splVsxgb.train = rmse.splVsxgb.train,
                rmse.splVsxgb = rmse.splVsxgb))
    
  }
  
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
  preds.s = as.vector(preds.s)
  preds.full = as.vector(preds.full)
  return(evalueLoss(preds.s, dat)/evalueLoss(preds.full, dat))
  # return(evalueLoss(preds.s, dat))
}

evalueLoss <- function(preds, dat) {
  N <- nrow(dat)
  labels <- rep(NA, N)
  labels[dat$trt01p == 1 & dat$evnt == 1] <- -1000 - dat$aval[dat$trt01p == 1 & dat$evnt == 1]
  labels[dat$trt01p == 1 & dat$evnt == 0] <- -1 - dat$aval[dat$trt01p == 1 & dat$evnt == 0]
  labels[dat$trt01p == 0 & dat$evnt == 1] <- 1000 + dat$aval[dat$trt01p == 0 & dat$evnt == 1]
  labels[dat$trt01p == 0 & dat$evnt == 0] <- dat$aval[dat$trt01p == 0 & dat$evnt == 0]
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
  aval[abs(labels) >= (1000)] <- abs(labels[abs(labels) >= (1000)]) - 1000
  aval[labels < 0 & labels > -1000] <- -labels[labels < 0 & labels > -1000] - 1
  aval[labels >= 0 & labels < 1000] <- labels[labels >=  0 & labels < 1000]
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
      gH1.r1.denom <- sum((denom$trt01p == arm.val[1]) * denom$pred)
      gH0.r1.denom <- sum((denom$trt01p == arm.val[2]) * denom$pred)
      gH1.r2.denom <- sum((denom$trt01p == arm.val[1]) * (1 - denom$pred))
      gH0.r2.denom <- sum((denom$trt01p == arm.val[2]) * (1 - denom$pred))
      if (gH1.r1.denom > 0) {
        H1.r1 <- H1.r1 + sum((nume$trt01p == arm.val[1]) * (nume$evnt == 1) * nume$pred) / gH1.r1.denom
      }
      if (gH0.r1.denom > 0) {
        H0.r1 <- H0.r1 + sum((nume$trt01p == arm.val[2]) * (nume$evnt == 1) * nume$pred) / gH0.r1.denom
      }
      if (gH1.r2.denom > 0) {
        H1.r2 <- H1.r2 + sum((nume$trt01p == arm.val[1]) * (nume$evnt == 1) * (1 - nume$pred)) / gH1.r2.denom
      }
      if (gH0.r2.denom > 0) {
        H0.r2 <- H0.r2 + sum((nume$trt01p == arm.val[2]) * (nume$evnt == 1) * (1 - nume$pred)) / gH0.r2.denom
      }
    }
    rmst.diff.r1 <- rmst.diff.r1 + (exp(-H1.r1) - exp(-H0.r1)) * dt[i + 1]
    rmst.diff.r2 <- rmst.diff.r2 + (exp(-H1.r2) - exp(-H0.r2)) * dt[i + 1]
  }
  # err <- (-1) * (exp(mean(log(km.dat$pred))) * rmst.diff.r1 - exp(mean(log((1 - km.dat$pred)))) * rmst.diff.r2)
  err <- (-1) * ((mean((km.dat$pred))) * rmst.diff.r1 - (mean(((1 - km.dat$pred)))) * rmst.diff.r2)
  return(err)
}
