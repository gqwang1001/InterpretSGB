prepareData <- 
  function(dat){
  N <- nrow(dat)
  labels <- rep(NA, N)
  labels[dat$trt01p == 1 & dat$evnt == 1] <- -1000 - dat$aval[dat$trt01p == 
                                                                1 & dat$evnt == 1]
  labels[dat$trt01p == 1 & dat$evnt == 0] <- -1 - dat$aval[dat$trt01p == 
                                                             1 & dat$evnt == 0]
  labels[dat$trt01p == 0 & dat$evnt == 1] <- 1000 + dat$aval[dat$trt01p == 
                                                               0 & dat$evnt == 1]
  labels[dat$trt01p == 0 & dat$evnt == 0] <- dat$aval[dat$trt01p == 
                                                        0 & dat$evnt == 0]
  dat$evnt <- NULL
  dat$aval <- NULL
  dat$trt01p <- NULL
  
  MatData <- xgb.DMatrix(as.matrix(dat), label = labels)
  return(list(dat=dat, labels = labels, DatMat = MatData))
}