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
