##---------------------------##
#     Input Parameters       #
##---------------------------##

## (1) dat: data frame with predictors ONLY, not include trt01p, aval and evnt
## (2) labels: embedded label
##     Use below mapping function 
## embed trt01p, aval and evnt into y, which serves as the 'label' for xgboost ##
#  labels <- rep(NA,N)
#  labels[dat$trt01p==1 & dat$evnt==1] <- -1000-dat$aval[dat$trt01p==1 & dat$evnt==1] 
#  labels[dat$trt01p==1 & dat$evnt==0] <- -1-dat$aval[dat$trt01p==1 & dat$evnt==0] 
#  labels[dat$trt01p==0 & dat$evnt==1] <- 1000+dat$aval[dat$trt01p==0 & dat$evnt==1] 
#  labels[dat$trt01p==0 & dat$evnt==0] <- dat$aval[dat$trt01p==0 & dat$evnt==0] 

##---------------------------##
#     Return Values          #
##---------------------------##
# Return a vector with 2 elements, first one for whether signal has top rank of variable importance
# second is a classification error in terms of signal group 


library(survminer)
library(survival)
library(dplyr)
#library(survRM2)
library(xgboost)
#library(DiagrammeR)
#library(purrr)
library(survRM2)

MyBoost <- function(dat,otr,dat.test,tc){
  
  if('s2' %in% colnames(dat)){
    p<-52
  } else {
    p<-51
  }
  
  if("Pembro" %in% unique(dat$trt01p)){
    dat$trt01p <- as.numeric(dat$trt01p=="Pembro")
  } else {
    dat$trt01p <- as.numeric(dat$trt01p==1)
  }
  
  N <- nrow(dat)
  labels <- rep(NA,N)
  labels[dat$trt01p==1 & dat$evnt==1] <- -1000-dat$aval[dat$trt01p==1 & dat$evnt==1] 
  labels[dat$trt01p==1 & dat$evnt==0] <- -1-dat$aval[dat$trt01p==1 & dat$evnt==0] 
  labels[dat$trt01p==0 & dat$evnt==1] <- 1000+dat$aval[dat$trt01p==0 & dat$evnt==1] 
  labels[dat$trt01p==0 & dat$evnt==0] <- dat$aval[dat$trt01p==0 & dat$evnt==0] 
  dat2 <- dat.test
  dat$evnt <- NULL; dat$aval <- NULL; dat$trt01p <- NULL;
  dat.test$evnt <- NULL; dat.test$aval<-NULL; dat.test$trt01p <- NULL;
  #----------------------------------------------------------------------#
  #                     Step 2: Customized loss and error function                  #
  #----------------------------------------------------------------------#
  
  rmst_loss <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")
    
    ## (0) Decode y to trt01p, aval and evnt ##
    trt01p<-rep(NA,length(labels))
    evnt<-rep(NA,length(labels))
    aval<-rep(NA,length(labels))
    trt01p[labels< 0] <- 1
    trt01p[labels>= 0] <- 0
    evnt[abs(labels)>= (1000)] <- 1
    evnt[abs(labels)< (1000)] <- 0
    aval[abs(labels)>= (1000)] <- abs(labels[abs(labels)>= (1000)])-1000
    aval[labels<0 & labels>-1000] <- -labels[labels<0 & labels>-1000]-1
    aval[labels>=0 & labels<1000] <- labels[labels>=0 & labels<1000]
    
    arm.val <- c(1,0)
    ## (1) Get Time to event Data Ready ##
    km.dat <- data.frame(trt01p,evnt,aval)
    n<-dim(km.dat)[1]
    km.dat$id<-c(1:n)
    km.dat$f <- preds
    km.dat$pred <- 1/(1+exp(-preds))
    km.dat$predg <- exp(preds)/(1+exp(preds))^2
    km.dat$predh <- exp(preds)*(1-exp(preds))/(1+exp(preds))^3
    km.dat<-km.dat[order(km.dat$aval),]
    
    
    ## (2) Set up gradient and Hessian ##
    utime <- unique(km.dat$aval[km.dat$evnt==1])
    #utime <-   utime[utime<=tc]
    dt <- utime-c(0,utime[1:length(utime)-1])
    
    rmst.diff.r1 <- 0
    rmst.diff.r2 <- 0
    rmst.diff.r1.g <- 0
    rmst.diff.r2.g <- 0
    rmst.diff.r1.h <- 0
    rmst.diff.r2.h <- 0
    
    
    for(i in 0:(length(utime)-1)){
      if(i==0){
        H1.r1 <- 0
        gH1.r1 <- 0
        hH1.r1 <- 0
        H0.r1 <- 0
        gH0.r1 <- 0
        hH0.r1 <- 0
        
        H1.r2 <- 0
        gH1.r2 <- 0
        hH1.r2 <- 0
        H0.r2 <- 0
        gH0.r2 <- 0
        hH0.r2 <- 0
        
      } else {
        denom <- subset(km.dat,aval>=utime[i])
        nume <- subset(km.dat,aval==utime[i] & evnt==1)
        
        gH1.r1.denom <- sum((denom$trt01p==arm.val[1])*denom$pred)
        gH0.r1.denom <- sum((denom$trt01p==arm.val[2])*denom$pred)
        
        gH1.r2.denom <- sum((denom$trt01p==arm.val[1])*(1-denom$pred))
        gH0.r2.denom <- sum((denom$trt01p==arm.val[2])*(1-denom$pred)) 
        
        ## H1 r1 ##
        if(gH1.r1.denom > 0){
          ## H1 ##
          H1.r1 <- H1.r1 + sum((nume$trt01p==arm.val[1])*(nume$evnt==1)*nume$pred) / gH1.r1.denom
          ## dH1/dp ##
          gH1.r1 <- gH1.r1 + (((km.dat$aval==utime[i])*(km.dat$trt01p==arm.val[1])*(km.dat$evnt==1)*sum((denom$trt01p==arm.val[1])*denom$pred)) - ( sum((nume$trt01p==arm.val[1])*(nume$evnt==1)*nume$pred)*(km.dat$aval>=utime[i])*(km.dat$trt01p==arm.val[1])  )) / gH1.r1.denom^2
          ## d2H1/dp2 ##
          #hH1.r1 <- hH1.r1 + (-2)*(((km.dat$aval==utime[i])*(km.dat$trt01p==arm.val[1])*(km.dat$evnt==1)*sum((denom$trt01p==arm.val[1])*denom$pred)) - ( sum((nume$trt01p==arm.val[1])*(nume$evnt==1)*nume$pred)*(km.dat$aval>=utime[i])*(km.dat$trt01p==arm.val[1])  ))*(km.dat$trt01p==arm.val[1])*(km.dat$aval>=utime[i])/gH1.r1.denom^3
          
          
        }
        
        ## H0 r1##
        if(gH0.r1.denom > 0){
          ## H0 ##
          H0.r1 <- H0.r1 + sum((nume$trt01p==arm.val[2])*(nume$evnt==1)*nume$pred) / gH0.r1.denom
          ## dH1/dp ##
          gH0.r1 <- gH0.r1 + (((km.dat$aval==utime[i])*(km.dat$trt01p==arm.val[2])*(km.dat$evnt==1)*sum((denom$trt01p==arm.val[2])*denom$pred)) - ( sum((nume$trt01p==arm.val[2])*(nume$evnt==1)*nume$pred)*(km.dat$aval>=utime[i])*(km.dat$trt01p==arm.val[2])  )) / gH0.r1.denom^2
          ## d2H1/dp2 ##
          #hH0.r1 <- hH0.r1 + (-2)*(((km.dat$aval==utime[i])*(km.dat$trt01p==arm.val[2])*(km.dat$evnt==1)*sum((denom$trt01p==arm.val[2])*denom$pred)) - ( sum((nume$trt01p==arm.val[2])*(nume$evnt==1)*nume$pred)*(km.dat$aval>=utime[i])*(km.dat$trt01p==arm.val[2])  ))*(km.dat$trt01p==arm.val[2])*(km.dat$aval>=utime[i])/gH0.r1.denom^3
          
        }
        
        #rmst.diff.r1 <- rmst.diff.r1 + (exp(-ch.trt.r1)-exp(-ch.cntl.r1))*dt[i]
        ## H1 r2 ##
        if(gH1.r2.denom > 0){
          ## H1 ##
          H1.r2 <- H1.r2 + sum((nume$trt01p==arm.val[1])*(nume$evnt==1)*(1-nume$pred)) / gH1.r2.denom
          ## dH1/dp ##
          gH1.r2 <- gH1.r2 + ((-1*(km.dat$aval==utime[i])*(km.dat$trt01p==arm.val[1])*(km.dat$evnt==1)*sum((denom$trt01p==arm.val[1])*(1-denom$pred))) - ( sum((nume$trt01p==arm.val[1])*(nume$evnt==1)*(1-nume$pred))*(-1)*(km.dat$aval>=utime[i])*(km.dat$trt01p==arm.val[1])  )) / gH1.r2.denom^2
          ## d2H1/dp2 ##
          #hH1.r2 <- hH1.r2 + (2)*((-1*(km.dat$aval==utime[i])*(km.dat$trt01p==arm.val[1])*(km.dat$evnt==1)*sum((denom$trt01p==arm.val[1])*(1-denom$pred))) - ( sum((nume$trt01p==arm.val[1])*(nume$evnt==1)*(1-nume$pred))*(-1)*(km.dat$aval>=utime[i])*(km.dat$trt01p==arm.val[1])  ))*(km.dat$trt01p==arm.val[1])*(km.dat$aval>=utime[i])/gH1.r2.denom^3
          
        }
        
        ## H0 r2 ##
        if(gH0.r2.denom > 0){
          ## H0 ##
          H0.r2 <- H0.r2 + sum((nume$trt01p==arm.val[2])*(nume$evnt==1)*(1-nume$pred)) / gH0.r2.denom
          ## dH1/dp ##
          gH0.r2 <- gH0.r2 + ((-1*(km.dat$aval==utime[i])*(km.dat$trt01p==arm.val[2])*(km.dat$evnt==1)*sum((denom$trt01p==arm.val[2])*(1-denom$pred))) - ( sum((nume$trt01p==arm.val[2])*(nume$evnt==1)*(1-nume$pred))*(-1)*(km.dat$aval>=utime[i])*(km.dat$trt01p==arm.val[2])  )) / gH0.r2.denom^2
          ## d2H1/dp2 ##
          #hH0.r2 <- hH0.r2 + 2*((-1*(km.dat$aval==utime[i])*(km.dat$trt01p==arm.val[2])*(km.dat$evnt==1)*sum((denom$trt01p==arm.val[2])*(1-denom$pred))) - ( sum((nume$trt01p==arm.val[2])*(nume$evnt==1)*(1-nume$pred))*(-1)*(km.dat$aval>=utime[i])*(km.dat$trt01p==arm.val[2])  ))*(km.dat$trt01p==arm.val[2])*(km.dat$aval>=utime[i])/gH0.r2.denom^3
          
        }
        
      }
      
      
      rmst.diff.r1 <- rmst.diff.r1 + (exp(-H1.r1)-exp(-H0.r1))*dt[i+1]
      rmst.diff.r2 <- rmst.diff.r2 + (exp(-H1.r2)-exp(-H0.r2))*dt[i+1]
      ## Gradient ##
      rmst.diff.r1.g <- rmst.diff.r1.g + (-(exp(-H1.r1)*gH1.r1)+(exp(-H0.r1)*gH0.r1))*dt[i+1]
      rmst.diff.r2.g <- rmst.diff.r2.g + (-(exp(-H1.r2)*gH1.r2)+(exp(-H0.r2)*gH0.r2))*dt[i+1]
      
      ## Hessian ##
      #rmst.diff.r1.h <- rmst.diff.r1.h + (exp(-H1.r1)*gH1.r1^2 - exp(-H1.r1)*hH1.r1 - exp(-H0.r1)*gH0.r1^2 + exp(-H0.r1)*hH0.r1)*dt[i+1] 
      #rmst.diff.r2.h <- rmst.diff.r2.h + (exp(-H1.r2)*gH1.r2^2 - exp(-H1.r2)*hH1.r2 - exp(-H0.r2)*gH0.r2^2 + exp(-H0.r2)*hH0.r2)*dt[i+1]
      
    }
    
    
    g.p <- (sum(km.dat$pred)*rmst.diff.r1.g + rmst.diff.r1 - sum(1-km.dat$pred)*rmst.diff.r2.g + rmst.diff.r2)
    #h.p <- (2*rmst.diff.r1.g + sum(km.dat$pred)*rmst.diff.r1.h + 2*rmst.diff.r2.g - sum(1-km.dat$pred)*rmst.diff.r2.h)
    g <-  km.dat$predg*(-1)*g.p
    #h <- (-1)*( (km.dat$predg)^2 * h.p  + g.p*km.dat$predh)
    g <- g[order(km.dat$id)]
    #h <- h[order(km.dat$id)]
    h<-rep(0.00001,n)
    
    
    return(list(grad = g, hess = h))
    #return(list(grad = g, hess = rep(1,n)))
  }
  
  
  
  evalerror <- function(preds, dtrain) {
    ## (0) Decode y to trt01p, aval and evnt ##
    labels <- getinfo(dtrain, "label")
    trt01p<-rep(NA,length(labels))
    evnt<-rep(NA,length(labels))
    aval<-rep(NA,length(labels))
    trt01p[labels< 0] <- 1
    trt01p[labels>= 0] <- 0
    evnt[abs(labels)>= (1000)] <- 1
    evnt[abs(labels)< (1000)] <- 0
    aval[abs(labels)>= (1000)] <- abs(labels[abs(labels)>= (1000)])-1000
    aval[labels<0 & labels>-1000] <- -labels[labels<0 & labels>-1000]-1
    aval[labels>=0 & labels<1000] <- labels[labels>=0 & labels<1000]
    
    arm.val <- c(1,0)
    
    ## (1) Get Time to event Data Ready ##
    km.dat <- data.frame(trt01p,evnt,aval)
    km.dat$pred <- 1/(1+exp(-preds))
    km.dat<-km.dat[order(km.dat$aval),]
    n<-dim(km.dat)[1]
    
    ## (2) Set up gradient and Hessian ##
    utime <- unique(km.dat$aval[km.dat$evnt==1])
    #utime <-   utime[utime<=tc]
    dt <- utime-c(0,utime[1:length(utime)-1])
    
    rmst.diff.r1 <- 0
    rmst.diff.r2 <- 0
    
    for(i in 0:(length(utime)-1)){
      if(i==0){
        H1.r1 <- 0
        H0.r1 <- 0
        H1.r2 <- 0
        H0.r2 <- 0
        
      } else {
        denom <- subset(km.dat,aval>=utime[i])
        nume <- subset(km.dat,aval==utime[i] & evnt==1)
        
        gH1.r1.denom <- sum((denom$trt01p==arm.val[1])*denom$pred)
        gH0.r1.denom <- sum((denom$trt01p==arm.val[2])*denom$pred)
        
        gH1.r2.denom <- sum((denom$trt01p==arm.val[1])*(1-denom$pred))
        gH0.r2.denom <- sum((denom$trt01p==arm.val[2])*(1-denom$pred)) 
        
        ## H1 r1 ##
        if(gH1.r1.denom > 0){
          ## H1 ##
          H1.r1 <- H1.r1 + sum((nume$trt01p==arm.val[1])*(nume$evnt==1)*nume$pred) / gH1.r1.denom
          
        }
        
        ## H0 r1##
        if(gH0.r1.denom > 0){
          ## H0 ##
          H0.r1 <- H0.r1 + sum((nume$trt01p==arm.val[2])*(nume$evnt==1)*nume$pred) / gH0.r1.denom
          
        }
        
        #rmst.diff.r1 <- rmst.diff.r1 + (exp(-ch.trt.r1)-exp(-ch.cntl.r1))*dt[i]
        ## H1 r2 ##
        if(gH1.r2.denom > 0){
          ## H1 ##
          H1.r2 <- H1.r2 + sum((nume$trt01p==arm.val[1])*(nume$evnt==1)*(1-nume$pred)) / gH1.r2.denom
          
        }
        
        ## H0 r2 ##
        if(gH0.r2.denom > 0){
          ## H0 ##
          H0.r2 <- H0.r2 + sum((nume$trt01p==arm.val[2])*(nume$evnt==1)*(1-nume$pred)) / gH0.r2.denom
          
        }
        
      }
      
      
      rmst.diff.r1 <- rmst.diff.r1 + (exp(-H1.r1)-exp(-H0.r1))*dt[i+1]
      rmst.diff.r2 <- rmst.diff.r2 + (exp(-H1.r2)-exp(-H0.r2))*dt[i+1]
      
    }
    
    err <- (-1)*( sum(km.dat$pred)*rmst.diff.r1 - sum(1-km.dat$pred)*rmst.diff.r2    )
    
    return(list(metric = "RMST_error", value = err))
  }
  
  
  # evalerror3 <- function(preds, dtrain) {
  #   ## (0) Decode y to trt01p, aval and evnt ##
  #   labels <- getinfo(dtrain, "label")
  #   trt01p<-rep(NA,length(labels))
  #   evnt<-rep(NA,length(labels))
  #   aval<-rep(NA,length(labels))
  #   trt01p[labels< 0] <- 1
  #   trt01p[labels>= 0] <- 0
  #   evnt[abs(labels)>= (1000)] <- 1
  #   evnt[abs(labels)< (1000)] <- 0
  #   aval[abs(labels)>= (1000)] <- abs(labels[abs(labels)>= (1000)])-1000
  #   aval[labels<0 & labels>-1000] <- -labels[labels<0 & labels>-1000]-1
  #   aval[labels>=0 & labels<1000] <- labels[labels>=0 & labels<1000]
  #   
  #   arm.val <- c(1,0)
  #   
  #   ## (1) Get Time to event Data Ready ##
  #   km.dat <- data.frame(trt01p,evnt,aval)
  #   km.dat$pred <- 1/(1+exp(-preds))
  #   km.dat$pred <-1*(km.dat$pred>=0.5)
  #   km.dat<-km.dat[order(km.dat$aval),]
  #   n<-dim(km.dat)[1]
  #   
  #   ## (2) Set up gradient and Hessian ##
  #   utime <- unique(km.dat$aval[km.dat$evnt==1])
  #   dt <- utime-c(0,utime[1:length(utime)-1])
  #   
  #   rmst.diff.r1 <- 0
  #   rmst.diff.r2 <- 0
  #   
  #   for(i in 0:(length(utime)-1)){
  #     if(i==0){
  #       H1.r1 <- 0
  #       H0.r1 <- 0
  #       H1.r2 <- 0
  #       H0.r2 <- 0
  #       
  #     } else {
  #       denom <- subset(km.dat,aval>=utime[i])
  #       nume <- subset(km.dat,aval==utime[i] & evnt==1)
  #       
  #       gH1.r1.denom <- sum((denom$trt01p==arm.val[1])*denom$pred)
  #       gH0.r1.denom <- sum((denom$trt01p==arm.val[2])*denom$pred)
  #       
  #       gH1.r2.denom <- sum((denom$trt01p==arm.val[1])*(1-denom$pred))
  #       gH0.r2.denom <- sum((denom$trt01p==arm.val[2])*(1-denom$pred)) 
  #       
  #       ## H1 r1 ##
  #       if(gH1.r1.denom > 0){
  #         ## H1 ##
  #         H1.r1 <- H1.r1 + sum((nume$trt01p==arm.val[1])*(nume$evnt==1)*nume$pred) / gH1.r1.denom
  #         
  #       }
  #       
  #       ## H0 r1##
  #       if(gH0.r1.denom > 0){
  #         ## H0 ##
  #         H0.r1 <- H0.r1 + sum((nume$trt01p==arm.val[2])*(nume$evnt==1)*nume$pred) / gH0.r1.denom
  #         
  #       }
  #       
  #       #rmst.diff.r1 <- rmst.diff.r1 + (exp(-ch.trt.r1)-exp(-ch.cntl.r1))*dt[i]
  #       ## H1 r2 ##
  #       if(gH1.r2.denom > 0){
  #         ## H1 ##
  #         H1.r2 <- H1.r2 + sum((nume$trt01p==arm.val[1])*(nume$evnt==1)*(1-nume$pred)) / gH1.r2.denom
  #         
  #       }
  #       
  #       ## H0 r2 ##
  #       if(gH0.r2.denom > 0){
  #         ## H0 ##
  #         H0.r2 <- H0.r2 + sum((nume$trt01p==arm.val[2])*(nume$evnt==1)*(1-nume$pred)) / gH0.r2.denom
  #         
  #       }
  #       
  #     }
  #     
  #     
  #     rmst.diff.r1 <- rmst.diff.r1 + (exp(-H1.r1)-exp(-H0.r1))*dt[i+1]
  #     rmst.diff.r2 <- rmst.diff.r2 + (exp(-H1.r2)-exp(-H0.r2))*dt[i+1]
  #     
  #   }
  #   
  # 
  #  err <- (-1)*( sum(km.dat$pred)*rmst.diff.r1 - sum(1-km.dat$pred)*rmst.diff.r2    )
  #   
  #   return(list(metric = "RMST_error", value = err))
  # }
  # 
  # 
  # ## Loss functon based on partial log likelihood ##
  # rmst_loss_pl3 <- function(preds, dtrain) {
  #   labels <- getinfo(dtrain, "label")
  #   
  #   ## (0) Decode y to trt01p, aval and evnt ##
  #   trt01p<-rep(NA,length(labels))
  #   evnt<-rep(NA,length(labels))
  #   aval<-rep(NA,length(labels))
  #   trt01p[labels< 0] <- 1
  #   trt01p[labels>= 0] <- -1
  #   evnt[abs(labels)>= (1000)] <- 1
  #   evnt[abs(labels)< (1000)] <- 0
  #   aval[abs(labels)>= (1000)] <- abs(labels[abs(labels)>= (1000)])-1000
  #   aval[labels<0 & labels>-1000] <- -labels[labels<0 & labels>-1000]-1
  #   aval[labels>=0 & labels<1000] <- labels[labels>=0 & labels<1000]
  #   
  #   
  #   ## (1) Get Time to event Data Ready ##
  #   km.dat <- data.frame(trt01p,evnt,aval)
  #   n<-dim(km.dat)[1]
  #   km.dat$id<-c(1:n)
  #   km.dat$f <- preds
  #   
  #   logsum <- rep(0,n)
  #   Thesum <- rep(0,n)
  #   
  #   for(i in 1:n){
  #     id <- which(km.dat$aval>=km.dat$aval[i])
  #     tmp <- km.dat[id,]
  #     Thesum[i] <- sum(  exp(tmp$f*tmp$trt01p)  )
  #     #logsum[i] <- log(Thesum[i])
  #   }
  #   
  #   #l<-sum( 1*(km.dat$evnt==1)*(  km.dat$f*km.dat$trt01p  ) - logsum  )
  #   
  #   g.p <- rep(NA,n)
  #   for(i in 1:n){
  #     tmp.aval <- km.dat$aval[i]
  #     tmp.evnt <- km.dat$evnt[i]
  #     tmp.trt01p <- km.dat$trt01p[i]
  #     tmp.pred <- km.dat$f[i]
  #     g.p[i] <- 1*(tmp.evnt==1)*tmp.trt01p - sum(  (1*(km.dat$evnt==1)*(km.dat$aval<=tmp.aval)) *exp(tmp.pred*tmp.trt01p)*tmp.trt01p  / Thesum   )
  #   }
  #   
  #   
  #   g <-  -g.p
  #   
  #   h<- rep(0.000001,n)
  #   
  #   
  #   return(list(grad = g, hess = h))
  #   
  # }
  # 
  # 
  # ## Error functon based on partial log likelihood ##
  # evalerrorPL3 <- function(preds, dtrain) {
  #   labels <- getinfo(dtrain, "label")
  #   
  #   ## (0) Decode y to trt01p, aval and evnt ##
  #   trt01p<-rep(NA,length(labels))
  #   evnt<-rep(NA,length(labels))
  #   aval<-rep(NA,length(labels))
  #   trt01p[labels< 0] <- 1
  #   trt01p[labels>= 0] <- -1
  #   evnt[abs(labels)>= (1000)] <- 1
  #   evnt[abs(labels)< (1000)] <- 0
  #   aval[abs(labels)>= (1000)] <- abs(labels[abs(labels)>= (1000)])-1000
  #   aval[labels<0 & labels>-1000] <- -labels[labels<0 & labels>-1000]-1
  #   aval[labels>=0 & labels<1000] <- labels[labels>=0 & labels<1000]
  #   
  #   
  #   ## (1) Get Time to event Data Ready ##
  #   km.dat <- data.frame(trt01p,evnt,aval)
  #   n<-dim(km.dat)[1]
  #   km.dat$f <- preds
  #   
  #   logsum <- rep(0,n)
  #   Thesum <- rep(0,n)
  #   
  #   for(i in 1:n){
  #     id <- which(km.dat$aval>=km.dat$aval[i])
  #     tmp <- km.dat[id,]
  #     Thesum[i] <- sum(  exp(tmp$f*tmp$trt01p)  )
  #     logsum[i] <- log(Thesum[i])
  #   }
  #   
  #   l<-sum( 1*(km.dat$evnt==1)*(  km.dat$f*km.dat$trt01p   - logsum  ) )
  #   err <- -l
  #   
  #   
  #   
  #   
  #   return(list(metric = "PL_error", value = err))
  # }
  
  
  #----------------------------------------------------------------------#
  #                     Step 3: Let's boost                  #
  #----------------------------------------------------------------------#
  
  #----------------------------------------------------------------------#
  # (3.1)                    Grid Search                 #
  #----------------------------------------------------------------------#
  hyper_grid <- expand.grid(
    eta =c(.001, .005, .01, .05,0.1), #c(.01, .025, .05, .075),
    max_depth = c(2,3),  #c(2, 4, 6),
    # min_child_weight = c(1, 3, 5, 7),
    subsample = 1,
    #colsample_bytree = c(.8, 1),
    #lambda =c(1,3,5),
    optimal_trees = 0,
    min_error = 0
  )
  
  start.time <- Sys.time()
  for(i in 1:nrow(hyper_grid)) {
    print(i)
    
    # create parameter list
    params <- list(
      eta = hyper_grid$eta[i],
      max_depth = hyper_grid$max_depth[i],
      lambda = 1,
      min_child_weight = 0,
      subsample = 1,
      colsample_bytree = 1
    )
    
    # reproducibility
    set.seed(123)
    
    # train model
    xgb.tune <- xgb.cv(
      params = params,
      data = as.matrix(dat),
      label = labels,
      nrounds = 500,
      nfold = 5,
      objective = rmst_loss,
      eval_metric = evalerror,
      maximize=F,
      verbose = 0,               # silent,
      early_stopping_rounds = 10 # stop if no improvement for 5 consecutive trees
    )
    
    # add min training error and trees to grid
    hyper_grid$optimal_trees[i] <- which.min(xgb.tune$evaluation_log$test_RMST_error_mean)
    hyper_grid$min_error[i] <- min(xgb.tune$evaluation_log$test_RMST_error_mean)
  }
  end.time <- Sys.time()
  end.time-start.time
  
  # hyper_grid %>%
  #   dplyr::arrange(min_error) %>%
  #   head(10)
  # print('Grid Search Finished')
  # #----------------------------------------------------------------------#
  # #  (3.2) CV to find the 'optimal' # of trees              #
  # #----------------------------------------------------------------------#
  # hyper_grid <- hyper_grid[order(hyper_grid$min_error),]
  # dtrain <- xgb.DMatrix(as.matrix(dat),label = labels)
  # param <- list(max_depth = hyper_grid$max_depth[1], eta = hyper_grid$eta[1], silent = 1, nthread = 2, 
  #               objective = rmst_loss, eval_metric = evalerror3,verbose = 2,lambda=hyper_grid$lambda[1],base_score=0,min_child_weight=0,colsample_bytree=hyper_grid$colsample_bytree[1])
  # model.cv <- xgb.cv(param, dtrain, nrounds = 500,nfold = 5,maximize=F,early_stopping_rounds = 10)
  # 
  # sum.cv <- model.cv$evaluation_log %>%
  # dplyr::summarise(
  # ntrees.train = which(train_RMST_error_mean == min(train_RMST_error_mean))[1],
  # rmse.train   = min(train_RMST_error_mean),
  # ntrees.test  = which(test_RMST_error_mean == min(test_RMST_error_mean))[1],
  # rmse.test   = min(test_RMST_error_mean)
  #    )
  # 
  # print('CV for ntree Finished')
  # print(sum.cv)
  # print(sum.cv[[3]])
  # ggplot(model.cvl$evaluation_log) +
  #   geom_line(aes(iter, train_RMST_error_mean), color = "red") +
  #   geom_line(aes(iter, test_RMST_error_mean), color = "blue")
  
  #----------------------------------------------------------------------#
  #  (3.3) Train a model based on the best CV parameter set                #
  #----------------------------------------------------------------------#
  
  hyper_grid <- hyper_grid[order(hyper_grid$min_error),]
  print (hyper_grid[1:5,])
  dtrain <- xgb.DMatrix(as.matrix(dat),label = labels)
  param <- list(max_depth = hyper_grid$max_depth[1], eta = hyper_grid$eta[1], silent = 1, 
                objective = rmst_loss, eval_metric = evalerror,verbose = 1,lambda=1,base_score=0,colsample_bytree=1,min_child_weight=0,subsample = 1)
  watchlist <- list(train = dtrain)
  #param <- list(max_depth = 4, eta = 0.01, silent = 1, nthread = 2, 
  #objective = rmst_loss, eval_metric = evalerror,verbose = 2,lambda=1,base_score=0,min_child_weight=0,colsample_bytree=0.7)
  #param <- list(max_depth = hyper_grid$max_depth[1], eta = hyper_grid$eta[1], silent = 1, nthread = 2, 
  #objective = rmst_loss, eval_metric = evalerror3,verbose = 2,lambda=hyper_grid$lambda[1],base_score=0,min_child_weight=0,colsample_bytree=lambda=hyper_grid$colsample_bytree[1])
  model <- xgb.train(param, dtrain, nrounds = hyper_grid$optimal_trees[1],watchlist)
  #xgb.dump(model = model, with_stats = TRUE)
  print('model fitting finished')
  print(model)
  print('Predict using new X')
  f <- predict(model, as.matrix(dat.test))
  boost.pred<-1/(1+exp(-f))
  #error <- model$evaluation_log
  #plot(error$iter,-1*error$train_RMST_error,xlab='iter',ylab='obj')
  
  print('Extract importance ranking')
  splt.dtl <- xgb.dump(model = model, with_stats = TRUE)
  if(length(splt.dtl)<=2){
    s1.rank <- median(1:p)
    s2.rank <- median(1:p)
  } else {
    importance_matrix <- xgb.importance(colnames(dat), model = model)
    #importance_matrix <- importance_matrix[order(importance_matrix$Cover,decreasing = T),]
    #xgb.plot.importance(importance_matrix)
    #xgb.plot.tree(feature_names = colnames(dat), model = model)
    #xgb.plot.multi.trees(model = model, feature_names = sparse_matrix@Dimnames[[2]], features.keep = 3)
    
    if(sum(grepl('s1',importance_matrix$Feature))>0){
      s1.rank  <- which(grepl('s1',importance_matrix$Feature))
    } else {
      s1.rank  <- length(importance_matrix$Feature)+ median(1:(p-length(importance_matrix$Feature)))
    }
    
    if(sum(grepl('s2',importance_matrix$Feature))>0){
      s2.rank  <- which(grepl('s2',importance_matrix$Feature))
    } else {
      s2.rank  <- length(importance_matrix$Feature)+ median(1:(p-length(importance_matrix$Feature)))
    }
  }
  
  print('Compute classification performance')
  boost.acc <- mean(1*( boost.pred>=0.5)==otr)
  boost.sen <- sum(boost.pred>=0.5 & otr==1)/sum(otr==1)
  boost.spe <- sum(boost.pred<0.5 & otr==0)/sum(otr==0)
  
  print(boost.acc)
  
  ## Compute value function using training data ##
  #f.train <- predict(model, as.matrix(dat))
  #pred.train<-1/(1+exp(-f.train))
  dat2$potr <- 1*( boost.pred>=0.5)
  
  print("Calculate Value Function")
  print(sum(dat2$potr==1))
  print(sum(dat2$potr==0))
  
  if( sum(dat2$potr==1)>25 & sum(dat2$potr==0)>25 ){
   
    dat.perform <- subset(dat2,potr==1)
    #print(dim(dat.perform))
    trteff.perform <- rmst2(dat.perform$aval, dat.perform$evnt, dat.perform$trt01p, tau = NULL, covariates = NULL, alpha = 0.05)
    
    dat.nperform <- subset(dat2,potr==0)
    #print(dim(dat.nperform))
    cntleff.nperform <- rmst2(dat.nperform$aval, dat.nperform$evnt, dat.nperform$trt01p, tau = NULL, covariates = NULL, alpha = 0.05)
    
    diff.boost <- (nrow(dat.perform)/nrow(dat2)) * trteff.perform[[5]][1,1] - (nrow(dat.nperform)/nrow(dat2)) * cntleff.nperform[[5]][1,1] 
  } else {
    diff.boost <- 0
  }

  
  return( c(boost.acc,s1.rank,s2.rank,boost.sen,boost.spe,diff.boost)   )
  
}




