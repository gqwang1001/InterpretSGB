library(MASS)
library(dplyr)
library(survminer)
library(survival)
#library(rattle)
#library(purrr)
#----------------------------------------------------------------------#
#                     Step 0: Some functions                   #
#----------------------------------------------------------------------#

# km_subgroup <- function(survdata,subgroup)
# {
#   N<- dim(survdata)[1]
#   cutoff<-ceiling(N*0.05)
#   group_os<- survdata  %>% filter(!!sym(subgroup)!="",!is.na(!!sym(subgroup)),!!sym(subgroup)!="Not Reported",!!sym(subgroup)!="Missing") %>% group_by_(subgroup) %>% filter(n()>cutoff)
#   group_os %>% group_split -> split.os
#   group_os %>% group_keys -> key.os
#   
#   
#   split.os %>% map(~surv_fit(Surv(aval, evnt) ~ factor(trt01p), data = .x)) -> fit.os
#   names(fit.os)<-pull(key.os,subgroup)
#   
#   plotlist.os<- ggsurvplot_list(fit.os,split.os,risk.table = TRUE, palette = c("red", "black"), break.time.by = 5, 
#                                 surv.media.line = "hv", legend.title = "", xlab = "Month", 
#                                 ylab = "", legend.labs = NULL, tables.theme = theme_cleantable(), 
#                                 fontsize = 3) 
#   
#   arrange_ggsurvplots(plotlist.os,ncol=length(split.os))  
# 
#   
# }
# 
# 
# # plot km by treatment group 
# km.trt01p <- function(dat, maintxt=NULL, print.sum = T){
#   os.rf <-surv_fit(Surv(aval, evnt) ~ factor(trt01p), data = dat)
#   cox<-coxph(Surv(aval,evnt)~factor(trt01p),dat)
#   hr=round(summary(cox)$conf.int[-2],2)
#   summstat<-surv_median(os.rf)
#   if(print.sum)
#     
#     os.all<-surv_fit(Surv(aval, evnt) ~ factor(trt01p),data = dat)
#   p <- ggsurvplot(os.all, data = ,risk.table=TRUE,palette=c("red","black"),break.time.by=5,surv.media.line="hv",tables.colr=FALSE,legend.title="",xlab="Month",legend.labs=NULL,tables.theme=theme_cleantable(),tables.height=0.25,fontsize=3.5)+guides(col = guide_legend(keywidth=unit(2,"cm"))) 
#   return(p)
# }

#----------------------------------------------------------------------#
#                     Step 1: Generate Survival Data                   #
#----------------------------------------------------------------------#


MySimData <-function(N,M,enrollment,follow,lambda0,p,rho,alpha0,eta,mu0){
  require(MASS)
  trt01p <-rbinom(n=N, 1, prob=0.5)
  #trt01p[trt01p==0]<- -1
  dat <- data.frame(id=1:N,trt01p)
  dat$stdt <- runif(n=N, min=0, max=enrollment)
  
  ## Simulate Z ##
  cov.z  <- (1-rho)*diag(rep(1,(p+2) )) + rho*rep(1,  (p+2) )%*%t(rep(1,  (p+2) ))
  Z <- mvrnorm(n = N, mu=rep(0, (p+2)  ), Sigma=cov.z, tol = 1e-6)
  
  s1 <- Z[,1]
  s2 <- Z[,2]
  Z <- Z[,-c(1,2)]
  
  ## Pattern 3 ##
  sub1 <- (s1<=s2)
  sub2 <- (s1>s2)
  
  otr <- 1*(sub2)

  dat$tte<-NA

  dat$tte <- exp( mu0 + dat$trt01p*(  s1-s2    ) + alpha0*rnorm(N))


  analysistime <- enrollment+follow
  
  ## Simulate Censoring ##
  dat$cnsr <- rexp(N,eta)
  dat$edt <- dat$stdt + dat$tte
  dat$evnt <- as.integer(dat$tte<=dat$cnsr & dat$edt<=analysistime)
  #dat$evnt <- as.integer(dat$tte<=analysistime)
  dat$aval <- NA
  dat$aval[dat$evnt==1]<- dat$tte[dat$evnt==1]
  dat$aval[dat$evnt==0 & dat$cnsr<=(analysistime-dat$stdt)] <- dat$cnsr[dat$evnt==0 & dat$cnsr<=(analysistime-dat$stdt)]
  dat$aval[dat$evnt==0 & dat$cnsr>(analysistime-dat$stdt)] <- analysistime-dat$stdt[dat$evnt==0 & dat$cnsr>(analysistime-dat$stdt)]
  #dat$aval <- ifelse(dat$evnt==1,dat$tte,analysistime)
  dat.old <- dat
  dat.old$otr <- as.factor(otr)
  dat <- dat[,c('trt01p','evnt','aval')]

  
  if(exists("s2")){
    dat <- cbind(dat,s1,s2,Z)
    colnames(dat)[6:(p+5)] <- c(paste0("z", 1:p))
  } else {
    dat <- cbind(dat,s1,Z)
    colnames(dat)[5:(p+4)] <- c(paste0("z", 1:p))
    
  }
  
  ### Check Results ###
  #fit.surv <- coxph( Surv(time=aval, event=evnt) ~ factor(trt01p), data=dat)
  #summary(fit.surv)
  #km.trt01p(dat)
  #km_subgroup(dat.old,'otr')
  dat$aval[dat$aval==0] <- 0.001
  ## Split dat to training and test ##
  train.dat <- dat[(1:M),]
  test.dat <- dat[(M+1):N,]
  
  ## compute OTR for each subject ##
  #otr.ind <- c(1*(otr<0))
  otr.ind <- otr
  otr.train <- otr.ind[1:M]
  otr.test <- otr.ind[(M+1):N]
  
  
  return( list(train.dat,test.dat,otr.train,otr.test)   )
}


#----------------------------------------------------------------------#
#                     Step 2: Run a test                   #
#----------------------------------------------------------------------#
N<-5500
M <- 500
enrollment<-12
follow <- 18
lambda0<-0.1
p <- 20                                                # p is the # of covariates interacting with Treatment
rho <- 1/3                                             # rho is the correlation parameter
alpha0 <- sqrt(2)
eta <- -log(.9)/12
mu1 <- 1
mu2 <- -1
mu3 <- -1
mu4 <- 1
mu0 <- sqrt(6)
# 
sim62 <-MySimData(N,M,enrollment,follow,lambda0,p,rho,alpha0,eta,mu0)
  

