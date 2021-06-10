rm(list=ls())
setwd("/SFS/scratch/zhapingy/sim32/code/")
source("MySimData_sim32.R")
source("GBTE.R")
source("Myboost.R")
source("OneStep.R")
source("ROWSi.R")
source("method_2stepRF.1.0_ez.R")

### Case 1: Interaction: NO; Prognostic coefficient: 0.2;
MySim <- function(rep){
  print(paste('Iteration: ',rep,sep=''))
  data <- MySimData(N=5500,M=500,enrollment=12,follow=18,lambda0=0.1,p=50,rho=1/3,alpha0=0.4,eta=-log(.9)/12,mu0=sqrt(6))
  #data <- SimData_Lu(N=5500,m=500,enrollment=12,follow=18,lambda0=0.1,p=50,rho=1/3,beta=c(1/sqrt(3),rep(  1/(2*sqrt(3)),(10) ),rep(0,40)),gamma=c(0.4,rep(c(0.8,-0.8),1),rep(0,48)),alpha12=0.8,alpha0=sqrt(2),eta=-log(.9)/12)
 #data <- SimData_Lu(N=5500,m=500,enrollment=12,follow=18,lambda0=0.1,p=100,rho=1/3,beta=c(1/sqrt(3),rep(0,10),rep(  1/(2*sqrt(3)),10 ),rep(0,80)),gamma=c(0.4,rep(c(0.8,-0.8),1),rep(0,98)),alpha12=0.8,alpha0=sqrt(2),eta=-log(.9)/12)
  dat.train <- data[[1]]
  dat.test <- data[[2]]
  otr.train <- data[[3]]
  otr.test <- data[[4]]
  # mean(otr.train==1)
  #datt <- cbind(dat.train,otr=otr.train)
  #datt <- cbind(dat.test,otr=otr.test)
  #km_subgroup(datt,'otr')

  dat.rf.train <- dat.train
  dat.rf.test <- dat.test
  dat.rf.train$trt01p <- as.factor(ifelse(dat.rf.train$trt01p==1,'Pembro','Control'))
  dat.rf.test$trt01p <- as.factor(ifelse(dat.rf.test$trt01p==1,'Pembro','Control'))
  
  ## GBTE ##
  
  res.GBT <- GBTE(y=dat.train$aval,  X=dat.train[,c(4:54)],  Tr=dat.train$trt01p,  newX=dat.test[,c(4:54)],  ntr=1000,response="S",Event=dat.train$evnt,th=20,otr=otr.test,dat.test=dat.test)
  
  ## Boosting with proposed value function##
  res.boost.rmst1 <-MyBoost(dat=dat.train,otr=otr.test,dat.test=dat.test,tc=20)
  
  ## Boosting with ITR value function##
  #res.boost.rmst2 <-MyBoost2(dat=dat.train,otr=otr.test,dat.test=dat.test)
  
  trainX<-dat.train[,-c(1,2,3)]
  testX<-dat.test[,-c(1,2,3)]
  treatment <- dat.train$trt01p
  Y<- dat.train$aval
  Y.test<- dat.test$aval
  treatment.test <- dat.test$trt01p
  evnt.test <- dat.test$evnt
  print('ROWSi')
  ROWSi.res <- ROWSi(trainX=trainX,treatment=treatment,Y=Y, grp.lasso.lbd=NULL, testX=testX,testOTR=otr.test,Y.test=Y.test,treatment.test=treatment.test,evnt.test=evnt.test)
  
  
  if('s2' %in% colnames(dat.train)){
    #res.1step <- OneStep(dat=dat.train,p=52,otr=otr.test,dat.test=dat.test,type='regular')
    res.1step.lasso <- OneStep(dat=dat.train,p=52,otr=otr.test,dat.test=dat.test,type='lasso')
    #res.2step.rf.vt0 <- TwoStepRF.wd(dat=dat.rf.train,n.topVar=NULL,var.label=NULL, step1method="vt",mtry=3, vt.version="vt.0", xvar.wt=c( 200,rep(1,52)   ), otr=otr.test ,dat.test=dat.rf.test,wt2nd=TRUE)
    #res.2step.rf.vt1 <- TwoStepRF.wd(dat=dat.rf.train,n.topVar=NULL,var.label=NULL, step1method="vt",mtry=3, vt.version="vt.1", xvar.wt=c( 200,rep(1,52)   ), otr=otr.test ,dat.test=dat.rf.test,wt2nd=TRUE)
    res.2step.rf.vt2 <- TwoStepRF.wd(dat=dat.rf.train,n.topVar=NULL,var.label=NULL, step1method="vt",mtry=3, vt.version="vt.1", xvar.wt=c( 200,rep(1,52)   ), otr=otr.test ,dat.test=dat.rf.test,wt2nd=FALSE)
    
    
    #print(c(res.boost.rmst1,res.boost.rmst2,res.1step.lasso,res.2step.rf.vt0,res.2step.rf.vt1,res.2step.rf.vt2,ROWSi.res))     
    #return( c(res.boost.rmst1,res.boost.rmst2,res.1step.lasso,res.2step.rf.vt0,res.2step.rf.vt1,res.2step.rf.vt2,ROWSi.res)  )
    
    print(c(res.boost.rmst1,res.GBT,res.1step.lasso,res.2step.rf.vt2,ROWSi.res))
    return( c(res.boost.rmst1,res.GBT,res.1step.lasso,res.2step.rf.vt2,ROWSi.res) )
  } else {
    #res.1step <- OneStep(dat=dat.train,p=51,otr=otr.test,dat.test=dat.test,type='regular')
    print('one step')
    res.1step.lasso <- OneStep(dat=dat.train,p=51,otr=otr.test,dat.test=dat.test,type='lasso')
    # print('VT0')
    # res.2step.rf.vt0 <- TwoStepRF.wd(dat=dat.rf.train,n.topVar=NULL,var.label=NULL, step1method="vt",mtry=3, vt.version="vt.0", xvar.wt=c( 200,rep(1,51)   ), otr=otr.test ,dat.test=dat.rf.test,wt2nd=TRUE)
    # print('VT1')
    # res.2step.rf.vt1 <- TwoStepRF.wd(dat=dat.rf.train,n.topVar=NULL,var.label=NULL, step1method="vt",mtry=3, vt.version="vt.1", xvar.wt=c( 200,rep(1,51)   ), otr=otr.test ,dat.test=dat.rf.test,wt2nd=TRUE)
    #print('VT2')
    res.2step.rf.vt2 <- TwoStepRF.wd(dat=dat.rf.train,n.topVar=NULL,var.label=NULL, step1method="vt",mtry=3, vt.version="vt.1", xvar.wt=c( 200,rep(1,51)   ), otr=otr.test ,dat.test=dat.rf.test,wt2nd=FALSE)
    
    
    
    #print(c(res.boost.rmst1,res.boost.rmst2,res.1step.lasso,res.2step.rf.vt0,res.2step.rf.vt1,res.2step.rf.vt2,ROWSi.res))     
    #return( c(res.boost.rmst1,res.boost.rmst2, res.1step.lasso,res.2step.rf.vt0,res.2step.rf.vt1,res.2step.rf.vt2,ROWSi.res)  )
    
    print(c(res.boost.rmst1,res.GBT,res.1step.lasso,res.2step.rf.vt2,ROWSi.res))
    return( c(res.boost.rmst1,res.GBT,res.1step.lasso,res.2step.rf.vt2,ROWSi.res) )
    
  }
  
}


task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
set.seed(task_id)

res<-sapply(1:2,MySim )

setwd("/SFS/scratch/zhapingy/sim32/res")
save.image(file = paste0(task_id,".RData"))


#save(res, file = "C:/Users/zhapingy/OneDrive - Merck Sharp & Dohme, Corp/Documents/tmp/study/intern/res/res_my_sim3.Rdata")
#load("C:/Users/zhapingy/OneDrive - Merck Sharp & Dohme, Corp/Documents/tmp/study/intern/res/res_Lu_multiply_3way.Rdata")

