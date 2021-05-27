library(SubgroupBoost)

data("simdata")
names(simdata)
str(simdata$train)
str(simdata$test)

#----- data manipulation-----#
#training data
dat.train = simdata[[1]][,-c(1:5)]
info.train = simdata[[1]][, 1:5]
dtrain <-
  xgb.DMatrix(as.matrix(dat.train), label = info.train$trt01p)

#testing data
dat.test = simdata[[2]][,-c(1:5)]
info.test = simdata[[2]][, 1:5]
dtest <- xgb.DMatrix(as.matrix(dat.test), label = info.test$trt01p)

#----- RMST death -----#

dat1 = data.frame(
  dat.train,
  trt01p = info.train$trt01p,
  evnt = info.train$evnt1,
  aval = info.train$aval1
)
set.seed(123)
model_RMST <- SubgroupBoost.RMST(dat1)

#check variables selected
importance = xgb.importance(model_RMST$feature_names, model_RMST)
importance

xgb.plot.importance(importance)

# prediction --------------------------------------------------------------

# prediction on trainset
preds.train = predict(model_RMST, as.matrix(dat.train))
dat1$preds.train = preds.train
dat1$subGroup_train = ifelse(preds.train>0, "+","-")
library(ggplot2)

p1 <- 
  ggplot(dat1, aes(s1, preds.train))+
  geom_point()
p1
ggsave("Results/simulation_trainSet_top1.png", p1)
p2 <- ggplot(dat1, aes(s1, z18, color = subGroup_train))+
  geom_point()
p2
ggsave("Results/simulation_trainSet_top2.png", p2)


# prediction on test set

preds.test = predict(model_RMST, as.matrix(dat.test))
dat2 = data.frame(dat.test,
                  logOdds=preds.test,
                  subGroup = ifelse(preds.test>0, "+","-"))
library(ggplot2)

p1.test <- 
  ggplot(dat2, aes(s1, logOdds))+
  geom_point()
p1.test
ggsave("Results/simulation_testSet_top1.png", p1.test)

p2.test <- ggplot(dat2, aes(s1, z18, color = subGroup))+
  geom_point(alpha=0.5)
p2.test
ggsave("Results/simulation_testSet_top2.png", p2.test)



