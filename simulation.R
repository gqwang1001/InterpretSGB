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
                  subGroup = ifelse(preds.test>0, "+","-"),
                  subGroups = ifelse(preds.test>0, 1,0))
xnew <- expand.grid(s1 = seq(range(dat2$s1)[1], range(dat2$s1)[2], by=.1),
                    z18 = seq(range(dat2$z18)[1], range(dat2$z18)[2], by=0.1))
xmeans <- matrix(rep(t(colMeans(dat.test)), nrow(xnew)), nrow = nrow(xnew), byrow = T)
colnames(xmeans) = colnames(dat.test)
dat.new = as.data.frame(xmeans)#

dat.new[,c("s1", "z18")] = xnew

preds.new = predict(model_RMST, as.matrix(dat.new))
dat3 = data.frame(dat.new,
                  logOdds=preds.new,
                  subGroup = ifelse(preds.new>0, "+","-"),
                  subGroups = ifelse(preds.new>0, 1,0))
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


p3.test <-
  ggplot(dat3, aes(s1, z18, z=subGroups))+
  geom_contour(breaks = .5, size=2)+
  geom_point(data = dat2, aes(s1, z18, color = subGroup),alpha=0.5)
p3.test
ggsave("Results/simulation_testSet_top2_withDecisionBoundary.png", p3.test)

p4 <- ggplot(dat3, aes(s1, z18)) +
  geom_raster(aes(fill=logOdds)) +
  scale_fill_distiller(palette = "Spectral", direction = -1)+
  geom_point(data = dat2, aes(s1, z18, color = subGroup),alpha=0.5)+
geom_contour(aes(z=subGroups),color='black',breaks = .5, size=2)
  
p4
ggsave("Results/simulation_testSet_top2_raster.png", p4)
