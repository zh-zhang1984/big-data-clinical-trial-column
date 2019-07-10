#toy exmple for understanding 
library(rpart);
library(rpart.plot);
dd <- data.frame(BP=c(123,111,98,154,199,
                      101,91,133,116,121),
                 gender=c('m','f','f','m',
                          'm','f','m','f','f','m'),
                 age=c(54,66,23,59,76,33,35,54,21,26),
                 obesity=c(1,0,0,1,1,0,1,0,0,0))
dd
dd$F0 <- mean(dd$BP)
dd$PseudoResid1 <- dd$BP-dd$F0
tree1 <- rpart(PseudoResid1 ~ age+gender+obesity,
               data=dd,method='anova',
               control = rpart.control(maxdepth=2,minsplit = 5))
dd$h1 <- predict(tree1)
dd$gamma1 <- rep(sum(dd$h1*(dd$BP-dd$F0))/sum(dd$h1*dd$h1),10)
dd$F1 <- dd$F0+dd$gamma1*dd$h1
dd$PseudoResid2 <- dd$BP - dd$F1
tree2 <- rpart(PseudoResid2 ~ age+gender+obesity,
               data=dd,method='anova',
               control = rpart.control(maxdepth=2,minsplit = 5))
dd$h2 <- predict(tree2)
dd$gamma2<- rep(sum(dd$h2*(dd$BP-dd$F1))/sum(dd$h2*dd$h2),10);
dd$F2 <- dd$F1+dd$gamma2*dd$h2
dd$PseudoResid3 <- dd$BP - dd$F2
par(mfrow=c(1,2)) 
rpart.plot(tree1,main='h1')
rpart.plot(tree2,main='h2')
#working example
set.seed(123);
sampleSize=1000;
Ncat=3; Ncont=4
ContMean=sample(1:10,Ncont)
ContSD=sample(1:10,Ncont)
CatProb=runif(Ncat,min = 0.3,max = 0.7)
for (ii in 1:Ncont) {
  assign(paste('Xcont',ii,sep = '.'),
          rnorm(sampleSize,mean = ContMean[ii],
                sd = ContSD[ii]));
}
for (ii in 1:Ncat) {
  assign(paste('Xcat',ii,sep = '.'),
         sample(c(0,1),size = sampleSize,
                prob = c(1-CatProb[ii],CatProb[ii]),
                replace = T));
}
linpred <- 5*Xcat.1*Xcont.1-6*Xcat.2*Xcont.2-
  Xcont.3*2-Xcont.4*5*Xcat.3;
summary(linpred)
prob <- exp(linpred)/(1 + exp(linpred))
runis <- runif(sampleSize,0,1)
ytest <- ifelse(runis < prob,1,0)
df <- data.frame(Xcat.1=Xcat.1,Xcat.2=Xcat.2,
                    Xcat.3=Xcat.3,Xcont.1=Xcont.1,
                    Xcont.2=Xcont.2,Xcont.3=Xcont.3,
                    Xcont.4=Xcont.4,Y=as.factor(ytest))
rm(list=setdiff(ls(), "df"))
##logistic regression model
library(caret)
set.seed(123)
index <- createDataPartition(df$Y, p=0.75, list=FALSE)
trainSet <- df[ index,]
testSet <- df[-index,]
xnam <- c(paste("Xcat", 1:3,sep='.'),paste("Xcont", 1:4,sep='.'))
Logitformula <- as.formula(paste("Y ~ ", paste(xnam, collapse= "+")))
LogitMod <- glm(Logitformula,trainSet,family = 'binomial')
logitPred <- predict(LogitMod,newdata = testSet,type = 'response')
LogitROC <- roc(testSet$Y,logitPred)
ci.auc(LogitROC)
##gradient boosting mechine
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)
gbmGrid <-  expand.grid(interaction.depth = c(1,3,5,7), 
                        n.trees = (1:30)*10, 
                        shrinkage = 0.1,
                        n.minobsinnode = 50)
gbmMod <- train(Logitformula,
                 trainSet,method='gbm',verbose = FALSE,
                 trControl = fitControl,
                 tuneGrid = gbmGrid)
trellis.par.set(caretTheme())
plot(gbmMod);
gbmPred <- predict(gbmMod,newdata=testSet,
                   type='prob')
gbmROC <- roc(testSet[,'Y'],gbmPred[,2])
auc(gbmROC)
ci.auc(gbmROC)
##explain how the GBM works
library(gbm)
varImp(gbmMod)
summary(gbmMod$finalModel)
#Marginal plots of fitted gbm objects
plot(gbmMod$finalModel,i.var = c(5,4,1),type = 'response')
#Estimate The Strength Of Interaction Effects
interact.gbm(gbmMod$finalModel,
             data=testSet,
             i.var = c(2,5))
#GBM Performance\
gbm.perf(gbmMod$finalModel, 
         plot.it = TRUE, 
         oobag.curve = FALSE, 
         overlay = TRUE, 
         method='OOB')
#Calibration Plot
calibrate.plot(testSet$Y,gbmPred[,2],distribution="bernoulli")
gbm.roc.area(testSet$Y,gbmPred[,2])
