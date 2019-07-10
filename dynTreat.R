
##simulation study for described by Baqun Zhang
#####Baqun Zhang: Estimating optimal treatment regimes from a classification perspective
set.seed(123)
n=500
for (i in 1:3) {
  assign(paste('x1',i,sep = ''),rnorm(n))
}
A1linpred <- exp(-0.1+0.5*x11+0.5*x12)
A1pro <- A1linpred/(1+A1linpred)
A1 <- rbinom(n,1,prob = A1pro)
A1opt <- (x11 > -0.54)*(x12 < 0.54)
x21 <- 0.8*x11 + 0.6*A1 + rnorm(n)
x22 <- 0.8*x12 - 0.6*A1 + rnorm(n)
x23 <- 0.8*x13 + 0.7*A1 + rnorm(n)
A2linpred <- exp(-0.1+0.5*x21+0.5*x22)
A2pro <- A2linpred/(1+A2linpred)
A2 <- rbinom(n,1,prob = A2pro)
A2opt <- (x21 > 0.3)*(x23 < 0.46)
y <- 2 + 0.25*x11 + 0.25*x12 - 
           0.25*x13 - 0.5*(A1-A1opt)^2 -
           (A2-A2opt)^2 + rnorm(n)
dt <- data.frame(x11=x11,x12=x12,x13=x13,
                 x21=x21,x22=x22,x23=x23,
                 A1=A1, A2=A2,
                 y=y)
head(round(dt,2))
#Qlearning
library(DynTxRegime)
# outcome model
moMain <- buildModelObj(model = ~x11+x12+x13+
                          x21+x22+x23,
                        solver.method = 'lm')
moCont <- buildModelObj(model = ~x21+x22+x23,
                        solver.method = 'lm')
#### Second-Stage Analysis
fitSS <- qLearn(moMain = moMain, 
                moCont = moCont,
                data = dt, response = y, 
                txName = 'A2')
#confusion matrix
A2Qlearn <- optTx(fitSS)$optimalTx
table(A2,A2opt)
table(A2Qlearn,A2opt)
#### First-stage Analysis
# outcome model
moMain <- buildModelObj(model = ~x11+x12+x13,
                        solver.method = 'lm')
moCont <- buildModelObj(model = ~x11+x12+x13,
                        solver.method = 'lm')
fitFS <- qLearn(moMain = moMain, moCont = moCont,
                data = dt, response = fitSS, 
                txName = 'A1')
coef(fitFS)
#confusion matrix
A1Qlearn <- optTx(fitFS)$optimalTx
table(A1,A1opt)
table(A1Qlearn,A1opt)
##optimal class
# Define the propensity for treatment model and methods.
moPropen <- buildModelObj(model =  ~ x21+x22, 
                          solver.method = 'glm', 
                          solver.args = list('family'='binomial'),
                          predict.method = 'predict.glm',
                          predict.args = list(type='response'))
library(rpart)
moClass <- buildModelObj(model = ~x21+x22+x23,
                         solver.method = 'rpart',
                         solver.args = list(method="class"),
                         predict.args = list(type='class'))
#### Second-Stage Analysis using IPW
fitSS_IPW <- optimalClass(moPropen = moPropen, 
                          moClass = moClass,
                          data = dt, response = y,  
                          txName = 'A2')
rpart.plot(classif(object = fitSS_IPW))
#confusion matrix
A2IPW <- optTx(fitSS_IPW)$optimalTx
table(A2,A2opt)
table(A2IPW,A2opt)
# outcome model
moMain <- buildModelObj(model = ~x21+x22+x23+x11+x12+x13,
                        solver.method = 'lm')

moCont <- buildModelObj(model = ~x21+x22+x23,
                        solver.method = 'lm')

#### Second-Stage Analysis using AIPW
fitSS_AIPW <- optimalClass(moPropen = moPropen, 
                           moMain = moMain, moCont = moCont,
                           moClass = moClass,
                           data = dt, response = y,  
                           txName = 'A2')
#confusion matrix
A2AIPW <- optTx(fitSS_AIPW)$optimalTx
table(A2,A2opt)
table(A2AIPW,A2opt)
rpart.plot(classif(object = fitSS_AIPW))
#### First-stage Analysis using AIPW
# Define the propensity for treatment model and methods.
moPropen <- buildModelObj(model =  ~ x11+x12, 
                          solver.method = 'glm', 
                          solver.args = list('family'='binomial'),
                          predict.method = 'predict.glm',
                          predict.args = list(type='response'))

# classification model
moClass <- buildModelObj(model = ~x11 +x12+x13,
                         solver.method = 'rpart',
                         solver.args = list(method="class",
                                            control = rpart.control(minsplit = 50)),
                         predict.args = list(type='class'))
# outcome model
moMain <- buildModelObj(model = ~x11+x12+x13,
                        solver.method = 'lm')

moCont <- buildModelObj(model = ~x11+x12+x13,
                        solver.method = 'lm')

fitFS_AIPW <- optimalClass(moPropen = moPropen, 
                           moMain = moMain, moCont = moCont,
                           moClass = moClass,
                           data = dt, response = fitSS_AIPW,  
                           txName = 'A1')
rpart.plot(classif(object = fitFS_AIPW))
#confusion matrix
A1AIPW <- optTx(fitFS_AIPW)$optimalTx
table(A1,A1opt)
table(A1AIPW,A1opt)


