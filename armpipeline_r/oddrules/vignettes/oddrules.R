## ------------------------------------------------------------------------
library(oddrules)


## ------------------------------------------------------------------------
data(Pollution)


## ------------------------------------------------------------------------
binned = bin_quantile(Pollution, col=2:31, top=0.95)


## ------------------------------------------------------------------------
binned2 = bin_quantile(Pollution, col=2:31, top=0.95, group=5)


## ------------------------------------------------------------------------
binned3 = bin_quantile(Pollution, col=2:31, top=0.95, group=5, groupthreshold=1)

## ------------------------------------------------------------------------
binned <- as.data.frame(binned$data)
sample <- binned
sampleLog <- data.frame(lapply(sample, as.logical))


## ------------------------------------------------------------------------
splitdata <- random_split(sampleLog, seed=1, split=5000) 
whole_training <- as.data.frame(splitdata$whole_training)
testing <- as.data.frame(splitdata$testing)


## ------------------------------------------------------------------------
trainingTrans <- as(whole_training, "transactions")
rules <- apriori(trainingTrans, parameter = list(support = 0.005, confidence = 0.0, maxlen = 6))
rules <- subset(rules, subset = rhs %in% "event")

## ------------------------------------------------------------------------
sampleNot <- as.data.frame(!whole_training)
sampleNot$event <- whole_training$event
trainingNot <- as(sampleNot, "transactions")


## ------------------------------------------------------------------------
myOddsRatioExp(rules, trainingTrans, trainingNot, CI=TRUE, t=0.03)



## ------------------------------------------------------------------------
customCI(rules, trainingTrans, trainingNot, alpha = 1.598)


## ------------------------------------------------------------------------
improvement_customCI_par(rules)


## ------------------------------------------------------------------------
sampleNot <- as.data.frame(!testing)
sampleNot$event <- testing$event

testingTrans <- as(testing, "transactions")
testingNot <- as(sampleNot, "transactions")

myOddsRatioExp(rules, testingTrans, testingNot, CI=TRUE, t=0.03)
customCI(rules, testingTrans, testingNot, alpha = 1.598)


## ------------------------------------------------------------------------
overfit_test(whole_training, testing, seed=1,trainingsize=5000)


## ------------------------------------------------------------------------
Results <- overfit_test(whole_training, testing, seed=1:3, trainingsize=c(5000,2000,1000))
Results


## ------------------------------------------------------------------------
variables <- colnames(binned)[-1] 
nonredundant(Results$TestFrame,variables)


