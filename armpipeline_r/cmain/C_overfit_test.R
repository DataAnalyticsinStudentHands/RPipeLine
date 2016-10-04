overfit_test <- function (whole_training, testing, seed=1:3, trainingsize=c(5000,2000,1000), alpha = 1.598, support = 0.005) {
  
  #binned, splitseed=1, whole_training_size=5000
  # This function mines for rules and test for overfitting and sample sixe
  # whole_training: must already be binned, logical and split
  # testing: must already be binned, logical, and split
  # seed: what seeds to use for randomizing testing and training set
  # trainingsize: vector with all training sizes
  
  library("arules")
  
  A <- matrix(, ncol = length(seed), nrow = length(trainingsize),dimnames = list(trainingsize, seed))
  B <- matrix(, ncol = length(seed), nrow = length(trainingsize),dimnames = list(trainingsize, seed))
  C <- matrix(, ncol = length(seed), nrow = length(trainingsize),dimnames = list(trainingsize, seed))  
  TestFrame=data.frame()

  for (i in 1:length(seed)) {
    for (j in 1:length(trainingsize)){
      
      #sampleLog <- data.frame(lapply(sample, as.logical))
      #splitdata <- random_split(sampleLog,seed=splitseed,split=whole_training_size)
      #whole_training <- data.frame(lapply(splitdata$whole_training, as.logical))
      #testing <- data.frame(lapply(splitdata$testing, as.logical))
      
      # Training set can be further reduced in size, if necessary
      set.seed(seed[i])
      tr_size = trainingsize[j]
      training <- whole_training[sample(nrow(whole_training), tr_size, replace=TRUE), ]
      
      # # Expanding data with negated features
      sampleNot <- as.data.frame(!training)
      # newcols <- colnames(sampleNot[,-1])
      # newcols <- paste("Not",newcols,sep="_")
      # colnames(sampleNot) <- c("event",newcols)
      # sampleLog <- cbind(sampleLog,sampleNot[,-1])
      # sampleNot <- cbind(sampleNot,sampleNot[,-1])
      sampleNot$event <- training$event
      # colnames(sampleNot) <- colnames(sampleLog)
      if (sum(training$event) < 1){
        A[j,i]=0
        B[j,i]=0
        C[j,i]=0
      }
      else{
        trainingTrans <- as(training, "transactions")
        trainingNot <- as(sampleNot, "transactions")
        
        rules <- apriori(trainingTrans, parameter = list(support = support, confidence = 0.0, maxlen = 6))
        
        # subset rules that generate an event
        rulesCases <- subset(rules, subset=rhs %in% "event")
        
        # # If rules is too large, use
        # rr <- rhs(rules) %in% "event"
        # rulesCases <- rules[rr]
        
        # as data frame
        dataRules <- as(rulesCases, "data.frame")
        #create vector of number of rules
        A[j,i]<-nrow(dataRules)
        
        # add odds ratio to quality measures
        quality(rulesCases) <- cbind(quality(rulesCases), myOddsRatioExp(rulesCases, trainingTrans, trainingNot, CI=TRUE, t=0.03))
        
        # subset only significant rules
        significant <- subset(rulesCases, subset = oddsRatio != Inf)
        significant <- subset(significant, subset = significance == TRUE)
        signRules <- as(significant, "data.frame")
        
        if (nrow(signRules)<1) {
          B[j,i]=0
          C[j,i]=0
          
        }
        
        else {
          B[j,i]<-nrow(signRules)
          
          # add custom CI for method 4 (10%)
          quality(significant) <- cbind(quality(significant), customCI(significant, trainingTrans, trainingNot, alpha = alpha))
          
          quality(significant) <- cbind(quality(significant), improvement_customCI_par(significant)) #####YOU! YOU CAUSE MY CONNECTIONS ERRORS!
          
          significant <- subset(significant, subset = impCIcustom_par == TRUE)
          
          signRules <- as(significant, "data.frame")
          
          C[j,i]<-nrow(signRules)
          
          
          
          
          # What has been found in the training set, now has to be evaluated in the 
          # testing set
          if (C[j,i]>0){
            sampleNot <- as.data.frame(!testing)
            sampleNot$event <- testing$event
            
            testingTrans <- as(testing, "transactions")
            testingNot <- as(sampleNot, "transactions")
            
            ruleTest <- significant
            ruleTest@quality$support <- NULL
            ruleTest@quality$confidence <- NULL
            ruleTest@quality$lift <- NULL
            ruleTest@quality$oddsRatio <- NULL
            ruleTest@quality$standardError <- NULL
            ruleTest@quality$lowCI95 <- NULL
            ruleTest@quality$highCI95 <- NULL
            ruleTest@quality$significance <- NULL
            ruleTest@quality$chi_squared <- NULL
            ruleTest@quality$p_value <- NULL
            ruleTest@quality$lowCIcustom <- NULL
            ruleTest@quality$highCIcustom <- NULL
            ruleTest@quality$impCIcustom_par <- NULL
            if (C[j,i]==1){
              quality(ruleTest) <- cbind(quality(ruleTest), t(interestMeasure(ruleTest, c("support", "confidence", "lift"), transactions = testingTrans)))
            }
            if (C[j,i]>1){
              quality(ruleTest) <- cbind(quality(ruleTest), interestMeasure(ruleTest, c("support", "confidence", "lift"), transactions = testingTrans))
            }
            
            quality(ruleTest) <- cbind(quality(ruleTest), myOddsRatioExp(ruleTest, testingTrans, testingNot, CI=TRUE, t=0.03))
            quality(ruleTest) <- cbind(quality(ruleTest), customCI(ruleTest, testingTrans, testingNot, alpha = alpha))
            quality(ruleTest) <- cbind(quality(ruleTest), improvement_customCI_par(ruleTest)) ###and you you evil improvement customCI fiend!

            ruleTest=as(ruleTest,"data.frame")
            TestFrame=rbind(TestFrame,ruleTest)
            
            
          }
          
          
        }      
        
      }    
    }    
  }

  return(list(Support_Only=A,With_OR_Significance_Filter=B,With_Non_Redundancy_Filter=C,TestFrame=TestFrame))
  
}
