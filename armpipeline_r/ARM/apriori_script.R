# script to run apriori rule mining from a tranasction database

# from binned file, select, if needed, only subjects included in a particular
# stratus, and drop non-binary columns

S <- logical(nrow(binned))

for (i in 1:nrow(binned)) {
  S[i] <- !is.na(match(binned$UNQID[i],strata$Mild_Int))  #change strata if necessary
}

col <- seq(from = 7, to =65, by = 2)
sample <- binned[S,c(2,col)]

# sample <- binned[,c(2,col)]

sampleLog <- data.frame(lapply(sample, as.logical))

asthmaTrans <- as(sampleLog, "transactions")

rules <- apriori(asthmaTrans, parameter = list(support = 0.01, confidence = 0.01))

# subset rules that generate an event
rulesCases <- subset(rules, subset = rhs %in% "event")

# add odds ratio to quality measures
source('~/Desktop/armpipeline_r/ARM/myOddsRatio.R')
quality(rulesCases) <- cbind(quality(rulesCases), myOddsRatio(rulesCases, asthmaTrans, CI=TRUE, t=0.03))

# add improvement
source('~/Desktop/armpipeline_r/ARM/improvement.R')
quality(rulesCases) <- cbind(quality(rulesCases), improvement(rulesCases))

# subset only significant rules
significant <- subset(rulesCases, subset = significance == TRUE)
significant <- subset(significant, subset = imp == TRUE)

# as data frame
dataRules <- as(significant, "data.frame")

write.csv(dataRules,"002rules_03tolerance_Severe.csv")
