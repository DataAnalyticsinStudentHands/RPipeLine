##### GA ALGORITHM FOR RULE MINING - OR ORIENTED #####

# This script implements a genetic algorithm for the automatic
# mining of association rules of the form {exposures} -> event.
# The design is based on the algorithm proposed in Alatas-05.

# In the GAmulticlass version, features can be of class numeric,
# integer and boolean.

# Import utility functions
source('~/Desktop/my papers/7 - GA for ARM/R files/GAfunctions.R')

# First we import the data frame of data that we are going to
# use to test the quality of the rules. The first column of this 
# data frame must contain the event information, while the 
# remaining columns are exposures (continuous)

load("~/Desktop/my papers/7 - GA for ARM/embedded_V1V10_6TRUE_05.RData")
test <- X

# Choose gene pool size (number of rules)
upop_size <- 30
genes_num <- ncol(test) - 1

# We need to know the upper and lower boundaries for all the attributes
# Also computing median and distances for later use
LB <- sapply(test[-1], min)
UB <- sapply(test[-1], max)
Emedian <- sapply(test[-1], median)
Udist <- UB - Emedian
Ldist <- Emedian - LB
maxdist <- sapply(as.data.frame(rbind(Ldist, Udist)), max)

for (i in 1:ncol(test[-1])) {
  if (is.logical(test[,i+1])) {
    LB[i] <- NA
    UB[i] <- NA
    Emedian[i] <- NA
    maxdist[i] <- NA
  }
}

# Create matrix ant, which determines if an attribute is included in 
# each rule (0/1)
Pin <- 0.3 # initial probability of any attribute to be included in a rule
set.seed(3)

ant <- matrix(data = NA, nrow = upop_size, ncol = genes_num)
for (i in 1:genes_num) {
  ant[,i] <- as.numeric(runif(upop_size) < Pin)
}

# Create matrix threshold, which stores the exposure threshold for each attribute
threshold <- matrix(data = NA, nrow = upop_size, ncol = genes_num)

for (i in 1:genes_num) {
  if (!is.na(LB[i])) {
    threshold[,i] <- runif(upop_size, LB[i], UB[i])
  }
  else threshold[,i] <- NaN
}

# Create fitnes array
fitness <- array(data = NA, upop_size)

# fitness function coefficients
a1 <- 1 #supp
a2 <- 1/0.2 #conf
a3 <- 1/0.5 #length
a4 <- 1/5 #OR
a5 <- 1 #repetitions
a6 <- 1 # distance from extreme values

generations <- 200
pc <- 0.6 # probability of crossover
pc2 <- 0.5 # probability of crossover 2: percentage of exchanged genes
pm_range <- c(0.001, 0.01) # probability of mutation (range)
k <- 4
pt_range <- c(0.9, 0.6) # probability of a strong chromosome to win a tournament (range)
test_metrics <- data.frame(pm=array(), pc=array(), sim=array(), bestfit=array())

# ant, threshold and fitness are grouped in a list
Upop <- list(ant = ant, threshold = threshold, fitness = fitness)

# Testing fitness of all rules and storing the value in the list
for (i in 1:upop_size) {
  # first we binarize a copy of test according to the thresholds
  # of the rule under evaluation
  testbin <- test
  for (j in 1:genes_num) {
    if (!is.logical(testbin[,j+1])) { # logical features do not need to be converted
      testbin[,j+1] <- testbin[,j+1] > Upop$threshold[i,j]
    }
  }
  # then we preserve only the true attributes and the event column
  A <- which(Upop$ant[i,] %in% 1)
  A <- A + 1
  testbin <- testbin[,c(1,A)]
  
  # testbin needs to contain at least one attribute and the event column
  if (is.null(ncol(testbin))) 
    Upop$fitness[i] = NaN
  else
  {
    testbin[,1] <- as.logical(testbin[,1])
    
    # We compute the required rule quality measures
    support <- sum(rowSums(testbin)==length(A)+1)/nrow(testbin)
    coverage <- sum(rowSums(testbin[-1])==(length(A)))/nrow(testbin)
    confidence <- support/coverage
    length <- length(A)
    OR <- OR.fitness(testbin)
    extreme <- median.distance(Upop$threshold[i,(A-1)], Emedian[A-1], maxdist[A-1])
    
    Upop$fitness[i] <- a1*support + a2*confidence - a3*length + a4*OR - a6*extreme
    # some rules can end up having a "non numerical" fitness, such as
    # NaN or Inf. NaN are harmless, but Inf could bias the algorithm
    # and should be eliminated
    if (!is.na(Upop$fitness[i]) && !is.nan(Upop$fitness[i]))
      if (Upop$fitness[i] == Inf)
        Upop$fitness[i] = NaN
  }
}

# Upop needs to be sorted
Upop_rank <- order(Upop$fitness, decreasing=TRUE)
Upop$fitness <- Upop$fitness[Upop_rank]
Upop$ant <- Upop$ant[Upop_rank,]
Upop$threshold <- Upop$threshold[Upop_rank,]

# Adjust fitness to avoid rule repetitions
for (i in 2:upop_size) {
  Upop$fitness[i] <- Upop$fitness[i] - a5*rep.penalty(Upop,i,upop_size)
}

# sort again to account for adjusted fitness
Upop_rank <- order(Upop$fitness, decreasing=TRUE)
Upop$fitness <- Upop$fitness[Upop_rank]
Upop$ant <- Upop$ant[Upop_rank,]
Upop$threshold <- Upop$threshold[Upop_rank,]

# Loop to create and evaluate new chromosomes
for (n in 1:generations) {
  # Duplicate Upop (Upop2). This only serves the purpose of creating a list
  # of similar shape. Its content is going to be changed to host the new
  # chromosomes generated from Upop
  Upop2 <- Upop
  
  # Compute rules similarity to adjust pm and pt
  # currently based only on antecedents, not thresholds
  sim <- pop.similarity(Upop$ant)
  pt <- sim.adjust(pt_range, sim)
  pm <- sim.adjust(pm_range, sim)
  
  # Crossover and mutation
  # Each crossover generates 2 children. We want upop_size new children,
  # so we iterate upop_size/2 times
  for (i in 1:(upop_size/2)) {
    # pick two parents
    ip1 <- tournament(k, upop_size, pt)
    ip2 <- tournament(k, upop_size, pt)
    
    cross <- uni.crossover(Upop, ip1, ip2, pc, pc2)
    
    cross <- mutate(cross, pm, Emedian) # calling mutation function on 2 children
    
    Upop2$ant[i,] <- cross$c1$ant
    Upop2$threshold[i,] <- cross$c1$t
    
    Upop2$ant[i+(upop_size/2),] <- cross$c2$ant
    Upop2$threshold[i+(upop_size/2),] <- cross$c2$t  
  }
  
  
  # Measuring fitness of new population
  for (i in 1:upop_size) {
    # first we binarize a copy of test according to the thresholds
    # of the rule under evaluation
    testbin <- test
    for (j in 1:genes_num) {
      if (!is.logical(testbin[,j+1])) { # logical features do not need to be converted
        testbin[,j+1] <- testbin[,j+1] > Upop2$threshold[i,j]
      }
    }
    # then we preserve only the true attributes and the event column
    A <- which(Upop2$ant[i,] %in% 1)
    A <- A + 1
    testbin <- testbin[,c(1,A)]
    
    # testbin needs to contain at least one attribute and the event column
    if (is.null(ncol(testbin))) 
    {
      Upop2$fitness[i] = NaN}
    else
    {
      testbin[,1] <- as.logical(testbin[,1])
      
      # We compute the required rule quality measures
      support <- sum(rowSums(testbin)==length(A)+1)/nrow(testbin)
      coverage <- sum(rowSums(testbin[-1])==(length(A)))/nrow(testbin)
      confidence <- support/coverage
      length <- length(A)
      OR <- OR.fitness(testbin)
      extreme <- median.distance(Upop2$threshold[i,(A-1)], Emedian[A-1], maxdist[A-1])
      
      Upop2$fitness[i] <- a1*support + a2*confidence - a3*length + a4*OR - a6*extreme
      # some rules can end up having a "non numerical" fitness, such as
      # NaN or Inf. NaN are harmless, but Inf could bias the algorithm
      # and should be eliminated
      if (!is.na(Upop2$fitness[i]) && !is.nan(Upop2$fitness[i]))
        if (Upop2$fitness[i] == Inf)
          Upop2$fitness[i] = NaN
    }
  }
  
  # Merge, sort and preserve best half
  temp <- mapply(rbind, Upop, Upop2, SIMPLIFY=FALSE)
  temp$fitness <- c(Upop$fitness, Upop2$fitness)
  
  temp_rank <- order(temp$fitness, decreasing=TRUE)
  temp$fitness <- temp$fitness[temp_rank]
  temp$ant <- temp$ant[temp_rank,]
  temp$threshold <- temp$threshold[temp_rank,]
  
  # Adjust fitness to avoid rule repetitions
  for (i in 2:(upop_size*2)) 
    temp$fitness[i] <- temp$fitness[i] - a5*rep.penalty(temp,i,upop_size*2)
  
  #if(sum(temp$ant[,1])>0)
  #  browser()
  
  # Sort adjusted fitness
  temp_rank <- order(temp$fitness, decreasing=TRUE)
  temp$fitness <- temp$fitness[temp_rank]
  temp$ant <- temp$ant[temp_rank,]
  temp$threshold <- temp$threshold[temp_rank,]
  
  Upop$ant <- temp$ant[1:upop_size,]
  Upop$threshold <- temp$threshold[1:upop_size,]
  Upop$fitness <- temp$fitness[1:upop_size]
  
  test_metrics <- rbind(test_metrics, c(pm, pt, sim, Upop$fitness[1]))
}

#### OUTPUT REPORT ####
Output <- Upop$ant*Upop$threshold
# convert logical features to true/false
Output[is.na(Output)] <- Upop$ant[is.na(Output)] == 1
Output <- cbind(Output, fitness=Upop$fitness)
penalty <- array(NA, upop_size)
stats <- data.frame(OR = array(), low95CI = array(), high95CI = array(), p = array())

for (i in 1:upop_size) {
  # first we binarize a copy of test according to the thresholds
  # of the rule under evaluation
  testbin <- test
  for (j in 1:genes_num) {
    if (!is.logical(testbin[,j+1])) { # logical features do not need to be converted
      testbin[,j+1] <- testbin[,j+1] > Upop$threshold[i,j]
    }
  }
  # then we preserve only the true attributes and the event column
  A <- which(Upop$ant[i,] %in% 1)
  A <- A + 1
  testbin <- testbin[,c(1,A)]
  
  testbin[,1] <- as.logical(testbin[,1])
  
  # We compute the required rule quality measures
  support[i] <- sum(rowSums(testbin)==length(A)+1)/nrow(testbin)
  coverage[i] <- sum(rowSums(testbin[-1])==(length(A)))/nrow(testbin)
  confidence[i] <- support[i]/coverage[i]
  length[i] <- length(A)
  OR[i] <- OR.fitness(testbin)
  extreme[i] <- median.distance(Upop$threshold[i,(A-1)], Emedian[A-1], maxdist[A-1])
  if (identical(rep.penalty(Upop,i,upop_size), numeric(0)))
    penalty[i] <- 0
  else penalty[i] <- rep.penalty(Upop,i,upop_size)
  
  stats[i,] <- confidence.int(testbin)
}

Output <- cbind(Rule = (1:upop_size), Output, a1*support, a2*confidence, a3*length, a4*OR, a6*extreme, a5*penalty, stats)

require(ggplot2)
require(reshape)

col <- c(1,12:18)
df <- melt(Output[,col], id.vars = 'Rule', variable.name = 'metrics')

# plot on same grid, each series colored differently 
ggplot(df, aes(x=Rule,y=value,group=variable, colour=variable)) + geom_line() + geom_point()

write.csv(Output, file = "GAresults_normalWeight_newPenalty.csv")
