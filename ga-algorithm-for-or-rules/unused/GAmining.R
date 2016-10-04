##### GA ALGORITHM FOR RULE MINING - OR ORIENTED #####

# This script implements a genetic algorithm for the automatic
# mining of association rules of the form {exposures} -> event.
# The design is based on the algorithm proposed in Alatas-05.

#### UTILITY FUNCTIONS ####
odds.ratio <- function(testbin) {
  cases <- testbin[testbin[,1],]
  controls <- testbin[!testbin[,1],]
  x11 <- sum(rowSums(cases)==ncol(cases)) # exposed cases
  x10 <- nrow(cases) - x11 # unexposed cases (general definition)
  x01 <- sum(rowSums(controls)==ncol(controls)-1) # exposed controls
  x00 <- nrow(controls) - x01 # unexposed controls
  return((x11*x00)/(x01*x10))
}

# This function generates random indexes to pick chromosomes to
# use to generate new ones. Low indexes are favored, because they 
# represent higher ranking chromosomes. In this implementation, 
# chromosomes from the top half of the population have 68% of chances
# of being picked, against 32% of the bottom half.
pick.parent <- function(upop_size){
  pick <- rnorm(1, sd = upop_size/3)
  pick <- ceiling(abs(pick))
  if (pick > upop_size)
    pick <- upop_size
  return(pick)
}

# Function to perform uniform crossover of two parents
uni.crossover <- function(Upop, ip1, ip2, pc, pc2){
  c1 <- list(ant = Upop$ant[ip1,], t = Upop$threshold[ip1,])
  c2 <- list(ant = Upop$ant[ip2,], t = Upop$threshold[ip2,])
  
  # First it is decided if crossover is performed or not
  if(runif(1) < pc){
    # If crossover is performed, we go bit by bit
    for(i in 1:length(c1$ant)){
      if(runif(1) < pc2) {
        #swap ant
        temp <- c1$ant[i]
        c1$ant[i] <- c2$ant[i]
        c2$ant[i] <- temp
        #swap thresholds
        temp <- c1$t[i]
        c1$t[i] <- c2$t[i]
        c2$t[i] <- temp
      }
    }  
  } 
  return(children <- list(c1 = c1, c2 = c2))
}

# Mutation function
mutate <- function(cross, pm) {
  for(n in 1:2){
    for(i in 1:length(cross[[n]]$ant)){
      if(runif(1) < pm) {
        if (cross[[n]]$ant[i]==1)
          cross[[n]]$ant[i] <- 0
        else cross[[n]]$ant[i] <- 1
      }
    }
  }
  return(cross)
}

rep.penalty <- function(Upop, i, upop_size) {
  # Define array of penalties (monotonically decreasing)
  P <- seq(from = 1, to = 0, by = -1/upop_size) 
  Ptot <- 0
  
  for (j in 1:(i-1)){
    # check if chromosomes are equal to those to a previous rule
    equal <- sum(abs(Upop$ant[i,] - Upop$ant[j,]))  # equal = 0 if the rules share the same chromosomes
    if (equal == 0) {
      Ptot = Ptot + P[j]
    }
  } 
  # return total penalty
  return(Ptot)
}

#### END UTILITY FUNCTIONS ####


# First we import the data frame of data that we are going to
# use to test the quality of the rules. The first column of this 
# data frame must contain the event information, while the 
# remaining columns are exposures (continuous)

load("~/Desktop/my papers/7 - GA for ARM/embedded_V1_6_05.RData")
test <- X

# Choose gene pool size (number of rules)
upop_size <- 20
genes_num <- ncol(test) - 1

# We need to know the upper and lower boundaries for all the attributes
LB <- sapply(test[-1], min)
UB <- sapply(test[-1], max)

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
  threshold[,i] <- runif(upop_size, LB[i], UB[i])
}

# Create fitnes array
fitness <- array(data = NA, upop_size)

# fitness function coefficients
a1 <- 0.1 #supp
a2 <- 0.1 #conf
a3 <- 0.4 #length
a4 <- 0.1 #OR
a5 <- 0.1 #repetitions

generations <- 200
pc <- 0.6 # probability of crossover
pc2 <- 0.5 # probability of crossover 2: percentage of exchanged genes
pm <- 0.01 # probability of mutation

# ant, threshold and fitness are grouped in a list
Upop <- list(ant = ant, threshold = threshold, fitness = fitness)

# Testing fitness of all rules and storing the value in the list
for (i in 1:upop_size) {
  # first we binarize a copy of test according to the thresholds
  # of the rule under evaluation
  testbin <- test
  for (j in 1:genes_num) {
    testbin[,j+1] <- testbin[,j+1] > Upop$threshold[i,j]
  }
  # then we preserve only the true attributes and the event column
  A <- which(Upop$ant[i,] %in% 1)
  A <- A + 1
  testbin <- testbin[,c(1,A)]
  
  # testbin needs to contain at least one attribute and the event column
  if (!is.numeric(testbin)) {
    testbin[,1] <- as.logical(testbin[,1])
    
    # We compute the required rule quality measures
    support <- sum(rowSums(testbin)==length(A)+1)/nrow(testbin)
    coverage <- sum(rowSums(testbin[-1])==(length(A)))/nrow(testbin)
    confidence <- support/coverage
    length <- length(A)
    OR <- abs(1-odds.ratio(testbin))
    
    Upop$fitness[i] <- a1*support + a2*confidence - a3*length + a4*OR
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

# Crossover and mutation
# Each crossover generates 2 children. We want upop_size new children,
# so we iterate upop_size/2 times
for (i in 1:(upop_size/2)) {
  # pick two parents
  ip1 <- pick.parent(upop_size)
  ip2 <- pick.parent(upop_size)
  
  cross <- uni.crossover(Upop, ip1, ip2, pc, pc2)
  cross <- mutate(cross, pm) # calling mutation function on 2 children
  
  if(is.null(cross))
    browser()
  
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
    testbin[,j+1] <- testbin[,j+1] > Upop2$threshold[i,j]
  }
  # then we preserve only the true attributes and the event column
  A <- which(Upop2$ant[i,] %in% 1)
  A <- A + 1
  testbin <- testbin[,c(1,A)]
  
  # testbin needs to contain at least one attribute and the event column
  if (!is.numeric(testbin)) {
    testbin[,1] <- as.logical(testbin[,1])
    
    # We compute the required rule quality measures
    support <- sum(rowSums(testbin)==length(A)+1)/nrow(testbin)
    coverage <- sum(rowSums(testbin[-1])==(length(A)))/nrow(testbin)
    confidence <- support/coverage
    length <- length(A)
    OR <- abs(1-odds.ratio(testbin))
    
    Upop2$fitness[i] <- a1*support + a2*confidence - a3*length + a4*OR
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
for (i in 2:(upop_size*2)) {
  if (rep.penalty(temp,i,upop_size*2) > 0)
  temp$fitness[i] <- temp$fitness[i] - a5*rep.penalty(temp,i,upop_size*2)
}

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
}