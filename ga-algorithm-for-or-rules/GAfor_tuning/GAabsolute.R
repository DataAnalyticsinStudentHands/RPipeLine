##### GA ALGORITHM FOR RULE MINING - OR ORIENTED #####

# This script implements a genetic algorithm for the automatic
# mining of association rules of the form {exposures} -> event.
# The design is based on the algorithm proposed in Alatas-05.

GA.absolute <- function(data, generations = 200, upop_size = 30, Pin = 0.25, seed = NULL,
                      weights = rep(1,7), pc = 0.6, pc2 = 0.5, k = 4, pm_range = c(0.001, 0.01),
                      pt_range = c(0.9, 0.6)) {

  # ARGUMENTS: 
  # data: data frame to mine (must contain the event information in 
  #       the first column)
  # generations: number of generations for GA
  # upop_size: gene pool size
  # Pin: initial probability of any attribute to be included in a rule
  # seed: seed for random number generation
  # weights: weights for fitness function, in order support, confidence,
  #          OR, repetitions, redundancy
  # pc: crossover probability
  # pc2: proportion of genes from 2 parents (default 50-50)
  # k = tournament size
  # pm_range = probability of mutation (range)
  # pt_range = probability of win in tournament (range)
  # fitness: function used to account for different objective functions
  #          (aritmetic or geometric mean)
  
  
  # Import utility functions
  source('GAfunctions.R')
  source('GAfitness.R')
  
  # First we rename the data frame of data that we are going to
  # use to test the quality of the rules. The first column of this 
  # data frame must contain the event information, while the 
  # remaining columns are exposures.
  test <- data
  
  # Compute number of genes
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
  # If the user defines a seed for random number generation, it is set here:
  if(!is.null(seed))
    set.seed(seed)
  
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
  
  # Create data frame for stats (OR, CI and p-value)
  stats = data.frame(OR = rep(NA,upop_size), low95CI = rep(NA,upop_size), high95CI = rep(NA,upop_size), p = rep(NA,upop_size))
  
  # ant, threshold and fitness are grouped in a list
  Upop <- list(ant = ant, threshold = threshold, fitness = fitness, stats = stats)
  
  # compute fitness of Upop
  Upop <- fitness.fun(Upop, upop_size, genes_num, test, weights, Emedian, maxdist)
  
  # Upop needs to be sorted
  Upop <- pop.sort(Upop)

  # Adjust fitness to avoid rule repetitions and redundancy
  penalty<-1
  red<-1
  for (i in 1:upop_size) {
    Upop$fitness[i] <- Upop$fitness[i] - weights[6]*rep.penalty(Upop,i,upop_size)
    penalty[i]<-weights[6]*rep.penalty(Upop,i,upop_size)
    Upop$fitness[i] <- Upop$fitness[i] - weights[7]*redundancy(Upop,i,upop_size)
    red[i] <-weights[7]*redundancy(Upop,i,upop_size)
  }
  
  # sort again to account for adjusted fitness
  Upop <- pop.sort(Upop)
  
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
    Upop2 <- fitness.fun(Upop2, upop_size, genes_num, test, weights, Emedian, maxdist)
    
    # Merge, sort and preserve best half
    temp <- mapply(rbind, Upop, Upop2, SIMPLIFY=FALSE)
    temp$fitness <- c(Upop$fitness, Upop2$fitness)
    temp <- pop.sort(temp)
    
    # Adjust fitness to avoid rule repetitions
    for (i in 1:(upop_size*2)) {
      temp$fitness[i] <- temp$fitness[i] - weights[6]*rep.penalty(temp,i,upop_size*2)
      temp$fitness[i] <- temp$fitness[i] - weights[7]*redundancy(temp,i,upop_size*2)
    }
    
    # Sort adjusted fitness
    temp <- pop.sort(temp)
    
    Upop$ant <- temp$ant[1:upop_size,]
    Upop$threshold <- temp$threshold[1:upop_size,]
    Upop$fitness <- temp$fitness[1:upop_size]
    Upop$stats <- temp$stats[1:upop_size,]
  }
  return(Upop)
}
  
#### OUTPUT REPORT ####
GAreport <- function(test, Upop, upop_size = 30, weights = rep(1,7), filename = "GAreport.csv") {
  Output <- Upop$ant*Upop$threshold
  # convert logical features to true/false
  Output[is.na(Output)] <- Upop$ant[is.na(Output)] == 1
  Output <- cbind(Output, fitness=Upop$fitness)
  stats <- data.frame(OR = array(), low95CI = array(), high95CI = array(), p = array())
  genes_num <- ncol(Upop$ant)
  
  support <- array()
  coverage <- array()
  confidence <- array()
  length <- array()
  OR <- array()
  extreme <- array()
  penalty <- array()
  red <- array()
  
  LB <- sapply(test[-1], min)
  UB <- sapply(test[-1], max)
  Emedian <- sapply(test[-1], median)
  Udist <- UB - Emedian
  Ldist <- Emedian - LB
  maxdist <- sapply(as.data.frame(rbind(Ldist, Udist)), max)
  
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
    if (identical(redundancy(Upop,i,upop_size), numeric(0)))
      red[i] <- 0
    else red[i] <- redundancy(Upop,i,upop_size)
    
    stats[i,] <- confidence.int(testbin)
  }
  
  Output <- cbind(Rule = (1:upop_size), Output, supp = support, conf = confidence, length = length, OR = OR, distance = extreme, rep = penalty, red = red, stats)
  write.csv(Output, file = filename)
  return(Output)
}