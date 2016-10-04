# Script of some functions that did not make the final cut

# rep.penalty with partial overlap
rep.penalty.part <- function(Upop, i, upop_size) {
  # Define array of penalties (monotonically decreasing)
  P <- seq(from = 1, to = 0, by = -1/(upop_size-1)) 
  Ptot <- 0
  
  for (j in 1:(i-1)){
    # check if some genes are equal to those to a previous rule:
    # add chromosomes together
    equal <- Upop$ant[i,] + Upop$ant[j,]  
    # if a gene appears twice it adds up to two. Count 2s in equal
    reps <- sum(match(equal,2), na.rm=TRUE)
    # count total true genes in tested rule
    genetot <- sum(Upop$ant[i,])
    # fraction of genes repeated from a previous rule (range 0-1)
    if(genetot > 0){
      repf <- reps/genetot
      Ptot <- Ptot + P[j] * repf   
    }
  } 
  # return total penalty
  return(Ptot)
}

# fitness fun with relative control of weights
fitness.funR <- function(Upop, upop_size, genes_num, test, weights, Emedian, maxdist) {
  # Testing fitness of all rules and storing the value in the list
  max_weights <- array(0, length(weights))
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
      if(!is.na(support))
        if(support > max_weights[1])
          max_weights[1] <- support
      coverage <- sum(rowSums(testbin[-1])==(length(A)))/nrow(testbin)
      confidence <- support/coverage
      if(!is.na(confidence))
        if(confidence > max_weights[2])
          max_weights[2] <- confidence
      length <- length(A)
      if(!is.na(length))
        if(length > max_weights[3])
          max_weights[3] <- length
      OR <- OR.fitness(testbin)
      if(!is.na(OR))
        if(support > max_weights[4])
          max_weights[4] <- OR
      extreme <- median.distance(Upop$threshold[i,(A-1)], Emedian[A-1], maxdist[A-1])
      if(!is.na(extreme))
        if(support > max_weights[5])
          max_weights[5] <- extreme
      
      Upop$fitness[i] <- weights[1]*support + weights[2]*confidence - weights[3]*length + weights[4]*OR - weights[5]*extreme
      # some rules can end up having a "non numerical" fitness, such as
      # NaN or Inf. NaN are harmless, but Inf could bias the algorithm
      # and should be eliminated
      if (!is.na(Upop$fitness[i]) && !is.nan(Upop$fitness[i]))
        if (Upop$fitness[i] == Inf)
          Upop$fitness[i] = NaN
    }    
  }
  return(list(Upop = Upop, max_weights = max_weights))
}

