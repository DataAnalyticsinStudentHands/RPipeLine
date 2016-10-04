# Fitness function for GA rule mining based on aritmetic average

fitness.fun <- function(Upop, upop_size, genes_num, test, weights) {
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
      length <- length(A)/genes_num
      if(is.na(length))
        length <- 0
      OR <- OR.fitness(testbin)
      if(is.na(OR))
        OR <- 0
      
      Upop$fitness[i] <- (weights[1]*length + weights[2]*OR)/2
      Upop$stats[i,] <- confidence.int(testbin)
      # some rules can end up having a "non numerical" fitness, such as
      # NaN or Inf. NaN are harmless, but Inf could bias the algorithm
      # and should be eliminated
      if (!is.na(Upop$fitness[i]) && !is.nan(Upop$fitness[i]))
        if (Upop$fitness[i] == Inf)
          Upop$fitness[i] = NaN
    }    
  }
  return(Upop)
}