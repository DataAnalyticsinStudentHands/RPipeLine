#### UTILITY FUNCTIONS ####
odds.ratio <- function(testbin) {
  cases <- testbin[testbin[,1],]
  controls <- testbin[!testbin[,1],]
  x11 <- sum(rowSums(cases)==ncol(cases)) # exposed cases
  x01 <- nrow(cases) - x11 # unexposed cases (general definition)
  x10 <- sum(rowSums(controls)==ncol(controls)-1) # exposed controls
  x00 <- nrow(controls) - x10 # unexposed controls
  return((x11*x00)/(x01*x10))
}

OR.fitness <- function(testbin) {
  cases <- testbin[testbin[,1],]
  controls <- testbin[!testbin[,1],]
  x11 <- sum(rowSums(cases)==ncol(cases)) # exposed cases
  x01 <- nrow(cases) - x11 # unexposed cases (general definition)
  x10 <- sum(rowSums(controls)==ncol(controls)-1) # exposed controls
  x00 <- nrow(controls) - x10 # unexposed controls
  
  OR <- (x11*x00)/(x01*x10)
  logOR <- log(OR)
  logSE <- sqrt(1/x11 + 1/x01 + 1/x10 + 1/x00)
  lowCI <- exp(logOR - 1.96*logSE)
  highCI <- exp(logOR + 1.96*logSE)
  
  if (!is.na(lowCI) && !is.na(highCI)) {
    if (lowCI <= 1 && highCI >= 1)
      return (0)
    else
      if (OR < 1)
        return (1-OR)
    else
      return ((1/(1+exp(-(OR-1)))-1/2)*2) # OR fitness limited between 0 and 1
    # May add coefficient to increase 
    # steepness between 1 and 2 (i.e. -3*(OR-1))
  }    
  else return(NaN) # NaN is outputted when one slot of the contingency
  # table is empty (0)
}

confidence.int <- function(testbin) {
  cases <- testbin[testbin[,1],]
  controls <- testbin[!testbin[,1],]
  x11 <- sum(rowSums(cases)==ncol(cases)) # exposed cases
  x01 <- nrow(cases) - x11 # unexposed cases (general definition)
  x10 <- sum(rowSums(controls)==ncol(controls)-1) # exposed controls
  x00 <- nrow(controls) - x10 # unexposed controls
  
  logOR <- log((x11*x00)/(x01*x10))
  logSE <- sqrt(1/x11 + 1/x01 + 1/x10 + 1/x00)
  log_lowCI <- logOR - 1.96*logSE
  log_highCI <- logOR + 1.96*logSE
  
  perc_exp <- (x11 + x10)/(x11 + x01 + x10 + x00)
  C <- perc_exp*(x11+x01) # expected cases if independent
  K <- perc_exp*(x10+x00) # expected controls if independent
  chi_squared <- (x11-C)^2/C + (x01-(x11+x01-C))^2/(x11+x01-C) + (x10-K)^2/K + (x00-(x10+x00-K))^2/(x10+x00-K)
  p <- pchisq(chi_squared,1,lower.tail=FALSE) 
  
  data <- c((x11*x00)/(x01*x10), exp(log_lowCI), exp(log_highCI), p)
  return(data)
}

# This function generates random indexes to pick chromosomes to
# use to generate new ones. Low indexes are favored, because they 
# represent higher ranking chromosomes. In this implementation, 
# the strongest rule has a probability p to win the tournament. 
# The second strongest rule can win with a probability p*(1-p), 
# the third with a probability p(1-p)^2, and so on.
tournament <- function(k, upop_size, p){
  index <- sort(sample(1:upop_size, k))
  win_prob <- runif(1)
  t <- 0 # threshold
  for (i in 1:k) {
    t <- t + p*((1-p)^(i-1))
    if (win_prob < t)
      return(index[i])
  }
  # In the eventuality no gene is selected:
  return(index[k])
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
        #thresholds = arithmetic average
        c1$t[i] <- (c1$t[i]+c2$t[i])/2
        c2$t[i] <- (c1$t[i]+c2$t[i])/2
      }
    }  
  } 
  return(children <- list(c1 = c1, c2 = c2))
}

# Mutation function
mutate <- function(cross, pm, Emedian) {
  for(n in 1:2){
    for(i in 1:length(cross[[n]]$ant)){
      if(runif(1) < pm) {
        if (cross[[n]]$ant[i]==1)
          cross[[n]]$ant[i] <- 0
        else cross[[n]]$ant[i] <- 1
        # mutating threshold
        delta <- 10*pm*Emedian[i] # percentage of median dependent on pm, scaled by 10 to produce larger changes
        # coin toss to decide if subtract or add
        if(runif(1) < 0.5) 
          cross[[n]]$t[i] <- cross[[n]]$t[i] - delta
        else cross[[n]]$t[i] <- cross[[n]]$t[i] + delta       
      }
    }
  }
  return(cross)
}

rep.penalty <- function(Upop, i, upop_size) {
  # Define array of penalties (monotonically decreasing)
  P <- seq(from = 1, to = 0, by = -1/(upop_size-1)) 
  Ptot <- 0
  
  if (i != 1) {
    for (j in 1:(i-1)){
      # check if chromosomes are equal to those to a previous rule
      equal <- sum(abs(Upop$ant[i,] - Upop$ant[j,]))  # equal = 0 if the rules share the same chromosomes
      if (equal == 0) {
        Ptot = Ptot + P[j]
      }
    } 
  }
  # return total penalty
  return(Ptot)
}


pop.similarity <- function(ant) {
  sim <- abs(2*colSums(ant)-nrow(ant))
  sim <- sim/nrow(ant)
  return(sum(sim)/length(sim))
}

sim.adjust <- function(range, sim) {
  return(range[1] + (range[2]-range[1])*sim)
}

median.distance <- function(thresholds, medians, maxdist) {
  dist <- abs(thresholds - medians)/maxdist # distance from median normalized by maximum distance possible
  dist <- dist[!is.na(dist)] # drop NaN caused by logical features
  return(sum(dist)/length(dist))   
}


pop.sort <- function(Upop) {
  Upop_rank <- order(Upop$fitness, decreasing=TRUE)
  Upop$fitness <- Upop$fitness[Upop_rank]
  Upop$ant <- Upop$ant[Upop_rank,]
  Upop$threshold <- Upop$threshold[Upop_rank,]
  Upop$stats <- Upop$stats[Upop_rank,]
  return(Upop)
}


redundancy <- function(Upop, i, upop_size) {
  child <- Upop$ant[i,]
  if (sum(child) < 2)
    return(0) # rules of length 1 can not be redundant
  if(!complete.cases(Upop$stats[i,]))
    return(1) # maximal penalty for "impossible" rules
  
  penalty <- 0
  count <- 0

  for (j in 1:upop_size) {
    if (j!=i && complete.cases(Upop$stats[j,])) {
      testP <- child - Upop$ant[j,]
      if (sum(testP < 0) == 0 && sum(Upop$ant[j,])>0) {
        a <- Upop$stats$low95CI[j]
        b <- Upop$stats$high95CI[j]
        c <- Upop$stats$low95CI[i]
        d <- Upop$stats$high95CI[i]
        if (sum(testP) == 0) # Equal rules are not penalized by this function
          penalty <- penalty + 0
        else if (d<a || c>b)
          penalty <- penalty + 0
        else if (d>b && c<a) { # child covers parent
            penalty <- penalty + 1
            count <- count + 1
          }
        else if (c>a && d<b) { # parent covers child (not sure if it is possible)
          penalty <- penalty + (d-c)/(b-a)
          count <- count + 1
        }
        else { # partial overlap
          if (c>a) {
            penalty <- penalty + (b-c)/(b-a)
            count <- count + 1
          }
          else {
            penalty <- penalty + (d-a)/(b-a)
            count <- count + 1
          }
        }
        # no penalty is added in case of no overlap
      }
    }
  }
  if (count == 0)
    return(0)
  else return(penalty/count)
}
#### END UTILITY FUNCTIONS ####