customCI <- function(x, transactions, transNot, reuse = TRUE, alpha = 1.96, t=0) {
  
  # This function is based on a function included in the package arules:
  # .basicRuleMeasure
  # https://r-forge.r-project.org/scm/viewvc.php/pkg/R/interestMeasures.R?view=markup&root=arules&sortby=rev&pathrev=1537
  # originall designed to compute various measures of interest for association rules.
  #
  # I modified the function in order to create confidence interval based on 
  # different probabilities (default 95% CI)
  #
  # INPUT
  # x = a set of rules.
  # transactions = the transaction data set used to mine the associations or a set of 
  # different transactions to calculate interest measures from (Note: you need to set 
  # reuse=FALSE in the later case).
  # reuse = logical indicating if information in quality slot should be reuse for 
  # calculating the measures.
  # alpha = critical value from normal distribution
  
  counts <- .getCountsExp(x, transactions, transNot)
  N   <- counts$N #total entries
  f1x <- counts$f1x #total exposed
  fx1 <- counts$fx1 #total cases
  f11 <- counts$f11 #exposed cases
  f0x <- counts$f0x #total non-exposed
  fx0 <- counts$fx0 #total controls
  f10 <- counts$f10 #exposed controls
  f01 <- counts$f01 #non-exposed cases
  f00 <- counts$f00 #non-exposed controls
  
  OR <- (f11*f00)/(f10*f01)
  logSE <- sqrt(1/f11 + 1/f01 + 1/f10 + 1/f00)
  
  # the class rules does not handle well NAs and NaNs in quality measures, so they are
  # substituted with 0.
  ind <- complete.cases(OR)
  OR[!ind] <- 0
  logSE[!ind] <- 0

  logOR <- log(OR)
  log_lowCI <- logOR - alpha*logSE
  log_highCI <- logOR + alpha*logSE   
  data <- data.frame(lowCIcustom=exp(log_lowCI), highCIcustom=exp(log_highCI))
  
  return(data)  
}


is.significant <- function(low, high, t = 0) {
  
  if (is.na(low) | is.na(high))
  {return(NA)}
  else if (low >= 1)
  {return(TRUE)}
  else if (high <= 1)
  {return(TRUE)}
  else if (high <= 1+t & low <= 1-t)
  {return(TRUE)}
  else if (high >= 1+t & low >= 1-t)
  {return(TRUE)}
  else {return(FALSE)}
  
}

## count helpers
.getCountsExp <- function(x, transactions, transNot){
  N <- length(transactions)
  f11 <- interestMeasure(x, "support", transactions, reuse = TRUE) * N
  f1x <- interestMeasure(x, "coverage", transactions, reuse = TRUE) * N
  f01 <- interestMeasure(x, "support", transNot, reuse=FALSE) * N
  f0x <- interestMeasure(x, "coverage", transNot, reuse=FALSE) * N
  f10 <- f1x - f11
  f00 <- f0x - f01
  fx0 <- f00 + f10
  fx1 <- f11 + f01
  list(f11 = f11, f1x = f1x, fx1 = fx1, 
       f0x = f0x, fx0= fx0, 
       f10 = f10, f01 = f01, f00=f00, 
       N = N)
}


# .rhsSupport <- function(x, transactions, reuse = TRUE){
#   
#   if(is.null(transactions)) stop("transactions missing.")
#   N <- length(transactions)
#   
#   q <- quality(x)
#   if(reuse && !is.null(q$confidence) && !is.null(q$lift)) 
#     rhsSupport <- q$confidence / q$lift
#   else rhsSupport <- support(rhs(x), transactions)
#   
#   ## for consequents with only one item this might be faster
#   ## cons <- unlist(LIST(rhs(x), decode = FALSE))
#   ## that's low-level but way faster!
#   #cons <- x@rhs@data@i+1
#   ## check that the consequents are all singletons
#   #if (length(cons) != length(x)) stop("this implementation only works for
#   #    rules with one item in the rhs.")
#   #c_Y <- itemFrequency(transactions, type = "absolute")[cons]
#   #names(c_Y) <- NULL
#   
#   rhsSupport
# }
