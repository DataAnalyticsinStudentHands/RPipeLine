myOddsRatio <- function(x, transactions, reuse = TRUE, CI = FALSE, t=0) {
  
  # This function is based on a function included in the package arules:
  # .basicRuleMeasure
  # https://r-forge.r-project.org/scm/viewvc.php/pkg/R/interestMeasures.R?view=markup&root=arules&sortby=rev&pathrev=1537
  # originall designed to compute various measures of interest for association rules.
  #
  # I modified the function to make it dedicated to odds ratio only, with some new 
  # functionalities (standard error, confidence interval, statistical significance)
  
  # INPUT
  # x = a set of rules.
  # transactions = the transaction data set used to mine the associations or a set of 
  # different transactions to calculate interest measures from (Note: you need to set 
  # reuse=FALSE in the later case).
  # reuse = logical indicating if information in quality slot should be reuse for 
  # calculating the measures.
  # CI = include 95% confidence interval and statistical significance in the output 
  # data frame (logical). 
  # t = tolerance for statistical significance (0<=t<=1).
  
  counts <- .getCounts(x, transactions, reuse)
  N   <- counts$N
  f1x <- counts$f1x
  fx1 <- counts$fx1
  f11 <- counts$f11
  f0x <- counts$f0x 
  fx0 <- counts$fx0
  f10 <- counts$f10
  f01 <- counts$f01
  f00 <- counts$f00
  
  OR <- f11*f00/(f10*f01)
  logSE <- sqrt(1/f11 + 1/f01 + 1/f10 + 1/f00)
  
  # the class rules does not handle well NAs and NaNs in quality measures, so they are
  # substituted with 0.
  ind <- complete.cases(OR)
  OR[!ind] <- 0
  logSE[!ind] <- 0
  
  if (!CI)
    data <- data.frame(oddsRatio=OR, standardError=exp(logSE))
  else {
    logOR <- log(OR)
    log_lowCI <- logOR - 1.96*logSE
    log_highCI <- logOR + 1.96*logSE   
    S<-mapply(is.significant, exp(log_lowCI), exp(log_highCI), t)
    S[!ind] <- FALSE
    data <- data.frame(oddsRatio=OR, standardError=exp(logSE), lowCI95=exp(log_lowCI), highCI95=exp(log_highCI), significance = S)
  }
  
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
.getCounts <- function(x, transactions, reuse = TRUE){
  N <- length(transactions)
  f11 <- interestMeasure(x, "support", transactions, reuse) * N
  f1x <- interestMeasure(x, "coverage", transactions, reuse) * N
  fx1 <- .rhsSupport(x, transactions, reuse) * N
  f0x <- N - f1x
  fx0 <- N - fx1
  f10 <- f1x - f11
  f01 <- fx1 - f11
  f00 <- f0x - f01
  list(f11 = f11, f1x = f1x, fx1 = fx1, 
       f0x = f0x, fx0= fx0, 
       f10 = f10, f01 = f01, f00=f00, 
       N = N)
}


.rhsSupport <- function(x, transactions, reuse = TRUE){
  
  if(is.null(transactions)) stop("transactions missing.")
  N <- length(transactions)
  
  q <- quality(x)
  if(reuse && !is.null(q$confidence) && !is.null(q$lift)) 
    rhsSupport <- q$confidence / q$lift
  else rhsSupport <- support(rhs(x), transactions)
  
  ## for consequents with only one item this might be faster
  ## cons <- unlist(LIST(rhs(x), decode = FALSE))
  ## that's low-level but way faster!
  #cons <- x@rhs@data@i+1
  ## check that the consequents are all singletons
  #if (length(cons) != length(x)) stop("this implementation only works for
  #    rules with one item in the rhs.")
  #c_Y <- itemFrequency(transactions, type = "absolute")[cons]
  #names(c_Y) <- NULL
  
  rhsSupport
}
