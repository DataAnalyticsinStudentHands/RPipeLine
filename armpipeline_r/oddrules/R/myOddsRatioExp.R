#' My Odds Ratio
#' 
#' This function is based on a function included in the package arules basicRuleMeasure originally designed to compute various measures of interest for association rules
#' The function has been modified to make it dedicated to odds ratio only, with some new functionalities such as standard error, confidence interval, and statistical significance.
#' @usage myOddsRatioExp(x, transactions, transNot, reuse = TRUE, CI = FALSE, t=0)
#' @param x A set of rules.
#' @param transactions The transaction data set used to mine the associations or a set of different transactions to calculate interest measures from (Note: you need to set reuse=FALSE in the later case).
#' @param transNot The negated transaction data for calculations
#' @param reuse A logical indicating if information in quality slot should be reuse for calculating the measures.
#' @param CI Include 95 percent confidence interval and statistical significance in the output data frame (logical). Default = FALSE.
#' @param t Tolerance for statistical significance (0<=t<=1). Default = 0 (no tolerance)
#' @return return: a data frame containing at least 2 columns, oddsRatio and standardError. If CI = TRUE, 3 additional columns are included in the data frame: lowCI95, highCI95, and significance. lowCI95 and highCI95 store the margins of the 95 percent confidence interval. significance is a logical column that states if the rule is statistically significant given the specified tolerance. When t = 0, a rule is statistically significant if the confidence interval does not include 1. When t > 0, a rule is statistically significant in one of these cases:the confidence interval does not include 1, highCI95 <= 1+t and lowCI95 <= 1-t, highCI95 >= 1+t and lowCI95 >= 1-t
#' @export
myOddsRatioExp <- function(x, transactions, transNot, reuse = TRUE, CI = FALSE, t=0) {
  
  # This function is based on a function included in the package arules:
  # .basicRuleMeasure
  # https://r-forge.r-project.org/scm/viewvc.php/pkg/R/interestMeasures.R?view=markup&root=arules&sortby=rev&pathrev=1537
  # originall designed to compute various measures of interest for association rules.
  #
  # I modified the function to make it dedicated to odds ratio only, with some new 
  # functionalities (standard error, confidence interval, statistical significance,
  # chi-squared value and p-value)
  
  # WARNING: THE VALUES OUTPUTED BY THIS FUNCTION ARE ONLY VALID FOR RULES THAT
  # HAVE ONLY OUTPUT OF INTEREST IN RHS, FOR DEFINITION OF TRANSNOT
  
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
  
  # computing rules chi-squared value
  perc_exp <- f1x/(f1x+f0x)
  C <- perc_exp*fx1 # expected cases if independent
  K <- perc_exp*fx0 # expected controls if independent
  chi_squared <- (f11-C)^2/C + (f01-(fx1-C))^2/(fx1-C) + (f10-K)^2/K + (f00-(fx0-K))^2/(fx0-K)
  p <- pchisq(chi_squared,1,lower.tail=FALSE) 
  
  # the class rules does not handle well NAs and NaNs in quality measures, so they are
  # substituted with 0.
  ind <- complete.cases(OR)
  OR[!ind] <- 0
  logSE[!ind] <- 0
  
  if (!CI)
    data <- data.frame(oddsRatio=OR, standardError=exp(logSE), chi_squared=chi_squared, p_value=p)
  else {
    logOR <- log(OR)
    log_lowCI <- logOR - 1.96*logSE
    log_highCI <- logOR + 1.96*logSE   
    S<-mapply(is.significant, exp(log_lowCI), exp(log_highCI), t)
    S[!ind] <- FALSE
    data <- data.frame(oddsRatio=OR, standardError=exp(logSE), lowCI95=exp(log_lowCI), highCI95=exp(log_highCI), significance = S, chi_squared=chi_squared, p_value=p)
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
