#https://cran.r-project.org/web/packages/arules/arules.pdf
library("arules")

QuARM.perform_associative_rule_mining <- function(dataframe, event_column='event', seeds=1:3, training_sample_sizes=c(1, .4, .2), alpha=1.598, support=0.005) {
	training_size <- nrow(dataframe$training)
	
  support_matrix <- matrix(NA, ncol=length(seeds), nrow=length(training_sample_size), dimnames=list(training_sample_size, seeds))
  support_matrix_with_odds_ratio_significance_filter <- matrix(NA, ncol=length(seeds), nrow=length(training_sample_size), dimnames=list(training_sample_size, seeds))
  support_matrix_with_redundancy_filter <- matrix(NA, ncol=length(seeds), nrow=length(training_sample_size), dimnames=list(training_sample_size, seeds))  
	
	testing_dataframe <- data.frame()
	
	for (seed_index in 1:length(seeds)) {
		for (training_index in 1:length(training_sample_sizes)) {
		  set.seed(seed[seed_index])
			training_sample_data <- dataframe$training[sample(training_size, training_sample_sizes[training_index], replace=TRUE), ]
			
			if (sum(dataframe$training[, event_column]) < 1) {
				support_matrix[training_index, seed_index] <- 0
				support_matrix_with_odds_ratio_significance_filter[training_index, seed_index] <- 0
				support_matrix_with_redundancy_filter[training_index, seed_index] <- 0
			}
			else {
				training_transactions <- as(training_sample_data, "transactions")
				training_transactions_negated <- as(!training_sample_data, "transactions")
				rules <- apriori(training_transactions, parameter=list(support=support, confidence=0.0, maxlen=6))
				rules_cases <- subset(rules, subset=rhs %in% event_column)
				rules_cases_dataframe <- as.data.frame(rules_cases)
				support_matrix[j, i] <- nrow(rules_cases_dataframe)
				
				quality(rules_cases) <- cbind(quality(rules_cases, myOddsRatioExp(rules_cases, training_transactions, training_transactions_negated, confidence_interval=TRUE, t=0.03)))
		        significant <- subset(rules_cases, subset=(odds_ratio != Inf))
		        significant <- subset(significant, subset=(significance == TRUE))
		        signRules <- as(significant, "data.frame")     
		
		        if (nrow(signRules) < 1) {
		          support_matrix_with_odds_ratio_significance_filter[j, i] <- 0
		          support_matrix_with_redundancy_filter[j, i] <- 0        
		        }
				
		        else {
		          support_matrix_with_odds_ratio_significance_filter[j,i] <- nrow(signRules)
          
		          quality(significant) <- cbind(quality(significant), custom_confidence_interval(significant, training_transactions, training_transactions_negated, alpha=alpha))          
		          quality(significant) <- cbind(quality(significant), improvement_custom_confidence_interval_par(significant))
          
		          significant <- subset(significant, subset=(imp_confidence_interval_custom_par == TRUE))
		          signRules <- as(significant, "data.frame")
          
		          support_matrix_with_redundancy_filter[j, i]<-nrow(signRules)
				  
		          if (support_matrix_with_redundancy_filter[j,i]>0){
		            testingNotSample <- as.data.frame(!dataframe$testing)
		            testingNotSample[,event_column] <- dataframe$testing[,event_column]
            
		            testingTransactions <- as(dataframe$testing, "transactions")
		            testingNotTransactions <- as(testingNotSample, "transactions")
            
		            ruleTest <- significant
		            ruleTest@quality$support <- NULL
		            ruleTest@quality$confidence <- NULL
		            ruleTest@quality$lift <- NULL
		            ruleTest@quality$odds_ratio <- NULL
		            ruleTest@quality$standard_error <- NULL
		            ruleTest@quality$lowconfidence_interval95 <- NULL
		            ruleTest@quality$highconfidence_interval95 <- NULL
		            ruleTest@quality$significance <- NULL
		            ruleTest@quality$chi_square <- NULL
		            ruleTest@quality$p_value <- NULL
		            ruleTest@quality$lowconfidence_intervalcustom <- NULL
		            ruleTest@quality$highconfidence_intervalcustom <- NULL
		            ruleTest@quality$impconfidence_intervalcustom_par <- NULL
					
		            if (support_matrix_with_redundancy_filter[j, i] == 1) {
		              quality(ruleTest) <- cbind(quality(ruleTest), t(interestMeasure(ruleTest, c("support", "confidence", "lift"), transactions=testingTransactions)))
		            }
					
		            if (support_matrix_with_redundancy_filter[j, i] > 1) {
		              quality(ruleTest) <- cbind(quality(ruleTest), interestMeasure(ruleTest, c("support", "confidence", "lift"), transactions=testingTransactions))
		            }
            
		            quality(ruleTest) <- cbind(quality(ruleTest), myOddsRatioExp(ruleTest, testingTransactions, testingNotTransactions, confidence_interval=TRUE, t=0.03))
		            quality(ruleTest) <- cbind(quality(ruleTest), customconfidence_interval(ruleTest, testingTransactions, testingNotTransactions, alpha=alpha))
		            quality(ruleTest) <- cbind(quality(ruleTest), improvement_customconfidence_interval_par(ruleTest))

		            ruleTest <- as(ruleTest, "data.frame")
		            testing_dataframe <- rbind(testing_dataframe, ruleTest)            
		          }
			   }
			}
		}
	}
	
	return(list(rules=testing_dataframe, statistics=list(support=support_matrix, signficant=support_matrix_with_odds_ratio_significance_filter, nonredundant=support_matrix_with_redundancy_filter)))
}

#' My Odds Ratio
#' 
#' This function is based on a function included in the package arules basicRuleMeasure originally designed to compute various measures of interest for association rules
#' The function has been modified to make it dedicated to odds ratio only, with some new functionalities such as standard error, confidence interval, and statistical significance.
#' WARNING: THE VALUES OUTPUTED BY THIS FUNCTION ARE ONLY VALID FOR RULES THAT 
#' HAVE ONLY OUTPUT OF INTEREST IN RHS, FOR DEFINITION OF TRANSNOT
#' @usage myOddsRatioExp(x, transactions, negated_transactions, reuse=TRUE, confidence_interval=FALSE, t=0)
#' @param x A set of rules.
#' @param transactions The transaction data set used to mine the associations or a set of different transactions to calculate interest measures from (Note: you need to set reuse=FALSE in the later case).
#' @param negated_transactions The negated transaction data for calculations
#' @param reuse A logical indicating if information in quality slot should be reuse for calculating the measures.
#' @param confidence_interval Include 95 percent confidence interval and statistical significance in the output data frame (logical). Default=FALSE.
#' @param t Tolerance for statistical significance (0<=t<=1). Default=0 (no tolerance)
#' @return return: a data frame containing at least 2 columns, odds_ratio and standard_error. If confidence_interval=TRUE, 3 additional columns are included in the data frame: lowconfidence_interval95, highconfidence_interval95, and significance. lowconfidence_interval95 and highconfidence_interval95 store the margins of the 95 percent confidence interval. significance is a logical column that states if the rule is statistically significant given the specified tolerance. When t=0, a rule is statistically significant if the confidence interval does not include 1. When t > 0, a rule is statistically significant in one of these cases:the confidence interval does not include 1, highconfidence_interval95 <= 1+t and lowconfidence_interval95 <= 1-t, highconfidence_interval95 >= 1+t and lowconfidence_interval95 >= 1-t
#' @export
help.QuARM.measure_rules_with_odds_ratio <- function(rules, transactions, transactions_negated, reuse_quality_data=TRUE, filter_by_confidence_interval=FALSE, significance_tolerance=0) {
  samples_total <- length(transactions)
  
  exposed_cases <- interestMeasure(x, "support", transactions, reuse=TRUE) * samples_total
  exposed_total <- interestMeasure(x, "coverage", transactions, reuse=TRUE) * samples_total
  exposed_controls <- exposed_total - exposed_cases
  
  unexposed_cases <- interestMeasure(x, "support", negated_transactions, reuse=FALSE) * samples_total
  unexposed_total <- interestMeasure(x, "coverage", negated_transactions, reuse=FALSE) * samples_total
  unexposed_controls <- unexposed_total - unexposed_cases
  
  controls_total <- unexposed_controls + exposed_controls
  cases_total <- exposed_cases + unexposed_cases
  
  proportion_exposed <- (exposed_cases + exposed_controls)/nrow(binarize_dataframe)
  expected_exposed_cases <- proportion_exposed * number_of_cases # assuming independence
  expected_unexposed_cases <- number_of_cases - expected_exposed_cases
  expected_exposed_controls <- proportion_exposed * number_of_controls # assuming independence
  expected_unexposed_controls <- number_of_control - expected_exposed_controls
  
  chi_square <- (exposed_cases - expected_exposed_cases)^2/expected_exposed_cases + (unexposed_cases - expected_unexposed_cases)^2/expected_unexposed_cases + (exposed_controls - expected_exposed_controls)^2/expected_exposed_controls + (unexposed_controls - expected_unexposed_controls)^2/expected_unexposed_controls
  p_value <- pchisq(chi_square, 1, lower.tail=FALSE) 
  
  odds_ratio <- (exposed_cases*unexposed_controls)/(exposed_controls*unexposed_cases)
  log_odds_ratio_standard_error <- sqrt(1/exposed_cases + 1/unexposed_cases + 1/exposed_controls + 1/unexposed_controls)
  
  incomplete_odds_ratios <- !complete.cases(odds_ratio)
  odds_ratio[incomplete_odds_ratios] <- 0
  log_odds_ratio_standard_error[incomplete_odds_ratios] <- 0
  
  if (filter_by_confidence_interval) {
    log_odds_ratio <- log(odds_ratio)
    log_low_confidence_interval <- log_odds_ratio - 1.96 * log_odds_ratio_standard_error
    log_high_confidence_interval <- log_odds_ratio + 1.96 * log_odds_ratio_standard_error 
    low_confidence_interval <- exp(log_low_confidence_interval)
    high_confidence_interval <- exp(log_high_confidence_interval)  
    
    signficance_tests <- mapply(help.QuARM.test_significance, low_confidence_interval, high_confidence_interval, significance_tolerance)
    significance_tests[incomplete_odds_ratios] <- FALSE
    return_data <- data.frame(odds_ratio=odds_ratio, standard_error=exp(log_odds_ratio_standard_error), low_confidence_interval95=low_confidence_interval, high_confidence_interval95=high_confidence_interval, significance=signficicance, chi_square=chi_square, p_value=p_value)
  }
  else
    return_data <- data.frame(odds_ratio=odds_ratio, standard_error=exp(log_odds_ratio_standard_error), chi_square=chi_square, p_value=p_value)

  
  return(return_data)  
}


help.QuARM.test_significance <- function(low, high, t=0) {
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

help.improvement <- function(rules, t = 0.1) {
  
  # This function determines whether the OR of a rule is interesting or if it just 
  # a consequence of a parent rule. The goal is to eliminate redundant rules with a 
  # Occam's Razor strategy
  
  # Exctract interesting data from rules structure
  i <- rules@lhs@data@i   #tags
  p <- rules@lhs@data@p   #indexes
  length <- length(p)-1   #total rules
  
  Rules <- list()
  
  # First, rules are organized in a list structure
  for (ind in 1:length) {
    size <- p[ind+1] - p[ind]
    start<-p[ind]+1
    end <- p[ind]+size
    s <- start:end
    temp <- i[s]
    Rules <- c(Rules,list(temp))
  }
  
  # initialize logical array "imp" -> tag improved rules
  imp <- logical(length(Rules))
  imp[1:length(imp)] <- TRUE
  
  for (i in 1:length(Rules)) {
    if (length(Rules[[i]])>1) { # rules of length 1 are all significant 
      Sets <- list()
      # produce Sets: list of all subset of current rule
      for (j in 1:length(Rules[[i]])-1)
        Sets <- c(Sets, combn(Rules[[i]],j,simplify=FALSE))
      # look up elements of Sets in Rules and evaluate improvement
      for (j in Sets) {
        if (length(j) > 0) {
          pos <- match(list(j),Rules)
          if (!is.na(pos)) {
            if (abs(rules@quality$oddsRatio[i]-rules@quality$oddsRatio[pos])<t)
              imp[i] <- FALSE
          }
        }
      }
    }
  }
  
  return(data.frame(imp=imp))
}

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

or.table <- function(rules, poll, lag, sortby = 0, file="ORplot.png") {
  
  # This function organizes rules (from a data frame) in tables according to the
  # specified pollutant, such as:
  #
  # O3 table
  # 
  # single  <rule1> <rule2> ... <rulen>
  # O3_0    val0    val0    ... val0
  # O3_1    val1    val1    ... val1
  # O3_2    val2    val2    ... val2
  #
  # These tables are meant to be used to produce readable graphs of the rules
  
  # convert rules from factor to character if necessary
  if (class(rules$rules)=="factor")
    rules$rules <- as.character(rules$rules)
  
  # clean data frame "rules" from useless characters
  for (i in 1:nrow(rules)) {
    rules$rules[i] <- gsub("event", "", rules$rules[i])
    rules$rules[i] <- gsub("val", "", rules$rules[i])
    rules$rules[i] <- gsub("([{=>}])", "", rules$rules[i])
  }
  
  # separate single rules from the rest
  ind <- grep(",", rules$rules)
  singles <- rules[-ind,]
  rules <- rules[ind,]
  
  # order single rules and find single OR for the required pollutant
  singles <- singles[order(singles$rules),]
  singleOR <- rep(0,lag)
  
  for (j in 0:lag) {
    temp <- singles$oddsRatio[grep(paste(j,"_",poll,sep=""), singles$rules)]
    if (length(temp)>0)
      singleOR[j+1] <- temp
    else singleOR[j+1] <- NA
  }
  
  # start building the data frame
  ORtable <- data.frame(singleOR)
  
  r <- grep(poll,rules$rules)
  
  for (i in r) {
    # each rule is unpacked and rebuilt to have the rule for 0_poll, 1_poll and 2_poll
    str <- rules$rules[i]
    str <- gsub(" ","",str) #removing spaces
    str <- unlist(strsplit(str,","))
    or <- rep(0,lag+1)
    for (j in 0:lag) {
      str0 <- str
      str0[grep(poll,str)[1]] <- paste(j,"_",poll,sep="")
      if (length(unique(str0))==length(str0)) {
        tab <- rules
        for (x in str0) {
          ind <- grep(x,tab$rules)
          tab <- tab[ind,]
        } 
        if (length(tab$oddsRatio)>0) {
          # count "," in rule to be sure to have the one of right length
          s2 <- gsub(",","",tab$rules)
          if (nchar(tab$rules) - nchar(s2) < length(str0))
            or[j+1] <- tab$oddsRatio[1]
          else or[j+1] <- NA
        }
        else or[j+1] <- NA
      }
      else or[j+1] <- NA 
    }
    ORtable <- cbind(ORtable,or)
    name<-str0
    name<-name[-grep(poll,name)[1]]    
    colnames(ORtable)[ncol(ORtable)] <- paste(name,collapse=",")
  }
  
  ORtable <- ORtable[,unique(colnames(ORtable))]
  ORtable <- t(ORtable)
  ORtable <- ORtable[order(ORtable[,sortby+1
                                   ]),] # sort table
  colnames(ORtable) <- c(paste("0_",poll,sep=""),paste("1_",poll,sep=""),paste("2_",poll,sep=""),paste("3_",poll,sep=""),paste("4_",poll,sep=""))
  
  png(filename=file)
  op <- par(mar = c(10,4,3,2) + 0.1)
  plot(ORtable[,1],type="b",lwd=2, ann=FALSE, xaxt="n",ylim=c(0.5,1.5),col="black",ylab="OR",main=paste("OR for combinations with",poll))
  axis(1,at=1:length(rownames(ORtable)),labels=rownames(ORtable),las=3)
  lines(ORtable[,2],col="red",type="b",lwd=2)
  lines(ORtable[,3],col="green",type="b",lwd=2)
  lines(ORtable[,4],col="purple",type="b",lwd=2)
  lines(ORtable[,5],col="blue",type="b",lwd=2)
  legend("bottomright",legend=c(paste("0_",poll,sep=""),paste("1_",poll,sep=""),paste("2_",poll,sep=""),paste("3_",poll,sep=""),paste("4_",poll,sep="")),
         border = "black",lty=1,lwd=2,pch=21,col=c("black","red","green","purple","blue"),horiz=TRUE,bty="n",cex=0.8,
         text.col=c("black","red","green","purple","blue"),
         inset=0.01)
  grid()
  par(op)
  dev.off()
  
  return(ORtable)
}

#'Nonredundant
#' 
#' This function determines whether the OR of a rule is interesting or if it just a consequence of a parent rule. The goal is to eliminate redundant rules with a Occam's Razor strategy. 
#' It is very similar to improvementCI() but accepts a data frame input.It is assumed all rules in the data frame will be unique.
#' @usage nonredundant(rules, variables, lowCI, highCI)
#' @param rules A data frame of rules. The data frame must have a column labeled rules with the variables used in each rule, and columns for both the lower and higher bounds of the odds ratio confidence interval. The columns for the the confidence interval should be labeled lowCIcustom and highCIcustom. The TestFrame from the list output of overfit_test() is already formated in this way.
#' @param variables This should be a vector of the names of each variable used in the making of the rules as they will appear in the rules. They should be in string form.
#' @param lowCI the name of the column with the lower bound of the confidence interval
#' @param highCI the name of the column with the higher bound of the confidence interval
#' @return A logical vector of the length rules with non-redundant rules labeled TRUE and redundant rules labeled FALSE.
#' @export
nonredundant <- function(rules, variables, lowCI="lowCIcustom", highCI="highCIcustom") {
  
  rules=unique(rules)
  rules$rules=as.character(rules$rules)
  rrules=rules
  rrules$lowCIcustom=rules[,eval(lowCI)]
  rrules$highCIcustom=rules[,eval(highCI)]
  nonredundant <- logical(length(rules$rules))
  nonredundant[1:length(nonredundant)] <- TRUE
  
  for (i in 1:length(variables)){
    b=vector()
    a=c(grep(paste(variables[i],"}",sep=""),rules$rules),grep(paste(variables[i],",",sep=""),rules$rules))
    # if a variable is used in multiple rules check odds ratio with odds ratio of shortest rule, then change the shortest rule to NA and repeat 
    if (length(a)>1){
      for (j in 1:length(a)){
        c=grep(",",rules$rules[a[j]])
        b[j]=length(c)+1
      }
      for (s in 1:length(a-2)){
        for (k in 1: length(b)){
          if (is.na(b[k]) == FALSE){
            if (b[k] != min(b,na.rm=TRUE)){
              if (rrules$lowCIcustom[a[which.min(b)]] <= rrules$lowCIcustom[a[k]] & rrules$lowCIcustom[a[k]] <= rrules$highCIcustom[a[which.min(b)]]){
                nonredundant[a[k]]=FALSE
              }
              if (rrules$highCIcustom[a[which.min(b)]] >= rrules$highCIcustom[a[k]] & rrules$highCIcustom[a[k]] >= rrules$lowCIcustom[a[which.min(b)]]){
                nonredundant[a[k]]=FALSE
              }
              if (rrules$lowCIcustom[a[which.min(b)]] >= rrules$lowCIcustom[a[k]] & rrules$highCIcustom[a[k]] >= rrules$highCIcustom[a[which.min(b)]]){
                nonredundant[a[k]]=FALSE
              }
            }
            
          }
          
          
        }
        b[which.min(b)]=NaN
        
      }
      
    }
    
  }
  rules$nonredundant=nonredundant
  return(rules)
}

#' Improvement Custom CI
#' 
#' This function determines whether the OR of a rule is interesting or if its importance depends only on one of its subsets. The goal is to identify and possibly eliminate redundant rules with a Occam's Razor strategy
#' @usage improvement_customCI_par(rules, t = 0.1)
#' @param rules An object of class rules. This object must include the column oddsRatio, added using the function myOddsRatio().
#' @param t The threshold of minimum improvement in odds ratio to declare a rule more important than all its subsets. Default = 0.1
#' @return d: data frame containing one logical column "imp". Rules that do not show significant improvement are labeled as FALSE. This single column data frame can be added to the other quality measures of the object of class rules.
#' @export
improvement_customCI_par <- function(rules, t=0.1) {
  
  # This function determines whether the OR of a rule is interesting or if it just 
  # a consequence of a parent rule. The goal is to eliminate redundant rules with a 
  # Occam's Razor strategy
  
  # The parameter used to avoid redundancy is OR minimum difference (defined by t)
  
  # Exctract interesting data from rules structure
  i <- rules@lhs@data@i   #tags
  p <- rules@lhs@data@p   #indexes
  length <- length(p)-1   #total rules
  
  Rules <- list()
  
  # First, rules are organized in a list structure
  for (ind in 1:length) {
    size <- p[ind+1] - p[ind]
    start<-p[ind]+1
    end <- p[ind]+size
    s <- start:end
    temp <- i[s]
    Rules <- c(Rules,list(temp))
  }
  
  # uncomment if not testing only significant rules
  # Rules[1] <- 0
  
  # initialize logical array "imp95CI" -> tag improved rules
  impCIcustom <- logical(length(Rules))
  impCIcustom[1:length(impCIcustom)] <- TRUE
  
  # parallelization of improvement function
  library(doParallel)
  library(foreach)
  registerDoParallel(4)
  
  imptag <- function(i, Rules, rules) {
    
    if (length(Rules[[i]])==1) { return(TRUE) } # rules of length 1 are all significant 
    else {
      Iflag = TRUE
      Sets <- list()
      # produce Sets: list of all subset of current rule
      for (j in 1:length(Rules[[i]])-1)
        Sets <- c(Sets, combn(Rules[[i]],j,simplify=FALSE))
      # look up elements of Sets in Rules and evaluate improvement
      for (j in Sets) {
        if (length(j) > 0) {
          pos <- match(list(j),Rules)
          if (!is.na(pos)) {
            if (rules@quality$oddsRatio[i]>rules@quality$oddsRatio[pos]) {
              if (rules@quality$lowCIcustom[i]<rules@quality$highCIcustom[pos])
                Iflag = FALSE
            }
            if (rules@quality$oddsRatio[i]<rules@quality$oddsRatio[pos]) {
              if (rules@quality$highCIcustom[i]>rules@quality$lowCIcustom[pos])
                Iflag = FALSE
            }
          }
        }
      }
    }
    
    return(Iflag)
  }
  
  strt<-Sys.time()
  
  impCIcustom <- foreach(i=1:length(Rules)) %dopar% {
    imptag(i, Rules, rules)
  }
  stopImplicitCluster()
  print(Sys.time()-strt)
  
  A <- array(unlist(impCIcustom), dim = c(nrow(impCIcustom[[1]]), ncol(impCIcustom[[1]]), length(impCIcustom)))
  
  return(data.frame(impCIcustom_par=A))
  
}


