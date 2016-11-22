ARM.perform_associative_rule_mining <- function(dataframe, event_column='event', seeds=1:3, training_sample_sizes=c(1, .4, .2), alpha=1.598, support=0.005) {
	library("arules")
	
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
		            ruleTest@quality$poisson_fitness <- NULL
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
help.ARM.measure_rules_with_odds_ratio <- function(rules, transactions, transactions_negated, reuse_quality_data=TRUE, filter_by_confidence_interval=FALSE, significance_tolerance=0) {
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
  poisson_fitness <- pchisq(chi_square, 1, lower.tail=FALSE) 
  
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
    
    signficance_tests <- mapply(help.ARM.test_significance, low_confidence_interval, high_confidence_interval, significance_tolerance)
    significance_tests[incomplete_odds_ratios] <- FALSE
    return_data <- data.frame(odds_ratio=odds_ratio, standard_error=exp(log_odds_ratio_standard_error), low_confidence_interval95=low_confidence_interval, high_confidence_interval95=high_confidence_interval, significance=signficicance, chi_square=chi_square, poisson_fitness=poisson_fitness)
  }
  else
    return_data <- data.frame(odds_ratio=odds_ratio, standard_error=exp(log_odds_ratio_standard_error), chi_square=chi_square, poisson_fitness=poisson_fitness)

  
  return(data)  
}


help.ARM.test_significance <- function(low, high, t=0) {
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
