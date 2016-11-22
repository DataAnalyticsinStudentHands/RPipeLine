process.using_ARM <- function(dataframe, event_column='event', seeds=1:3, training_sample_sizes=c(1, .4, .2), alpha=1.598, support=0.005, remove_insignificant_results=TRUE, remove_redundant_results=TRUE) {
	source('libs/ARM.R')
	results <- ARM.perform_associative_rule_mining(dataframe, event_column, seeds, training_sample_sizes, alpha, support)
	filtered_results <- results
	if (remove_insignificant_results) {
	  filtered_results = unique(subset(filtered_results$dataframe, significance == TRUE))
	}
	if (remove_redundant_results) {
	  filtered_results$nonredundant = nonredundant(filtered_results, colnames(dataframe)[-1])
	  filtered_results = subset(filtered_results, nonredundant == TRUE)
	}
	
	return(filtered_results)
}

process.using_GA <- function(dataframe, number_of_generations=200, gene_pool_size=30, initial_attribute_probability=0.25, seed=NULL,
                                      weights=rep(1,4), crossover_probability=0.6, parental_inheritance_proportion=0.5, tournament_size=4, mutation_probability_range=c(0.001, 0.01),
                                      tournament_win_probability_range=c(0.9, 0.6), fitness_mode=c("sum", "arithmetic", "geometric"), with_statistics=TRUE) {
  source('libs/GA.R')
  rules <- GA.perform_genetic_algorithm_rule_mining(dataframe, number_of_generations, gene_pool_size, initial_attribute_probability, seed,
             weights, crossover_probability, parental_inheritance_proportion, tournament_size, mutation_probability_range,
             tournament_win_probability_range, fitness_mode)
  
  statistics <- NULL
  if (with_statistics)
    statistics <- GA.statistics_report(dataframe, rules, gene_pool_size, weights)
    
  return(list(rules=rules, statistics=statistics))
}