##### GenARM ALGORITHM FOR RULE MINING - OR ORIENTED #####

# This script implements a genetic algorithm for the automatic
# mining of association rules of the form {exposures} -> event.
# The design is based on the algorithm proposed in Alatas-05.

GenARM.perform_genetic_algorithm_rule_mining <- function(dataframe, number_of_generations=200, gene_pool_size=30, initial_attribute_probability=0.25, seed=NULL,
                                                     weights=rep(1,4), crossover_probability=0.6, parental_inheritance_proportion=0.5, tournament_size=4, mutation_probability_range=c(0.001, 0.01),
                                                     tournament_win_probability_range=c(0.9, 0.6), fitness_mode=c("sum", "arithmetic", "geometric")) {
  
  # ARGUMENTS: 
  # data: data frame to mine (must contain the event information in 
  #       the first column) The first column of this 
  # data frame must contain the event information, while the 
  # remaining columns are exposures.
  # number_of_generations: number of number_of_generations for GenARM
  # gene_pool_size: gene pool size
  # initial_attribute_probability: initial probability of any attribute to be included in a rule
  # seed: seed for random number generation
  # weights: weights for fitness function, in order support, confidence,
  #          odds_ratio, repetitions, redundancy
  # crossover_probability: crossover probability
  # parental_inheritance_proportion: proportion of genes from 2 parents (default 50-50)
  # k = tournament size
  # mutation_probability_range = probability of mutation (range)
  # tournament_win_probability_range = probability of win in tournament (range)
  # fitness_mode: function used to account for different objective functions
  #          (simple sum, aritmetic mean, or geometric mean)
  
  fitness_mode <- match.arg(fitness_mode)
  if(!is.null(seed))
    set.seed(seed)
  
  number_of_genes <- ncol(dataframe) - 1
  
  minimum_attribute_values <- sapply(dataframe[-1], min)
  maximum_attribute_values <- sapply(dataframe[-1], max)
  median_attribute_values <- sapply(dataframe[-1], median)
  
  #antecedents matrix determines if an attribute is included in each rule 
  antecedents <- matrix(data=NA, nrow=gene_pool_size, ncol=number_of_genes)
  
  #threshold matrix threshold stores the exposure threshold for each attribute
  thresholds <- matrix(data=NA, nrow=gene_pool_size, ncol=number_of_genes)
  
  #initialize/normalize values 
  for (gene_index in 1:number_of_genes) {
    if (is.logical(dataframe[,gene_index+1])) {
      minimum_attribute_values[gene_index] <- NA
      maximum_attribute_values[gene_index] <- NA
      median_attribute_values[gene_index] <- NA
    }
    
    antecedents[, gene_index] <- as.numeric(runif(gene_pool_size) < initial_attribute_probability)
    
    if (!is.na(minimum_attribute_values[gene_index])) {
      thresholds[, gene_index] <- runif(gene_pool_size, minimum_attribute_values[gene_index], maximum_attribute_values[gene_index])
    }
    else thresholds[, gene_index] <- NaN
  }
  
  fitness_values <- array(data = NA, gene_pool_size)  
  statistics = data.frame(odds_ratio=rep(NA, gene_pool_size), odds_ratio_fitness=rep(NA, gene_pool_size), low_confidence_interval=rep(NA, gene_pool_size), high_confidence_interval=rep(NA, gene_pool_size), p_value=rep(NA, gene_pool_size))
  population <- list(antecedents=antecedents, thresholds=thresholds, fitness_values=fitness_values, statistics=statistics)
  
  #calculate initial population fitness and sort
  population <- help.GenARM.calculate_population_fitness(population, gene_pool_size, number_of_genes, dataframe, weights, mode=fitness_mode)
  population <- help.GenARM.sort_population(population)
  
  help.GenARM.adjust_population_fitness(population, gene_pool_size, weights, sort=TRUE)
  
  for (n in 1:number_of_generations) {
    # Duplicate population (new_population). This only serves the purpose of creating a list
    # of similar shape. Its content is going to be changed to host the new
    # chromosomes generated from population
    new_population <- population
    
    # Compute rules (population) similarity to adjust mutation_probability and tournament_win_probability
    # currently based only on antecedents, not thresholds
    similarity <- help.GenARM.calculate_population_similarity(population$antecedents)
    tournament_win_probability <- help.GenARM.adjust_probability_to_similarity(tournament_win_probability_range, similarity)
    mutation_probability <- help.GenARM.adjust_probability_to_similarity(mutation_probability_range, similarity)
    
    # Crossover and mutation
    # Each crossover generates 2 children. We want gene_pool_size new children,
    # so we iterate gene_pool_size/2 times
    for (population_index in 1:(gene_pool_size/2)) {
      # pick two parents
      parent_one <- help.GenARM.hold_tournament(tournament_size, gene_pool_size, tournament_win_probability)
      parent_two <- help.GenARM.hold_tournament(tournament_size, gene_pool_size, tournament_win_probability)
      
      children <- help.GenARM.generate_children_from_uniform_crossover(population, parent_one, parent_two, crossover_probability, parental_inheritance_proportion)
      children <- help.GenARM.mutate_children(children, mutation_probability, median_attribute_values)
      
      new_population$antecedents[population_index,] <- children[[1]]$antecedents
      new_population$thresholds[population_index,] <- children[[1]]$thresholds
      
      new_population$antecedents[population_index + (gene_pool_size / 2),] <- children[[2]]$antecedents
      new_population$thresholds[population_index + (gene_pool_size / 2),] <- children[[2]]$thresholds
    }
    
    # Measuring fitness of new population
    new_population <- help.GenARM.calculate_population_fitness(new_population, gene_pool_size, number_of_genes, dataframe, weights, mode=fitness_mode)
    
    # Merge, sort and preserve best half
    merged_populations <- mapply(rbind, population, new_population, SIMPLIFY=FALSE)
    merged_populations$fitness_values <- c(population$fitness_values, new_population$fitness_values)
    merged_populations <- help.GenARM.sort_population(merged_populations)
    help.GenARM.adjust_population_fitness(merged_populations, gene_pool_size*2, weights, sort=TRUE)
    
    population$antecedents <- merged_populations$antecedents[1:gene_pool_size,]
    population$thresholds <- merged_populations$thresholds[1:gene_pool_size,]
    population$fitness_values <- merged_populations$fitness_values[1:gene_pool_size]
    population$statistics <- merged_populations$statistics[1:gene_pool_size,]
  }
  
  return(population)
}

GenARM.statistics_report <- function(dataframe, population, gene_pool_size = 30, weights = rep(1,5)) {
  return_values <- population$antecedents * population$threshold
  # convert logical features to true/false
  return_values[is.na(return_values)] <- (population$antecedents[is.na(return_values)] == 1)
  return_values <- cbind(return_values, fitness_values=population$fitness_values)
  number_of_genes <- ncol(population$antecedents)
  
  supports <- array()
  lengths <- array()
  repetition_penalties <- array(NA, gene_pool_size)
  redundancy_penalties <- array()
  statistics <- data.frame(odds_ratio=array(), odds_ratio_fitness=array(), low_confidence_interval=array(), log_high_confidence_interval=array(), p_value=array())
  
  for (population_index in 1:gene_pool_size) {
    # first we binarize a copy of dataframe according to the thresholds
    # of the rule under evaluation
    binarized_dataframe <- dataframe
    for (gene_index in 1:number_of_genes) {
      if (!is.logical(binarized_dataframe[, gene_index+1])) { # logical features do not need to be converted
        binarized_dataframe[, gene_index + 1] <- binarized_dataframe[, gene_index + 1] > population$thresholds[population_index, gene_index]
      }
    }
    
    # then we preserve only the true attributes and the event column
    preserve_columns <- 1 + which(population$antecedents[population_index, ] %in% 1)
    binarized_dataframe <- binarized_dataframe[, c(1, preserve_columns)]
    binarized_dataframe[, 1] <- as.logical(binarized_dataframe[, 1])
    
    # We compute the required rule quality measures
    supports[population_index] <- sum(rowSums(binarized_dataframe) == ncol(binarized_dataframe)) / nrow(binarized_dataframe) # support is no o.f., but it's useful
    lengths[population_index] <- ncol(binarized_dataframe - 1) / number_of_genes
    
    repetition_penalty = help.GenARM.calculate_repetition_penalty(population, population_index, gene_pool_size)
    if (identical(repetition_penalty, numeric(0)))
      repetition_penalties[i] <- 0
    else 
      repetition_penalties[i] <- repetition_penalty
    
    redundancy_penalty = help.GenARM.calculate_redundancy_penalty(population, population_index, gene_pool_size)
    if (identical(redundancy_penalty, numeric(0)))
      redundancy_penalties[population_statistics] <- 0
    else 
      redundancy_penalties[population_statistics] <- redundancy_penalty
    
    statistics[population_index, ] <- help.GenARM.calculate_population_statistics(binarized_dataframe)
  }
  
  return_values <- cbind(rule_index=(1:gene_pool_size), return_values, length=lengths, repetition=repetition_penalties, redundancy=redundancy_penalties, support=supports, statistics)
  
  return(return_values)
}

#' help.GenARM functions
#' 
#' This section of file contains auxillary helper functions for taking care of various calculation
#' functions as well as the modular components of the genetic algorithm itself. 
#' 
#' These functions are used for code simplicity and will likely never be called outside of this file. 
#' All helper function names should begin help.GenARM. to aid contributors to this file in locating them.

help.GenARM.adjust_population_fitness <- function(population, gene_pool_size, weights, sort=TRUE) {
  # Adjust fitness to avoid rule repetitions and redundancy and sort again
  for (population_index in 1:gene_pool_size) {
    repetition_penalty <- help.GenARM.calculate_repetition_penalty(population, population_index, gene_pool_size)
    redundancy_penalty <- help.GenARM.calculate_redundancy_penalty(population, population_index, gene_pool_size)
    
    if(!is.na(repetition_penalty))
      population$fitness_values[ population_index] <- population$fitness_values[ population_index] - weights[3]*repetition_penalty
    if(!is.na(redundancy_penalty))
      population$fitness_values[ population_index] <- population$fitness_values[ population_index] - weights[4]*redundancy_penalty
  }
  if (sort) 
    population <- help.GenARM.sort_population(population)
  return(population)
}

help.GenARM.adjust_probability_to_similarity <- function(probability_range, similarity) {
  return(probability_range[1] + (probability_range[2] - probability_range[1]) * similarity)
}

help.GenARM.calculate_median_distance <- function(thresholds, medians, maximum_possible_distance) {
  distance <- abs(thresholds - medians)/maximum_possible_distance # distance from median normalized by maximum distance possible
  distance <- distance[!is.na(distance)] # drop NaN caused by logical features
  return(sum(distance)/length(distance))   
}

help.GenARM.calculate_population_fitness <- function(population, population_size, number_of_genes, dataframe, weights, mode) {
  for (current_sample in 1:population_size) {
    
    binarized_dataframe <- dataframe
    for (current_gene in 1:number_of_genes) {
      if (!is.logical(binarized_dataframe[,current_gene+1])) {
        binarized_dataframe[,current_gene+1] <- (binarized_dataframe[,current_gene+1] > population$thresholds[current_sample, current_gene])
      }
    }
    
    binarized_dataframe <- binarized_dataframe[,c(1, which(population$antecedents[current_sample, ] %in% 1) + 1)]
    if (is.null(ncol(binarized_dataframe)))
      population$fitness_values[ current_sample] = NaN
    else {
      binarized_dataframe[,1] <- as.logical(binarized_dataframe[,1])
      rule_length <- (ncol(binarized_dataframe) - 1)/number_of_genes
      if (is.na(rule_length))
        rule_length <- 0
      odds_ratio <- odds_ratio.fitness(binarized_dataframe)
      if (is.na(odds_ratio))
        odds_ratio <- 0
      
      if (mode == 'sum') {
        #simple sum
        population$fitness_values[ current_sample] <- (weights[2]*odds_ratio - weights[1]*rule_length)
      }
      else if (mode == 'arithmetic') {
        #arithmetic average
        population$fitness_values[ current_sample] <- (weights[1]*rule_length + weights[2]*odds_ratio)/2
      }
      else if (mode == 'geometric') {
        #geometric average
        population$fitness_values[ current_sample] <- exp((log(weights[1]*rule_length) + log(weights[2]*odds_ratio))) #SHOULD BE /2
      }
      else
        stop("Improper mode parameter given to fitness function")
      
      population$statistics[current_sample,] <- confidence.int(binarized_dataframe)
      # some rules can end up having a "non numerical" fitness, such as
      # NaN or Inf. NaN are harmless, but Inf could bias the algorithm
      # and should be eliminated
      if (!is.na(population$fitness_values[ current_sample]) && !is.nan(populationp$fitness[current_sample]))
        if (population$fitness_values[ current_sample] == Inf)
          population$fitness_values[ current_sample] = NaN
    }
  }
  
  return(population)
}

help.GenARM.calculate_population_similarity <- function(antecedents) {
  similarity <- 2 * abs(colSums(antecedents) / nrow(antecedents) - .5)
  return(sum(similarity) / length(similarity))
}

help.GenARM.calculate_population_statistics <- function(binarized_dataframe) {
  cases <- binarized_dataframe[binarized_dataframe[,1], 2:ncol(binarize_dataframe)]
  controls <- binarized_dataframe[!binarized_dataframe[,1], 2:ncol(binarized_dataframe)]
  
  number_of_cases <- nrow(cases)
  number_of_controls <- nrow(controls)
  
  exposed_cases <- sum(rowSums(cases) == ncol(cases))
  unexposed_cases <- nrow(cases) - exposed_cases
  
  exposed_controls <- sum(rowSums(controls) == ncol(controls))
  unexposed_controls <- nrow(controls) - exposed_controls
  
  odds_ratio <- ((exposed_cases * unexposed_controls) / (unexposed_cases * exposed_controls))
  log_odds_ratio <- log(odds_ratio)
  log_odds_ratio_standard_error <- sqrt(1/x11 + 1/x01 + 1/x10 + 1/x00)
  
  log_low_confidence_interval <- log_odds_ratio - 1.96 * log_odds_ratio_standard_error
  log_high_confidence_interval <- log_odds_ratio + 1.96 * log_odds_ratio_standard_error
  
  low_confidence_interval <- exp(log_low_confidence_interval)
  high_confidence_interval <- exp(log_high_confidence_interval)
  
  odds_ratio_fitness <- NaN
  if (!is.na(low_confidence_interval) && !is.na(high_confidence_interval)) {
    if (low_confidence_interval <= 1 && high_confidence_interval >= 1)
      odds_ratio_fitness <- 0
    else if (odds_ratio < 1)
      odds_ratio_fitness <- 1 - odds_ratio
    else
      odds_ratio_fitness <- (2 * (1 / (1 + exp(1 - odds_ratio)) - .5)) # odds_ratio fitness limited between 0 and 1
    # May add coefficient to increase steepness between 1 and 2 (i.e. -3*(odds_ratio-1))
  }    
  
  proportion_exposed <- (exposed_cases + exposed_controls)/nrow(binarize_dataframe)
  expected_exposed_cases <- proportion_exposed * number_of_cases # assuming independence
  expected_unexposed_cases <- number_of_cases - expected_exposed_cases
  expected_exposed_controls <- proportion_exposed * number_of_controls # assuming independence
  expected_unexposed_controls <- number_of_control - expected_exposed_controls
  
  chi_square <- (exposed_cases - expected_exposed_cases)^2/expected_exposed_cases + (unexposed_cases - expected_unexposed_cases)^2/expected_unexposed_cases + (exposed_controls - expected_exposed_controls)^2/expected_exposed_controls + (unexposed_controls - expected_unexposed_controls)^2/expected_unexposed_controls
  p_value <- pchisq(chi_square, 1, lower.tail=FALSE) 
  
  return(c(odds_ratio, odds_ratio_fitness, low_confidence_interval, log_high_confidence_interval, p_value))
}

help.GenARM.calculate_redundancy_penalty <- function(population, sample_index, gene_pool_size) {
  if (sum(population$antecedents[sample_index,]) < 2)
    return(0) # rules of length 1 can not be redundant
  
  if(!complete.cases(population$statistics[i,]))
    return(1) # maximal penalty for "impossible" rules
  
  cumulative_penalty <- 0
  count <- 0
  
  for (parent_index in 1:gene_pool_size) {
    if (parent_index != sample_index && complete.cases(population$statistics[parent_index,])) {
      difference_vector <- population$antecedents[sample_index,] - population$antecedents[parent_index,]
      if (sum(difference_vector < 0) == 0 && sum(population$antecedents[parent_index, ]) > 0) {
        parent_low_CI <- population$statistics$low_confidence_interval[parent_index]
        parent_high_CI <- population$statistics$high_confidence_interval[parent_index]
        child_low_CI <- population$statistics$low_confidence_interval[sample_index]
        child_high_CI <- population$statistics$high_confidence_interval[sample_index]
        if (sum(difference_vector) == 0) 
          # equal rules are not penalized by this function
          cumulative_penalty <- cumulative_penalty + 0
        else if (child_high_CI < parent_low_CI || child_low_CI > parent_high_CI) 
          # no penalty is added in case of no overlap
          cumulative_penalty <- cumulative_penalty + 0
        else if (child_high_CI > parent_high_CI && child_low_CI < parent_high_CI) { 
          # child covers parent
          cumulative_penalty <- cumulative_penalty + 1
          count <- count + 1
        }
        else if (child_low_CI > parent_low_CI && child_high_CI < parent_high_CI) { 
          # parent covers child (not sure if it is possible)
          cumulative_penalty <- cumulative_penalty + ((parent_high_CI - parent_low_CI) / (child_high_CI - child_low_CI))
          count <- count + 1
        }
        else { 
          # partial overlaps
          if (child_low_CI > parent_low_CI) {
            cumulative_penalty <- cumulative_penalty + ((parent_high_CI - child_low_CI) / (child_high_CI - child_low_CI))
            count <- count + 1
          }
          else {
            cumulative_penalty <- cumulative_penalty + ((child_high_CI - parent_low_CI) / (child_high_CI - child_low_CI))
            count <- count + 1
          }
        }
      }
    }
  }
  if (count == 0)
    return(0)
  else return(cumulative_penalty/count)
}

help.GenARM.calculate_repetition_penalty <- function(population, sample_index, gene_pool_size) {
  penalty_weights <- seq(from=1, to=0, by=(-1 / (gene_pool_size - 1))) # monotonically decreasing
  cumulative_penalty <- 0
  
  for (parent_index in 1:(sample_index - 1)) {
    # check if the chromosomes are equal to those in a previous rule, if so increase the penalty
    if (sum(abs(population$antecedents[sample_index,] - population$antecedents[parent_index,])) == 0) {
      cumulative_penalty <- cumulative_penalty + penalty_weights[parent_index]
    }
  } 
  
  return(cumulative_penalty)
}

help.GenARM.chance_event_occurs <- function(event_probability) {
  return(runif(1) < event_probability)
}

help.GenARM.coin_toss <- function() {
  return(help.GenARM.chance_event_occurs(.5))
}

# Function to perform uniform crossover of two parents
help.GenARM.generate_children_from_uniform_crossover <- function(population, parent_one, parent_two, crossover_probability, parental_inheritance_proportion=.5){
  child_one <- list(antecedents=population$antecedents[parent_one,], thresholds=population$thresholds[parent_one,])
  child_two <- list(antecedents=population$antecedents[parent_two,], thresholds=population$thresholds[parent_two,])
  
  # First it is decided if crossover is performed or not
  if(runif(1) < crossover_probability){
    # If crossover is performed, we go bit by bit
    for(gene_index in 1:length(child_one$antecedents)){
      if(runif(1) < parental_inheritance_proportion) {
        #swap antecedents
        temp <- child_one$antecedents[gene_index]
        child_one$antecedents[gene_index] <- child_two$antecedents[gene_index]
        child_two$antecedents[gene_index] <- temp
        #thresholds = arithmetic average
        child_one$thresholds[gene_index] <- (child_one$thresholds[gene_index]+child_two$thresholds[gene_index])/2
        child_two$thresholds[gene_index] <- (child_one$thresholds[gene_index]+child_two$thresholds[gene_index])/2
      }
    }  
  } 
  return(list(child_one, child_two))
}

# This function generates random indexes to pick chromosomes to
# use to generate new ones. Low indexes are favored, because they 
# represent higher ranking chromosomes. In this implementation, 
# the strongest rule has a probability p to win the tournament. 
# The second strongest rule can win with a probability p*(1-p), 
# the third with a probability p(1-p)^2, and so on.
help.GenARM.hold_tournament <- function(tournament_size, gene_pool_size, tournament_win_probability){
  tournament_lose_probability <- 1 - tournament_win_probability
  competitors <- sort(sample(1:gene_pool_size, tournament_size))
  competitor_win_probability <- 0
  for (competitor_index in 1:tournament_size) {
    competitor_win_probability <- competitor_win_probability + tournament_win_probability*(tournament_lose_probability^(competitor_index - 1))
    if (help.GenARM.chance_event_occurs(competitor_win_probability))
      return(competitors[competitor_index])
  }
  # In the eventuality no gene is selected:
  return(competitors[tournament_size])
}

help.GenARM.mutate_children <- function(children, mutation_probability, median_attribute_values) {
  for(child_index in 1:2){
    for(gene_index in 1:length(children[[child_index]]$antecedents)){
      if(GenARM.chance_event_occurs(mutation_probability)) {
        if (children[[child_index]]$antecedents[gene_index]==1)
          children[[child_index]]$antecedents[gene_index] <- 0
        else children[[child_index]]$antecedents[gene_index] <- 1
        
        change_in_mutation_threshold <- 10 * mutation_probability * median_antecedent_values[gene_index] # percentage of median dependent on pm, scaled by 10 to produce larger changes
        
        if(help.GenARM.coin_toss()) 
          children[[child_index]]$thresholds[i] <- children[[child_index]]$thresholds[i] - change_in_mutation_threshold
        
        else children[[child_index]]$thresholds[i] <- children[[child_index]]$thresholds[i] + change_in_mutation_threshold       
      }
    }
  }
  return(children)
}

help.GenARM.sort_population <- function(population) {
  population_rank <- order(population$fitness_values, decreasing=TRUE)
  population$fitness_values <- population$fitness_values[population_rank]
  population$antecedents <- population$antecedents[population_rank, ]
  population$thresholds <- population$thresholds[population_rank, ]
  population$statistics <- population$statistics[population_rank, ]
  return(population)
}