#' dataload.R
#' 
#' This file contains all necessary functions for preparing data stored in data frame objects. 
#' These are functions that alter data stored in the data frames in preparation for intended
#' analysis, since many types of analysis will share common needs of this type; These functions 
#' are not themselves the functions that are specific to one particular type of analysis.
#' General categories of these functions include: data subselection/pruning, data normalization, 
#' data format conversions, and data splitting. In general these functions will require a dataframe 
#' object as their first parameter and will return a modified dataframe (or list of dataframes) as 
#' their first return value.
 
#' prepare.binarize_dataframe(dataframe, quantile_threshold, column_list, group_size, group_columns_to_use, invert_binarization) {
#'
#' Single-purpose function to convert quantitative data into binary data. This function works by 
#' analyzing the quantitative data in quantiles, using a set quantile threshhold to binarize the
#' original data into TRUE (above threshhold) and FALSE (below threshhold) values. These values can
#' be inverted by setting invert_binarization=TRUE. This function can be applied to a specified 
#' subset of the columns in the dataframe by passing a list of columns (by name or index) in the 
#' column_list parameter. Additionally, columns can be grouped to be analyzed together to produce a 
#' quantile threshhold shared among them. If certain columns are to be ignored in the calculation
#' of the group threshhold, the group_columns_to_use parameter takes a list of which columns to use.
#' The columns not specified will still be binarized according to the group threshhold that was
#' calculated without them. By default, the column group_size is 1, meaning that every column's 
#' threshhold is analyzed individually and binarized accordingly.
#'    
#' @param dataframe             Parameter containing dataframe object
#' @param quantile_threshold     Optional quantile threshhold [0.75]
#' @param column_list            Optional list of columns to binarize [1:ncol(dataframe)]
#' @param group_size             Optional size of column groups [1]
#' @param group_columns_to_use   Optional masking of column groups [1:group_size]
#' @param invert_binarization    Optional parameter to invert TRUE and FALSE binarization [FALSE]
#'
#' @return data     dataframe with quantitative data binarized
#' @return output   quantitative value for each column (or group) at the given quantile threshhold
#'
#' @examples
prepare.binarize_dataframe <- function(dataframe, quantile_threshold = 0.75, column_list = 1:ncol(dataframe), group_size = 1, group_columns_to_use = 1:group_size, invert_binarization = FALSE) {
  normalized_column_list <- help.prepare.normalize_and_validate_column_list(dataframe, column_list)
	
  if (length(normalized_column_list)%%group_size != 0) {
    stop("Impossible to group columns as requested.")
  }
  
  if (sum(group_columns_to_use > group_size) > 0) {
    stop("Impossible to calculate threshold as requested")
  }
  
  group_thresholds <- array(dim = length(normalized_column_list)/group_size)
  index_of_current_group <- 1
  
  for (group_column_offset in seq(from=0, to=length(normalized_column_list)-1, by=group_size)) {
	  values_of_current_group <- NULL
	  for (column_iterator in (group_columns_to_use + group_column_offset)) {
		  index_of_current_column <- normalized_column_list[column_iterator]
		  values_of_current_group <- cbind(values_of_current_group, dataframe[,index_of_current_column])
	  }
	  group_thresholds[index_of_current_group] <- quantile(values_of_current_group, quantile_threshold, na.rm = TRUE)
	  for (column_iterator in seq(from=group_column_offset, to=group_column_offset+group_size-1)) {
		  index_of_current_column <- normalized_column_list[j]
		  binarized_data <- (dataframe[,index_of_current_column] >= group_thresholds[index_of_current_group])
		  dataframe[which(binarized_data == TRUE), index_of_current_column] <- (TRUE != invert_binarization)
		  dataframe[which(binarized_data == FALSE), index_of_current_column] <- (FALSE != invert_binarization)
	  }
	  index_of_current_group <- index_of_current_group + 1
  }
  
  return(list(dataframe = prepare.convert_dataframe_to_type(dataframe, 'logical'), output = group_thresholds))
}

prepare.convert_dataframe_to_type <- function(dataframe, to="logical") {
  if (to == "logical")
    return(data.frame(lapply(dataframe, as.logical)))
  else 
    return(dataframe)
}

#' prepare.prune_dataframe(dataframe, omit_columns, omit_rows, omit_NA_cases)
#'
#' General pruning function that can remove rows and columns from a dataframe. Once the 
#' rows and columns have been removed, the function can then optionally remove any rows that
#' include any NA values. 
#'
#' @param dataframe 
#'
#' @return dataframe with rows omitted
#'
#' @examples
prepare.prune_dataframe <- function(dataframe, omit_rows=c(), omit_columns=c(), omit_NA_rows=FALSE) {
	normalized_column_list <- help.prepare.normalize_and_validate_column_list(dataframe, omit_columns)
	if (is.numeric(omit_rows)) {
	  row_mask <- array(TRUE, dim=nrow(dataframe))
	  row_mask[omit_rows] <- FALSE
		dataframe <- dataframe[row_mask,]
	}	
	if (is.numeric(normalized_column_list)) {
	  col_mask <- array(TRUE, dim=ncol(dataframe))
	  col_mask[normalized_column_list] <- FALSE
	  dataframe <- dataframe[,col_mask]
	}
	if (omit_NA_rows) {
		dataframe <- na.omit(dataframe)
	}
	return(dataframe)
}

prepare.split_dataframe_for_training_and_testing <- function(dataframe, sampling_method='random', training_size=.5, seed=1) {
	numTraining <- training_size
	numTotal <- nrow(dataframe)
	if (training_size <= 1 && training_size >= 0) {
		numTraining <- as.int(training_size * numTotal)
	}
	
	if (numTotal <= numTraining) {
		stop("Not enough data to split into training and testing samples")
	}
	
	if (sampling_method == 'random') {
		set.seed(seed)
		randomRowOrder <- sample(nrow(dataframe))
		dataframe <- dataframe[randomRowOrder, ]
	}
	
	return (list(training=dataframe[1:numTraining], testing=dataframe[(numTraining+1):numTotal,]))
}

#' help.prepare functions
#' 
#' This section of file contains auxillary helper functions for taking care of basic tasks needed at various 
#' points in data preparation functions, mostly for parameter normalization to aid in the versatility
#' of the code. 
#' 
#' These functions are used for code simplicity and will likely never be called outside of this file. 
#' All helper function names should begin help.prepare. to aid contributors to this file in locating them.


#' help.prepare.normalize_and_validate_column_list(dataframe, column_list)
#' 
#' A helper function that accepts a dataframe and a list of columns either specified by names or by
#' indexes, and returns a list of indexes. Additionally, the function validates that the indexes and
#' names exist in the given dataframe. This helper function can be used in any data preparation
#' function to allow it to accept different column specifiers.
#'
#' @param dataframe    a data frame object
#' @param column_list   a list of strings or integers
#'
#' @return a list of integers corresponding to the indexes of the columns specified in column_list
#'
#' @examples
help.prepare.normalize_and_validate_column_list <- function(dataframe, column_list) {
  normalized_column_list <- column_list
  if (is.character(normalized_column_list)) {
    normalized_column_list <- match(column_list, colnames(dataframe))
    if (any(is.na(column_list))) 
      stop("One or more column tag not found.")
  }
  if (is.numeric(normalized_column_list)) {
    if (any(normalized_column_list > ncol(dataframe))) {    
      stop("Index exceeds data frame dimensions.")
    }
  } else if (length(normalized_column_list) > 0) {
    stop("Column parameter not appropriately designated by indexes or tag names")
  }
  return(normalized_column_list)
}

