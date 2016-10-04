#' Quantile binarization
#' 
#' The function bin_quantile() allows to transform numeric columns of data frames in logical values based on a threshold. The threshold is computed from the distribution of the numeric values in the column(s), divided in quantiles as specified by the user.
#' @usage bin_quantile(D, cols = 1:ncol(D), top = 0.75, group = 1, groupthreshold=1:group, reverse = FALSE, dropNA = TRUE)
#' @param D A data frame
#' @param  cols The indexes or names of the columns the users wants to binarize. Default = all the columns
#' @param  top The probability of the quantile (in the interval [0,1]). For example, when top = 0.75, the bottom 75 percent of the values in the column will be tagged as FALSE, and the remaining 25 percent will be tagged as TRUE
#' @param  group The interval of columns that should be processed together. Default = 1
#' @param  groupthreshold What columns of each group to use to find the threshold. For example, when groupthreshold=3, the third column of each group will be used to determine the threshold of each group. Default = use all columns of a group
#' @param  reverse A logical parameter; when reverse=TRUE, the bottom part of the value is considered TRUE, and the top is FALSE. Default = FALSE
#' @param dropNA  A logical parameter: when dropNA=TRUE output rows containing at least one NA are deleted. Default = TRUE
#' @return test$data: database with binarized listed columns
#' @return test$thresholds: thresholds used for the binarization (they follow the column order)
#' @export
bin_quantile <- function (D, cols = 1:ncol(D), top = 0.75, group = 1, groupthreshold=1:group, reverse = FALSE, dropNA = TRUE) {
  
  # This function performs binarization by column.
  # Inputs:
  # D: input data frame
  # cols: column indexes or names to binarize
  # top: threshold above which data should be considered TRUE (default = 0.75, for top quartile)
  # group: number of columns that should be evaluated together (same threshold). Default = 1
  # groupthreshold: what columns of each group to use to find the threshold. Default = all columns of a group
  # reverse: option to label as TRUE values below the threshold (default = FALSE)
  
  # If cols has been passed as column names instead of indexes, the indexes are found by the function
  if (is.character(cols)) {
    cols <- match(cols, colnames(D))
    if (any(is.na(cols))) 
      stop("One or more column tag not found.")
  }
  
  if (is.numeric(cols)) {
    if (any(cols>ncol(D))) 
      
      
      stop("Index exceeds data frame dimensions.")
  }
  
  if (length(cols)%%group != 0) {
    stop("Impossible to group columns as requested.")
  }
  
  if (sum(groupthreshold>group)>0) {
    stop("Impossible to calculate threshold as requested")
  }
  
  # index for threshold array
  n <- 1
  threshold <- array(dim = length(cols)/group)
  
  # Binarization cycle (by columns listed in cols)
  for (i in seq(from=1, to=length(cols), by=group)) {
    col_group <- NULL
    for (j in groupthreshold+i-1) {
      index <- cols[j]
      col_group <- cbind(col_group,D[,index])
    }
    threshold[n] <- quantile(col_group,top,na.rm = TRUE)  
    for (j in seq(from=i, to=i+group-1)) {
      index <- cols[j]
      ind <- (D[,index]>=threshold[n])
      q <- which(ind)
      w <- which(!ind)
      D[q,index] <- TRUE
      D[w,index] <- FALSE
      if (reverse) {
        D[,index]=!D[,index]
      }
    }
    n <- n + 1
  }
  
  if (dropNA) {
    D <- D[complete.cases(D),]
  }
  
  return(list(data=D,thresholds=threshold))
}