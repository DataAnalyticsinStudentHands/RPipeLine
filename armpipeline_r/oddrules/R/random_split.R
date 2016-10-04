#' Random Split
#' 
#' Splits data into testing and training sets
#' @usage random_split(sampleLog, seed=1, split=5000)
#' @param sampleLog The data to be split into training and testing
#' @param split The number of rows to be put into the training set. The rest of the data will be put into the testing set.
#' @param seed The seed used for the random split. Default is 1.
#' @return return$whole_training: data sorted into training
#' @return return$testing: data sorted into testing group
#' @export
random_split <- function (sampleLog, seed=1, split=5000) {
  
  if (nrow(sampleLog) <= split){
    stop("Not enough data for split")
  }
  
  #First the data must be reordered in a random order
  set.seed(seed)
  rand <- sample(nrow(sampleLog))
  sampleLog <- sampleLog[rand, ]
  
  #Then it must be split in testing or training
  # First Rows go to Training
  whole_training <- sampleLog[1:split, ]
  # Remaining rows go in the testing set
  testing <- sampleLog[(split+1):nrow(sampleLog), ]
  return(list(whole_training=whole_training,testing=testing))
}
