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
  