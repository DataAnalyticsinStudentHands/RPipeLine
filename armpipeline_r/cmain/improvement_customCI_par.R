improvement_customCI_par <- function(rules) {
  
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
  registerDoParallel(4)# causes connections to open

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


