improvement <- function(rules, t = 0.1) {
  
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