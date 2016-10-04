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
