nonredundant <- function(rules, columnnames) {
  
  rules=unique(rules)
  rules$rules=as.character(rules$rules)
  nonredundant <- logical(length(rules$rules))
  nonredundant[1:length(nonredundant)] <- TRUE
  
  #Check to see if a variable is used in a rule
  for (i in 1:length(columnnames)){
    b=vector()
    a=c(grep(paste(columnnames[i],"}",sep=""),rules$rules),grep(paste(columnnames[i],",",sep=""),rules$rules))
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
              if (rules$lowCIcustom[a[which.min(b)]] <= rules$lowCIcustom[a[k]] & rules$lowCIcustom[a[k]] <= rules$highCIcustom[a[which.min(b)]]){
                nonredundant[a[k]]=FALSE
              }
              if (rules$highCIcustom[a[which.min(b)]] >= rules$highCIcustom[a[k]] & rules$highCIcustom[a[k]] >= rules$lowCIcustom[a[which.min(b)]]){
                nonredundant[a[k]]=FALSE
              }
              if (rules$lowCIcustom[a[which.min(b)]] >= rules$lowCIcustom[a[k]] & rules$highCIcustom[a[k]] >= rules$highCIcustom[a[which.min(b)]]){
                nonredundant[a[k]]=FALSE
              }
            }
            
          }
          
          
        }
        b[which.min(b)]=NaN
        
      }
      
    }
  
  }
return(nonredundant)
}
