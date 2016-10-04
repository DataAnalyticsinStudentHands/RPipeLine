weekday.stat <-function(subjects, tags) {
  
  # This function finds the group of cases connected to the rule specified by "tags"
  # and shows the distribution of weekdays on which they visited the ER. This way it
  # is possible to see if some days are more "responsible" for the outcome than others
  
  # INPUT
  # subjects: binned dataset of complete cases used to extract the rules
  # tags: pollutants included in the rules, in the form  n_XX where n (0:4) 
  # indicates the day lag and XX is the pollutant chemical formula (O3, PM...). 
  # An example of a valid entry is c("0_O3","1_CO").
  
  # Although the input should have only complete cases, the function starts by filtering 
  # out non-complete cases, in case of user mistakes
  ind <- complete.cases(subjects)
  C <- subjects[ind,]
  
  # Preserving only cases
  ind <- C$event == 1
  C <- C[ind,]
  
  t <- rep(0, length(tags))
  n <- 1
  for (name in tags) {
    t[n] <- match(paste("val",name,sep=""), colnames(C))
    n <- n + 1
  }
  
  c <- as.data.frame(C[,t]) # subdatabase of values (binary) 
  
  # Preserving only exposed cases
  ind <- apply(c, 1, function(row) all(row !=0 ))
  
  # Subset of exposed cases DOS
  dates <- C$DOS[ind]
  
  # Convert Matlab format to Dates format
  dates <- dates - 1
  dates <- as.Date(dates, "0000-01-01")
  
  # weekdays barplot
  labels <- format(seq(ISOdate(2010,7,11), ISOdate(2010,7,17), "24 hours"),"%a")
  barplot(table(format((dates), "%w")), names.arg=labels, main=paste(c(tags, "DOS"), sep=" ", collapse=" "))
  
}