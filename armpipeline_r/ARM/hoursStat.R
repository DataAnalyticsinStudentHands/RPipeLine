hours.stat <- function(subjects,tags) {
  
  # Function to compare peak hour distribution of a group of cases to the global
  # distribution for that pollutant.
  # For a more generalize information, the data frame "subjects" should be the complete one
  # (binned or not) and not the one including only the complete cases
  
  # Global distribution of single pollutants
  CO <- cbind(c(subjects$hour0_CO, subjects$hour1_CO, subjects$hour2_CO),c(subjects$val0_CO, subjects$val1_CO, subjects$val2_CO))
  SO2 <- cbind(c(subjects$hour0_SO2, subjects$hour1_SO2, subjects$hour2_SO2),c(subjects$val0_SO2, subjects$val1_SO2, subjects$val2_SO2))
  NO <- cbind(c(subjects$hour0_NO, subjects$hour1_NO, subjects$hour2_NO),c(subjects$val0_NO, subjects$val1_NO, subjects$val2_NO))
  NO2 <- cbind(c(subjects$hour0_NO2, subjects$hour1_NO2, subjects$hour2_NO2),c(subjects$val0_NO2, subjects$val1_NO2, subjects$val2_NO2))
  O3 <- cbind(c(subjects$hour0_O3, subjects$hour1_O3, subjects$hour2_O3),c(subjects$val0_O3, subjects$val1_O3, subjects$val2_O3))
  PM <- cbind(c(subjects$hour0_PM, subjects$hour1_PM, subjects$hour2_PM),c(subjects$val0_PM, subjects$val1_PM, subjects$val2_PM))
  
  # Preserving only high pollution days
  ind <- (CO[,2]==1)
  CO <- CO[ind,]
  CO <- CO[complete.cases(CO),]
  ind <- (SO2[,2]==1)
  SO2 <- SO2[ind,]
  SO2 <- SO2[complete.cases(SO2),]
  ind <- (NO[,2]==1)
  NO <- NO[ind,]
  NO <- NO[complete.cases(NO),]
  ind <- (NO2[,2]==1)
  NO2 <- NO2[ind,]
  NO2 <- NO2[complete.cases(NO2),]
  ind <- (O3[,2]==1)
  O3 <- O3[ind,]
  O3 <- O3[complete.cases(O3),]
  ind <- (PM[,2]==1)
  PM <- PM[ind,]
  PM <- PM[complete.cases(PM),]
  
  H <- list(CO=CO[,1],SO2=SO2[,1],NO=NO[,1],NO2=NO2[,1],O3=O3[,1],PM=PM[,1])

  ind <- complete.cases(subjects)
  C <- subjects[ind,]
  
  # Preserving only cases
  ind <- C$event == 1
  C <- C[ind,]
  
  s <- rep(0, length(tags))
  t <- rep(0, length(tags))
  n <- 1
  for (name in tags) {
    s[n] <- match(paste("hour",name,sep=""), colnames(C))
    t[n] <- match(paste("val",name,sep=""), colnames(C))
    n <- n + 1
  }
  
  c <- C[,t] # subdatabase of values (binary) 
  h <- C[,s] # subdatabase of hours 
  
  ind <- apply(c, 1, function(row) all(row !=0 ))
  c <- c[ind,]
  h <- h[ind,] # database including only peak hours related to the rule
  
  # removing *_ from tags
  tags <- gsub("^.*?_","",colnames(h))
  H <- H[match(tags, names(H))]
  
  boxplot(c(h,H))
  
  t <- rep(0, ncol(h))
  for (i in 1:ncol(h))
    t[i] <- wilcox.test(h[,i],H[[i]])[["p.value"]]
  
  return(list(rule_hours = h, original_hours = H, p_value = t))  
}