or.table <- function(rules, poll, lag, sortby = 0, file="ORplot.png") {
  
  # This function organizes rules (from a data frame) in tables according to the
  # specified pollutant, such as:
  #
  # O3 table
  # 
  # single  <rule1> <rule2> ... <rulen>
  # O3_0    val0    val0    ... val0
  # O3_1    val1    val1    ... val1
  # O3_2    val2    val2    ... val2
  #
  # These tables are meant to be used to produce readable graphs of the rules
  
  # convert rules from factor to character if necessary
  if (class(rules$rules)=="factor")
    rules$rules <- as.character(rules$rules)
  
  # clean data frame "rules" from useless characters
  for (i in 1:nrow(rules)) {
    rules$rules[i] <- gsub("event", "", rules$rules[i])
    rules$rules[i] <- gsub("val", "", rules$rules[i])
    rules$rules[i] <- gsub("([{=>}])", "", rules$rules[i])
  }
  
  # separate single rules from the rest
  ind <- grep(",", rules$rules)
  singles <- rules[-ind,]
  rules <- rules[ind,]
   
  # order single rules and find single OR for the required pollutant
  singles <- singles[order(singles$rules),]
  singleOR <- rep(0,lag)
  
  for (j in 0:lag) {
  temp <- singles$oddsRatio[grep(paste(j,"_",poll,sep=""), singles$rules)]
  if (length(temp)>0)
    singleOR[j+1] <- temp
  else singleOR[j+1] <- NA
  }
  
  # start building the data frame
  ORtable <- data.frame(singleOR)
  
  r <- grep(poll,rules$rules)
  
  for (i in r) {
    # each rule is unpacked and rebuilt to have the rule for 0_poll, 1_poll and 2_poll
    str <- rules$rules[i]
    str <- gsub(" ","",str) #removing spaces
    str <- unlist(strsplit(str,","))
    or <- rep(0,lag+1)
    for (j in 0:lag) {
      str0 <- str
      str0[grep(poll,str)[1]] <- paste(j,"_",poll,sep="")
      if (length(unique(str0))==length(str0)) {
      tab <- rules
      for (x in str0) {
        ind <- grep(x,tab$rules)
        tab <- tab[ind,]
      } 
      if (length(tab$oddsRatio)>0) {
        # count "," in rule to be sure to have the one of right length
        s2 <- gsub(",","",tab$rules)
      if (nchar(tab$rules) - nchar(s2) < length(str0))
        or[j+1] <- tab$oddsRatio[1]
      else or[j+1] <- NA
      }
      else or[j+1] <- NA
      }
      else or[j+1] <- NA 
    }
    ORtable <- cbind(ORtable,or)
    name<-str0
    name<-name[-grep(poll,name)[1]]    
    colnames(ORtable)[ncol(ORtable)] <- paste(name,collapse=",")
  }

  ORtable <- ORtable[,unique(colnames(ORtable))]
  ORtable <- t(ORtable)
  ORtable <- ORtable[order(ORtable[,sortby+1
                                   ]),] # sort table
  colnames(ORtable) <- c(paste("0_",poll,sep=""),paste("1_",poll,sep=""),paste("2_",poll,sep=""),paste("3_",poll,sep=""),paste("4_",poll,sep=""))
  
  png(filename=file)
  op <- par(mar = c(10,4,3,2) + 0.1)
  plot(ORtable[,1],type="b",lwd=2, ann=FALSE, xaxt="n",ylim=c(0.5,1.5),col="black",ylab="OR",main=paste("OR for combinations with",poll))
  axis(1,at=1:length(rownames(ORtable)),labels=rownames(ORtable),las=3)
  lines(ORtable[,2],col="red",type="b",lwd=2)
  lines(ORtable[,3],col="green",type="b",lwd=2)
  lines(ORtable[,4],col="purple",type="b",lwd=2)
  lines(ORtable[,5],col="blue",type="b",lwd=2)
  legend("bottomright",legend=c(paste("0_",poll,sep=""),paste("1_",poll,sep=""),paste("2_",poll,sep=""),paste("3_",poll,sep=""),paste("4_",poll,sep="")),
         border = "black",lty=1,lwd=2,pch=21,col=c("black","red","green","purple","blue"),horiz=TRUE,bty="n",cex=0.8,
         text.col=c("black","red","green","purple","blue"),
         inset=0.01)
  grid()
  par(op)
  dev.off()
  
  return(ORtable)
}