merge_poll <- function(C, E, labels, dist = 20) {
	# This function is a particular form of "merge" built for subjects and pollutants. It is more complexe than regular merge because one tag is not sufficiant to link a subjects to exposures, and exposures are listed differently.
	
	# C in input should have 4 columns: event, DOS, latitude, longitude.
	# E in input should have 7 columns: poll_code, site, latitude, longitude, epoch, hour, value.
# labels: poll_code, name. A name should be available for every pollutant in E.

	# First, an output database (Out) is created based on the column in C and with 7 new columns for every pollutant: site_name, hour0_name, val0_name, hour1_name, val1_name, hour2_name, val2_name, where name is the label for that pollutant and 0-1-2 indicates how many days before the event
	
  allcodes <- (unique(E$poll_code))
  
	Out <- C
	for (code in allcodes) {
		ind <- match(code, labels$poll_code)
		label <- labels[ind, 2]
		# site <- paste("site", label, sep = "_")
		h0 <- paste("hour0", label, sep = "_")
		v0 <- paste("val0", label, sep = "_")
		h1 <- paste("hour1", label, sep = "_")
		v1 <- paste("val1", label, sep = "_")
		h2 <- paste("hour2", label, sep = "_")
		v2 <- paste("val2", label, sep = "_")
		h3 <- paste("hour3", label, sep = "_")
		v3 <- paste("val3", label, sep = "_")
		h4 <- paste("hour4", label, sep = "_")
		v4 <- paste("val4", label, sep = "_")
		newcolumns <- c(h0, v0, h1, v1, h2, v2, h3, v3, h4, v4)
		Out[, newcolumns] <- NA
	}
	
# 	# A table including sites numbers and coordinates is extracted from E for later use
# 	n <- 1
# 	u <- unique(E$site)
# 	u <- sort(u)
# 	sitetable <- data.frame(NA, NA, NA)
# 	colnames(sitetable) <- c('site', 'latitude', 'longitude')
# 	for (i in u) {
# 		ind <- match(i,E$site)foreach(i=1:3) %dopar% sqrt(i)

# 		sitetable[n,] <- E[ind,2:4]
# 		n <- n + 1
# 	}

	# New cycle over the poll_codes to fill the columns 
  strt<-Sys.time()	

  f <- factor(E[,1])
  splitted <- split(E,f,drop=FALSE)

  lag0 <- foreach (n=1:nrow(Out)) %dopar% {
    foreach (code=allcodes)  %dopar%{
      exposure.index(C[n,],code,allcodes,splitted,dist,0)
    }
  }

  lag1 <- foreach (n=1:nrow(Out)) %dopar% {
    foreach (code=allcodes)  %dopar%{
      exposure.index(C[n,],code,allcodes,splitted,dist,1)
    }
  }

  lag2 <- foreach (n=1:nrow(Out)) %dopar% {
    foreach (code=allcodes)  %dopar%{
      exposure.index(C[n,],code,allcodes,splitted,dist,2)
    }
  }

  lag3 <- foreach (n=1:nrow(Out)) %dopar% {
    foreach (code=allcodes)  %dopar%{
      exposure.index(C[n,],code,allcodes,splitted,dist,3)
    }
  }

  lag4 <- foreach (n=1:nrow(Out)) %dopar% {
    foreach (code=allcodes)  %dopar%{
      exposure.index(C[n,],code,allcodes,splitted,dist,4)
    }
  }
#   foreach (n=1:nrow(Out))  %dopar% {
# 		out_col <- 5	#index for output column
# 			foreach (code=allcodes)  %dopar%{
# 
# 			flag <- TRUE	#marks the beginning of a new pollutant
# 			index <- 0
# 
# 			#find pollutant in E
# 			dd <- match(42401,allcodes)
# 			if(!is.na(dd)){
# 				for (day in 0:2) {
# 					#find date
# 					d <- (splitted[[dd]]$epoch==(C[n,2]-day))
# 					reduced2 <- splitted[[dd]][d,]
#           localdist <- dist
#           if (nrow(reduced2)>0){
#             point1 <- c(C[n,4],C[n,3])
#             for (j in 1:nrow(reduced2)) {            
#               point2 <- c(reduced2[j,4],reduced2[j,3])
#               temp <- geodetic.distance(point1,point2)
#               if(!is.nan(temp)){
#                 if(temp<localdist){
#                   localdist <- temp
#                   index <- j
#                 }
#               }
#             }
#           }
# 					
#           if(index!=0){
# 					#write values in Out
# 					if(flag){
# 						Out[n,out_col] <- reduced2[index,2]
# 						flag <- FALSE
# 						out_col <- out_col + 1
# 					}
# 					Out[n,out_col] <- reduced2[index,6]
# 					out_col <- out_col + 1
# 					Out[n,out_col] <- reduced2[index,7]
# 					out_col <- out_col + 1
# 				}
# 				}
# 			}
# 		}
# 	}
  
  print(Sys.time()-strt)

  # Use indexes lag0-2 to copy from E to Out
  for (i in 1:nrow(Out)) {
    for (j in 1:length(allcodes)) {
      exprow0 <- as.integer(lag0[[i]][j])
      exprow1 <- as.integer(lag1[[i]][j])
      exprow2 <- as.integer(lag2[[i]][j])
      exprow3 <- as.integer(lag3[[i]][j])
      exprow4 <- as.integer(lag4[[i]][j])
      if(!is.na(exprow0)) {
      Out[i,j*9-4] <- E$hour[exprow0]
      Out[i,j*9-3] <- E$value[exprow0]
      }
      if(!is.na(exprow1)) {
      Out[i,j*9-2] <- E$hour[exprow1]
      Out[i,j*9-1] <- E$value[exprow1]
      }
      if(!is.na(exprow2)) {
      Out[i,j*9] <- E$hour[exprow2]
      Out[i,j*9+1] <- E$value[exprow2]
      }
      if(!is.na(exprow3)) {
        Out[i,j*9+2] <- E$hour[exprow3]
        Out[i,j*9+3] <- E$value[exprow3]
      }
      if(!is.na(exprow4)) {
        Out[i,j*9+4] <- E$hour[exprow4]
        Out[i,j*9+5] <- E$value[exprow4]
      }
    }
  }

  print(Sys.time()-strt)
	return(Out)

} # end function

exposure.index <- function(Crow, code, allcodes, E, dist, day) {
  index <- 0
  # match code to go in right branch of E
  dd <- match(code,allcodes)
  if(!is.na(dd)) {
    # find correct days
    d <- (E[[dd]]$epoch==Crow$DOS-day)
    reduced2 <- E[[dd]][d,]
    localdist <- dist
    
    if (nrow(reduced2)>0){
      point1 <- c(Crow$longitude,Crow$latitude)
      for (j in 1:nrow(reduced2)) {            
        point2 <- c(reduced2$longitude[j],reduced2$latitude[j])
        #compute and find minimum distance
        temp <- geodetic.distance(point1,point2)
        if(!is.nan(temp)){
          if(temp<localdist){
            localdist <- temp
            index <- j
          }
        }
      }
    }
  }
  return(as.integer(row.names(reduced2)[index]))
}

#The following program computes the distance on the surface of the earth between two points point1 and point2. Both the points are of the form (Longitude, Latitude)
#Output in km
geodetic.distance <- function(point1, point2)
{
  R <- 6371
  p1rad <- point1 * pi/180
  p2rad <- point2 * pi/180
  d <- sin(p1rad[2])*sin(p2rad[2])+cos(p1rad[2])*cos(p2rad[2])*cos(abs(p1rad[1]-p2rad[1]))  
  d <- acos(d)
  R*d
}

