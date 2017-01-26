TCEQ.daily_aggregate <- function(dataframe, aggregate_method="max", minimum_number_of_measurements=18) {
  aggregate_by <- list(day=cut(dataframe$time, breaks="day"), site=dataframe$site, region=dataframe$region, param_id=dataframe$param_id)
  
  aggregated_dataframe <- aggregate(dataframe$value, aggregate_by, aggregate_method)
  aggregated_dataframe$number_of_measurements <- aggregate(dataframe$value, aggregate_by, "length")$x
  
  # aggregated_dataframe$hour <- dataframe$time[match(aggregated_dataframe$x, dataframe$value)]
  aggregated_dataframe$value <- help.TCEQ.normalize_numeric_values(aggregated_dataframe$x)
  
  aggregated_dataframe[which(aggregated_dataframe$number_of_measurements <= minimum_number_of_measurements),"value"] <- NA
  return(aggregated_dataframe)  
}

TCEQ.hourly_aggregate <- function(dataframe, aggregate_method="mean", minimum_number_of_measurements=8) {
  aggregate_by <- list(hourofday=cut(dataframe$time, breaks="hour"), site=dataframe$site, region=dataframe$region, param_id=dataframe$param_id)
  
  aggregated_dataframe <- aggregate(dataframe$value, aggregate_by, aggregate_method)
  aggregated_dataframe$number_of_measurements <- aggregate(dataframe$value, aggregate_by, "length")$x
  
  aggregated_dataframe$time <- help.TCEQ.normalize_datetime_values(aggregated_dataframe$hourofday)
  aggregated_dataframe$hourofday <- as.numeric(format(aggregated_dataframe$time, format="%H"))
  aggregated_dataframe$value <- help.TCEQ.normalize_numeric_values(aggregated_dataframe$x)
  
  aggregated_dataframe[which(aggregated_dataframe$number_of_measurements <= minimum_number_of_measurements),"value"] <- NA
  return(aggregated_dataframe)
}

TCEQ.rolling_average <- function(value_vector, rolling_length=8, maximum_allowed_missing=2) {
  average_vector <- rep(NA, length(value_vector))
  
  for (value_index in 1:length(value_vector)) {
    number_missing <- 0
    value_total <- 0
    for (rolling_offset in 0:(rolling_length-1)) {
      if (is.na(value_vector[value_index + rolling_offset]))
        number_missing <- number_missing + 1
      else {
        value_total <- value_total + value_vector[value_index + rolling_offset]
      }
    }
    if (number_missing <= maximum_allowed_missing)
      average_vector[value_index] <- (value_total + .1 * number_missing) / rolling_length
  }
  
  return(average_vector)
}

TCEQ.load_aqs_sites <- function(filename) {
  sites <- read.csv(filename)
  return (sites)
}
TCEQ._sites <- TCEQ.load_aqs_sites("./prepare_libs/TCEQ_AQS_sites.csv")

TCEQ.sites <- function(region = NA) {
  sites <- TCEQ._sites
  if (!is.na(region)) {
    sites <- TCEQ._sites[which(sites$Reg. == region),]
  }
  return(sites)
}

#converter
help.TCEQ.normalize_source_data <- function(dataframe, omit_NA_rows=TRUE, omit_flagged_rows=TRUE) {
  dataframe$time <- help.TCEQ.normalize_datetime_values(dataframe$epoch)
  dataframe$value <- help.TCEQ.normalize_numeric_values(dataframe$value)
  if (omit_NA_rows)
    dataframe <- dataframe[complete.cases(dataframe), ]
  if (omit_flagged_rows)
    dataframe <- dataframe[which(dataframe$flag == ""), ]
  return(dataframe)
}

help.TCEQ.normalize_datetime_values <- function(datetime_values) {
  return(as.POSIXct(datetime_values, origin="1970-01-01", tz="UTC"))
}

help.TCEQ.normalize_numeric_values <- function(numeric_values) {
  return(as.numeric(as.character(numeric_values)))
}

# This function is the first block of the ARM pipeline. 
# It expands the matrix of cases (C) adding exposures from the matrix (E) according to a common parameter (for example, a column "date").

# C in input should have 2 columns, the first one being cases/controls.
# E can have an arbitrary number of columns. The column of the common parameter that the user wishes to use for merging must contain unique values.

# IMPORTANT: the class of vectors in the database are important for the correct functioning of the function mergeCE. To be sure to preserve the class when importing the data from a csv file, use the command read.csv("my_exposures.csv",stringsAsFactors=FALSE). 

# If m and/or n have been passed as column names instead of index, the index is found by the function
help.TCEQ.merge_cases_and_exposures <- function(case_dataframe, exposure_dataframe, case_match_column=2, exposure_match_column=1, drop_match_column=FALSE) {
  case_match_column <- help.prepare.normalize_and_validate_column_list(case_dataframe, case_match_column)
  exposure_match_column <- help.prepare.normalize_and_validate_column_list(exposure_dataframe, exposure_match_column)
  
  if (length(unique(exposure_dataframe[, exposure_match_column])) < nrow(exposure_dataframe))
    stop("parameter column in exposure matrix contains non-unique values.")
  
  if (class(exposure_dataframe[, exposure_match_column]) != class(case_dataframe[, case_match_column])) 
    stop("selected parameter columns are of different data type.")
  
  number_of_case_features <- ncol(case_dataframe)
  number_of_cases <- nrow(case_dataframe)
  exposure_feature_labels <- colnames(exposure_dataframe)
  number_of_exposure_features <- length(exposure_feature_labels)
  
  for (exposure_feature_index in 1:number_of_exposure_features) {
    case_dataframe[, exposure_feature_labels[exposure_feature_index]] <- vector(class(exposure_dataframe[, exposure_feature_index]), number_of_cases)
  }
  
  for (case_index in 1:number_of_cases) {
    matched_exposure_index <- match(case_dataframe[case_index, case_match_column], exposure_dataframe[, exposure_match_column])
    if (is.na(matched_exposure_index)) {
      warning("Match not found.")
      cat(sprintf('%s \t %s\n', as.character(case_dataframe[case_index, case_match_column]), as.character(case_index)))
    }
    else {
      case_dataframe[case_index, (number_of_case_features + 1):(number_of_case_features + number_of_exposure_features - 1)] <- exposure_dataframe[matched_exposure_index, -exposure_match_column]
    }
  }
  
  if (drop_match_column) {
    case_dataframe[, case_match_column] <- NULL
  }
  
  return(case_dataframe)
}

TCEQ.generate_air_pollution_statistics <- function (dataframe, region=NA, pollutants=NA, AQS_codes=NA, stats=NA) {
  dataframe <- help.TCEQ.normalize_source_data(dataframe)
  if (is.na(AQS_codes)) {
    if (is.na(region)) {
      AQS_codes <- unique(dataframe$site)
    }
    else {
      AQS_codes <- TCEQ.sites(region = region)$AQS.Code
      # if (is.element(10, region))
      #   AQS_codes <- c(AQS_codes, 480055013, 482450009, 482450011, 482450014, 482450017, 482450018, 482450019, 482450021, 482450022, 482450101, 482450102, 482450628, 482451035, 483611001, 483611100)
      # if (is.element(11, region))
      #   AQS_codes <- c(AQS_codes, 480210684, 480535009, 480551604, 481490001, 482090614, 482091675, 484530014, 484530020, 484530021, 484530326, 484531068, 484531603, 484531605, 484535001, 484535003, 484910690, 484916602)
      # if (is.element(9, region))
      #   AQS_codes <- c(AQS_codes, 480271045, 480271047, 480415011, 483091037, 483095010)
      # if (is.element(13, region))
      #   AQS_codes <- c(AQS_codes, 480290032, 480290051, 480290052, 480290053, 480290055, 480290059, 480290060, 480290501, 480290502, 480290622, 480290623, 480290625, 480290626, 480290676, 480290677, 480291069, 480910503, 480910505, 481870504, 481870506, 481875004, 482551070, 484931038)
      # if (is.element(5, region))
      #   AQS_codes <- c(AQS_codes, 480370004, 480371031, 481830001, 482030002, 484230007)
      # if (is.element(12, region))
      #   AQS_codes <- c(AQS_codes, 480390618, 480391003, 480391004, 480391012, 480391016, 480710013, 480711606, 481570696, 481670004, 481670005, 481670056, 481670571, 481670615, 481670616, 481670621, 481670683, 481670697, 481671034, 481675005, 482010024, 482010026, 482010029, 482010036, 482010046, 482010047, 482010051, 482010055, 482010057, 482010058, 482010060, 482010061, 482010062, 482010066, 482010069, 482010071, 482010075, 482010307, 482010416, 482010551, 482010552, 482010553, 482010554, 482010556, 482010557, 482010558, 482010559, 482010560, 482010561, 482010562, 482010563, 482010570, 482010572, 482010617, 482010669, 482010670, 482010671, 482010673, 482010695, 482010803, 482011015, 482011017, 482011034, 482011035, 482011039, 482011042, 482011043, 482011049, 482011050, 48201105, 482011066, 48201600, 48291069, 483390078, 483390698, 483395006, 484715012)
      # if (is.element(6, region))
      #   AQS_codes <- c(AQS_codes, 480430101, 481095018, 481410029, 481410037, 481410038, 481410044, 481410047, 481410054, 481410055, 481410057, 481410058, 481410693, 481411021)
      # if (is.element(15, region))
      #   AQS_codes <- c(AQS_codes, 480610006, 480611023, 480612004, 482150043, 482151046)
      # if (is.element(1, region))
      #   AQS_codes <- c(AQS_codes, 480650004, 480650005, 480650007, 483750024, 483750320, 483751025)
      # if (is.element(4, region))
      #   AQS_codes <- c(AQS_codes, 480850003, 480850005, 480850007, 480850009, 480850029, 480971504, 481130018, 481130050, 481130061, 481130069, 481130075, 481130087, 481131067, 481131500, 481131505, 481210034, 481211007, 481211013, 481211032, 481215008, 481390016, 481391044, 482210001, 482311006, 482510003, 482511008, 482511063, 482511501, 482570005, 482570020, 483491051, 483495014, 483631502, 483670081, 483671506, 483970001, 484390075, 484391002, 484391006, 484391009, 484391018, 484391053, 484391062, 484391065, 484391503, 484392003, 484393009, 484393010, 484393011, 484395007, 484970088, 484971064)
      # if (is.element(14, region))
      #   AQS_codes <- c(AQS_codes, 481231602, 481750624, 482730314, 483550025, 483550026, 483550029, 483550032, 483550034, 483550041, 483550083, 483550660, 483550664, 483551024, 484090659, 484090685, 484090686, 484090687, 484690003, 484690609)
      # if (is.element(7, region))
      #   AQS_codes <- c(AQS_codes, 481350003, 481351014)
      # if (is.element(16, region))
      #   AQS_codes <- c(AQS_codes, 483230004, 484655017, 484790016, 484790017, 484790313)
      # if (is.element(3, region))
      #   AQS_codes <- c(AQS_codes, 483371507, 484411509, 484415015, 484851508)
      # if (is.element(8, region))
      #   AQS_codes <- c(AQS_codes, 484515016)
    }
  }
  if (!is.na(pollutants)) {
    pollutants <- unique(dataframe[which(dataframe$param_name %in% tolower(pollutants)), "param_id"])
    # pollutant_codes <- vector(mode="numeric", length=0)
    # if (is.element("CO", pollutants))
    #   pollutant_codes <- c(pollutant_codes, 42101)
    # if (is.element("SO2", pollutants))
    #   pollutant_codes <- c(pollutant_codes, 42401)
    # if (is.element("NO", pollutants))
    #   pollutant_codes <- c(pollutant_codes, 42601)
    # if (is.element("NO2", pollutants))
    #   pollutant_codes <- c(pollutant_codes, 42602)
    # if (is.element("NOx", pollutants))
    #   pollutant_codes <- c(pollutant_codes, 42603)
    # if (is.element("O3", pollutants))
    #   pollutant_codes <- c(pollutant_codes, 44291)
    # if (is.element("PM", pollutants))
    #   pollutant_codes <- c(pollutant_codes, 88502, 88101)
    # if (is.element("Temperature", pollutants))
    #   pollutant_codes <- c(pollutant_codes, 62101)
    # if (is.element("Wind Speed", pollutants))
    #   pollutant_codes <- c(pollutant_codes, 61103)
    # if (is.element("Wind Direction", pollutants))
    #   pollutant_codes <- c(pollutant_codes, 61104)
    # if (is.element("Humidity", pollutants))
    #   pollutant_codes <- c(pollutant_codes, 62201)
    # if (is.element("Solar Radiation", pollutants))
    #   pollutant_codes <- c(pollutant_codes, 63301)
    # pollutants <- pollutant_codes
  } 
  else {
    pollutants <- unique(dataframe$param_id)
  }
  if (is.na(stats)){
    stats <- c('maxhr','eighthrrolling','eighthrmax')
  }
  if(!"AQS_codes" %in% colnames(dataframe)){
    AQS_codes[AQS_codes == 481490001] <- 1
    AQS_codes[AQS_codes == 800060011] <- 10
    AQS_codes[AQS_codes == 484970088] <- 100
    AQS_codes[AQS_codes == 480290626] <- 1001
    AQS_codes[AQS_codes == 481410055] <- 1015
    AQS_codes[AQS_codes == 482010060] <- 1016
    AQS_codes[AQS_codes == 480290059] <- 1020
    AQS_codes[AQS_codes == 482010553] <- 1022
    AQS_codes[AQS_codes == 483550660] <- 1029
    AQS_codes[AQS_codes == 800060007] <- 1034
    AQS_codes[AQS_codes == 800060007] <- 1036
    AQS_codes[AQS_codes == 483550660] <- 1049
    AQS_codes[AQS_codes == 481211013] <- 108
    AQS_codes[AQS_codes == 481231602] <- 109
    AQS_codes[AQS_codes == 481231602] <- 110
    AQS_codes[AQS_codes == 484415015] <- 114
    AQS_codes[AQS_codes == 800060007] <- 139
    AQS_codes[AQS_codes == 800060011] <- 145
    AQS_codes[AQS_codes == 800060007] <- 146
    AQS_codes[AQS_codes == 800060011] <- 147
    AQS_codes[AQS_codes == 480535009] <- 15
    AQS_codes[AQS_codes == 800060007] <- 150
    AQS_codes[AQS_codes == 481231602] <- 152
    AQS_codes[AQS_codes == 800060011] <- 154
    AQS_codes[AQS_codes == 800060011] <- 1615
    AQS_codes[AQS_codes == 800060011] <- 1616
    AQS_codes[AQS_codes == 800060011] <- 1621
    AQS_codes[AQS_codes == 484970088] <- 165
    AQS_codes[AQS_codes == 484971064] <- 166
    AQS_codes[AQS_codes == 482010058] <- 167
    AQS_codes[AQS_codes == 484655017] <- 169
    AQS_codes[AQS_codes == 481410055] <- 18
    AQS_codes[AQS_codes == 481231602] <- 181
    AQS_codes[AQS_codes == 481490001] <- 2003
    AQS_codes[AQS_codes == 484655017] <- 22
    AQS_codes[AQS_codes == 482010058] <- 235
    AQS_codes[AQS_codes == 800060011] <- 240
    AQS_codes[AQS_codes == 800060007] <- 243
    AQS_codes[AQS_codes == 800060011] <- 245
    AQS_codes[AQS_codes == 484393009] <- 26
    AQS_codes[AQS_codes == 481490001] <- 309
    AQS_codes[AQS_codes == 481131067] <- 34
    AQS_codes[AQS_codes == 481211013] <- 35
    AQS_codes[AQS_codes == 484655017] <- 404
    AQS_codes[AQS_codes == 800060007] <- 405
    AQS_codes[AQS_codes == 800060007] <- 406
    AQS_codes[AQS_codes == 482450014] <- 407
    AQS_codes[AQS_codes == 800060007] <- 408
    AQS_codes[AQS_codes == 483550660] <- 409
    AQS_codes[AQS_codes == 484531603] <- 410
    AQS_codes[AQS_codes == 800060007] <- 411
    AQS_codes[AQS_codes == 800060011] <- 416
    AQS_codes[AQS_codes == 800060011] <- 45
    AQS_codes[AQS_codes == 800060007] <- 48
  }
  
  hourly_averages <- data.frame()
  hourly_maximums <- data.frame()
  eight_hour_rolling_averages <- data.frame()
  eight_hour_maximums <- data.frame()
  
  for (AQS_code_index in (1:length(AQS_codes))) {
    matching_rows <- dataframe[which(dataframe$site == AQS_codes[AQS_code_index]),]
    for (pollutant_index in (1:length(pollutants))) {
      matching_rows <- matching_rows[which(matching_rows$param_id == pollutants[pollutant_index]),]
      if (nrow(matching_rows) > 0) {
        hourly_averages <- TCEQ.hourly_aggregate(matching_rows, "mean")
        if (is.element('maxhr', stats))
          hourly_maximums <- TCEQ.hourly_aggregate(matching_rows, "max")
        if (is.element('eighthrrolling', stats) || is.element('eighthrmax', stats)){
          eight_hour_rolling_averages <- hourly_averages
          eight_hour_rolling_averages$value <- TCEQ.rolling_average(hourly_averages$value)
        }
        if (is.element('eighthrmax', stats)) {
          before_4pm <- eight_hour_rolling_averages[which(eight_hour_rolling_averages$hourofday <= 16), ]
          eight_hour_maximums <- TCEQ.daily_aggregate(before_4pm)
        }
      }
    }
  }
  return(list(hourly_averages=hourly_averages, hourly_maximums=hourly_maximums, eight_hour_rolling_averages=eight_hour_rolling_averages, eight_hour_maximums=eight_hour_maximums))
}

TCEQ.merge_pollutants <- function(case_dataframe, exposure_dataframe, pollutant_labels, dist = 20) {
  
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
  #           if (nrow(reduced2) > 0){
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
        Out[i,j*10-4] <- E$hour[exprow0]
        Out[i,j*10-3] <- E$value[exprow0]
      }
      if(!is.na(exprow1)) {
        Out[i,j*10-2] <- E$hour[exprow1]
        Out[i,j*10-1] <- E$value[exprow1]
      }
      if(!is.na(exprow2)) {
        Out[i,j*10] <- E$hour[exprow2]
        Out[i,j*10+1] <- E$value[exprow2]
      }
      if(!is.na(exprow3)) {
        Out[i,j*10+2] <- E$hour[exprow3]
        Out[i,j*10+3] <- E$value[exprow3]
      }
      if(!is.na(exprow4)) {
        Out[i,j*10+4] <- E$hour[exprow4]
        Out[i,j*10+5] <- E$value[exprow4]
      }
    }
  }
  
  print(Sys.time()-strt)
  return(Out)
  
}

TCEQ.exposure_index <- function(Crow, code, allcodes, E, dist, day) {
  index <- 0
  # match code to go in right branch of E
  dd <- match(code,allcodes)
  if(!is.na(dd)) {
    # find correct days
    d <- (E[[dd]]$epoch==Crow$DOS-day)
    reduced2 <- E[[dd]][d,]
    localdist <- dist
    
    if (nrow(reduced2) > 0){
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
TCEQ.geodetic_distance <- function(point1, point2)
{
  R <- 6371
  p1rad <- point1 * pi/180
  p2rad <- point2 * pi/180
  d <- sin(p1rad[2])*sin(p2rad[2])+cos(p1rad[2])*cos(p2rad[2])*cos(abs(p1rad[1]-p2rad[1]))  
  d <- acos(d)
  R*d
}

help.TCEQ.get_chloropleth_color <- function(chor) {
  if (chor$val1_O3 > 55)
    return(2)
  else if ((chor$val0_SO2 > 7) & (chor$val0_NO > 37) & (chor$val0_NO2 > 34) & (chor$val1_PM > 24))
    return(2)
  else if ((chor$val3_NO > 37) & (chor$val4_NO2 > 34) & (chor$val1_PM > 24))
    return(2)
  else if ((chor$val0_NO > 37) & (chor$val1_NO > 37) & (chor$val1_NO2 > 34) & (chor$val0_PM > 24) & (chor$val1_PM > 24))
    return(2)
  else if ((chor$val2_NO > 37) & (chor$val4_NO > 37) & (chor$val1_PM > 24))
    return(2)
	else if ((chor$val0_NO > 37) & (chor$val0_NO2 > 34) & (chor$val1_NO2) & (chor$val0_PM) & (chor$val1_PM))
		return(2)
	else if ((chor$val0_SO2 > 7) & (chor$val3_O3) & (chor$val1_PM > 24))
		return(2)
	else if ((chor$val0_SO2 > 7) & (chor$val0_PM > 24) & (chor$val2_PM > 24))
		return(2)
	else if ((chor$val1_SO2 > 7) & (chor$val3_NO2 > 34) & (chor$val2_O3 > 55))
		return(2)
	else if ((chor$val1_NO > 37) & (chor$val2_NO2 > 34) & (chor$val0_PM > 24))
		return(2)
	else if ((chor$val1_SO2 > 7) & (chor$val1_NO > 37) & (chor$val0_O3 > 55))
		return(2)
	else if ((chor$val3_NO > 37) & (chor$val4_NO > 37) & (chor$val1_NO2 > 34))
		return(2)
	else if ((chor$val3_NO > 37) & (chor$val0_NO2 > 34) & (chor$val2_NO2 > 34))
		return(2)
	else if ((chor$val0_NO2 > 34) & (chor$val2_O3 > 55) & (chor$val0_PM > 24))
		return(2)
	else if ((chor$val3_NO2 > 34) & (chor$val1_PM > 24))
		return(2)
	else if ((chor$val1_NO2 > 34) & (chor$val2_O3 > 55) & (chor$val0_PM > 24))
		return(2)
	else if ((chor$val2_O3 > 55) & (chor$val1_PM > 24) & (chor$val4_PM > 24))
		return(2)
	else if ((chor$val3_NO > 37) & (chor$val4_NO > 37) & (chor$val2_NO2 > 34))
		return(2)
	else if ((chor$val2_O3 > 55) & (chor$val0_PM > 24) & (chor$val2_PM > 24))
		return(2)
	else if ((chor$val4_SO2 > 7) & (chor$val0_PM > 24))
		return(2)
	else if ((chor$val1_SO2 > 7) & (chor$val2_O3 > 55) & (chor$val3_O3 > 55))
		return(2)
	else if ((chor$val1_NO > 37) & (chor$val4_O3 > 55))
		return(2)
	else if ((chor$val4_NO > 37) & (chor$val0_O3))
		return(2)
	else if ((chor$val0_SO2 > 7) & (chor$val0_O3 > 55))
		return(2)
	else if ((chor$val0_O3 > 55) & (chor$val1_PM))
		return(2)
	else if ((chor$val0_O3 > 55) & (chor$val4_O3 > 55))
		return(2)
	else if ((chor$val0_O3 > 55) & (chor$val0_PM > 24))
		return(2)
  return()
}
