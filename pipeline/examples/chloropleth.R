library(doParallel)
library(foreach)

#Load hourly averages
houravg=read.csv('1houravg_allsites_0209_7P.csv')
colnames(houravg) <- c("param_id", "region", "site","year","month","day","hour","value")
houravg$date <- as.Date(with(houravg, paste(year, month, day,sep="-")), "%Y-%m-%d")
houravg$time <- paste(houravg$hour,"00","00",sep=":")
houravg$hourofday <- as.POSIXct(paste(houravg$date, houravg$time), format="%Y-%m-%d %H:%M:%S")

#Get Daily averages and prepare for merge
source('maxhour.R')
daily=maxhr(houravg)
colnames(daily)[colnames(daily) == 'x'] <- 'value'
colnames(daily)[colnames(daily) == 'param_id'] <- 'poll_code'
daily$epoch=as.numeric(daily$hour)
daily$hour=as.numeric(format(daily$hour, format="%H"))
daily$date=as.Date(format(daily$day, format="%Y-%m-%d"))
latlong=read.csv('latandlongofsites.csv')
daily <- (merge(latlong, daily, by = 'site'))
daily$epoch=as.numeric(daily$date)
day<- daily[c("poll_code", "site", "latitude","longitude","epoch","hour","value")]

#Get Zipcodes and prepare for merge
zipcode=read.csv('latandlongofzipcodecentroids.csv')
zipcode$date=as.Date("2003-2-6",format="%Y-%m-%d")
zipcode$DOS=as.numeric(zipcode$date)
zipcode$event=1
colnames(zipcode)[colnames(zipcode) == 'Lat'] <- 'latitude'
colnames(zipcode)[colnames(zipcode) == 'Long'] <- 'longitude'
zipcodes=zipcode[c("Zipcode","event","DOS","latitude","longitude")]

#Get pollcodes
pollcodes=read.csv('pollcodes.csv')
pollcodes$poll_code=as.factor(pollcodes$poll_code)

#Merge
source('mergepoll5.R')
chor=merge_poll(zipcodes,day,pollcodes)

#Get Colors
source('getColor.R')
chor$color=getColor(chor)

#load choropleth maker
install.packages("devtools")
library(devtools)
install_github('arilamstein/choroplethrZip@v1.5.0')
library(choroplethr)
library(choroplethrZip)
library(ggplot2)

#Current Texas Map
colnames(chor)[colnames(chor) == 'Zipcode'] <- 'region'
colnames(chor)[colnames(chor) == 'color'] <- 'value'
chor$region=as.character(chor$region)
chor$value[is.na(chor$value)] <- 0
zip_choropleth(chor, state_zoom="texas")

#Messing with stuff
choro = ZipChoropleth$new(df_pop_zip)
choro$title = "2012 ZCTA Population Estimates"
choro$ggplot_scale = scale_fill_brewer(name="Population", palette=2, drop=FALSE)
choro$set_zoom_zip(state_zoom="new york", county_zoom=NULL, msa_zoom=NULL, zip_zoom=NULL)
choro$render()

colnames(chor)[colnames(chor) == 'Zipcode'] <- 'region'
colnames(chor)[colnames(chor) == 'color'] <- 'value'
chor$region=as.character(chor$region)
chor$value[is.na(chor$value)] <- 0
zip_choropleth(chor, state_zoom="texas")
a <- subset(zip.regions,zip.regions$state.name=="texas") 

#Messing with choropleths from map-grod
library(foreign)	# To read .dbf file
library(maps)		# To draw map
library(plyr)		# Data formatting
library(mapproj)

a <- subset(zip.regions,zip.regions$state.name=="texas") 
a <- (merge(a, chor, by = 'region'))
a$value[a$value==0] <- "#13373e"
a$value[a$value==2] <- "#246571"
map(a, regions=a$region, proj="albers", param=c(39,45), fill=TRUE, col=a$value, border=NA, resolution=0)
