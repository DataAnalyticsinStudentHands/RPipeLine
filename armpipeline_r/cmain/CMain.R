###Main###

source('Cbin_quantile.R')
source('random_split.R')
source('C_overfit_test.R')

##Load Data
data <-read.csv("C:/Users/owner/Desktop/pollutant stats/Carol - ARM/002subjects_poll_merged_5dlag.csv")
data$X<-NULL

##Bin Data
#It is recomended the User clean the data before binning, but the User can only bin certain columns by specifying the parameter "col"
#The User may also choose to group columns in order to use the same threshold for multiple columns this can be done with the parameter "group", the User can also chose a column to create the threshold for the group using "groupthreshold"
col <- seq(from = 6, to =64, by = 2)
binneddata <- bin_quantile(data, cols=col, top=0.95)
binned <- as.data.frame(binneddata$data)

##Make Data Logical
#if the user has not cleaned the data, the user will have to specify which columns they want to use for rule mining
col <- seq(from = 6, to =64, by = 2) 
sample <- binned[c(1,col)] #index which colss here
sampleLog <- data.frame(lapply(sample, as.logical))

##Randomize and split data into training and testing
splitdata <- random_split(sampleLog) # can specify seed with parameter "seed" or how many to put in training with "split"
whole_training <- as.data.frame(splitdata$whole_training)
testing <- as.data.frame(splitdata$testing)

##Mine and test rules
#Before mining the data frame must have a column named event that the User is interested in finding rules to predict for.
Results <- overfit_test(whole_training,testing, seed=1:10, trainingsize=c(5000,2000,1000,500,200,100,50,20,10))

##start filtering

Significant=subset(Results$TestFrame,significance==TRUE)
Significant=unique(Significant)
Significant$nonredundant=nonredundant(Significant,colnames(sample)[-1])
Significant=subset(Significant,nonredundant==TRUE)

##For G4 and G5
data <-read.csv("C:/Users/owner/Desktop/pollutant stats/Carol - ARM/GRules5.csv")
var=colnames(sample)[-1]
nonredundant(data,var)
data$nonredundant=nonredundant(data,var)
data=subset(data,nonredundant==TRUE)
write.csv(data,"sheet5rules")
