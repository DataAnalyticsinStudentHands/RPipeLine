#### MAIN FOR GA MINING ####

# Load data to mine
load("~/Desktop/my papers/7 - GA for ARM/synthetic datasets/GAdataset4.RData")
source('GAmodules.R')
source('GAfunctions.R')

testPop <- GA.absolute(X, seed=3, weights=c(1,1,0,1,0,1,1))
testPop <- GA.modules(X, seed=3, weights=c(1,1,1,1,1), fitness = "arithmetic")
report <- GAreport(X, testPop, weights=c(1,1,1,1,1), filename = "report2")

# # Plot of fitness and other metrics
# require(ggplot2)
# require(reshape)
# 
# col <- c(1,12:18)
# df <- melt(report[,col], id.vars = 'Rule', variable.name = 'metrics')
# # plot on same grid, each series colored differently 
# ggplot(df, aes(x=Rule,y=value,group=variable, colour=variable)) + geom_line() + geom_point()
