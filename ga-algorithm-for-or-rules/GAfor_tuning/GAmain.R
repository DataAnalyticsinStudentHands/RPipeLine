#### MAIN FOR GA MINING ####

# Load data to mine
load("GAdataset5.RData")
source('GAabsolute.R')

testPop <- GA.absolute(X, seed=3, weights=c(0,0,0,1,0,0,0))
report <- GAreport(X, testPop, weights=c(0,0,0,1,0,0,0))

# # Plot of fitness and other metrics
# require(ggplot2)
# require(reshape)
# 
# col <- c(1,12:18)
# df <- melt(report[,col], id.vars = 'Rule', variable.name = 'metrics')
# # plot on same grid, each series colored differently 
# ggplot(df, aes(x=Rule,y=value,group=variable, colour=variable)) + geom_line() + geom_point()
