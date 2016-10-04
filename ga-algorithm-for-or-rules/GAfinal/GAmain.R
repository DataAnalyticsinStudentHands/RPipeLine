#### MAIN FOR GA MINING ####

# Load data to mine
load("/Users/gtoti/Desktop/my papers/1 - asthma analysis with log reg and ARM/002_CompleteCases_prebin.RData")
cols <- seq(7,65,by=2)
X <- S[,c(2,cols)]
source('GAmodules.R')

testPop <- GA.modules(X, upop_size = 30, seed=3, weights=c(1,1,1,1), fitness="sum")
report <- GAreport(X, upop_size = 30, testPop, weights=c(1,1,1,1))

# # Plot of fitness and other metrics
# require(ggplot2)
# require(reshape)
# 
# col <- c(1,12:18)
# df <- melt(report[,col], id.vars = 'Rule', variable.name = 'metrics')
# # plot on same grid, each series colored differently 
# ggplot(df, aes(x=Rule,y=value,group=variable, colour=variable)) + geom_line() + geom_point()
