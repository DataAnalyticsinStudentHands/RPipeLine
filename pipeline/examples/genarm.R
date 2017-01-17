source('../pipeline.R')

# data loading functions
load.saved_rdata("/Users/gtoti/Desktop/my papers/1 - asthma analysis with log reg and ARM/002_CompleteCases_prebin.RData")

#user-define data parameters
columns_to_process <- seq(from=7, to=65, by=2)

#data preparing functions
dataframe <- S[, c(2, columns_to_process)]

#data processing functions
results <- process.using_GenARM(dataframe, gene_pool_size=30, seed=3, weights=c(1,1,1,1), fitness_mode="sum", with_statistics=TRUE)

#data output functions
output.dataframe_to_file(results$statistics, "GenARMreport.csv")
output.plot_dataframe(results$statistics, x_axis_column="rule_index", columns_to_plot=c(1,12:18))