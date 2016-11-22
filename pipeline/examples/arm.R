source('../pipeline.R')

#data loading functions
dataframe <- load.file_as_dataframe('sample_data/ARM.csv')

#user-define data parameters
columns_to_process <- seq(from=6, to=64, by=2)

#data preparing functions
dataframe <- prepare.prune_dataframe(dataframe, omit_columns=c('X'), omit_NA_rows=TRUE)
binarized_dataframe <- prepare.binarize_data_frame(dataframe, quantile_threshold=0.95, column_list=columns_to_process)$dataframe
split_data <- prepare.split_dataframe_for_training_and_testing(binarized_dataframe, sampling_method='random', training_size=5000)

#data processing functions
results <- process.using_ARM(split_data, seeds=1:10, training_sample_sizes=c(5000,2000,1000,500,200,100,50,20,10), remove_insignificant_results=TRUE, remove_redundant_results=TRUE)
