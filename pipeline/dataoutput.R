output.dataframe_to_file <- function(dataframe, filename, datatype=FALSE) {
  if (is.null(dataframe)) {
    stop("No output data to write to file")
  }
  
  if (typeof(datatype) != "str") {
    library(tools)
    datatype <- file_ext(filename) 
  }
  
  if (datatype == "csv") {
    write.csv(dataframe, file = filename)
  }
}

output.plot_dataframe <- function(dataframe, x_axis_column, columns_to_plot) {
  require(ggplot2)
  require(reshape)

  plot_dataframe <- melt(dataframe[, columns_to_plot], id.vars=x_axis_column, variable.name='metrics')
  ggplot(plot_dataframe, aes(x=x_axis_column, y=value, group=variable, colour=variable)) + geom_line() + geom_point()
}

output.save_all_rdata <- function(filename) {
  save(file=filename)
}