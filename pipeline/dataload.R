#' dataload.R
#' 
#' This file contains all necessary functions for loading data from files of certain
#' types and returning the data as a 'data frame' object. These are strictly file loading functions.
#' They do not investigate, analyze, or modify any aspect of the data itself. Whatever the content of 
#' the datafiles is should not be considered in these functions (i.e. there do not need to be different 
#' versions of loading csv files). If, for example, certain columns or rows need to be ignored, or 
#' certain values need to be changed or removed, these types of activities should be defined in the 
#' "dataprepare.R" file.
#' 
#' When support needs to be added for a new datatype, please create a stand-alone function first, 
#' and then add it as a new case in the "if" conditions of the generic load.file_as_dataframe function. 


#' load.file_as_dataframe(datafile, datatype = FALSE)
#' 
#' Generic function to load a file of data that will be sent through the rest of the pipeline
#'
#' All pipeline functions assume that data is loaded as an R 'data frame' object. Actual loading
#' is delegated out to a stand-alone datatype-specific loading function. This generic function should
#' be used so that code is consistent, readable, and because if the datafiles themselves are properly
#' named, no actual datatype parameter needs to be specified and code that uses this function will not
#' have to be changed to support different datafiles. 
#' 
#' @param datafile  The path to the datafile given as a string
#' @param datatype  Optional override of automatic data format detection
#'
#' @return a dataframe object containing the contents of the file
#' @export
#'
#' @examples
#+ example_dataframe <- load.file_as_dataframe("./examples/sample_data/test.csv")
#+ example_dataframe <- load.file_as_dataframe("./examples/sample_data/test.dat", datatype = "csv")
load.file_as_dataframe <- function(datafile, datatype = FALSE, maximum_rows = NA) {
  dataframe <- NULL
  
  if (typeof(datatype) != "str") {
    library(tools)
    datatype <- file_ext(datafile) 
  }
  if (datatype == "csv") {
    dataframe <- load.csv_as_dataframe(datafile, maximum_rows)
  }
  
  if (is.null(dataframe)) {
    stop("File generated an empty data frame")
  }
  
  return(dataframe)
}

#' load.csv_as_dataframe(datafile)
#' 
#' Single-purpose function to load data from file in comma-seperated-value format
#' 
#' @param datafile  The path to the datafile given as a string
#'
#' @return a dataframe object containing the contents of the file
#'
#' @examples
#+ example_dataframe <- load.csv_as_dataframe("./examples/sample_data/test.csv")
load.csv_as_dataframe <- function(datafile, maximum_rows = NA) {
  dataframe <- read.csv(datafile, nrows = maximum_rows)
  return(dataframe)
}

load.saved_rdata <- function(datafile) {
  load(datafile)
}