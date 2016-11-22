#' DASH R PipeLine
#' 
#' The R Pipeline is a set of modular functions that can be used for quick data analysis tasks.
#' It is designed to make performing multiple types of analysis on multiple types of data 
#' straight-forward and easily readable.
#' 
#' The pipeline is divided into four distinct steps:
#' Loading the data
#' Preparing the data
#' Processing the data
#' Outputting the data
#' 
#' Each of these steps can be extended by adding new functions that follow the same style.
#' Please follow the R coding style guide here: \url{http://adv-r.had.co.nz/Style.html}
#' Please follow the naming conventions of the Pipeline itself, which generally uses verbose names
#' to increase human readability and intelligability to non-experts. 

source('dataload.R')
source('dataoutput.R')
source('dataprepare.R')
source('dataprocess.R')