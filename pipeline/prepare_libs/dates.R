dates.normalize <- function(dates, origin = "0000-01-01", adjust = 0) {
  # INPUT
  # dates = array of dates to test. Accepted formats are Date or numerical (date vector).
  # origin = origin date to use for conversion from date vector to Date.
  # adjust = shift to apply to the array of date vectors before conversion to Date. 
  # For example, conversion from Matlab date vectors requires adjust = -1.
  if ("Date" != class(dates) & !is.numeric(dates)) {
    stop("Invalid date format")
  }
  
  if ("Date" == class(dates) & adjust!=0) {
    stop("Input in Date format can not be adjusted")
  }
  
  if (is.numeric(dates)) {
    dates <- dates + adjust
    dates <- as.Date(dates, origin)
  }
  
  M <- as.numeric(format(dates, "%m"))
  
  return(M)
}

dates.in_months <- function(dates, months) {
  # This function returns a logical array of the same length of "dates", where TRUE means
  # that the date falls within the indicated months.
  #
  # INPUT
  # dates = array of dates to test. Accepted formats are Date
  # months = array of months (values 1-12) that form the season of interest.
  
  if (min(months)<1 | max(months)>12) {
    stop("Invalid month(s) entered")
  }
  
  ind <- match(M,months)
  ind <- !is.na(ind)
  
  return(ind)
}

