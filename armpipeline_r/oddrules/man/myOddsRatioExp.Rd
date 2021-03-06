% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myOddsRatioExp.R
\name{myOddsRatioExp}
\alias{myOddsRatioExp}
\title{My Odds Ratio}
\usage{
myOddsRatioExp(x, transactions, transNot, reuse = TRUE, CI = FALSE, t=0)
}
\arguments{
\item{x}{A set of rules.}

\item{transactions}{The transaction data set used to mine the associations or a set of different transactions to calculate interest measures from (Note: you need to set reuse=FALSE in the later case).}

\item{transNot}{The negated transaction data for calculations}

\item{reuse}{A logical indicating if information in quality slot should be reuse for calculating the measures.}

\item{CI}{Include 95 percent confidence interval and statistical significance in the output data frame (logical). Default = FALSE.}

\item{t}{Tolerance for statistical significance (0<=t<=1). Default = 0 (no tolerance)}
}
\value{
return: a data frame containing at least 2 columns, oddsRatio and standardError. If CI = TRUE, 3 additional columns are included in the data frame: lowCI95, highCI95, and significance. lowCI95 and highCI95 store the margins of the 95 percent confidence interval. significance is a logical column that states if the rule is statistically significant given the specified tolerance. When t = 0, a rule is statistically significant if the confidence interval does not include 1. When t > 0, a rule is statistically significant in one of these cases:the confidence interval does not include 1, highCI95 <= 1+t and lowCI95 <= 1-t, highCI95 >= 1+t and lowCI95 >= 1-t
}
\description{
This function is based on a function included in the package arules basicRuleMeasure originally designed to compute various measures of interest for association rules
The function has been modified to make it dedicated to odds ratio only, with some new functionalities such as standard error, confidence interval, and statistical significance.
}

