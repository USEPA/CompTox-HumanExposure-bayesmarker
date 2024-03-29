% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{lnlike}
\alias{lnlike}
\title{lnlike}
\usage{
lnlike(x, mean, lsd)
}
\arguments{
\item{x}{log of measurement value ??? what is in column 2?}

\item{mean}{same as x}

\item{lsd}{log of standard deviation}
}
\value{
lognormal distribution with left censoring fit to the inputs
}
\description{
Fit lognormal distribution with left censoring to a dataset
Assume we are passing the value in x
cens is 0 for non-censored, 1 for left-censored (below LOD, eg)
llod is the log-limit of detection
}
