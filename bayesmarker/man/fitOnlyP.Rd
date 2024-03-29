% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitOnlyP.R
\name{fitOnlyP}
\alias{fitOnlyP}
\title{fitOnlyP}
\usage{
fitOnlyP(
  SUBPOP,
  Measured,
  mapping,
  pred.data,
  quick = FALSE,
  cores = NULL,
  save_directory = "."
)
}
\arguments{
\item{SUBPOP}{The subpopulation name for which we run the inference}

\item{Measured}{Output from the examine_error function (saved as
NewMeasured), with the last column being PosteriorShape}

\item{mapping}{Output from the readNHANES function, saved in MCMCdata file}

\item{pred.data}{Output from the readNHANES function, saved in MCMCdata file}

\item{quick}{Default is FALSE.  If TRUE, the calculations will be run
with less iterations, resulting in a shorter overall runtime, however,
model convergence may not be reached.  Check convergence using the fromMolar function.}

\item{cores}{Number of cores to use during computation}

\item{save_directory}{String providing the directory in which to save the results.  If
left as the default, ".", it will save to ./OnlyPparms3_SUBPOP.RData
with a date stamp.  Otherwise, it will save to save_directory.}
}
\value{
out.samps3R
out.coda3R
nhanesdata: a list containing all the data needed to run the jags model.
Used as input to the fromMolar function to assign chemical
information to the final estimates.
}
\description{
Runs Bayesian inference for a provided subpopulation
}
