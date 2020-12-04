#' fromMolar
#'
#' This is generic code for checking the results of the MCMC analyses of
#' the NHANES-based modeling and for converting the estimates from units
#' of nmole / kg bodyweight / day to units appropriate for exposure
#' mg / kg bodyweight / day.
#' We check convergence first, using heidel.diag() and gelman.diag
#' on the three chain results.  Check both lP and lU (for OnlyP).
#' Then, compare the lU estimates to the values of Measured$loggm.  Just a plot
#' should do. This is not a particularly rigorous test, but if it fails,
#' certainly we have problems.
#' Finally, do the conversion, and save.
#'
#' @param SUBPOP The subpopulation name for which we want to look at convergence
#' @param Measured Output from the examine_error function (saved as
#'                 NewMeasured), with the last column being PosteriorShape.
#'                 Same as input for fitOnlyP
#' @param pred.data Output from the readNHANES function, saved in MCMCdata file
#' @param nhanesdata Output from fitOnlyP. A list containing all the data needed to run the jags model.
#'                   Used as input to the fromMolar function to assign chemical
#'                   information to the final estimates.
#' @param out.coda3R Output from fitOnlyP. A MCMC list object.
#' @param doplot Logical input to indicate whether a plot should be generated and saved.
#'               Default is TRUE.
#'
#'
#' @import ggplot2
#' @importFrom coda varnames gelman.diag heidel.diag
#'
#'
#' @return lPsampsgm: data frame whose rows are the parent chemicals, columns are
#'                    the different subpopulations, and the values are the exposure
#'                    estimates in mg/kg bodyweight/day for that cheimcal in that
#'                    subpopulation.
#'         Measured: same as the input to fitOnlyP, reduced to rows matching the
#'                   input subpopulation. Geometric means of the metabolites that
#'                   are created from the parents in lPsampsgm.
#'
#' @export
#'
fromMolar <- function(SUBPOP, Measured, pred.data, nhanesdata, out.coda3R, doplot = TRUE) {

  ## Take samples of lP from OnlyP, convert to a matrix, and convert from
  ## nmoles/kg bodyweight/day to mg/kg bodyweight/day
  ## nmoles/1e9 = moles
  ## moles * MW = grams
  ## grams * 1000 = mg
  # nmoles / 1e9 * MW * 1000 = nmoles * MW / 1e6 = mg
  ## so, add log(MW) - log(1e6) to each lP value to get log(mg / kg bodyweight / day)


  ## -----------------------------------------------------------------------
  ##                   Test Convergence
  ##   lU
  XX <- heidel.diag(out.coda3R[,grep("^lU\\[", varnames(out.coda3R), value=TRUE)])
  ## Summaries for each chain
  writeLines("Convergence of individual lU parameters")
  for (i in 1:3){
    writeLines(paste("Chain ",i,": ",sum(XX[[i]][,"stest"])," of ",nrow(XX[[i]])," passed",sep=""))
  }

  writeLines("Gelman - Rubin diagnostic")
  print(ZZ <- gelman.diag(out.coda3R[,grep("^lU\\[", varnames(out.coda3R), value=TRUE)]))

  XX <- heidel.diag(out.coda3R[,grep("^lP\\[", varnames(out.coda3R), value=TRUE)])
  writeLines("Convergence of individual lP parameters")
  for (i in 1:3){
    writeLines(paste("Chain ",i,": ",sum(XX[[i]][,"stest"])," of ",nrow(XX[[i]])," passed",sep=""))
  }

  writeLines("Gelman - Rubin diagnostic")
  print(ZZ <- gelman.diag(out.coda3R[,grep("^lP\\[", varnames(out.coda3R), value=TRUE)]))


  ## -----------------------------------------------------------------------
  ## Replicate Measured and pred.data as they existed when the inputs for this
  ## data set were created.
  Measured <- subset(Measured, subpop == SUBPOP)
  M <- nrow(Measured)
  Mn <- sum(Measured$PosteriorShape == "normal")
  ## Resort Measured so that the normal posteriors come first, then
  ## the uniform ones.  Break ties by sorting on CAS
  indx <- order(Measured$PosteriorShape, Measured$CAS)
  Measured <- Measured[indx,]


  ## ------------------------------------------------------------------------
  ## Now, merge all the lU (metabolite) estimates and compare to corresponding
  ## loggm values from the data.
  lU <- as.matrix(out.coda3R[,grep("^lU\\[", varnames(out.coda3R))])

  lUests <- t(apply(lU, 2, function(z) {
    x <- quantile(z, pr = c(0.025, 0.975))
    c(x[1], mean(z), x[2])
  }))

  if (doplots){
    tdta <- data.frame(y=lUests[,2], x=Measured$loggm, ylower=lUests[,1], yupper=lUests[,3],
                       xlower=Measured$loggm - Measured$loggm_se * 1.96,
                       xupper=Measured$loggm + Measured$loggm_se * 1.96,
                       PosteriorShape=Measured$PosteriorShape)
    p <- ggplot(data=tdta, aes(x=x, y=y, ymin=ylower, ymax=yupper,
                               xmin=xlower, xmax=xupper, color=PosteriorShape)) + geom_point() +
      geom_errorbar(width=0) +
      geom_errorbarh(width=0) +
      geom_abline(xintercept=0, yintercept=1, lty=3, color="darkgray") +
      scale_y_continuous("Estimated lU", limits=c(-7.5,2.5)) +
      scale_x_continuous("Measurement from NHANES", limits=c(-40,5)) +
      ggtitle("Excretion Rates for Urinary Metabolites")
    pdf(file=paste("MetabComparison_", SUBPOP, "_", format(Sys.time(), "%Y-%m-%d"), ".pdf", sep=""))
    print(p)
    dev.off()
  }

  ## Match up rows in pred.data to lP. They are in the same order, but when not all
  ## the parents have a prediction (e.g. limited data for 0 - 5 subpop), we need to
  ## use the right MWs and reduce rows of pred.data for output.
  if (MissingCols) {
    ind <- match(row.names(nhanesdata$Phi), pred.data$CAS)
    pred.data <- pred.data[ind,]
  }

  lPsampsmo <- as.matrix(out.coda3R[,grep("^lP\\[", varnames(out.coda3R))])

  adjust <- log(pred.data$MW) - log(1e6)

  lPsampsgm <- sweep(lPsampsmo, 2 ,adjust,"+")

  row.names(lPsampsgm) <- pred.data$CAS

  save(lPsampsgm, Measured, file=paste("lPsamps-gm_kg_day_", SUBPOP, "_",
                                                  format(Sys.time(), "%Y-%m-%d"), ".RData", sep = ""))

  return(lPsampsgm, Measured)

}


