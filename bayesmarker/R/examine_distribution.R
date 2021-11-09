#' examine_error
#'
#' Examines the standard errors of the maximum likelihood estimates and
#' the number of measurements below the LOD.  If either of those
#'
#' @param Measured A data frame containing a row for each metabolite with columns
#'                   providing the chemical identifiers, log geometric means, subpopulation,
#'                   and other data carried on from the codes file.
#' @param codes_file Manually created xls or xlsx file containing 3 sheets:  1. NHANES
#'                   chemicals to include (with identifier, code, file, demographic, and units),
#'                   2. Associated weights, filenames, and column names associated with each phase used,
#'                   and 3. Parent-metabolite map containing chemical identifiers and molecular weights.
#' @param data_path String providing the path to the raw data.  Default is ".",
#'                  the direct output from running get_NHANES_data would be "rawData", resulting
#'                  in e.g. ./rawData/1999-2000.
#' @param save_directory String giving the path of where to save the output table
#'                  and plot.
#'
#' @import survey
#' @import ggplot2
#' @import logspline
#' @import foreach
#' @importFrom parallel mclapply
#' @importFrom stats vcov
#' @importFrom gdata read.xls
#' @importFrom doMC registerDoMC
#'
#' @return Measured: same as the input, but with 4 additional columns.
#'                   1. Prabove: fraction of measurements above the LOD.
#'                   2. Rightseg: right interval of curve around expected
#'                   mean of distribution
#'                   3. Leftseg: left interval of curve around expected
#'                   mean of distribution
#'                   4. PosteriorShape: either uniform or normal. Decided
#'                   based on Prabove and the standard error distribution.
#'
#'
#'
examine_error <- function(Measured, codes_file, data_path = ".", save_directory = ".") {

  nms <- unique(Measured$subpop)

  ## Look at the distribution of loggm_se
  print("Generate distribution plots for the chemical geometric mean standard errors")
  pdf(file = file.path(save_directory, paste("dists_loggm_se_", format(Sys.time(), "%Y-%m-%d"), ".pdf", sep = "")))
  for (nm in nms) {
    print(nm)
    tdta <- Measured[Measured$subpop == nm,]
    dta2 <- data.frame(Lloggm_se = log10(tdta$loggm_se),
                       lsdlog_se = tdta$lsdlog_se)
    #print(dim(tdta))
    lgmse.fit <- logspline(dta2$Lloggm_se)

    x <- with(dta2, seq(min(Lloggm_se, na.rm=TRUE)-.5, max(Lloggm_se, na.rm=TRUE)+.5, length=300))
    y <- dlogspline(x, lgmse.fit)

    dta <- data.frame(x=x, y=y)

    p <- ggplot(data=dta, aes(x=x, y=y)) + geom_line()  +
      geom_rug( data=dta2, mapping=aes(x=Lloggm_se, y=rep(0, nrow(dta2))),sides="b") +
      scale_y_continuous("Density") + scale_x_continuous("log10(loggm_se)")+
      ggtitle(nm)

    print(p)
  }
  dev.off()

  ## So, look at loggm_se  < 0.1, and loggm_se > 10.
  tightse <- Measured[!is.na(Measured$loggm_se) & (Measured$loggm_se < 0.1),]
  loosese <- Measured[is.na(Measured$loggm_se) | (Measured$loggm_se > 10),]


  ## datapaths gives the paths for the NHANES individual .xpt data files.  This
  ## construct gives a vector of paths, one for each sample set.  We then
  ## name the files using the shorthand name for each sample set
  phases <- unique(unlist(sapply(Measured$recent_sample, function(x) strsplit(x, ","))))
  ind <- sort(substring(phases, 4, 5), index.return = TRUE)
  long <- sapply(phases[ind$ix], function(x) ifelse(as.numeric(substring(x, 1, 2)) < 50,
                                                    paste("20", substring(x, 1, 3), "20", substring(x, 4, 5), sep = ""),
                                                    paste("19", substring(x, 1, 3), "20", substring(x, 4, 5), sep = "")))
  datapaths <- structure(file.path(data_path, "rawData", long), names= names(long))


  ## EXCEL spreadsheet giving information about the variables and files used in
  ## this analysis
  NHANEScodes <- codes_file

  ## -------------------------------------------------------------------
  ##             Set up files and codes for main function

  ## convtbl was constructed manually, starting with an earlier list of NHANES urine
  ## products. Age cutoffs were based on the highest age group reported in the 4th report.
  convtbl <- read.xls(NHANEScodes,as.is=TRUE)
  wtvars <- read.xls(NHANEScodes, sheet=2, as.is=TRUE)

  # Reduce these 2 tables to include cohorts only in Measured
  convtbl <- convtbl[convtbl$recent_sample %in% phases,]
  wtvars <- wtvars[wtvars$sample %in% phases,]
                 
  demofiles <- wtvars$demofile
  names(demofiles) <- wtvars$sample

  if (length(grep(",", Measured$NHANESfile)) > 0){
    # Combined data, have multiple files
    tmp2 <- list()
    metabs <- unique(Measured$CAS)
    for (i in 1:length(metabs)){
      ind <- match(metabs[i], Measured$CAS)
      tmp2[[i]] <- cbind(strsplit(Measured$recent_sample[ind], split = ",")[[1]],
                   paste(strsplit(Measured$NHANESfile[ind], split = ",")[[1]], ".XPT", sep = ""))
    }
    names(tmp2) <- metabs

    datafiles <- lapply(tmp2, function(x) file.path(datapaths[as.character(x[,1])],
                                                    tolower(x[,2])))
    demofiles <- lapply(tmp2, function(x) file.path(datapaths[as.character(x[,1])],
                                                    demofiles[as.character(x[,1])]))
    bwtfiles <- lapply(tmp2, function(x) file.path(datapaths[as.character(x[,1])],
                                                   tolower(wtvars$BWfile[match(x[,2], wtvars$file)])))
    creatfiles <- lapply(tmp2, function(x) file.path(datapaths[as.character(x[,1])],
                                                     tolower(wtvars$creatfile[match(x[,2], wtvars$file)])))
    chemwt <- lapply(tmp2, function(x) wtvars$wtvariable[match(x[,2], wtvars$file)])

  } else {
    
    # Some variables are missing a weight
    if (any(wtvars$wtvariable == "")){
      print("Warning:  at least one data file is missing the 2 year weight. This file and its chemicals will be discarded.")
      convtbl <- convtbl[!paste(convtbl$NHANESfile, ".XPT", sep = "") == wtvars$file[wtvars$wtvariable == ""],]
      wtvars <- wtvars[!wtvars$wtvariable == "",]
    }
    tmp <- unique(convtbl[,c("recent_sample","NHANESfile")])
    # Only keep the phases in Measured
    tmp <- tmp[tmp$recent_sample %in% phases,]
    wtvars <- wtvars[wtvars$sample %in% phases, ]
    # Obtain list of files
    datafiles <- file.path(datapaths[as.character(tmp$recent_sample)], paste(tolower(tmp$NHANESfile), ".xpt", sep = ""))
    names(datafiles) <- tmp$NHANESfile
    demofiles <- file.path(datapaths[as.character(tmp$recent_sample)], demofiles[as.character(tmp$recent_sample)])
    #cat("demofiles: ", demofiles, "\n")
    names(demofiles) <- tmp$NHANESfile ## associate a demo file with each data file
    bwtfiles <- file.path(datapaths[as.character(tmp$recent_sample)],
                          tolower(wtvars$BWfile[match(paste(names(demofiles), ".XPT", sep = ""), wtvars$file)]))
    names(bwtfiles) <- tmp$NHANESfile
    creatfiles <- file.path(datapaths[as.character(tmp$recent_sample)],
                            tolower(wtvars$creatfile[match(paste(names(demofiles), ".XPT", sep = ""), wtvars$file)]))
    names(creatfiles) <- tmp$NHANESfile
    chemwt <- wtvars$wtvariable
    names(chemwt) <- tolower(gsub(".XPT", "", wtvars$file))
  }

  # Check to see if any chemicals use units we haven't accounted for
  lunitscale <- log(c("ug/L" = 1.0,"ng/mL" = 1.0,"ng/L" = 0.001, "pg/mL" = 0.001))
  tmp <- setdiff(unique(Measured$units), names(lunitscale))
  if (length(tmp) > 0) {
    print("Error:  incorrect units (must be either ug/L, ng/mL, ng/L, or pg/mL")
    print(paste("Problem units:  ", tmp, sep = ""))
    stop()
  }

  ### Run doplots
  Measured$Prabove <- (Measured$Sample_Size - Measured$BelowLOD)/Measured$Sample_Size
  indx <- order(Measured$Prabove)
  Measured <- Measured[indx,]

  J <- 1:nrow(Measured)

  #outplots <- mclapply(J, doplots, dsgn=Measured, mc.preschedule=TRUE, mc.cores=11)
  print("Run doplots function for each chemical")
  registerDoMC(cores = 10)
  outplots <- foreach(i = J) %dopar% {
    tmp <- doplots(i, dsgn = Measured, demofiles, datafiles, chemwt, bwtfiles, creatfiles, codes_table = convtbl)
    return(tmp)
  }

  Rightseg <- sapply(outplots, function(z) z$Rightseg)
  Measured$Rightseg <- c(Rightseg, numeric(nrow(Measured) - length(J)))
  Leftseg <- sapply(outplots, function(z) z$Leftseg)
  Measured$Leftseg <- c(Leftseg, numeric(nrow(Measured) - length(J)))


  ### Examine error
  Asymerrordata <- t(sapply(outplots, function(z) c(z$se.robust, z$se.model, z$Leftseg/z$Rightseg)))
  colnames(Asymerrordata) <- c("se.robust","se.model","LeftRightRatio")
  Asymerrordata <- as.data.frame(Asymerrordata)
  Asymerrordata$Prabove <- Measured$Prabove
  Asymerrordata$AboveLOD <- Measured$Sample_Size - Measured$BelowLOD

  print(paste("Plotting standard error and LOD probability in the file leftRightRatio_",
              format(Sys.time(), "%Y-%m-%d"), ".pdf", sep = ""))
  p <- ggplot(data=Asymerrordata, aes(x=LeftRightRatio, y=se.robust, size=Prabove, color=AboveLOD)) +
    geom_point() + scale_x_log10("Left:Right Ratio") + scale_y_log10("Robust SE") + geom_vline(xintercept=2)
  pdf(file.path(save_directory, paste("leftrightratioplot", format(Sys.time(), "%Y-%m-%d"), ".pdf", sep = "")),
                width = 10, height = 10)
  print(p)
  dev.off()

  Asymerrordata$asuniform <- factor(with(Asymerrordata, Prabove < 0.2 | se.robust > 0.5))

  print(paste("Plotting error distribution in the file AsymmetryvsSE_",
              format(Sys.time(), "%Y-%m-%d"), ".pdf", sep = ""))
  p <- ggplot(data=Asymerrordata, aes(x=LeftRightRatio, y=se.robust, color=asuniform)) + geom_point(alpha=0.5) +
    scale_x_log10("Left:Right Ratio") + scale_y_log10("Robust SE") + geom_vline(xintercept=2, lty=3) +
    annotate("text", label="Use Uniform if Pr{above LOD} < 0.2 or \nrobust SE > 0.5", x=0.1, y=1, hjust=0)
  pdf(file.path(save_directory, paste("AsymmetryvsSE_", format(Sys.time(), "%Y-%m-%d"), ".pdf", sep = "")),
                width=11, height=11)
  print(p)
  dev.off()

  Measured$PosteriorShape <- factor(Asymerrordata$asuniform, labels=c("normal","uniform"))

  print(paste("Saving outputs to NewMeasured_", format(Sys.time(), "%Y-%m-%d"), ".RData", sep = ""))
  save(Measured, file = file.path(save_directory, paste("NewMeasured_", format(Sys.time(), "%Y-%m-%d"), ".RData", sep = "")))
  return(Measured)

}


#' doplots
#'
#' Derivative of the loglike function.  Obtain the gradient
#' to help in calculating variance estimates.
#'
#' @param j Row number of dsgn
#' @param dsgn Measured (output from running the readNHANES function)
#'             with the Prabove column, sorted by this column.
#' @param demofiles Character vector of demographic files in the codes file
#' @param datafiles Character vector of NHANES chemical files in the codes file
#' @param chemwt Character vector of chemical weight files in the codes file
#' @param bwtfiles Character vector of bodyweight files in the codes file
#' @param creatfiles Character vector of creatinine files in the codes file
#' @param codes_table R dataframe of the first sheet of the codes_file in the 
#'                    readNHANES() function reduced to the cohorts of interest
#'                    (those contained in Measured or dsgn here)
#'
#' @import survey
#' @import ggplot2
#' @import stats
#'
#' @return outplots a list with length equal to the number of chemicals,
#'         each list element being a list with 5 components.
#'         1. plot: ???
#'         2. Leftseg: left interval of curve around expected
#'         mean of distribution
#'         3. Rightseg:right interval of curve around expected
#'         mean of distribution
#'         4. se.robust: ???
#'         5. se.model: ???
#'         Measured with 3 columns added.
#'         1. Prabove: fraction of measurements above the LOD.
#'         2. Rightseg: right interval of curve around expected
#'         mean of distribution
#'         3. Leftseg: left interval of curve around expected
#'         mean of distribution
#'
#'
#' @export
#'
#'
doplots <- function(j, dsgn, demofiles, datafiles, chemwt, bwtfiles, creatfiles, codes_table) {
  if (is.na(dsgn$loggm_se[j]))
    return(list(plot=paste("Row",j,"loggm_se is NA"),
                Rightseg=NA,
                Leftseg=NA,
                se.robust=NA,
                se.model=NA))

  lunitscale <- log(c("ug/L" = 1.0,"ng/mL" = 1.0,"ng/L" = 0.001, "pg/mL" = 0.001))
  measurehead <- "URX"
  lodindhd <- "URD"
  lodtail <- "LC"

  NM <- dsgn$NHANEScode[j]
  meascore <- sub(paste("^",measurehead,"(.+)","$",sep=""),"\\1",NM)
  nmlc <- paste(lodindhd,meascore,lodtail,sep="")

  if (class(demofiles) == "list"){
    # There is combined data (more than 1 file per chemical).
    # Do getDesign per chemical instead of per file
    chem <- dsgn$CAS[j]
    design0 <- getDesign(demof=demofiles[[chem]], chemdtaf=datafiles[[chem]], chem2yrwt=chemwt[[chem]],
                         nm=NM, CreatFun=creatinine, bodywtfile=bwtfiles[[chem]],
                         creatfile=creatfiles[[chem]], codes_table = codes_table)
  } else {
    # Run per file (multiple chemicals per file)
    datafile <- dsgn$NHANESfile[j]
    design0 <- getDesign(demof=demofiles[datafile], chemdtaf=datafiles[datafile], chem2yrwt=chemwt[tolower(datafile)],
                         nm=NM, CreatFun=creatinine, bodywtfile=bwtfiles[datafile],
                         creatfile=creatfiles[datafile], codes_table = codes_table)
  }

  keep <- switch(dsgn$subpop[j],
                 Total=rep(TRUE, nrow(design0$variables)),              # This and next line reduces design to the age group
                 Male=design0$variables$RIAGENDR == "Male",             # for this row in measured
                 Female=design0$variables$RIAGENDR == "Female",
                 `0 - 5` = design0$variables$AgeGroup == "0 - 5",
                 `6 - 11 years` = design0$variables$AgeGroup == "6 - 11 years",
                 `12 - 19 years` = design0$variables$AgeGroup == "12 - 19 years",
                 `20 - 65 years` = design0$variables$AgeGroup == "20 - 65 years",
                 `66 years and older` = design0$variables$AgeGroup == "66 years and older",
                 `BMI > 30` = design0$variables$Obesity == "BMI > 30",
                 `BMI <= 30` = design0$variables$Obesity == "BMI <= 30",
                 ReproAgeFemale = design0$variables$ReproAgeFemale)
  design <- subset(design0, keep)
  ## if we don't have any above LOD, then stop right here and move on
  if (dsgn$BelowLOD[j] == dsgn$Sample_Size[j])
    return(list(plot=paste("Row",j,"all below LOD"),
                Rightseg=NA,
                Leftseg=NA,
                se.robust=NA,
                se.model=NA))

  forms <- list(mean=as.formula(paste("I(cbind(Measure,",nmlc,")) ~ 0 + ones",sep="")),
                lsd=~ 0 + ones)
  ## First try gradient free
  xmean <- (dsgn$loggm[j] + log(dsgn$MW[j])) - lunitscale[dsgn$units[j]]
  Start <- list(mean=unname(xmean), lsd=unname(dsgn$lsdlog[j]))
  out1a <- try(svymle(lnlike, design=design, formulas=forms,
                      start=Start),
               silent=TRUE)
  if (is(out1a, "try-error"))
  {
    attr(out1a,"Sample") <- paste("out1a, row",j)
    return(list(plot=out1a,
                Rightseg=NA,
                Leftseg=NA,
                se.robust=NA,
                se.model=NA))
  }
  Start <- list(mean=unname(out1a$par[1]), lsd=unname(out1a$par[2]))
  out1b <- try(svymle(lnlike, gradient=gcens, design=design, formulas=forms, method="BFGS",
                      start=Start, control=list(maxit=500, reltol=1e-10)),
               silent=TRUE)
  if (is(out1b, "try-error"))
  {
    attr(out1b, "Sample") <- paste("out1b, row",j)
    return(list(plot=out1b,
                Rightseg=NA,
                Leftseg=NA,
                se.robust=NA,
                se.model=NA))
  }
  se <- sqrt(diag(vcov(out1b))[1])
  se.model <- sqrt(diag(vcov(out1b, stderr="model"))[1])
  ### Fit lognormal distribution with left censoring to a dataset
  ### assume we are passing the value in x
  ### cens is 0 for non-censored, 1 for left-censored (below LOD, eg)
  ### llod is the log-limit of detection

  weights <- 1/design$prob
  wtotal <- sum(weights) / nrow(design$variables)

  Start <- list(lsd=unname(out1b$par[2]))
  zz <- unname(out1b$par[1])

  formsR <- list(as.formula(paste("~ I(cbind(Measure,",nmlc,"))")), lsd=~ 0 + ones)

  tmpstart <- svymle(lnlike, gradient=gcens, design=design, formulas=formsR, method="BFGS",
                     start=Start, control=list(maxit=500, reltol=1e-10), mean=zz)
  valmax <- tmpstart$value / wtotal

  ## Step down 3 standard errors, taking 60 steps.

  StepsDn <- seq(zz, zz - 3*se, length=60)

  Downpll <- numeric(length(StepsDn))
  tmp <- tmpstart
  for (i in seq_along(Downpll)) {
    Start <- list(lsd=tmp$par[1])

    tmp <- svymle(lnlike, gradient=gcens, design=design, formulas=formsR, method="BFGS",
                  start=Start, control=list(maxit=500, reltol=1e-10), mean=StepsDn[i])
    Downpll[i] <- tmp$value / wtotal
  }

  ## Step up 3 standard errors, taking 200 steps
  zz <- unname(out1b$par[1])
  StepsUp <- seq(zz, zz + 3*se, length=200)

  Uppll <- numeric(length(StepsUp))
  tmp <- tmpstart
  for (i in seq_along(Uppll)) {
    Start <- list(lsd=tmp$par[1])

    tmp <- svymle(lnlike, gradient=gcens, design=design, formulas=formsR, method="BFGS",
                  start=Start, mean=StepsUp[i])
    Uppll[i] <- tmp$value / wtotal
  }

  pLL <- c(rev(Downpll), Uppll[-1])

  lgm <- c(rev(StepsDn), StepsUp[-1])
  Vlod <- dsgn$LOD[j]
  Vlod <- Vlod * dsgn$MW[j]
  Vlod <- Vlod / exp(lunitscale[dsgn$units[j]])

  pdta <- data.frame(x=unname(lgm), y=unname(pLL))
  ## Get left and right interval, where the curve defined in pdta crosses the horizontal
  ## line defined by qchisq(0.95,1)/2 down from valmax, which should be max(y)
  ## splitpoint is at
  Nsplit <- length(Downpll)
  Xcrit <- max(pdta$y) - qchisq(0.95, 1)/2
  ## right side
  tpdta <- pdta[Nsplit:(nrow(pdta)),]
  Rightseg <-
    if (any(tpdta$y < Xcrit))
    {
      tpdta$y <- tpdta$y - Xcrit
      fn <- splinefun(tpdta)
      ## Bracket the solution
      x0 <- min(tpdta$x)
      x1 <- tpdta$x[min(which(tpdta$y < 0))]
      zout <- uniroot(fn, c(x0,x1))
      zout$root
    }
  else Inf
  Rightseg <- Rightseg - out1b$par[1]

  tpdta <- pdta[1:Nsplit,]
  Leftseg <-
    if (any(tpdta$y < Xcrit))
    {
      tpdta$y <- tpdta$y - Xcrit
      fn <- splinefun(tpdta)
      ## Bracket the solution
      x1 <- max(tpdta$x)
      x0 <- tpdta$x[max(which(tpdta$y < 0))]
      zout <- uniroot(fn, c(x0,x1))
      zout$root
    }
  else Inf
  Leftseg <- out1b$par[1] - Leftseg

  p <- ggplot(data=pdta, aes(x=x, y=y)) +
    geom_line() + geom_hline(yintercept=valmax - qchisq(0.95,1)/2) +
    geom_vline(xintercept=out1b$par[1]) +
    geom_vline(xintercept = out1b$par[1] + se * c(-1.96,1.96), lty=3) +
    geom_vline(xintercept = out1b$par[1] + se.model * c(-1.96,1.96), lty=2, col="cyan") +
    geom_vline(xintercept = log(Vlod), col="red") +
    ggtitle(paste(j,
                  "se:",signif(dsgn$loggm_se[j], digits=2),
                  "LOD:",signif(log(Vlod), digits=2),
                  "Above LOD:", dsgn$Sample_Size[j] - dsgn$BelowLOD[j],
                  "Pr{Above LOD}:", signif(dsgn$Prabove[j], digits=2)))
  plots <- list(plot=p,
                Leftseg=Leftseg,
                Rightseg=Rightseg,
                se.robust=se,
                se.model=se.model)
  return(plots)
}


