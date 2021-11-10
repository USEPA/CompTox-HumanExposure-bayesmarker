#' getNhanesQuantiles
#'
#' Population Quantiles, MLEs of lognormal population parameters
#' w/standard errors, and comparisons of fitted lognormal curve to
#' empirical CDF for data.
#'
#' This code attempts to reproduce the Phthalate quantiles and
#' confidence limits reported in 'Fourth National Exposure Report,
#' Updated Tables, March 2013', estimate parameters for an assumed
#' lognormal population distribution of the urinary metabolites in
#' mg/day excretion rates, and compare the estimated lognormal
#' distribution to the estimated population CDF for the sample data,
#' all using proper adjustments for the NHANES sampling design.
#'
#' Function for reading relevant files, merging them, and estimating
#' quantiles for urine concentrations
#'
#'
#' @param demof String containing name of SAS transport file containing
#'           demographic data
#' @param chemdtaf String containing name of SAS transport file containing
#'           the urine concentrations for a group of chemicals.
#'           presumed to have cols SEQN, WTSB2YR, and pairs of column for each chemical: one
#'           for the measurement and one for an indicator for < LOD.  The value in the measurement
#'           when the measurement is < LOD is LOD/sqrt(2), rounded to about 2 digits.
#' @param measurehead String matching the initial part of the measurement variable
#'           names.  The defaults are correct for files for laboratory measurement of
#'           products in urine.
#' @param measuretail Strings matching the final part of the measurement variable
#'           names.  The defaults are correct for files for laboratory measurement of
#'           products in urine.
#' @param lodindhd Regular expression for matching the 'comment' fields (lod indicator), initial part
#' @param lodtail Regular expressions for matching the 'comment' fields (lod indicator), final part
#' @param seq ID, used to merge different files (assumed to be the same variable name across files)
#' @param demoageyr Name of the variable giving age in years in the demographic data file.
#' @param demogendr Name of the variable giving gender in the demographic data file.
#' @param demoeth Name of the variable giving ethnicity in the demographic data file.
#' @param PSU Name of the variable in the demographic data file giving the sampling unit
#' @param STRA Name of the variable in the demographic data file giving the stratum for each observation
#' @param creatinine Name of the variable in the file named by chemdataf that contains the creatinine concentration
#'             in each urine sample
#' @param CreatFun Function to compute median or random daily estimates of creatinine excretion.  Used to convert
#'           outputs scaled by creatinine concentration to daily excretion rates. The argument is a data frame with variables:
#'           RIAGENDR, RIDRETH1, BMXWT, and RIDAGEYR
#' @param creatfile Name of file with creatinine measurements
#' @param bodywt Variable giving bodyweight in kg
#' @param bodywtcomment Variable flagging special cases for bodyweight
#' @param MECwt Variable giving the weight to be used when analyzing a full MEC variable
#' @param chem2yrwt Variable in chemdataf giving the sampling weight to be used for those data
#' @param chemvars Variables in chemdataf that are to be analyzed
#' @param urinefile File containing basic information about urine rate and volume.  Only relevant to 2009 and
#'            later
#' @param urinerate Variable in urinefile giving the rate of urine production
#' @param bodywtfile File giving the body weight information
#' @param bodymassindex Code for BMI information in the NHANES file.  Default is "BMXBMI".
#' @param Q Vector of percentiles to estimate (e.g. 50 = 50th percentile)
#' @param plot If not NULL, the name to give a pdf file for plotting CDFs, empirical and lognormal
#'       from parameter estimates
#' @param codes_table R dataframe of the first sheet of the codes_file in the readNHANES() function
#' @param lognormfit Logical, fit a lognormal distribution to data (including the subpopulations)
#' @param code Information for adding chemical names and CAS numbers to the output.  Is a list with elements
#'         table: a character matrix or data frame with character (not
#'                factor) elements giving the names NHANES uses for
#'                chemical measures, the chemical names, and the CAS
#'                numbers as strings.
#'         codename: the name in table of the variable that contains the codes NHANES uses
#'         chemame: tha name in table of the variable containing the chemical names
#'         CAS: either NULL or the name in table of the variable containing the CAS numbers
#' @param LODfilter should the estimated quantiles be filtered so they are NA for values < LOD?  This
#'           is for matching NHANES summary reports.  Only relevant to quantile estimates and their
#'           confidence limits
#' @param MaximumAge What is the oldest age in years to consider?
#'
#'
#' @import survey
#' @importFrom methods is
#' @importFrom stats as.formula dnorm na.omit pnorm vcov
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom graphics lines points segments
#' @importFrom foreign read.xport
#'
#'
#' @return A dataframe with columns:
#'   NHANEScode: nhanes variable name
#'   Chem: chemical name (if !is.null(code) in arguments)
#'   CAS: CAS number (if !is.null(code) and !is.null(code$CAS))
#'   subpop: subpopulation identifier: Total, Male, Female, 0 - 5, 6 - 11
#'           years, 12 - 19 years, 20 - 65 years, 66 years and older,
#'           ReproAgeFemale, BMI <= 30, BMI > 30
#'   then triplets of variables for each quantile requested.  For example, if Q is c(50, 90),
#'   then Q_50, Q_50_lcl, Q_50_ucl, Q_90, Q_90_lcl, Q_90_ucl
#'   LOD: the limit of detection for each variable
#'   if lognormfit, then
#'   loggm: natural log of the geometric mean (R plnorm parameter meanlog)
#'   loggm_se: standard error of the estimate of loggm; NA if not
#'             computable (in which case, probably loggm and lsdlog are not unique).
#'   lsdlog: the natural log of the standard deviation of the log-transformed variable
#'           (the natural logarithm of the R plnorm parameter sdlog).
#'   lsdlog_se: the standard error of the estimate of lsdlog.  See loggm_se for meaning of NA.
#'
#'   In addition, if a value is given for plot, a pdf file is produced
#'   with estimated population cdfs (i.e., taking account of the NHANES
#'   sampling design), with superimposed quantile estimates and
#'   confidence limits.  If lognormfit is true, then the fitted curve is
#'   superimposed.
#'
#' @export
#'
#' @examples
#' # indata <- c(
#' # "Mono(carboxyoctyl) phthalate;URXCOP
#' # Mono (carboxynonyl) phthalate;URXCNP
#' # Mono-n-butyl phthalate;URXMBP
#' # Mono-isobutyl phthalate;URXMIB
#' # Mono-n-methyl phthalate;URXMNM
#' # Mono-benzyl phthalate;URXMZP")
#' # con <- textConnection(indata)
#' # tbl <- read.delim(con, sep = ";", header = FALSE, col.names = c("Chem", "varname"),
#' #                  stringsAsFactors = FALSE)
#' # close(con)
#'
#' # inpath <- file.path(find.package("bayesmarker"), "extdata")
#' # phthalatesQ <- getNhanesQuantiles(demof = file.path(inpath, "DEMO_F.XPT"),
#' # chemdtaf = file.path(inpath, "PHTHTE_F.XPPPPT"), lognormfit = TRUE, plot = "phthalates.pdf",
#' # code = list(table = tbl, chemname = "Chem", codename = "varname", CAS = NULL))
#'
getNhanesQuantiles <- function(demof="DEMO_F.XPT", chemdtaf, measurehead="URX", measuretail=NULL,
                               lodindhd="URD", lodtail="LC", seq="SEQN",
                               demoageyr="RIDAGEYR",demogendr="RIAGENDR",demoeth="RIDRETH1",
                               PSU="SDMVPSU",
                               STRA="SDMVSTRA", chem2yrwt=NULL,chemvars,
                               creatinine="URXUCR", CreatFun=NULL, creatfile=NULL,
                               urinerate=NULL,urinefile=NULL,bodywtfile=NULL,
                               bodywt="BMXWT",bodywtcomment="BMIWT",
                               bodymassindex="BMXBMI",
                               MECwt="WTMEC2YR",
                               codes_table = NULL,
                               Q=c(50,75,90,95), lognormfit=TRUE, plot=NULL, code=NULL,
                               LODfilter=TRUE, MaximumAge = 150) {
  ## ----------------------------------------------------------------

  # For scaling data before 2007
  oldmethod <- c("1999-2000", "2001-2002", "2003-2004", "2005-2006")
  scale_old <- function(z) {
    if (!is.na(z)) {
      if (z < 75) {
        z <- (1.02*sqrt(z) - 0.36)^2
      } else if (z >= 75 && z < 250) {
        z <- (1.05*sqrt(z) - 0.74)^2
      } else {
        z <- (1.01*sqrt(z) - 0.1)^2
      }
    } else {
      z <- NA
    }
    return(z)
  }

  unitscale <- c("ug/L" = 1.0,"ng/mL" = 1.0,"ng/L" = 0.001, "pg/mL" = 0.001)
  
  ## Sanity Checks
  if (length(demof) > 1){

    demo <- c()
    cdta <- c()
    bwt <- c()
    for (i in 1:length(demof)){
      if (!file.exists(demof[i])) stop(paste(demof[i], "not found"))
      miss_test <- chemdtaf[i]
      if (missing(miss_test) || !file.exists(chemdtaf[i])) {
        stop(paste("File for chemical data", chemdtaf[i], "not found"))
      }
      tmp <- read.xport(demof[i])
      demovars <- c(seq,PSU,STRA,demoageyr,demogendr,demoeth,MECwt)
      if (any(!(XX <- demovars %in% colnames(tmp)))) {
        stop(paste("Variables requested NOT in ",demof[i],": ",
                   paste(demovars[!XX],collapse=", "), sep=""))
      }
      tmp <- tmp[,colnames(tmp) %in% demovars]

      # Deal with creatinine: get file if column not there and scale old data
      tmp2 <- read.xport(chemdtaf[i])
      z <- grep(toupper(gsub(".xpt", "", basename(chemdtaf[i]))), codes_table$NHANESfile, fixed = TRUE)
      y <- unname(z[which(codes_table$NHANEScode[z] == chemvars)])
      tmp2[,chemvars] <- tmp2[,chemvars] * unitscale[codes_table$units[y]]
      chemvars2 <- c(seq, chem2yrwt[i], chemvars)
      if (creatinine %in% colnames(tmp2)) {
        chemvars2 <- c(chemvars2, creatinine)
        scalebycreatinine <- TRUE
      } else if (!(creatinine %in% colnames(tmp2)) & !is.null(creatinine)) {
        if (is.null(creatfile[i])) stop("No creatinine file provided")
        if (!file.exists(creatfile[i])) stop(paste(createfile[i], "not found"))
        creat <- read.xport(creatfile[i])
        tmp2 <- merge(tmp2, creat[,c(seq, creatinine)], all.x = TRUE)
        chemvars2 <- c(chemvars2, creatinine)
        scalebycreatinine <- TRUE
      } else {
        scalebycreatinine = FALSE
      }

      ind <- unlist(sapply(oldmethod, function(x) grep(x, demof[i])))
      # Apply if this phase is in oldmethod
      if (length(ind) > 0) {
        tmp2[,creatinine] <- sapply(tmp2[,creatinine], scale_old)
      }

      ##  cat("getNhanesQuantiles, chemvars2: \"", chemvars2, "\"\n")
      if (any(!(XX <- chemvars2 %in% colnames(tmp2)))) {
        stop(paste("Variables requested NOT in ",chemdtaf[i],": ",
                   paste(chemvars2[!XX], collapse=", "), sep=""))
      }
      tmp2 <- tmp2[,colnames(tmp2) %in% chemvars2]
      if (i > 1) {
        colnames(tmp2)[colnames(tmp2) == chem2yrwt[i]] <- chem2yrwt[1]
      }

      if (!is.null(urinefile[i])) {
        if (is.null(urinerate[i]))
          stop("File for urinerate specified without variable name for urine rate.")
        if (!file.exists(urinefile[i])) stop(paste(urinefile[i],"not found."))
        urine <- read.xport(urinefile[i])
        if (!urinerate %in% colnames(urine[i]))
          stop(paste(urinrate[i],"not in",urinefile[i]))
        scalebyurinerate <- TRUE
      } else {
        scalebyurinerate <- FALSE
      }
      ## scalebyurinerate dominates scalebycreatinine, so set the latter
      ## to false if we are scaling by urinerate
      if (scalebyurinerate) scalebycreatinine <- FALSE

      if (!is.null(bodywtfile[i])) {
        if (is.null(bodywt))
          stop("File for bodywt specified without variable name for body weight.")
        if (is.null(bodywtcomment))
          stop("Must specify name for bodywtcomment")
        if (is.null(MECwt))
          stop("Must specify name for Exam Weights")
        if (!file.exists(bodywtfile[i])) stop(paste(bodywtfile[i], "not found."))
        tmp3 <- read.xport(bodywtfile[i])
        if (any(!(XX <- c(bodywt, bodywtcomment,bodymassindex) %in% colnames(tmp3))))
          stop(paste(paste(c(bodywt,bodywtcomment,bodymassindex)[XX],collapse=", "),"not in",bodywtfile[i]))
        tmp3 <- tmp3[,colnames(tmp3) %in% c(seq, bodywt, bodywtcomment, bodymassindex)]
        scalebybodywt <- TRUE
      } else {
        scalebybodywt <- FALSE
        bodywt <- NULL
        bodywtcomment <- NULL
        MECwt <- NULL
      }

      demo <- rbind(tmp, demo)
      cdta <- rbind(tmp2, cdta)
      bwt <- rbind(tmp3, bwt)
    }
    chem2yrwt <- chem2yrwt[1]

    # Deal with combined weights
    demo[,MECwt] <- demo[,MECwt]/length(demof)


  } else {

    if (!file.exists(demof)) stop(paste(demof, "not found"))
    if (missing(chemdtaf) || !file.exists(chemdtaf)) {
      stop(paste("File for chemical data",chemdtaf, "not found"))
    }
    demo <- read.xport(demof)
    cdta <- read.xport(chemdtaf)
    
    # Scale units
    z <- sapply(codes_table$NHANESfile, function(x) grep(tolower(x), chemdtaf))
    z <- which(lengths(z) > 0)
    y <- unname(z[which(codes_table$NHANEScode[z] == chemvars)])
    cdta[,chemvars] <- cdta[,chemvars] * unitscale[codes_table$units[y]]
    
    demovars <- c(seq,PSU,STRA,demoageyr,demogendr,demoeth,MECwt)
    if (any(!(XX <- demovars %in% colnames(demo)))) {
      stop(paste("Variables requested NOT in ",demof,": ",
                 paste(demovars[!XX],collapse=", "), sep=""))
    }


    # For newer NHANES phases, need to get creatinine column from creatinine file.
    # Otherwise just use the creatinine column in the metabolite measurement file.
    chemvars2 <- c(seq, chem2yrwt, chemvars)
    if (creatinine %in% colnames(cdta)) {
      chemvars2 <- c(chemvars2, creatinine)
      scalebycreatinine <- TRUE
    } else if (!is.null(creatinine)) {
      creat <- read.xport(creatfile)
      cdta <- merge(cdta, creat[,c(seq, creatinine)], all.x = TRUE)
      chemvars2 <- c(chemvars2, creatinine)
      scalebycreatinine <- TRUE
    } else {
      scalebycreatinine <- FALSE
    }

    # NHANES suggests scaling creatinine measurements before 2007 using a pairwise
    # equation to be more comparable.  Do this here.
    oldmethod <- c("1999-2000", "2001-2002", "2003-2004", "2005-2006")
    ind <- unlist(sapply(oldmethod, function(x) grep(x, demof)))
    # Apply if this phase is in oldmethod
    if (length(ind) > 0) {
      cdta[,creatinine] <- sapply(cdta[,creatinine], scale_old)
    }

    ##  cat("getNhanesQuantiles, chemvars2: \"", chemvars2, "\"\n")
    if (any(!(XX <- chemvars2 %in% colnames(cdta)))) {
      stop(paste("Variables requested NOT in ",chemdtaf,": ",
                 paste(chemvars2[!XX], collapse=", "), sep=""))
    }

    if (!is.null(urinefile)) {
      if (is.null(urinerate))
        stop("File for urinerate specified without variable name for urine rate.")
      if (!file.exists(urinefile)) stop(paste(urinefile,"not found."))
      urine <- read.xport(urinefile)
      if (!urinerate %in% colnames(urine))
        stop(paste(urinrate,"not in",urinefile))
      scalebyurinerate <- TRUE
    } else {
      scalebyurinerate <- FALSE
    }
    ## scalebyurinerate dominates scalebycreatinine, so set the latter
    ## to false if we are scaling by urinerate
    if (scalebyurinerate) scalebycreatinine <- FALSE

    if (!is.null(bodywtfile)) {
      if (is.null(bodywt))
        stop("File for bodywt specified without variable name for body weight.")
      if (is.null(bodywtcomment))
        stop("Must specify name for bodywtcomment")
      if (is.null(MECwt))
        stop("Must specify name for Exam Weights")
      if (!file.exists(bodywtfile)) stop(paste(bodywtfile, "not found."))
      bwt <- read.xport(bodywtfile)
      if (any(!(XX <- c(bodywt, bodywtcomment,bodymassindex) %in% colnames(bwt))))
        stop(paste(paste(c(bodywt,bodywtcomment,bodymassindex)[XX],collapse=", "),"not in",bodywtfile))
      scalebybodywt <- TRUE
    } else {
      scalebybodywt <- FALSE
      bodywt <- NULL
      bodywtcomment <- NULL
      MECwt <- NULL
    }
  }

  if (length(MaximumAge) == 1) MaximumAge <- rep(MaximumAge, length(chemvars))
  ## End of sanity checks
  ## ----------------------------------------------------------------
  ## Scaling?
  ## Possibilities:
  ##   - no scaling.  chemvars from chemdataf processed as is: creatinine null and urinefile null.
  ##   - scale by creatinine mg/DL.  Convert creatinine to g/L (divide it by 100), then
  ##     convert chemvars from mass unit/ml -> mass unit/mg creatinine: chemvars/creatinine
  ##     creatinine Not null and urinefile null.  If CreatFun is not null, then the result of the
  ##     previous computation is multiplied by the estimate of daily creatinine excretion from CreatFun
  ##     to give mass units / day.
  ##   - scale by urine flow rate (in ml/min) by multiplying chemvars by urinerate.  This gives
  ##     mass unit / min.  Multiply by 24*60 min/day to get mass unit / day, which is
  ##     what is analyzed
  ##     and reported.  Happens when urinefile is not null, regardless of creatinine.
  ##   - bodywt scaling: in addition to the above, of bodywtfile is not null, the result
  ##     of the previous
  ##     scaling is divided by bodywt to give mass units / mg creatinine / kg bodywt Or
  ##     mass unit / day / kg bodywt or mass units / kg bodywt / day
  if (doPlot <- !is.null(plot)) {
    pdf(file=plot)
    on.exit(dev.off())
  }
  ## -----------------------------------------------------------------
  ##   Demographics

  ## Interesting demographic variables are:
  ## demoageyr: age in years
  ## demogendr: gender: male=1, female=2
  ## demoeth: Race/Ethnicity:
  ##        Mexican American=1
  ##        Other Hispanic=2
  ##        Non-Hispanic White=3
  ##        Non-Hispanic Black=4
  ##        Other Race-Including Multi-Racial=5

  ## Select out the variables we'll need going forward
  demo <- demo[,c(seq,PSU,STRA,demoageyr,demogendr,demoeth,MECwt)]
  ## Set up gender, age, and ethnicity as factors using the same levels as the NHANES reports
  demo[,demogendr] <- factor(demo[,demogendr], labels=c("Male","Female"))
  demo$AgeGroup <- cut(demo[,demoageyr], breaks=c(-1,5.5,11.5,19.5,65.5, 100.5),
                       labels=c("0 - 5","6 - 11 years","12 - 19 years", "20 - 65 years", "66 years and older"))
  demo$RaceEthn <- factor(demo[,demoeth],
                          labels=c("Mexican American","Other Hispanic","Non-Hispanic White",
                                   "Non-Hispanic Black","Other"))
  demo$ChildBearingAgeFemale <- factor(demo[,demogendr] == "Female" & (demo[,demoageyr] >= 16 & demo[,demoageyr] <= 49),
                                       labels=c("NotReproAgeFemale","ReproAgeFemale"))
  ## -----------------------------------------------------------------
  ## Chemical Data


  ## We keep from cdta seq, chem2yrwt, all the chemvars, and the
  ## comment variables for the chemvars, if they exist.  Not all
  ## chemvars have comment variables, so they need to be constructed.
  ## create the lodindicator variable names for all chemvars, then,
  ## construct the ones that don't exist.  We can identify missing
  ## values because they are the smallest variable in the variable,
  ## will probably have multiple instances, and will be approximately
  ## the second smallest in the variable / sqrt(2) (rounded to 2
  ## significant digits).
  meascore <- sub(paste("^",measurehead,"(.+)",measuretail,"$",sep=""),"\\1",chemvars)
  LODnames <- paste(lodindhd,meascore,lodtail,sep="")
  ## nms is a data frame with a record for each chemical.  It keeps
  ## track of the names of the measurement and LODind variables, the
  ## LOD, and space for the chemical name and CAS. LOD here is the
  ## maximum of the individual level LODs.
  nms <- data.frame(Measurement=chemvars,
                    LODind=LODnames,
                    LOD=numeric(length(chemvars)),
                    Chem = character(length(chemvars)),
                    CAS = character(length(chemvars)),
                    stringsAsFactors=FALSE)
  ## loop through the variables and find the lods.
  ## Need to confirm that there is only one value associated with an lod,
  ## then, the next to lowest value in the variable is the lod.
  for (i in seq_len(nrow(nms))) {
    nm <- nms$Measurement[i]
    nmlc <- nms$LODind[i]
    ## Do we already have an LOD indicator?
    if (!nmlc %in% colnames(cdta)) {
      ## We have to create the LOD indicator.  In this case, there is only 1 LOD
      ## and it will be sqrt(2) times the smallest value
      zlod <- sort(unique(cdta[,nm]))[1] * sqrt(2)
      cdta[,LODnames[i]] <- ifelse(cdta[,nm] >= zlod, 0, 1)
    }
  }
  ## ---------------------------------------------------------------
  ## Merge data sets
  ##
  ## Merge demo and cdta on seq variables, retaining all the variables still in demo,
  ## and in cdta: seq, chem2yrwt, chemvars, all LODind, and creatinine if it exists
  ## To be safe, if we are going to scale by urinerate, then reset creatinine to NULL.
  ## This protects us against unnecessary losses due to NAs in creatinine.
  if (scalebyurinerate) creatinine <- NULL
  alldata <- merge(demo, cdta[,c(seq, chem2yrwt, chemvars,LODnames, creatinine)],
                   by.x=seq, by.y=seq, all.y=TRUE)
  ## If scalebyurinerate, merge alldata and urine[,c(seq, urinerate)]
  if (scalebyurinerate) {
    alldata <- merge(alldata, urine[,c(seq, urinerate)],
                     by.x=seq, by.y=seq, all.x=TRUE)
  }
  ## if scalebybodywt, merge alldata and bwt[,c(seq, bodywt)]
  if (scalebybodywt) {
    ## fixup bodywt to include missings due to special codes in
    ## bodywtcomment.  Basically, set bodywt to NA for any non-missing bodywtcomment
    bwt[!is.na(bwt[,bodywtcomment]),bodywt] <- NA
    ## Create obesity factor from bodymassindex
    bwt$Obesity <- cut(bwt[,bodymassindex], breaks=c(-0.5,30,500),
                       labels=c("BMI <= 30", "BMI > 30"))
    alldata <- merge(alldata, bwt[,c(seq, bodywt, bodymassindex, "Obesity")],
                     by.x=seq, by.y=seq, all.x=TRUE)
    ## Impute missing values for bodywt based on age, gender, and ethnicity
    if (any(is.na(alldata[,bodywt]))) {
      ##writeLines(paste(sum(is.na(alldata[,bodywt])),"missing values in",bodywt))
      ##browser()
      dta <- merge(demo, bwt[,c(seq, bodywt)], by.x=seq, by.y=seq, all.x=TRUE, all.y=TRUE)
      dsg <- na.omit(svydesign(ids=make.formula(PSU), strata=make.formula(STRA),
                               weights=make.formula(MECwt), nest=TRUE, data=dta))
      ## Just fit separate means by age and gender
      selectmales <- dsg$variables[,demogendr] == "Male"
      selectfemales <- dsg$variables[,demogendr] == "Female"
      males <- svyby(make.formula(bodywt), by=make.formula(demoageyr),
                     subset(dsg, selectmales), svymean)
      females <- svyby(make.formula(bodywt), by=make.formula(demoageyr),
                       subset(dsg, selectfemales), svymean)
      imp <- data.matrix(cbind(males[,bodywt],females[,bodywt]))
      ## Now, do the imputation
      ##browser()
      isnabw <- is.na(alldata[,bodywt])
      alldata[isnabw,bodywt] <- imp[data.matrix(alldata[isnabw,c(demoageyr, demogendr)])]
      ##writeLines(">>> Imputed bodyweights")
      ##print(alldata[isnabw, c(demoageyr, demogendr, bodywt)])
    }
  }
  ## Fixup factors so there are no missing levels
  for (nm in names(alldata)) {
    if (is.factor(alldata[,nm])) {
      alldata[,nm] <- factor(alldata[,nm])
    }
  }
  ## if !is.null(CreatFun), add the variable DailyCreatinine
  if (!is.null(CreatFun)) {
    alldata$DailyCreatinine <-
      CreatFun(data.frame(RIAGENDR = unname(alldata[,demogendr]),
                          RIDRETH1 = unname(alldata[,"RaceEthn"]),
                          BMXWT = unname(alldata[,bodywt]),
                          RIDAGEYR = unname(alldata[,demoageyr])))
  }

  ## --------------------------------------------------------------
  ## Estimation: Quantiles and distribution parameters
  ## Qnames: names for the quantiles and their confidence limits
  Qnames <- c(sapply(Q,
                     function(z) c(paste("Q_",z,sep=""),
                                   paste("Q_",z,"_lcl",sep=""),
                                   paste("Q_",z,"_ucl",sep=""))))

  ## Now, loop through the variables characterized by the rows in nms
  ## Fill a dataframe with the estimates and confidence limits
  ## qouts has columns:
  qoutnames <- c("NHANEScode", "Chem", "CAS", "subpop",
                 Qnames,
                 "LOD","Sample_Size","BelowLOD","loggm","loggm_se","lsdlog","lsdlog_se")
  ## Some indices to make it easier to find particular values (mainly used down
  ## in LODfilter
  Qcenindx <- match(Qnames[seq(1,length(Qnames), by=3)], qoutnames)
  Qlowindx <- match(Qnames[seq(2,length(Qnames), by=3)], qoutnames)
  Qhighindx <- match(Qnames[seq(3,length(Qnames), by=3)], qoutnames)
  Q <- Q/100

  qout <- list()
  for ( i in 1:4) qout[[i]] <- character(6*nrow(nms))

  for (i in 5:length(qoutnames)) qout[[i]] <- numeric(6*nrow(nms))

  names(qout) <- qoutnames

  qout <- as.data.frame(qout,stringsAsFactors=FALSE)

  N <- 3*length(Q) - 1
  irow <- 0
  for (i in seq_len(nrow(nms))) {
    nmlc <- nms[i,"LODind"]
    nm <- nms[i,"Measurement"]
    ## Set up data frame for measurement variable nm.  We keep this
    ## data frame as small as possible, so we don't lose records because
    ## of NAs in variables we don't use (a problem for nlm, used in svymle).
    keepvars <- c(colnames(demo),chem2yrwt, nm, nmlc, creatinine, urinerate, bodywt,"Obesity")
    if (!is.null(CreatFun)) keepvars <- c(keepvars, "DailyCreatinine")
    ndta <- alldata[,keepvars]
    ndta$Measure <- ndta[,nm]
    ## implement scaling options
    if (scalebycreatinine) {
      ndta$Measure <- ndta$Measure * 100 / ndta[,creatinine]
      if (!is.null(CreatFun)) {
        ndta$Measure <- ndta$Measure * ndta$DailyCreatinine
      }
    }
    if (scalebyurinerate) {
      ndta$Measure <- ndta$Measure * ndta[,urinerate] * 24*60 ## gives mass unit / day
    }
    if (scalebybodywt) {
      ndta$Measure <- ndta$Measure / ndta[,bodywt]
    }
    ## replace CDC's below LOD value with the LOD
    ## First, save the original value: that gets used for computing quantiles
    ndta$Measure2 <- ndta$Measure
    formla <- make.formula("Measure2")
    isBelowLOD <- !is.na(ndta[,nmlc]) & !is.na(ndta$Measure) & ndta[,nms[i,"LODind"]] == 1
    ndta$Measure[isBelowLOD] <- ndta$Measure[isBelowLOD]*sqrt(2)

    ## The LOD reported is the unweighted geometric mean of the LODs.
    nms$LOD[i] <-
      if (any(isBelowLOD)) {
        exp(mean(log(ndta$Measure[isBelowLOD])))
      } else {
        0.0
      }

    ## lMeasure is used to get initial distribution parameter estimates, and uses
    ## the log(LOD)/sqrt(2) value for < LOD values
    ndta$lMeasure <- log(ndta$Measure)
    ndta$lMeasure[isBelowLOD] <- ndta$lMeasure[isBelowLOD]/sqrt(2)
    ## lMeasure2 is the variable fit by the distribution function.
    ndta$lMeasure2 <- log(ndta$Measure)
    ## Column of ones to trick the modeling function.
    ndta$ones <- rep(1, nrow(ndta))
    ## Set up age filter here
    KeepbyAge <- ndta[,demoageyr] <= MaximumAge[i]
    ## Set up indicator for ReproAgeFemale
    ndta$ReproAgeFemale <- ndta$ChildBearingAgeFemale == "ReproAgeFemale"
    #print(paste("Number of NAs in weights:  ", sum(is.na(ndta[, chem2yrwt])), sep = ""))
    na_ind <- !is.na(ndta[chem2yrwt])
    design <-
      na.omit(subset(svydesign(ids=make.formula(PSU), strata=make.formula(STRA),
                               weights=make.formula(chem2yrwt), nest=TRUE, data=ndta[na_ind,]),KeepbyAge[na_ind]))

    ## =======================================================================
    ## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ## =======================================================================
    ## Do the estimates, a separate section for each subpopulation:
    ## Total, Males, Females, the four age groups, BMI < 30 and BMI
    ## > 30, reproductive age females

    ## -----------------------------------------------------------------------
    ## --------------------------------- Total -------------------------------
    ## -----------------------------------------------------------------------
    irow <- irow + 1

    out <- svyquantile(formla, design, Q, ci=TRUE, se=TRUE,
                       na.rm=TRUE, interval.type="betaWald", ties="discrete")
    qout[irow,"NHANEScode"] <- nm
    qout[irow,"Chem"] <- nms[i,"Chem"]
    qout[irow,"subpop"] <- "Total"
    qout[irow,"LOD"] <- nms[i,"LOD"]
    qout[irow,qoutnames[4+1+seq(0,N, by=3)]] <- out$quantiles[1,]
    qout[irow,qoutnames[4+2+seq(0,N, by=3)]] <- out$CIs[1,,1]
    qout[irow,qoutnames[4+3+seq(0,N, by=3)]] <- out$CIs[2,,1]

    ## Sample Size
    qout[irow,"Sample_Size"] <- nrow(design$variables)
    qout[irow,"BelowLOD"] <- sum(!is.na(design$variables[,nmlc]) & !is.na(design$variables$Measure) & design$variables[,nms[i,"LODind"]] == 1)
    ## Fit lognormal
    if (lognormfit) {
      forms <- list(mean=as.formula(paste("I(cbind(Measure,",nmlc,")) ~ 0 + ones",sep="")),
                    lsd=~ 0 + ones)
      ## initial values
      xmean <- svymean(~lMeasure, design)
      xsd <- sqrt(svyvar(~lMeasure, design))
      ## First try gradient free
      out1a <- try(svymle(lnlike, design=design, formulas=forms,
                          start=list(mean=xmean, lsd=log(xsd))),
                   silent=TRUE)
      if (!is(out1a, "try-error")) {
        ## Use that as initial value
        out1b <- try(svymle(lnlike, gradient=gcens, design=design, formulas=forms,
                            start=list(mean=out1a$par[1], lsd=out1a$par[2])),
                     silent=TRUE)
        if (!is(out1b, "try-error") && out1b$ierr > 3) {
          vc <- try(vcov(out1b))
          qout[irow,"loggm"] <- out1b$par[1]
          qout[irow,"lsdlog"] <- out1b$par[2]
          qout[irow,"loggm_se"] <- sqrt(vc[1,1])
          qout[irow,"lsdlog_se"] <- sqrt(vc[2,2])
        } else {
          if (out1a$conv == 0) {
            vc <- vcov(out1a, stderr="model")
            qout[irow,"loggm"] <- out1a$par[1]
            qout[irow,"lsdlog"] <- out1a$par[2]
            qout[irow,"loggm_se"] <- sqrt(vc[1,1])
            qout[irow,"lsdlog_se"] <- sqrt(vc[2,2])
          } else {
            qout[irow,"loggm"] <- NA
            qout[irow,"lsdlog"] <- NA
            qout[irow,"loggm_se"] <- NA
            qout[irow,"lsdlog_se"] <- NA
          }
        }
      } else {
        qout[irow,"loggm"] <- NA
        qout[irow,"lsdlog"] <- NA
        qout[irow,"loggm_se"] <- NA
        qout[irow,"lsdlog_se"] <- NA
      }
    }
    if (doPlot) {
      plot(svycdf(~ lMeasure2, design), main=paste(nm, "Total", sep=":"), do.points=FALSE)
      ## plot requested quantiles and confidence limits
      ## set up plotting coordinates
      xlow <- unlist(qout[irow, qoutnames[4+2+seq(0,N,by=3)]])
      xlow[xlow < qout[irow,"LOD"]] <- NA
      x <- unlist(qout[irow, qoutnames[4+1+seq(0,N,by=3)]])
      x[is.na(xlow)] <- NA
      xhigh <- unlist(qout[irow, qoutnames[4+3+seq(0,N,by=3)]])
      xhigh[is.na(xlow)] <- NA
      points(log(x), Q, pch=20)
      segments(log(xlow),Q, log(xhigh),Q)

      if (lognormfit && !is.na(loggm <- qout[irow, "loggm"])
          && !is.na(lsdlog <- qout[irow, "lsdlog"])) {
        xx <- seq(min(ndta[,"lMeasure"], na.rm=TRUE),
                  max(ndta[,"lMeasure"], na.rm=TRUE), length=200)
        yy <- pnorm(xx, mean=loggm, sd=exp(lsdlog))
        lines(xx,yy,col="green")
      }
    }

    ## -----------------------------------------------------------------------
    ## ------------------------------ Gender ---------------------------------
    ## -----------------------------------------------------------------------

    out2 <- svyby(formula=formla, by=~RIAGENDR, design, FUN=svyquantile,
                  quantiles=Q, ci=TRUE,
                  na.rm=TRUE, interval.type="betaWald", ties="discrete", vartype="ci")
    if (lognormfit) {
      ## Fit lognormal
      ## xmean
      xmean <- svyby(~lMeasure, ~RIAGENDR, design, svymean)$lMeasure
      xsd <- svyby(~lMeasure, ~RIAGENDR, design, svyvar)
      xsd <- sqrt(xsd$lMeasure)

      forms <-
        list(mean=as.formula(paste("I(cbind(Measure,",nmlc,")) ~ 0 + RIAGENDR",sep="")),
             lsd=~ 0 + RIAGENDR)
      out2a <- try(svymle(lnlike, design=design, formulas=forms,
                          start=list(mean=xmean, lsd=log(xsd)),
                          method="BFGS"),
                   silent=TRUE)
      if (!is(out2a, "try-error")) {
        out2b <-  try(svymle(lnlike, gradient=gcens, design=design, formulas=forms,
                             start=list(mean=out2a$par[1:2], lsd=out2a$par[3:4]),
                             method="BFGS"),
                      silent=TRUE)
        if (!is(out2b, "try-error") && out2b$convergence == 0) {
          se <- sqrt(diag(vcov(out2b)))
          parms <- out2b$par
        } else {
          if (out2a$convergence == 0) {
            se <- sqrt(diag(vcov(out2a, stderr="model")))
            parms <- out2a$par
          } else {
            se <- rep(as.numeric(NA), 4)
            parms <- rep(as.numeric(NA),4)
          }
        }
      } else {
        se <- rep(as.numeric(NA), 4)
        parms <- rep(as.numeric(NA),4)
      }
    }
    for (j in 1:2) {
      irow <- irow + 1
      qout[irow,"NHANEScode"] <- nm
      qout[irow,"Chem"] <- nms[i,"Chem"]
      qout[irow,"LOD"] <- nms[i,"LOD"]
      qout[irow,"subpop"] <- levels(out2[,demogendr])[j]
      qout[irow,qoutnames[4+1+seq(0,N, by=3)]] <- out2[j,1+seq_len(length(Q))]
      qout[irow,qoutnames[4+2+seq(0,N, by=3)]] <- out2[j,1+length(Q)+seq_len(length(Q))]
      qout[irow,qoutnames[4+3+seq(0,N, by=3)]] <- out2[j,1+2*length(Q)+seq_len(length(Q))]
      ## Sample Size
      Indices <- which(design$variables[,demogendr] == levels(out2[,demogendr])[j])
      qout[irow,"Sample_Size"] <- length(Indices)
      qout[irow,"BelowLOD"] <- sum(!is.na(design$variables[Indices,nmlc]) & !is.na(design$variables$Measure[Indices]) &
                                     design$variables[Indices,nmlc] == 1)

      if (lognormfit) {
        qout[irow,"loggm"] <- parms[j]
        qout[irow, "lsdlog"] <- parms[2 + j]
        qout[irow,"loggm_se"] <- se[j]
        qout[irow,"lsdlog_se"] <- se[2 + j]
      }
    }
    if (doPlot) {
      gendr <- c("Male","Female")
      for (j in 1:2) {
        jrow <- irow - 2 + j
        plot(svycdf(~lMeasure2, subset(design, RIAGENDR == gendr[j])),
             main=paste(nm,paste(gendr[j],"s",sep=""),sep=":"), do.points=FALSE)
        ## plot requested quantiles and confidence limits
        ## set up plotting coordinates
        xlow <- unlist(qout[jrow, qoutnames[4+2+seq(0,N,by=3)]])
        xlow[xlow < qout[jrow,"LOD"]] <- NA
        x <- unlist(qout[jrow, qoutnames[4+1+seq(0,N,by=3)]])
        x[is.na(xlow)] <- NA
        xhigh <- unlist(qout[jrow, qoutnames[4+3+seq(0,N,by=3)]])
        xhigh[is.na(xlow)] <- NA
        points(log(x), Q, pch=20)
        segments(log(xlow),Q, log(xhigh),Q)
        if (lognormfit && !is.na(loggm <- qout[jrow, "loggm"]) && !is.na(lsdlog <- qout[jrow, "lsdlog"])) {
          xx <- seq(min(ndta[,"lMeasure"], na.rm=TRUE),
                    max(ndta[,"lMeasure"], na.rm=TRUE), length=200)
          yy <- pnorm(xx, mean=loggm, sd=exp(lsdlog))
          lines(xx,yy,col="green")
        }
      }
    }

    ## -----------------------------------------------------------------------
    ## ----------------------------- Age -------------------------------------
    ## -----------------------------------------------------------------------

    out3 <- svyby(formla, by=~AgeGroup, design, FUN=svyquantile, quantiles=Q, ci=TRUE,
                  na.rm=TRUE, interval.type="betaWald", ties="discrete", vartype="ci")
    NAges <- nlevels(factor(out3$AgeGroup))
    if (lognormfit) {
      ## Fit lognormal
      ## xmean
      xmean <- svyby(~lMeasure, ~AgeGroup, design, svymean)$lMeasure
      xsd <- svyby(~lMeasure, ~AgeGroup, design, svyvar)
      xsd <- sqrt(xsd$lMeasure)

      forms <-
        list(mean=as.formula(paste("I(cbind(Measure,",nmlc,")) ~ 0 + AgeGroup",sep="")),
             lsd=~ 0 + AgeGroup)
      out3a <- try(svymle(lnlike, design=design, formulas=forms,
                          start=list(mean=xmean, lsd=log(xsd)),
                          method="BFGS"),
                   silent=TRUE)
      ##browser()
      if (!is(out3a, "try-error")) {
        out3b <-  try(svymle(lnlike, gradient=gcens, design=design, formulas=forms,
                             start=list(mean=out3a$par[1:NAges], lsd=out3a$par[(NAges + 1):(2 * NAges)]),
                             method="BFGS"),
                      silent=TRUE)
        if (!is(out3b, "try-error") && out3b$convergence == 0) {
          se <- sqrt(diag(vcov(out3b)))
          parms <- out3b$par
        } else {
          if (out3a$convergence == 0) {
            se <- sqrt(diag(vcov(out3a, stderr="model")))
            parms <- if (out3a$convergence == 0) out3a$par else rep(as.numeric(NA), 2*NAges)
          } else {
            se <- rep(as.numeric(NA), 2*NAges)
            parms <- rep(as.numeric(NA),2*NAges)
          }
        }
      } else {
        se <- rep(as.numeric(NA), 2*NAges)
        parms <- rep(as.numeric(NA),2*NAges)
      }
    }
    out3$AgeGroup <- factor(out3$AgeGroup)
    for (j in seq_len(nlevels(out3$AgeGroup))) {
      irow <- irow + 1
      qout[irow,"NHANEScode"] <- nm
      qout[irow,"Chem"] <- nms[i,"Chem"]
      qout[irow,"LOD"] <- nms[i,"LOD"]
      qout[irow,"subpop"] <- levels(out3$AgeGroup)[j]
      qout[irow,qoutnames[4+1+seq(0,N, by=3)]] <- out3[j,1+seq_len(length(Q))]
      qout[irow,qoutnames[4+2+seq(0,N, by=3)]] <- out3[j,1+length(Q)+seq_len(length(Q))]
      qout[irow,qoutnames[4+3+seq(0,N, by=3)]] <- out3[j,1+2*length(Q)+seq_len(length(Q))]
      ## Sample Size
      Indices <- which(design$variables$AgeGroup == levels(out3$AgeGroup)[j])
      qout[irow,"Sample_Size"] <- length(Indices)
      qout[irow,"BelowLOD"] <- sum(!is.na(design$variables[Indices,nmlc]) & !is.na(design$variables$Measure[Indices]) &
                                     design$variables[Indices,nmlc] == 1)
      if (lognormfit) {
        qout[irow,"loggm"] <- parms[j]
        qout[irow, "lsdlog"] <- parms[NAges + j]
        qout[irow,"loggm_se"] <- se[j]
        qout[irow,"lsdlog_se"] <- se[NAges + j]
      }
    }
    if (doPlot) {
      gendr <- levels(out3$AgeGroup)
      for (j in seq_len(length(gendr))) {
        jrow <- irow - length(gendr) + j
        plot(svycdf(~lMeasure2, subset(design, AgeGroup == gendr[j])),
             main=paste(nm,paste(gendr[j],sep=""),sep=":"), do.points=FALSE)
        ## plot requested quantiles and confidence limits
        ## set up plotting coordinates
        xlow <- unlist(qout[jrow, qoutnames[4+2+seq(0,N,by=3)]])
        xlow[xlow < qout[jrow,"LOD"]] <- NA
        x <- unlist(qout[jrow, qoutnames[4+1+seq(0,N,by=3)]])
        x[is.na(xlow)] <- NA
        xhigh <- unlist(qout[jrow, qoutnames[4+3+seq(0,N,by=3)]])
        xhigh[is.na(xlow)] <- NA
        points(log(x), Q, pch=20)
        segments(log(xlow),Q, log(xhigh),Q)
        if (lognormfit && !is.na(loggm <- qout[jrow, "loggm"]) && !is.na(lsdlog <- qout[jrow, "lsdlog"])) {
          xx <- seq(min(ndta[,"lMeasure"], na.rm=TRUE),
                    max(ndta[,"lMeasure"], na.rm=TRUE), length=200)
          yy <- pnorm(xx, mean=loggm, sd=exp(lsdlog))
          lines(xx,yy,col="green")
        }
      }
    }

    ## -----------------------------------------------------------------------
    ## ---------------------------- Obesity ----------------------------------
    ## -----------------------------------------------------------------------

    out4 <- svyby(formla, by=~Obesity, design, FUN=svyquantile, quantiles=Q, ci=TRUE,
                  na.rm=TRUE, interval.type="betaWald", ties="discrete", vartype="ci")
    if (lognormfit) {
      ## Fit lognormal
      ## xmean
      xmean <- svyby(~lMeasure, ~Obesity, design, svymean)$lMeasure
      xsd <- svyby(~lMeasure, ~Obesity, design, svyvar)
      xsd <- sqrt(xsd$lMeasure)

      forms <-
        list(mean=as.formula(paste("I(cbind(Measure,",nmlc,")) ~ 0 + Obesity",sep="")),
             lsd=~ 0 + Obesity)
      out4a <- try(svymle(lnlike, design=design, formulas=forms,
                          start=list(mean=xmean, lsd=log(xsd))),
                   silent=TRUE)
      if (!is(out4a, "try-error")) {
        out4b <-  try(svymle(lnlike, gradient=gcens, design=design, formulas=forms,
                             start=list(mean=out4a$par[1:2], lsd=out4a$par[3:4])),
                      silent=TRUE)
        if (!is(out4b, "try-error") && out4b$ierr > 3) {
          se <- sqrt(diag(vcov(out4b)))
          parms <- out4b$par
        } else {
          if (out4a$conv == 0) {
            se <- sqrt(diag(vcov(out4a, stderr="model")))
            parms <- if (out4a$convergence == 0) out4a$par else rep(as.numeric(NA), 4)
          } else {
            se <- rep(as.numeric(NA), 4)
            parms <- rep(as.numeric(NA),4)
          }
        }
      } else {
        se <- rep(as.numeric(NA), 4)
        parms <- rep(as.numeric(NA),4)
      }
    }
    out4$Obesity <- factor(out4$Obesity)
    for (j in seq_len(nlevels(out4$Obesity))) {
      irow <- irow + 1
      qout[irow,"NHANEScode"] <- nm
      qout[irow,"Chem"] <- nms[i,"Chem"]
      qout[irow,"LOD"] <- nms[i,"LOD"]
      qout[irow,"subpop"] <- levels(out4$Obesity)[j]
      qout[irow,qoutnames[4+1+seq(0,N, by=3)]] <- out4[j,1+seq_len(length(Q))]
      qout[irow,qoutnames[4+2+seq(0,N, by=3)]] <- out4[j,1+length(Q)+seq_len(length(Q))]
      qout[irow,qoutnames[4+3+seq(0,N, by=3)]] <- out4[j,1+2*length(Q)+seq_len(length(Q))]
      ## Sample Size
      Indices <- which(design$variables$Obesity == levels(out4$Obesity)[j])
      qout[irow,"Sample_Size"] <- length(Indices)
      qout[irow,"BelowLOD"] <- sum(!is.na(design$variables[Indices,nmlc]) & !is.na(design$variables$Measure[Indices]) &
                                     design$variables[Indices,nmlc] == 1)
      if (lognormfit) {
        qout[irow,"loggm"] <- parms[j]
        qout[irow, "lsdlog"] <- parms[2 + j]
        qout[irow,"loggm_se"] <- se[j]
        qout[irow,"lsdlog_se"] <- se[2 + j]
      }
    }
    if (doPlot) {
      gendr <- levels(out4$Obesity)
      for (j in seq_len(length(gendr))) {
        jrow <- irow - length(gendr) + j
        plot(svycdf(~lMeasure2, subset(design, Obesity == gendr[j])),
             main=paste(nm,paste(gendr[j],sep=""),sep=":"), do.points=FALSE)
        ## plot requested quantiles and confidence limits
        ## set up plotting coordinates
        xlow <- unlist(qout[jrow, qoutnames[4+2+seq(0,N,by=3)]])
        xlow[xlow < qout[jrow,"LOD"]] <- NA
        x <- unlist(qout[jrow, qoutnames[4+1+seq(0,N,by=3)]])
        x[is.na(xlow)] <- NA
        xhigh <- unlist(qout[jrow, qoutnames[4+3+seq(0,N,by=3)]])
        xhigh[is.na(xlow)] <- NA
        points(log(x), Q, pch=20)
        segments(log(xlow),Q, log(xhigh),Q)
        if (lognormfit && !is.na(loggm <- qout[jrow, "loggm"]) && !is.na(lsdlog <- qout[jrow, "lsdlog"])) {
          xx <- seq(min(ndta[,"lMeasure"], na.rm=TRUE),
                    max(ndta[,"lMeasure"], na.rm=TRUE), length=200)
          yy <- pnorm(xx, mean=loggm, sd=exp(lsdlog))
          lines(xx,yy,col="green")
        }
      }
    }

    ## -----------------------------------------------------------------------
    ## ------------------- Women of Child-Bearing Age (16 - 49) --------------
    ## -----------------------------------------------------------------------

    dsgnsub <- subset(design, ReproAgeFemale)
    out5 <- svyquantile(formla, dsgnsub, Q, ci=TRUE, se=TRUE,
                        na.rm=TRUE, interval.type="betaWald", ties="discrete")
    irow <- irow + 1
    qout[irow,"NHANEScode"] <- nm
    qout[irow,"Chem"] <- nms[i,"Chem"]
    qout[irow,"LOD"] <- nms[i,"LOD"]
    qout[irow,"subpop"] <- "ReproAgeFemale"
    qout[irow,qoutnames[4+1+seq(0,N, by=3)]] <- out5$quantiles[1,]
    qout[irow,qoutnames[4+2+seq(0,N, by=3)]] <- out5$CIs[1,,1]
    qout[irow,qoutnames[4+3+seq(0,N, by=3)]] <- out5$CIs[2,,1]
    ## Sample Size
    qout[irow,"Sample_Size"] <- nrow(dsgnsub$variables)
    qout[irow,"BelowLOD"] <- sum(!is.na(dsgnsub$variables[,nmlc]) & !is.na(dsgnsub$variables$Measure) &
                                   dsgnsub$variables[,nmlc] == 1)

    ## Fit lognormal
    if (lognormfit) {
      forms <- list(mean=as.formula(paste("I(cbind(Measure,",nmlc,")) ~ 0 + ones",sep="")),
                    lsd=~ 0 + ones)
      ## initial values
      xmean <- svymean(~lMeasure, dsgnsub)
      xsd <- sqrt(svyvar(~lMeasure, dsgnsub))
      ## First try gradient free
      out1a <- try(svymle(lnlike, design=dsgnsub, formulas=forms,
                          start=list(mean=xmean, lsd=log(xsd))),
                   silent=TRUE)
      if (!is(out1a, "try-error")) {
        ## Use that as initial value
        out1b <- try(svymle(lnlike, gradient=gcens, design=dsgnsub, formulas=forms,
                            start=list(mean=out1a$par[1], lsd=out1a$par[2])),
                     silent=TRUE)
        if (!is(out1b, "try-error") && out1b$ierr > 3) {
          vc <- try(vcov(out1b))
          qout[irow,"loggm"] <- out1b$par[1]
          qout[irow,"lsdlog"] <- out1b$par[2]
          qout[irow,"loggm_se"] <- sqrt(vc[1,1])
          qout[irow,"lsdlog_se"] <- sqrt(vc[2,2])
        } else {
          if (out1a$conv == 0) {
            vc <- vcov(out1a, stderr="model")
            qout[irow,"loggm"] <- out1a$par[1]
            qout[irow,"lsdlog"] <- out1a$par[2]
            qout[irow,"loggm_se"] <- sqrt(vc[1,1])
            qout[irow,"lsdlog_se"] <- sqrt(vc[2,2])
          } else {
            qout[irow,"loggm"] <- NA
            qout[irow,"lsdlog"] <- NA
            qout[irow,"loggm_se"] <- NA
            qout[irow,"lsdlog_se"] <- NA
          }
        }
      } else {
        qout[irow,"loggm"] <- NA
        qout[irow,"lsdlog"] <- NA
        qout[irow,"loggm_se"] <- NA
        qout[irow,"lsdlog_se"] <- NA
      }
    }
    if (doPlot) {
      plot(svycdf(~ lMeasure2, dsgnsub), main=paste(nm, "Total", sep=":"), do.points=FALSE)
      ## plot requested quantiles and confidence limits
      ## set up plotting coordinates
      xlow <- unlist(qout[irow, qoutnames[4+2+seq(0,N,by=3)]])
      xlow[xlow < qout[irow,"LOD"]] <- NA
      x <- unlist(qout[irow, qoutnames[4+1+seq(0,N,by=3)]])
      x[is.na(xlow)] <- NA
      xhigh <- unlist(qout[irow, qoutnames[4+3+seq(0,N,by=3)]])
      xhigh[is.na(xlow)] <- NA
      points(log(x), Q, pch=20)
      segments(log(xlow),Q, log(xhigh),Q)

      if (lognormfit && !is.na(loggm <- qout[irow, "loggm"]) && !is.na(lsdlog <- qout[irow, "lsdlog"])) {
        xx <- seq(min(ndta[,"lMeasure"], na.rm=TRUE),
                  max(ndta[,"lMeasure"], na.rm=TRUE), length=200)
        yy <- pnorm(xx, mean=loggm, sd=exp(lsdlog))
        lines(xx,yy,col="green")
      }
    }

    ## =======================================================================
    ## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ## =======================================================================


    ## Replace all values that are <nms$LOD[i] with NA
    if (LODfilter) {
      ## This is an attempt to match the way the NHANES reports mask estimates.
      ## all estimates that are < LOD are replaced with NA
      ## if the central estimate is < LOD, then the UCL is replaced with NA, too.
      ## So, first, replace everything < LOD with NA
      indx <- qout[(i - 1)*11 + 1:11,5:(ncol(qout) - 4)] < nms$LOD[i]
      qout[(i - 1)*11 + 1:11,5:(ncol(qout) - 4)][indx] <- NA
      ## Now, set ucl to NA if central estimate is NA
      for (j in 1:11) {
        indx <- is.na(qout[(i - 1)*11 + j,Qcenindx])
        if (any(indx)) qout[(i - 1)*11 + j,Qhighindx[indx]] <- NA
      }
    }
  }
  ## add chemical names, and CAS numbers, if we've got them
  if (!is.null(code)) {
    #qout$Chem <- code$table[match(qout$NHANEScode, code$table[,code$table$NHANEScode]),code$chemname]  # new
    qout$Chem <- code$table[match(qout$NHANEScode, code$table[,code$codename]),code$chemname]  # old
    if (!is.null(code$CAS)) {
      #qout$CAS <- code$table[match(qout$NHANEScode, code$table[,code$table$NHANEScode]),code$CAS]  # new
      qout$CAS <- code$table[match(qout$NHANEScode, code$table[,code$codename]),code$CAS]   # old
    }
  }
  ## If we did not fit lognormal, drop those columns
  if (!lognormfit) {
    indx <- match(c("loggm","loggm_se","lsdlog","lsdlog_se"), colnames(qout))
    qout <- qout[,-indx]
  }
  return(qout)
}



#' getDesign
#'
#' Obtain the design object for a given name nm and subpop.  Same as the first part of
#' getNhanesQuantiles.
#'
#' @param demof String containing name of SAS transport file containing
#'           demographic data
#' @param chemdtaf String containing name of SAS transport file containing
#'           the urine concentrations for a group of chemicals.
#'           presumed to have cols SEQN, WTSB2YR, and pairs of column for each chemical: one
#'           for the measurement and one for an indicator for < LOD.  The value in the measurement
#'           when the measurement is < LOD is LOD/sqrt(2), rounded to about 2 digits.
#' @param measurehead String matching the initial part of the measurement variable
#'           names.  The defaults are correct for files for laboratory measurement of
#'           products in urine.
#' @param measuretail Strings matching the final part of the measurement variable
#'           names.  The defaults are correct for files for laboratory measurement of
#'           products in urine.
#' @param lodindhd Regular expression for matching the 'comment' fields (lod indicator), initial part
#' @param lodtail Regular expressions for matching the 'comment' fields (lod indicator), final part
#' @param seq ID, used to merge different files (assumed to be the same variable name across files)
#' @param demoageyr Name of the variable giving age in years in the demographic data file.
#' @param demogendr Name of the variable giving gender in the demographic data file.
#' @param demoeth Name of the variable giving ethnicity in the demographic data file.
#' @param PSU Name of the variable in the demographic data file giving the sampling unit
#' @param STRA Name of the variable in the demographic data file giving the stratum for each observation
#' @param nm Column name, also chemical code, for the chemical to get design for
#' @param creatinine Name of the variable in the file named by chemdataf that contains the creatinine concentration
#'             in each urine sample
#' @param CreatFun Function to compute median or random daily estimates of creatinine excretion.  Used to convert
#'           outputs scaled by creatinine concentration to daily excretion rates. The argument is a data frame with variables:
#'           RIAGENDR, RIDRETH1, BMXWT, and RIDAGEYR
#' @param creatfile Name of file with creatinine measurements
#' @param bodywt Variable giving bodyweight in kg
#' @param bodywtcomment Variable flagging special cases for bodyweight
#' @param MECwt Variable giving the weight to be used when analyzing a full MEC variable
#' @param chem2yrwt Variable in chemdataf giving the sampling weight to be used for those data
#' @param urinefile File containing basic information about urine rate and volume.  ONly relevant to 2009 and
#'            later
#' @param urinerate Variable in urinefile giving the rate of urine production
#' @param bodywtfile File giving the body weight information
#' @param bodymassindex Code for BMI information in the NHANES file.  Default is "BMXBMI".
#' @param codes_table R dataframe of the first sheet of the codes_file in the readNHANES() function
#' @param code Information for adding chemical names and CAS numbers to the output.  Is a list with elements
#'         table: a character matrix or data frame with character (not
#'                factor) elements giving the names NHANES uses for
#'                chemical measures, the chemical names, and the CAS
#'                numbers as strings.
#'         codename: the name in table of the variable that contains the codes NHANES uses
#'         chemame: tha name in table of the variable containing the chemical names
#'         CAS: either NULL or the name in table of the variable containing the CAS numbers
#'
#' @import survey
#' @importFrom stats as.formula dnorm na.omit pnorm
#' @importFrom foreign read.xport
#'
#' @return A design object
#' @export
#'
getDesign <- function(demof="DEMO_F.XPT", chemdtaf, measurehead="URX", measuretail=NULL,
                      lodindhd="URD", lodtail="LC", seq="SEQN",
                      demoageyr="RIDAGEYR",demogendr="RIAGENDR",demoeth="RIDRETH1",
                      PSU="SDMVPSU",
                      STRA="SDMVSTRA", chem2yrwt=NULL,nm=NULL,
                      creatinine="URXUCR", CreatFun=NULL, creatfile=NULL,
                      urinerate=NULL,urinefile=NULL,bodywtfile=NULL,
                      bodywt="BMXWT",bodywtcomment="BMIWT",
                      bodymassindex="BMXBMI",
                      MECwt="WTMEC2YR",
                      codes_table = NULL,
                      code=NULL) {
  ## ----------------------------------------------------------------
  ## Sanity Checks
  oldmethod <- c("1999-2000", "2001-2002", "2003-2004", "2005-2006")
  scale_old <- function(z) {
    if (!is.na(z)) {
      if (z < 75) {
        z <- (1.02*sqrt(z) - 0.36)^2
      } else if (z >= 75 && z < 250) {
        z <- (1.05*sqrt(z) - 0.74)^2
      } else {
        z <- (1.01*sqrt(z) - 0.1)^2
      }
    } else {
      z <- NA
    }
    return(z)
  }

  unitscale <- c("ug/L" = 1.0,"ng/mL" = 1.0,"ng/L" = 0.001, "pg/mL" = 0.001)
  
  ## Sanity Checks
  if (length(demof) > 1){
    demo <- c()
    cdta <- c()
    bwt <- c()
    for (i in 1:length(demof)){
      if (!file.exists(demof[i])) stop(paste(demof[i], "not found"))
      miss_test <- chemdtaf[i]
      if (missing(miss_test) || !file.exists(chemdtaf[i])) {
        stop(paste("File for chemical data", chemdtaf[i], "not found"))
      }
      tmp <- read.xport(demof[i])
      demovars <- c(seq,PSU,STRA,demoageyr,demogendr,demoeth,MECwt)
      if (any(!(XX <- demovars %in% colnames(tmp)))) {
        stop(paste("Variables requested NOT in ",demof[i],": ",
                   paste(demovars[!XX],collapse=", "), sep=""))
      }
      tmp <- tmp[,colnames(tmp) %in% demovars]

      # Deal with creatinine: get file if column not there and scale old data
      tmp2 <- read.xport(chemdtaf[i])
      z <- grep(toupper(gsub(".xpt", "", basename(chemdtaf[i]))), codes_table$NHANESfile, fixed = TRUE)
      y <- unname(z[which(codes_table$NHANEScode[z] == nm)])
      tmp2[,nm] <- tmp2[,nm] * unitscale[codes_table$units[y]]
      chemvars2 <- c(seq, chem2yrwt[i], nm)
      if (creatinine %in% colnames(tmp2)) {
        chemvars2 <- c(chemvars2, creatinine)
        scalebycreatinine <- TRUE
      } else if (!(creatinine %in% colnames(tmp2)) & !is.null(creatinine)) {
        if (is.null(creatfile[i])) stop("No creatinine file provided")
        if (!file.exists(creatfile[i])) stop(paste(createfile[i], "not found"))
        creat <- read.xport(creatfile[i])
        tmp2 <- merge(tmp2, creat[,c(seq, creatinine)], all.x = TRUE)
        chemvars2 <- c(chemvars2, creatinine)
        scalebycreatinine <- TRUE
      } else {
        scalebycreatinine = FALSE
      }

      ind <- unlist(sapply(oldmethod, function(x) grep(x, demof[i])))
      # Apply if this phase is in oldmethod
      if (length(ind) > 0) {
        tmp2[,creatinine] <- sapply(tmp2[,creatinine], scale_old)
      }

      ##  cat("getNhanesQuantiles, chemvars2: \"", chemvars2, "\"\n")
      if (any(!(XX <- chemvars2 %in% colnames(tmp2)))) {
        stop(paste("Variables requested NOT in ",chemdtaf[i],": ",
                   paste(chemvars2[!XX], collapse=", "), sep=""))
      }
      tmp2 <- tmp2[,colnames(tmp2) %in% chemvars2]
      if (i > 1) {
        colnames(tmp2)[colnames(tmp2) == chem2yrwt[i]] <- chem2yrwt[1]
      }

      if (!is.null(urinefile[i])) {
        if (is.null(urinerate[i]))
          stop("File for urinerate specified without variable name for urine rate.")
        if (!file.exists(urinefile[i])) stop(paste(urinefile[i],"not found."))
        urine <- read.xport(urinefile[i])
        if (!urinerate %in% colnames(urine[i]))
          stop(paste(urinrate[i],"not in",urinefile[i]))
        scalebyurinerate <- TRUE
      } else {
        scalebyurinerate <- FALSE
      }
      ## scalebyurinerate dominates scalebycreatinine, so set the latter
      ## to false if we are scaling by urinerate
      if (scalebyurinerate) scalebycreatinine <- FALSE

      if (!is.null(bodywtfile[i])) {
        if (is.null(bodywt))
          stop("File for bodywt specified without variable name for body weight.")
        if (is.null(bodywtcomment))
          stop("Must specify name for bodywtcomment")
        if (is.null(MECwt))
          stop("Must specify name for Exam Weights")
        if (!file.exists(bodywtfile[i])) stop(paste(bodywtfile[i], "not found."))
        tmp3 <- read.xport(bodywtfile[i])
        if (any(!(XX <- c(bodywt, bodywtcomment,bodymassindex) %in% colnames(tmp3))))
          stop(paste(paste(c(bodywt,bodywtcomment,bodymassindex)[XX],collapse=", "),"not in",bodywtfile[i]))
        tmp3 <- tmp3[,colnames(tmp3) %in% c(seq, bodywt, bodywtcomment, bodymassindex)]
        scalebybodywt <- TRUE
      } else {
        scalebybodywt <- FALSE
        bodywt <- NULL
        bodywtcomment <- NULL
        MECwt <- NULL
      }

      demo <- rbind(tmp, demo)
      cdta <- rbind(tmp2, cdta)
      bwt <- rbind(tmp3, bwt)
    }
    chem2yrwt <- chem2yrwt[1]

    # Deal with combined weights
    demo[,MECwt] <- demo[,MECwt]/length(demof)

  } else {

    if (!file.exists(demof)) stop(paste(demof, "not found"))
    if (missing(chemdtaf) || !file.exists(chemdtaf))
      stop(paste("File for chemical data",chemdtaf, "not found"))

    demo <- read.xport(demof)
    cdta <- read.xport(chemdtaf)
    
    # Scale units
    z <- sapply(codes_table$NHANESfile, function(x) grep(tolower(x), chemdtaf))
    z <- which(lengths(z) > 0)
    y <- unname(z[which(codes_table$NHANEScode[z] == nm)])
    cdta[,nm] <- cdta[,nm] * unitscale[codes_table$units[y]]
    
    demovars <- c(seq,PSU,STRA,demoageyr,demogendr,demoeth,MECwt)
    if (any(!(XX <- demovars %in% colnames(demo)))) {
      stop(paste("Variables requested NOT in ",demof,": ",
                 paste(demovars[!XX],collapse=", "), sep=""))
    }

    # For newer NHANES phases, need to get creatinine column from creatinine file.
    # Otherwise just use the creatinine column in the metabolite measurement file.
    chemvars2 <- c(seq, chem2yrwt, nm)
    if (creatinine %in% colnames(cdta)) {
      chemvars2 <- c(chemvars2, creatinine)
      scalebycreatinine <- TRUE
    } else if (!is.null(creatinine)) {
      creat <- read.xport(creatfile)
      cdta <- merge(cdta, creat[,c(seq, creatinine)], all.x = TRUE)
      chemvars2 <- c(chemvars2, creatinine)
      scalebycreatinine <- TRUE
    } else {
      scalebycreatinine <- FALSE
    }

    # NHANES suggests scaling creatinine measurements before 2007 using a pairwise
    # equation to be more comparable.  Do this here.
    ind <- unlist(sapply(oldmethod, function(x) grep(x, demof)))
    # Apply if this phase is in oldmethod
    if (length(ind) > 0) {
      cdta[,creatinine] <- sapply(cdta[,creatinine], scale_old)
    }

    if (any(!(XX <- chemvars2 %in% colnames(cdta))))
      stop(paste("Variables requested NOT in ",chemdtaf,": ",
                 paste(chemvars2[!XX], collapse=", "), sep=""))

    if (!is.null(urinefile)) {
      if (is.null(urinerate))
        stop("File for urinerate specified without variable name for urine rate.")
      if (!file.exists(urinefile)) stop(paste(urinefile,"not found."))
      urine <- read.xport(urinefile)
      if (!urinerate %in% colnames(urine))
        stop(paste(urinrate,"not in",urinefile))
      scalebyurinerate <- TRUE
    } else {
      scalebyurinerate <- FALSE
    }
    ## scalebyurinerate dominates scalebycreatinine, so set the latter
    ## to false if we are scaling by urinerate
    if (scalebyurinerate) scalebycreatinine <- FALSE

    if (!is.null(bodywtfile))
    {
      if (is.null(bodywt))
        stop("File for bodywt specified without variable name for body weight.")
      if (is.null(bodywtcomment))
        stop("Must specify name for bodywtcomment")
      if (is.null(MECwt))
        stop("Must specify name for Exam Weights")
      if (!file.exists(bodywtfile)) stop(paste(bodywtfile, "not found."))
      bwt <- read.xport(bodywtfile)
      if (any(!(XX <- c(bodywt, bodywtcomment,bodymassindex) %in% colnames(bwt))))
        stop(paste(paste(c(bodywt,bodywtcomment,bodymassindex)[XX],collapse=", "),"not in",bodywtfile))
      scalebybodywt <- TRUE
    }
    else
    {
      scalebybodywt <- FALSE
      bodywt <- NULL
      bodywtcomment <- NULL
      MECwt <- NULL
    }
  }

  ## End of sanity checks
  ## ----------------------------------------------------------------
  ## Scaling?
  ## Possibilities:
  ##   - no scaling.  chemvars from chemdataf processed as is: creatinine null and urinefile null.
  ##   - scale by creatinine mg/DL.  Convert creatinine to g/L (divide it by 100), then
  ##     convert chemvars from mass unit/ml -> mass unit/mg creatinine: chemvars/creatinine
  ##     creatinine Not null and urinefile null.  If CreatFun is not null, then the result of the
  ##     previous computation is multiplied by the estimate of daily creatinine excretion from CreatFun
  ##     to give mass units / day.
  ##   - scale by urine flow rate (in ml/min) by multiplying chemvars by urinerate.  This gives
  ##     mass unit / min.  Multiply by 24*60 min/day to get mass unit / day, which is
  ##     what is analyzed
  ##     and reported.  Happens when urinefile is not null, regardless of creatinine.
  ##   - bodywt scaling: in addition to the above, of bodywtfile is not null, the result
  ##     of the previous
  ##     scaling is divided by bodywt to give mass units / mg creatinine / kg bodywt Or
  ##     mass unit / day / kg bodywt or mass units / kg bodywt / day
  ## -----------------------------------------------------------------
  ##   Demographics

  ## Interesting demographic variables are:
  ## demoageyr: age in years
  ## demogendr: gender: male=1, female=2
  ## demoeth: Race/Ethnicity:
  ##        Mexican American=1
  ##        Other Hispanic=2
  ##        Non-Hispanic White=3
  ##        Non-Hispanic Black=4
  ##        Other Race-Including Multi-Racial=5

  ## Select out the variables we'll need going forward
  demo <- demo[,c(seq,PSU,STRA,demoageyr,demogendr,demoeth,MECwt)]
  ## Set up gender, age, and ethnicity as factors using the same levels as the NHANES reports
  demo[,demogendr] <- factor(demo[,demogendr], labels=c("Male","Female"))
  demo$AgeGroup <- cut(demo[,demoageyr], breaks=c(-1,5.5,11.5,19.5,65.5, 100.5),
                       labels=c("0 - 5","6 - 11 years","12 - 19 years", "20 - 65 years", "66 years and older"))
  demo$RaceEthn <- factor(demo[,demoeth],
                          labels=c("Mexican American","Other Hispanic","Non-Hispanic White",
                                   "Non-Hispanic Black","Other"))
  demo$ChildBearingAgeFemale <- factor(demo[,demogendr] == "Female" & (demo[,demoageyr] >= 16 & demo[,demoageyr] <= 49),
                                       labels=c("NotReproAgeFemale","ReproAgeFemale"))
  ## -----------------------------------------------------------------
  ## Chemical Data


  ## We keep from cdta seq, chem2yrwt, all the chemvars, and the
  ## comment variables for the chemvars, if they exist.  Not all
  ## chemvars have comment variables, so they need to be constructed.
  ## create the lodindicator variable names for all chemvars, then,
  ## construct the ones that don't exist.  We can identify missing
  ## values because they are the smallest variable in the variable,
  ## will probably have multiple instances, and will be approximately
  ## the second smallest in the variable / sqrt(2) (rounded to 2
  ## significant digits).
  meascore <- sub(paste("^",measurehead,"(.+)",measuretail,"$",sep=""),"\\1",nm)
  LODnames <- paste(lodindhd,meascore,lodtail,sep="")
  ## nms is a data frame with a record for each chemical.  It keeps
  ## track of the names of the measurement and LODind variables, the
  ## LOD, and space for the chemical name and CAS. LOD here is the
  ## maximum of the individual level LODs.
  nms <- list(Measurement=nm,
              LODind=LODnames,
              LOD=numeric(1),
              Chem = character(1),
              CAS = character(1))

  ## Need to confirm that there is only one value associated with an lod,
  ## then, the next to lowest value in the variable is the lod.
  nmlc <- nms$LODind
  ## Do we already have an LOD indicator?
  if (!nmlc %in% colnames(cdta))
  {
    ## We have to create the LOD indicator.  In this case, there is only 1 LOD
    ## and it will be sqrt(2) times the smallest value
    zlod <- sort(unique(cdta[,nm]))[1] * sqrt(2)
    cdta[,LODnames] <- ifelse(cdta[,nm] >= zlod, 0, 1)
  }

  ## ---------------------------------------------------------------
  ## Merge data sets
  ##
  ## Merge demo and cdta on seq variables, retaining all the variables still in demo,
  ## and in cdta: seq, chem2yrwt, chemvars, all LODind, and creatinine if it exists
  ## To be safe, if we are going to scale by urinerate, then reset creatinine to NULL.
  ## This protects us against unnecessary losses due to NAs in creatinine.
  if (scalebyurinerate) creatinine <- NULL
  alldata <- merge(demo, cdta[,c(seq, chem2yrwt, nm, LODnames, creatinine)],
                   by.x=seq, by.y=seq, all.y=TRUE)
  ## If scalebyurinerate, merge alldata and urine[,c(seq, urinerate)]
  if (scalebyurinerate)
    alldata <- merge(alldata, urine[,c(seq, urinerate)],
                     by.x=seq, by.y=seq, all.x=TRUE)
  ## if scalebybodywt, merge alldata and bwt[,c(seq, bodywt)]
  if (scalebybodywt)
  {
    ## fixup bodywt to include missings due to special codes in
    ## bodywtcomment.  Basically, set bodywt to NA for any non-missing bodywtcomment
    bwt[!is.na(bwt[,bodywtcomment]),bodywt] <- NA
    ## Create obesity factor from bodymassindex
    bwt$Obesity <- cut(bwt[,bodymassindex], breaks=c(-0.5,30,500),
                       labels=c("BMI <= 30", "BMI > 30"))
    alldata <- merge(alldata, bwt[,c(seq, bodywt, bodymassindex, "Obesity")],
                     by.x=seq, by.y=seq, all.x=TRUE)
    ## Impute missing values for bodywt based on age, gender, and ethnicity
    if (any(is.na(alldata[,bodywt])))
    {
      ##writeLines(paste(sum(is.na(alldata[,bodywt])),"missing values in",bodywt))
      ##browser()
      dta <- merge(demo, bwt[,c(seq, bodywt)], by.x=seq, by.y=seq, all.x=TRUE, all.y=TRUE)
      dsg <- na.omit(svydesign(ids=make.formula(PSU), strata=make.formula(STRA),
                               weights=make.formula(MECwt), nest=TRUE, data=dta))
      ## Just fit separate means by age and gender
      selectmales <- dsg$variables[,demogendr] == "Male"
      selectfemales <- dsg$variables[,demogendr] == "Female"
      males <- svyby(make.formula(bodywt), by=make.formula(demoageyr),
                     subset(dsg, selectmales), svymean)
      females <- svyby(make.formula(bodywt), by=make.formula(demoageyr),
                       subset(dsg, selectfemales), svymean)
      imp <- data.matrix(cbind(males[,bodywt],females[,bodywt]))
      ## Now, do the imputation
      ##browser()
      isnabw <- is.na(alldata[,bodywt])
      alldata[isnabw,bodywt] <- imp[data.matrix(alldata[isnabw,c(demoageyr, demogendr)])]
      ##writeLines(">>> Imputed bodyweights")
      ##print(alldata[isnabw, c(demoageyr, demogendr, bodywt)])
    }
  }
  ## Fixup factors so there are no missing levels
  for (colnm in names(alldata))
    if (is.factor(alldata[,colnm]))
      alldata[,colnm] <- factor(alldata[,colnm])
  ## if !is.null(CreatFun), add the variable DailyCreatinine
  if (!is.null(CreatFun))
    alldata$DailyCreatinine <-
    CreatFun(data.frame(RIAGENDR = unname(alldata[,demogendr]),
                        RIDRETH1 = unname(alldata[,"RaceEthn"]),
                        BMXWT = unname(alldata[,bodywt]),
                        RIDAGEYR = unname(alldata[,demoageyr])))

  ## Now, loop through the variables characterized by the rows in nms
  ## Fill a dataframe with the estimates and confidence limits
  ## qouts has columns:

  ## nm is an argument; no longer looping over i.
  ## Set up data frame for measurement variable nm.  We keep this
  ## data frame as small as possible, so we don't lose records because
  ## of NAs in variables we don't use (a problem for nlm, used in svymle).
  keepvars <- c(colnames(demo),chem2yrwt, nm, nmlc, creatinine, urinerate, bodywt,"Obesity")
  if (!is.null(CreatFun)) keepvars <- c(keepvars, "DailyCreatinine")
  ndta <- alldata[,keepvars]
  ndta$Measure <- ndta[,nm]
  ## implement scaling options
  if (scalebycreatinine)
  {
    ndta$Measure <- ndta$Measure * 100 / ndta[,creatinine]
    if (!is.null(CreatFun))
    {
      ndta$Measure <- ndta$Measure * ndta$DailyCreatinine
    }
  }
  if (scalebyurinerate)
    ndta$Measure <- ndta$Measure * ndta[,urinerate] * 24*60 ## gives mass unit / day
  if (scalebybodywt)
    ndta$Measure <- ndta$Measure / ndta[,bodywt]

  ## replace CDC's below LOD value with the LOD
  ## First, save the original value: that gets used for computing quantiles
  ndta$Measure2 <- ndta$Measure
  formla <- make.formula("Measure2")
  isBelowLOD <- !is.na(ndta[,nmlc]) & !is.na(ndta$Measure) & ndta[,nms$LODind] == 1
  ndta$Measure[isBelowLOD] <- ndta$Measure[isBelowLOD]*sqrt(2)

  ## The LOD reported is the unweighted geometric mean of the LODs.
  nms$LOD <-
    if (any(isBelowLOD))
      exp(mean(log(ndta$Measure[isBelowLOD]))) else 0.0

  ## lMeasure is used to get initial distribution parameter estimates, and uses
  ## the log(LOD)/sqrt(2) value for < LOD values
  ndta$lMeasure <- log(ndta$Measure)
  ndta$lMeasure[isBelowLOD] <- ndta$lMeasure[isBelowLOD]/sqrt(2)
  ## lMeasure2 is the variable fit by the distribution function.
  ndta$lMeasure2 <- log(ndta$Measure)
  ## Column of ones to trick the modeling function.
  ndta$ones <- rep(1, nrow(ndta))
  ## Set up indicator for ReproAgeFemale
  ndta$ReproAgeFemale <- ndta$ChildBearingAgeFemale == "ReproAgeFemale"
  #print(paste("Number of NAs in weights:  ", sum(is.na(ndta[, chem2yrwt])), sep = ""))

  na_ind <- !is.na(ndta[chem2yrwt])
  design <-
    na.omit(svydesign(ids=make.formula(PSU), strata=make.formula(STRA),
                      weights=make.formula(chem2yrwt), nest=TRUE, data=ndta[na_ind,]))

  #design <-
  # na.omit(svydesign(id=make.formula(PSU), strata=make.formula(STRA),
  #                           weights=make.formula(chem2yrwt), nest=TRUE, data=ndta))

  return(design)

}



#' lnlike
#'
#' Fit lognormal distribution with left censoring to a dataset
#' Assume we are passing the value in x
#' cens is 0 for non-censored, 1 for left-censored (below LOD, eg)
#' llod is the log-limit of detection
#'
#' @param x log of measurement value ??? what is in column 2?
#' @param mean same as x
#' @param lsd log of standard deviation
#'
#' @importFrom stats dnorm pnorm
#'
#' @return lognormal distribution with left censoring fit to the inputs
#' @export
#'
#'
lnlike <- function(x, mean, lsd)
{
  sd <- exp(lsd) + .Machine$double.eps
  ## Contribution of censored values to log-likelihood:
  ifelse(x[,2] == 1,
         pnorm(log(x[,1]), mean=mean, sd=sd, log.p=TRUE),
         dnorm(log(x[,1]), mean=mean, sd=sd, log=TRUE))
}



#' gcens
#'
#' Derivative of the loglike function.  Obtain the gradient
#' to help in calculating variance estimates.
#'
#' @param x Log (ln) of measurement value
#' @param mean Unnamed x
#' @param lsd Log (ln) of standard deviation
#'
#' @importFrom stats dnorm pnorm
#'
#' @return ???
#' @export
#'
#'
gcens<- function(x,mean,lsd)
{
  sd <- exp(lsd) + .Machine$double.eps
  dz<- exp(dnorm(log(x[,1]),mean,sd, log=TRUE) - pnorm(log(x[,1]), mean, sd, log.p=TRUE))

  dm<-ifelse(x[,2]==0,
             (log(x[,1]) - mean)/(sd^2),
             dz*(-1/sd))
  ds<-ifelse(x[,2]==0,
             (log(x[,1])-mean)^2 / (sd^3) - 1/sd,
             ds<- dz*(-(log(x[,1])-mean))*exp(lsd)/(sd^2))
  cbind(dm,ds)
}



#' geometric.mean.to.IUR
#'
#' Obtains IUR production volume code from geometric mean production
#' volume
#'
#' @param x Geometric mean production volumne
#'
#' @return IUR production volume
#' @export
#'
geometric.mean.to.IUR <-
  function (x)
  {
    if (abs(x - (exp((log(1 * 1e+06/2.204/365.25) + log(10 * 1e+06/2.204/365.25))/2))) < 1e-14)
    {
      string <- "1 to < 10 million lbs"
      lower <- 1 * 1e+06/2.204/365.25
      upper <- 10 * 1e+06/2.204/365.25
    }
    else
      if (abs(x - (exp((log(0.5 * 1e+06/2.204/365.25) + log(1 *
                                                            1e+06/2.204/365.25))/2))) < 1e-14)
      {
        string <- "500,000 to < 1 million lbs"
        lower <- 0.5 * 1e+06/2.204/365.25
        upper <- 1 * 1e+06/2.204/365.25
      }
    else
      if (abs(x - exp((log(25000/2.204/365.25) + log(0.5 * 1e+06/2.204/365.25))/2)) < 1e-14)
      {
        string <- "25,000 to < 500,000 lbs"
        lower <- 25000/2.204/365.25
        upper <- 0.5 * 1e+06/2.204/365.25
      }
    else
      if (abs(x - exp((log(6250/2.204/365.25) + log(25000/2.204/365.25))/2)) < 1e-14)
      {
        string <- "< 25,000 lbs"
        lower <- 6250/2.204/365.25
        upper <- 25000/2.204/365.25
      }
    else
      if (abs(x - exp((log(10 * 1e+06/2.204/365.25) + log(50 *
                                                          1e+06/2.204/365.25))/2)) < 1e-14)
      {
        string <- "10 to < 50 million lbs"
        lower <- 10 * 1e+06/2.204/365.25
        upper <- 50 * 1e+06/2.204/365.25
      }
    else
      if (abs(x - exp((log(50 * 1e+06/2.204/365.25) + log(100 *
                                                          1e+06/2.204/365.25))/2)) < 1e-14)
      {
        string <- "50 to < 100 million lbs"
        lower <- 50 * 1e+06/2.204/365.25
        upper <- 100 * 1e+06/2.204/365.25
      }
    else
      if (abs(x - exp((log(100 * 1e+06/2.204/365.25) +
                       log(500 * 1e+06/2.204/365.25))/2)) < 1e-14)
      {
        string <- "100 to < 500 million lbs"
        lower <- 100 * 1e+06/2.204/365.25
        upper <- 500 * 1e+06/2.204/365.25
      }
    else
      if (abs(x - exp((log(500 * 1e+06/2.204/365.25) +
                       log(1000 * 1e+06/2.204/365.25))/2)) < 1e-14)
      {
        string <- "500 million to < 1 billion lbs"
        lower <- 500 * 1e+06/2.204/365.25
        upper <- 1000 * 1e+06/2.204/365.25
      }
    else
      if (abs(x - 5000 * 1e+06/2.204/365.25) < 1e-14)
      {
        string <- "1 billion lbs and greater"
        lower <- 1000 * 1e+06/2.204/365.25
        upper <- Inf
      }
    else
    {
      string <- x
      lower <- NA
      upper <- NA
    }
    out <- data.frame(string)
    out <- cbind(out, data.frame(lower))
    out <- cbind(out, data.frame(upper))
    colnames(out) <- c("IUR Code", "Lower", "Upper")
    return(out)
  }



#' convert.urine.to.exposure
#'
#' Converts urine concentration in ug/g creatinine/kg bodyweight
#' to mg/kg bodyweight/day
#'
#' @param urine.conc Urine concentration in ug/g creatinine/kg bodyweight
#' @param urine.vol.per.day Urine volume per day in liters.  Default is 1.4 L.
#' @param average.creatinine Average creatinine concentration in mg/DL.
#' Default is 122.6 mg/DL.
#'
#' @return Exposure dose in mg/kg bodyweight/day
#' @export
#'
convert.urine.to.exposure <-
  function(urine.conc,urine.vol.per.day=1.4,average.creatinine=122.6)
  {
    ##      mg/g creatinine      g creatinine / L urine     L urine / day
    return((urine.conc/1000)*(average.creatinine*10/1000)*urine.vol.per.day)
  }



#' creatinine
#'
#' Estimates daily creatinine excretion in mg creatinine/day
#'
#' @param newdata A data.frame with individual-specific data including the
#' following data columns: RIAGENDR, RIDRETH1, BMXWT, RIDAGEYR
#'
#' @return Estimated daily creatinine excretion in mg creatinine/day
#'
#' @importFrom stats model.matrix predict
#'
#' @export
#'
#'
creatinine <- function(newdata){
  X <- cbind(model.matrix(~ 0 + RIAGENDR + RIDRETH1, data=newdata),
             predict(Wtns, newdata$BMXWT),
             predict(Agens, newdata$RIDAGEYR))
  return(10^(X %*% modparms[-length(modparms)]))
}



#' rcreatinine
#'
#' Similar to creatinine() but adds random noise to the result
#'
#' @param newdata A data.frame with individual-specific data including the
#' following data columns: RIAGENDR, RIDRETH1, BMXWT, RIDAGEYR
#'
#' @return Estimated daily creatinine excretion in mg creatinine/day
#'
#' @importFrom stats model.matrix predict rt
#'
#' @export
#'
#'
rcreatinine <- function(newdata){
  X <- cbind(model.matrix(~ 0 + RIAGENDR + RIDRETH1, data=newdata),
             predict(Wtns, newdata$BMXWT),
             predict(Agens, newdata$RIDAGEYR))
  Y <- X %*% modparms[-length(modparms)]
  z <- rt(nrow(newdata), df=3.5)
  return(10^(Y + z * exp(modparms[length(modparms)])))
}


