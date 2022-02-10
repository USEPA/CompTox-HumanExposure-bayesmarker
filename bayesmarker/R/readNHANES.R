
#' readNHANES
#'
#' Read the NHANES raw data and calls getNHANESquantile to obtain geometric mean
#' of concentration for each metabolite.  Matches these to those in the metabolite
#' map and some other objects needed for the next step.
#'
#' @param codes_file Manually created xls or xlsx file containing 3 sheets:  1. NHANES
#' chemicals to include (with identifier, code, file, demographic, and units),
#' 2. Associated weights, filenames, and column names associated with each phase used,
#' and 3. Parent-metabolite map containing chemical identifiers and molecular weights.
#' @param data_path String providing the path to the raw data.  Is "rawData",
#' the direct output from running get_NHANES_data would be "rawData", if left NULL, resulting
#' in e.g. ./rawData/1999-2000  If the directory is a subdirectory of your working one,
#' include "./" (e.g. "./test" will look in "./test/rawData")
#' @param cohort NHANES cycle/phase/cohort to analyze within the codes file. Default is to use the most recent
#' cycle. Other options include "oldest", "all"
#' @param save_directory String providing the directory in which to save the NHANES data files.  If left as the default,
#'                       NULL, it will save to ./rawData.  Otherwise, it will save to save_directory/rawData.
#'
#' @importFrom gdata read.xls
#' @importFrom parallel mclapply
#' @importFrom utils write.csv
#' @importFrom prodlim row.match
#'
#' @return Measured: a data frame containing a row for each metabolite with columns
#'                   providing the chemical identifiers, log geometric means, subpopulation,
#'                   and other data carried on from the codes file.
#'        pred.data: a data frame containing a row for each parent chemical, identifiers,
#'                   MW, and a near field exposure indicator (obtained from a use database).
#'        Uses:      a data frame with indicator columns for various chemical uses
#'        mapping:   the mapping information in the 3rd sheet of the codes file, now matched
#'                   on the chemicals in Measured
#'
#' @export
#'
#' @examples # readNHANES("NHANEScodes.xlsx", "rawData")
#'
readNHANES <- function(codes_file, data_path = NULL, cohort = "newest", save_directory = ".") {

  ## Curated EXCEL spreadsheet giving information about the variables and files used in
  ## this analysis
  NHANEScodes <- codes_file

  ## -------------------------------------------------------------------
  ##             Construct object 'measured', which contains all the
  ##             information we will use about the products measured
  ##             in urine from NHANES.

  ## convtbl was constructed manually, starting with an earlier list of NHANES urine
  ## products. Age cutoffs were based on the highest age group reported in the 4th report.
  convtbl <- read.xls(NHANEScodes, as.is = TRUE)
  wtvars <- read.xls(NHANEScodes, sheet = 2, as.is = TRUE)

  cycles <- c("99-00", "01-02", "03-04", "05-06", "07-08", "09-10", "11-12", "13-14", "15-16")
  if (length(cohort) == 1){
    if (cohort == "newest"){
      chems <- unique(convtbl$NHANEScode)
      ind <- sapply(chems, function(x) max(which(cycles %in% convtbl$recent_sample[convtbl$NHANEScode == x])))
      ind2 <- row.match(data.frame(chems, cycles[ind]), convtbl[,c("NHANEScode", "recent_sample")])
      convtbl <- convtbl[ind2,]
      wtvars <- wtvars[wtvars$file %in% convtbl$NHANESfile,]
    }
    if (cohort %in% cycles){
      convtbl <- convtbl[convtbl$recent_sample %in% cohort,]
      wtvars <- wtvars[wtvars$sample %in% cohort,]
    }
  } else if (length(cohort) > 1){
      convtbl <- convtbl[convtbl$recent_sample %in% cohort,]
      wtvars <- wtvars[wtvars$sample %in% cohort,]
  }


  # Which phases to combine
  byChem <- function(chem, ref){
    ind <- ref$CAS == chem
    res <- ref$recent_sample[ind]
  }
  phaseTbl <- lapply(unique(convtbl$CAS), byChem, ref = convtbl[,c("CAS", "recent_sample")])
  names(phaseTbl) <- unique(convtbl$CAS)


  #### Data checks
  print("Checking input codes file")
  if (any(convtbl$CAS == "")){
    missing <- sum(convtbl$CAS == "")
    convtbl <- convtbl[!(convtbl$CAS == ""),]
    phaseTbl[[which(names(phaseTbl) == "")]] <- NULL
    print(paste("Warning:  removed ", missing, " rows due to missing CAS"), sep = "")
  }
  tmp2 <- do.call(rbind, lapply(phaseTbl, function(x) c(length(unique(x)), length(x))))
  if (any(tmp2[,1] != tmp2[,2])) {
    issue <- which(tmp2[,1] != tmp2[,2])
    print("Duplicate chemical identifier - cohort pair(s) detected. Each chemical identifier must have a
          unique associated cohort.  Problem chemicals:")
    print(names(tmp2)[issue])
    stop()
  }
  if (any(tmp2[,2] > 1)){
    print("Chemical(s) spanning multiple cohorts. Data will be combined. See input argument 'group' for details")
    group <- TRUE
  } else {
    group <- FALSE
  }

  ## datapaths gives the paths for the NHANES individual .xpt data files.  This
  ## construct gives a vector of paths, one for each sample set.  We then
  ## name the files using the shorthand name for each sample set
  phases <- unique(convtbl$recent_sample)
  ind <- sort(substring(phases, 4, 5), index.return = TRUE)
  long <- sapply(phases[ind$ix], function(x) ifelse(as.numeric(substring(x, 1, 2)) < 50,
                                                    paste("20", substring(x, 1, 3), "20", substring(x, 4, 5), sep = ""),
                                                    paste("19", substring(x, 1, 3), "20", substring(x, 4, 5), sep = "")))
  if (is.null(data_path)){
    datapaths <- structure(file.path(".","rawData", long), names = names(long))
  } else {
    datapaths <- structure(file.path(data_path,"rawData", long), names = names(long))
  }

  ## First, run through all the input files and estimate quantiles.
  demofiles <- wtvars$demofile
  names(demofiles) <- wtvars$sample

  if (group){
    tmp2 <- list()
    for (i in 1:length(phaseTbl)){
      tmp3 <- c()
      for (j in 1:length(phaseTbl[[i]])){
        tmp4 <- c(phaseTbl[[i]][j],
                  paste(convtbl$NHANESfile[!is.na(row.match(convtbl[,c("CAS", "recent_sample")],
                                                      data.frame(names(phaseTbl)[i],
                                                                 phaseTbl[[i]][j])))], ".XPT", sep = ""))
        tmp3 <- rbind(tmp4, tmp3)
      }
      tmp2[[i]] <- tmp3
      colnames(tmp2[[i]]) <- c("recent_sample", "NHANESfile")
    }
    names(tmp2) <- names(phaseTbl)

    datafiles <- lapply(tmp2, function(x) file.path(datapaths[as.character(x[,1])],
                                                    tolower(x[,2])))
    demofiles <- lapply(tmp2, function(x) file.path(datapaths[as.character(x[,1])],
                                                    demofiles[as.character(x[,1])]))
    bwtfiles <- lapply(tmp2, function(x) file.path(datapaths[as.character(x[,1])],
                                                   tolower(wtvars$BWfile[match(x[,2], wtvars$file)])))
    creatfiles <- lapply(tmp2, function(x) file.path(datapaths[as.character(x[,1])],
                                                     tolower(wtvars$creatfile[match(x[,2], wtvars$file)])))

    scaledata <- vector("list", length(phaseTbl))
    # Need chemvars to be a single chemical and then pass the files as vectors
    chemvars <- convtbl$NHANEScode[match(names(tmp2), convtbl$CAS)]
    chemwt <- lapply(tmp2, function(x) wtvars$wtvariable[match(x[,2], wtvars$file)])

    print("Starting geometric mean estimations")
    scaledata <-
      mclapply(1:length(tmp2),
               FUN=function(i) {
                 getNhanesQuantiles(demof=demofiles[[i]], chemdtaf=datafiles[[i]], lognormfit=TRUE,
                                    chem2yrwt=chemwt[[i]],
                                    chemvars=chemvars[i],
                                    CreatFun=creatinine,
                                    bodywtfile=bwtfiles[[i]],
                                    creatfile=creatfiles[[i]],
                                    Q=c(50),
                                    codes_table = convtbl,
                                    code=list(table=convtbl[,c("Name","CAS","NHANEScode")],
                                              codename="NHANEScode",CAS="CAS",chemname="Name"),
                                    LODfilter=FALSE,
                                    MaximumAge=150)
               },
               mc.preschedule=FALSE, mc.cores=7)

  } else {

    # Some variables are missing a weight
    #if (any(wtvars$wtvariable == "")){
    #  print("Warning:  at least one data file is missing the 2 year weight. This file and its chemicals will be discarded.")
    #  convtbl <- convtbl[!paste(convtbl$NHANESfile, ".XPT", sep = "") == wtvars$file[wtvars$wtvariable == ""],]
    #  wtvars <- wtvars[!wtvars$wtvariable == "",]
    #}
    tmp <- unique(convtbl[,c("recent_sample","NHANESfile")])
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

    scaledata <- vector("list", length(datafiles))
    chemvars <- split(as.character(convtbl$NHANEScode), f=convtbl$NHANESfile)
    chemwt <- wtvars$wtvariable
    names(chemwt) <- wtvars$file


    print("Starting geometric mean estimations")
    scaledata <-
      mclapply(1:length(datafiles),
               FUN=function(i) {
                 getNhanesQuantiles(demof=demofiles[i], chemdtaf=datafiles[i],lognormfit=TRUE,
                                    chem2yrwt=chemwt[paste(names(datafiles)[i], ".XPT", sep = "")],
                                    chemvars=chemvars[[names(datafiles)[i]]],
                                    CreatFun=creatinine,
                                    bodywtfile=bwtfiles[names(datafiles)[i]],
                                    creatfile=creatfiles[names(datafiles)[i]],
                                    Q=c(50),
                                    codes_talbe = convtbl,
                                    code=list(table=convtbl[,c("Name","CAS","NHANEScode")],
                                              codename="NHANEScode",CAS="CAS",chemname="Name"),
                                    LODfilter=FALSE,
                                    MaximumAge=150)
               },
               mc.preschedule=FALSE, mc.cores = min(c(length(datafiles), 7)))
#    for (i in 1:length(datafiles)){
#      print(datafiles[i])
#      scaledata[[i]] <- getNhanesQuantiles(demof=demofiles[i], chemdtaf=datafiles[i],lognormfit=TRUE,
#                                           chem2yrwt=chemwt[paste(names(datafiles)[i], ".XPT", sep = "")],
#                                           chemvars=chemvars[[names(datafiles)[i]]],
#                                           CreatFun=creatinine,
#                                           bodywtfile=bwtfiles[names(datafiles)[i]],
#                                           creatfile=creatfiles[names(datafiles)[i]],
#                                           Q=c(50),
#                                           code=list(table=convtbl[,c("Name","CAS","NHANEScode")],
#                                                                   codename="NHANEScode",CAS="CAS",chemname="Name"),
#                                           LODfilter=FALSE,
#                                           MaximumAge=150)
#   }

  }

  scalequants <- do.call("rbind", scaledata)

  print("Complete.  Constructing Measured.")

  #######
  ## Build an object for the urine products with
  ## name, cas, sample, and quantiles, which should have names
  ## Value50 LCL50 UCL50
  NHANES.parent.metabolite.mapping <- read.xls(codes_file, sheet=3, as.is = TRUE)

  # Processing: Data check
  # Remove all rows without a MW
  if (any(is.na(NHANES.parent.metabolite.mapping$MW.1)) | any(is.na(NHANES.parent.metabolite.mapping$MW.1))) {
    print("Warning:  At least one parent or metabolite missing MW.  Removing these links.")
    NHANES.parent.metabolite.mapping <- NHANES.parent.metabolite.mapping[!is.na(NHANES.parent.metabolite.mapping$MW.1),]
    NHANES.parent.metabolite.mapping <- NHANES.parent.metabolite.mapping[!is.na(NHANES.parent.metabolite.mapping$MW),]
  }
  # Some CAS for metabolites have an extra space.  Remove these spaces by default.
  NHANES.parent.metabolite.mapping$CAS.1 <- gsub(" ", "", NHANES.parent.metabolite.mapping$CAS.1)
  NHANES.parent.metabolite.mapping$CAS <- gsub(" ", "", NHANES.parent.metabolite.mapping$CAS)

  chem.dfm <- data.frame(CAS=NHANES.parent.metabolite.mapping$CAS.1,
                         MW=NHANES.parent.metabolite.mapping$MW.1)
  chem.dfp <- data.frame(CAS=NHANES.parent.metabolite.mapping$CAS,
                         MW=NHANES.parent.metabolite.mapping$MW)
  chem.MW <-  unique(rbind(chem.dfm, chem.dfp))
  ################################################################################
  # Have DTXSIDs but only for some chemicals. Way to do by DTXSID then CAS?
  ################################################################################
  convtbl <- merge(convtbl, chem.MW, by.x="CAS", by.y="CAS")

  Measured <- merge(scalequants[,c("Chem","CAS","loggm","loggm_se","lsdlog","lsdlog_se",
                                   "Sample_Size","BelowLOD","LOD","subpop")],
                    convtbl[,c("CAS","recent_sample","MW","NHANESfile","NHANEScode","units")],
                    by.x="CAS",by.y="CAS",all=TRUE)

  # Remove rows for chemicals that didn't have a parent (NA for mult columns)
  z <- unique(Measured[is.na(Measured$recent_sample),c(1,2)])
  print(paste("Warning:  ", dim(z)[1], " chemicals have measurement data but are missing a parent.
              Chemical identifiers for these metabolites were saved in the file ChemsWithoutParents",
              format(Sys.time(), "%Y-%m-%d"), ".csv", sep = ""))
  if (!dir.exists(save_directory)){
    dir.create(save_directory)
  }
  write.csv(z, file = file.path(save_directory, paste("ChemsWithoutParent_", format(Sys.time(), "%Y-%m-%d"),
                                                      ".csv", sep = "")), row.names = FALSE)
  Measured <- Measured[!is.na(Measured$recent_sample),]
  #######
  nms <- names(Measured)
  nms <- sub("Chem","Name",nms)
  names(Measured) <- nms
  ## ng/L == pg/ml: z pg/ml = 0.001 ng/ml
  lunitscale <- log(c("ug/L" = 1.0,"ng/mL" = 1.0,"ng/L" = 0.001, "pg/mL" = 0.001))

  ## Need to change the current units.  The true units are in
  ## mass units expressed in Measured$units per kg bodyweight per day
  ## We got here by getting mass units per kg bodyweight per mg creatinine,
  ## times mg creatinine / day.

  Measured$loggm <- Measured$loggm + lunitscale[Measured$units]
  Measured$LOD <- Measured$LOD * exp(lunitscale[Measured$units])

  ## exp(loggm) is now ng chemical / kg bodyweight / day
  ## Convert to nmoles / kg bodyweight / day:
  ## x ng /kg /day / (MW g / mole) = (x * 1e-9 g / kg / day) / (MW g / mole) =
  ## (x * 1e-9 mole / kg / day) = x nmole / kg / day
  Measured$MW <- as.numeric(Measured$MW)
  Measured$loggm <- Measured$loggm - log(Measured[,"MW"])
  Measured$LOD <- Measured$LOD / Measured$MW


  ## ------------------------------------------------------------------------------
  ##              Mapping Parent to metabolite
  mapping <- NHANES.parent.metabolite.mapping

  ## Drop entries for metabolites we won't be following: Drop CAS.1 that are NOT in Measured$CAS
  indx <- which(!(mapping$CAS.1 %in% unique(Measured$CAS)))
  mapping <- mapping[-indx,]

  ## Sort this on CAS and Branch
  indx <- order(mapping$CAS, mapping$Branch)
  mapping <- mapping[indx,]

  ## -----------------------------------------------------------------------------
  ##                           pred.data
  pred.data <- unique(mapping[,c("Name","CAS","MW")])
  rownames(pred.data) <- as.character(pred.data$CAS)

  ## -----------------------------------------------------------------------------
  ##                          Xref
  ##
  ## Build a matrix with rows for parents and columns for metabolites
  ## that is, rows for the records in pred.data, and columns for the
  ## records in Measured.  use CAS numbers to do the linking.
  ## The entries will be NA for values that will be estimated,
  ## 1 for entries representing 1:1 molar relationships between parent and metabolite,
  ## and 0 everywhere else (the majority of entries).
  ## Set all the entries to 0, initially.
  ## Then, set the cells corresponding to entries in NHANES.parent.metabolite.mapping
  ## to NA.  Some of these will be reset to 1,
  ## and the directly observed parents.  These are the CAS numbers in pred.data that
  ## are also in NHANES.
  Xref <- matrix(0, nrow=nrow(pred.data), ncol=length(unique(Measured$CAS)),
                 dimnames=list(pred.data$CAS, unique(Measured$CAS)))
  indx <- cbind(match(mapping$CAS,rownames(Xref)),
                match(mapping$CAS.1,colnames(Xref)))

  Xref[indx] <- NA

  ## Make sure that the order of rows in Measured and
  ## pred.data match the order in Xref.
  pred.data <- pred.data[rownames(Xref),]


  ## ----------------------------------------------------------------------------------
  ## Add Use categories to pred.data
  ## Near field: Fragrance, food additive, consumer use, personal care product.
  indx <- match(pred.data$CAS, UseTable$CASRN)
  Uses <- UseTable[indx,]
  Uses$NearField <- as.numeric(with(Uses, FRAGRANCE | FOOD.ADDITIVE | CONSUMER.USE | PERSONAL.CARE.PRODUCT | PHARMACEUTICAL))
  pred.data$NearField <- Uses$NearField


  ## ----------------------------------------------------------------------------------
  ## Each chemical will have a row for each phase it was measured in.  So if data was merged,
  ## combine file and sample rows into one string.  Need this information moving forward.
  if (group){
    metabs <- unique(Measured[,c("CAS", "subpop")])
    newPhase <- c()
    newFile <- c()
    for (i in 1:dim(metabs)[1]){
      ind <- which(!is.na(row.match(Measured[,c("CAS", "subpop")], metabs[i,])))
      newPhase[i] <- paste(Measured$recent_sample[ind], collapse = ",")
      newFile[i] <- paste(Measured$NHANESfile[ind], collapse = ",")
    }
    Measured <- Measured[!duplicated(Measured[,c("CAS", "subpop")]),]
    Measured$recent_sample <- newPhase
    Measured$NHANESfile <- newFile
  }

  ## -------------------------------------------------------------------------------------------
  ## Outputs
  print(paste("Saving returned outputs as MCMCdata_", format(Sys.time(), "%Y-%m-%d"), ".RData", sep = ""))
  result <- list(Measured = Measured,
                 pred.data = pred.data,
                 Uses = Uses,
                 mapping = mapping)
  save(result, file = file.path(save_directory, paste("MCMCdata_", format(Sys.time(), "%Y-%m-%d"),
                                                      ".RData", sep = "")))

  return(result)

}
