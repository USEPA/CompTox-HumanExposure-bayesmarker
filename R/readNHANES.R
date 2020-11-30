
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
#' @param data_path String providing the path to the raw data.  Default is "rawData",
#' the direct output from running get_NHANES_data would be "rawData", resulting
#' in e.g. ./rawData/1999-2000
#'
#' @importFrom gdata read.xls
#' @importFrom parallel mclapply
#' @importFrom utils write.csv
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
readNHANES <- function(codes_file, data_path = "rawData") {


  datapaths <- structure(file.path(".",data_path,
                                   c("1999-2000","2001-2002","2003-2004","2007-2008","2009-2010","2011-2012","2013-2014","2015-2016")),
                         names= c("99-00","01-02","03-04","07-08","09-10","11-12","13-14","15-16"))

  ## Curated EXCEL spreadsheet giving information about the variables and files used in
  ## this analysis
  NHANEScodes <- codes_file

  ##########
  # Assumed that all files (demo, bmx, other) are lowercase. Search "tolower" when needing a more elegant solution.
  ##########


  ## -------------------------------------------------------------------
  ##             Construct object 'measured', which contains all the
  ##             information we will use about the products measured
  ##             in urine from NHANES.

  ## convtbl was constructed manually, starting with an earlier list of NHANES urine
  ## products. Age cutoffs were based on the highest age group reported in the 4th report.

  convtbl <- read.xls(NHANEScodes, as.is = TRUE)
  wtvars <- read.xls(NHANEScodes, sheet = 2, as.is = TRUE)

  #### Data checks
  print("Checking input codes file")
  if (any(convtbl$CAS == "")){
    convtbl <- convtbl[!(convtbl$CAS == ""),]
    print(paste("Warning:  removed ", as.character(sum(convtbl$CAS == "")), "rows due to missing CAS"), sep = "")
  }
  if (any(duplicated(convtbl$CAS))) {
    print("Error:  duplicate chemical identifiers, each row must have a unique identifier.")
    stop()
  }


  ## First, run through all the input files and estimate quantiles.
  demofiles <- wtvars$demofile
  names(demofiles) <- wtvars$sample

  tmp <- unique(convtbl[,c("recent_sample","NHANESfile")])
  datafiles <- file.path(datapaths[as.character(tmp$recent_sample)], tolower(tmp$NHANESfile))
  names(datafiles) <- tmp$NHANESfile
  demofiles <- file.path(datapaths[as.character(tmp$recent_sample)], demofiles[as.character(tmp$recent_sample)])
  #cat("demofiles: ", demofiles, "\n")
  names(demofiles) <- tmp$NHANESfile ## associate a demo file with each data file
  bwtfiles <- file.path(datapaths[as.character(tmp$recent_sample)],
                        tolower(wtvars$BWfile[match(names(demofiles), wtvars$file)]))
  names(bwtfiles) <- tmp$NHANESfile
  creatfiles <- file.path(datapaths[as.character(tmp$recent_sample)],
                          tolower(wtvars$creatfile[match(names(demofiles), wtvars$file)]))
  names(creatfiles) <- tmp$NHANESfile

  ## use getNHANESquantiles to estimate the 50th percentiles
  ## getNhanesQuantiles assumes a "comment" variable for each input
  ## variable, to contain the lod indicator.  Not all files have these.
  ## Create them on an ad-hoc basis.

  scaledata <- vector("list", length(datafiles))
  chemvars <- split(as.character(convtbl$NHANEScode), f=convtbl$NHANESfile)
  chemwt <- wtvars$wtvariable
  names(chemwt) <- wtvars$file


  scaledata <-
    mclapply(1:length(datafiles),
             FUN=function(i) {
               getNhanesQuantiles(demof=demofiles[i], chemdtaf=datafiles[i],lognormfit=TRUE,
                                  chem2yrwt=chemwt[names(datafiles)[i]],
                                  chemvars=chemvars[[names(datafiles)[i]]],
                                  CreatFun=creatinine,
                                  bodywtfile=bwtfiles[names(datafiles)[i]],
                                  creatfile=creatfiles[names(datafiles)[i]],
                                  Q=c(50),
                                  code=list(table=convtbl[,c("Name","CAS","NHANEScode")],
                                            codename="NHANEScode",CAS="CAS",chemname="Name"),
                                  LODfilter=FALSE,
                                  MaximumAge=150)
             },
             mc.preschedule=FALSE, mc.cores=7)

  scalequants <- do.call("rbind", scaledata)


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
  print(paste("Warning:  ", dim(z)[1], " have measurement data but are missing a parent.
              Chemical identifiers for these metabolites were saved in the file ChemsWithoutParents.csv", sep = ""))
  write.csv(z, file = paste("ChemsWithoutParent_", format(Sys.time(), "%Y-%m-%d"), ".csv", sep = ""), row.names = FALSE)
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

  ## -------------------------------------------------------------------------------------------
  ## Outputs
  print("Saving returned outputs as MCMCdata.RData")
  save(Measured, pred.data, Uses, mapping, file = "MCMCdata.RData")
  return(Measured, pred.data, Uses, mapping)

}
