#' get_NHANES_data
#'
#' Downloads NHANES data files needed for exposure inference based on a manually curated file describing which data to use.
#' Note that accuracy is very important for the input codes file. After "Starting Downloads" is printed to the console by
#' this function, no errors or warnings should be seen. Common issues when downloading the files are having incorrect
#' spelling of the file names (e.g. dem.xtpt instead of demo.xpt) and having the wrong phase and extension pair (e.g. 2015
#' - 2016 files end in "_I.XPT"; will fail if this is "_F.XPT").
#'
#' @param codes_file Manually created xls or xlsx file containing 3 sheets:  1. NHANES chemicals to include (with identifier,
#' code, file, demographic, and units), 2. Associated weights, filenames, and column names associated with each phase used,
#' and 3. Parent-metabolite map containing chemical identifiers and molecular weights.
#' @param cohort String or character vector indicating which cohorts to limit the data downloads to. Cohorts should be
#'               of length 5 representing the shortened dates for the cohort. For example "15-16" or c("01-02", "03-04").
#' @param save_directory String providing the directory in which to save the NHANES data files.  If left as the default,
#'                       NULL, it will save to ./rawData.  Otherwise, it will save to save_directory/rawData.
#'
#' @return
#' @importFrom readxl read_excel
#' @importFrom utils download.file
#' @export
#'
#' @examples
#'
#'
get_NHANES_data <- function(codes_file = NULL, cohort = NULL, save_directory = NULL) {

  if (is.null(codes_file)){
    print("Error: please provide a file name")
    stop()
  }

  codes <- as.data.frame(read_excel(codes_file, sheet = 1))
  weights <- as.data.frame(read_excel(codes_file, sheet = 2))


  # Generate list of needed laboratory files (has metabolite concentration measurements)
  samps <- unique(codes$recent_sample)

  # Limit by cohort
  if (!is.null(cohort)){
    samps <- samps[samps %in% cohort]
  }

  labFiles <- list()
  labFiles <- lapply(samps, function(x) unique(codes$NHANESfile[codes$recent_sample == x]))
  names(labFiles) <- samps


  ### Reference table for phases and file conventions
  ref <- data.frame("Full" = c("1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010",
                               "2011-2012", "2013-2014", "2015-2016", "2017-2018", "2019-2020"),
                    "Short" = c("99-00", "01-02", "03-04", "05-06", "07-08", "09-10", "11-12", "13-14", "15-16", "17-18", "19-20"),
                    "Ext" = c("", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"),
                    stringsAsFactors = FALSE)


  ### Add demographic, bodyweight, and creatinine files
  result <- list()
  for (i in 1:length(samps)) {
    ind <- match(samps[i], weights$sample)
    result[[i]] <- c(paste(labFiles[[i]], ".XPT", sep = ""),
                     paste(ifelse(samps[i] != "99-00", "DEMO_", "DEMO"), ref$Ext[ref$Short == samps[i]], ".XPT", sep = ""),
                     paste(ifelse(samps[i] != "99-00", "BMX_", "BMX"), ref$Ext[ref$Short == samps[i]], ".XPT", sep = ""),
                     toupper(weights$creatfile[ind]), toupper(weights$urineflow[ind]))
  }
  names(result) <- samps


  ### Download the needed files from the NHANES website and save in the correct directory
  print("Starting Downlodas")
  oldw <- getOption("warn")
  options(warn = -1)
  for (i in 1:length(samps)) {
    indP <- match(samps[i], ref$Short)

    # Create rawData directory
    if (is.null(save_directory)){
      if (!dir.exists("rawData")) {
        dir.create("./rawData")
      }
    } else {
      if (!dir.exists(save_directory)){
        dir.create(save_directory)
      }
      if (!dir.exists(file.path(save_directory, "rawData/"))) {
        dir.create(file.path(save_directory, "rawData/"))
      }
    }

    # Download the files
    for (j in 1:length(result[[i]])) {
      if (result[[i]][j] == ""){
        next
      }
      if (is.null(save_directory)){
        if (!dir.exists(paste("rawData/", ref$Full[indP], sep = ""))) {
          dir.create(paste("./rawData", ref$Full[indP], sep = ""))
        }
        print(paste("https://wwwn.cdc.gov/nchs/nhanes/", ref$Full[indP], "/", result[[i]][j], sep = ""))
        tryCatch(download.file(paste("https://wwwn.cdc.gov/nchs/nhanes/", ref$Full[indP], "/", result[[i]][j], sep = ""),
                      destfile = paste("./rawData/", ref$Full[indP], "/", tolower(result[[i]][j]), sep = ""),
                      mode = "wb"), error = function(e) print(paste(result[[i]][j], 'was not found', sep = " ")))
      } else {
        if (!dir.exists(file.path(save_directory, paste("rawData/", ref$Full[indP], sep = "")))) {
          dir.create(file.path(save_directory, paste("rawData/", ref$Full[indP], sep = "")))
        }
        print(paste("https://wwwn.cdc.gov/nchs/nhanes/", ref$Full[indP], "/", result[[i]][j], sep = ""))
        tryCatch(download.file(paste("https://wwwn.cdc.gov/nchs/nhanes/", ref$Full[indP], "/", result[[i]][j], sep = ""),
                      destfile = file.path(save_directory,
                                           paste("rawData/", ref$Full[indP], "/", tolower(result[[i]][j]), sep = "")),
                      mode = "wb"), error = function(e) print(paste(result[[i]][j], 'was not found', sep = " ")))
      }
    }
  }
  options(warn = oldw)

}




