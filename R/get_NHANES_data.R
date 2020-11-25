

### R script for downloading NHANES data files needed for exposure inference


#setwd("C:/Users/zstanfie/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/NHANES")
#library(plyr)
#library(readxl)

get_NHANES_data <- function(codes_file = NULL) {

  if (is.null(codes_file)){
    print("Error: please provide a file name")
    stop()
  }

  codes <- read.xls(codes_file, as.is = TRUE)
  weights <- read.xls(codes_file, sheet = 2, as.is = TRUE)
  mapping <- read.xls(codes_file, sheet = 3, as.is = TRUE)


  # Generate list of needed laboratory files (has metaboite concentration measurements)
  samps <- unique(codes$recent_sample)
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
  weights <- read.xls(codes_file, sheet = 2)
  result <- list()
  for (i in 1:length(samps)) {
    ind <- match(samps[i], weights$sample)
    result[[i]] <- c(labFiles[[i]], paste(ifelse(samps[i] != "99-00", "DEMO_", "DEMO"), ref$Ext[ref$Short == samps[i]], ".XPT", sep = ""),
                     paste(ifelse(samps[i] != "99-00", "BMX_", "BMX"), ref$Ext[ref$Short == samps[i]], ".XPT", sep = ""),
                     toupper(weights$creatfile[ind]))
  }
  names(result) <- samps


  ### Download the needed files from the NHANES website and save in the correct directory
  for (i in 1:length(samps)) {
    indP <- match(samps[i], ref$Short)
    if (!dir.exists(paste("rawData/", ref$Full[indP], sep = ""))) {
      dir.create(paste("rawData/", ref$Full[indP], sep = ""))
    }
    for (j in 1:length(result[[i]])) {
      #print(paste("file being downloaded:  ", "https://wwwn.cdc.gov/nchs/nhanes/", ref$Full[indP], "/", result[[i]][j], sep = ""))
      download.file(paste("https://wwwn.cdc.gov/nchs/nhanes/", ref$Full[indP], "/", result[[i]][j], sep = ""),
                    destfile = paste("rawData/", ref$Full[indP], "/", tolower(result[[i]][j]), sep = ""), mode = "wb")
    }
  }


}




