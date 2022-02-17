#' plot_subpops
#'
#' Plot the median exposure estimate for all subpopulations of each parent
#' chemical that have calculated exposure estimates.  Sorted from highest
#' to lowest predicted exposure.
#'
#' @param path_to_lPs String representing the location of the lPsamps outputs.
#'                    All subpopulations that are desired in the plot need to
#'                    be within the directory provided.
#' @param cohort A string indicating one of the NHANES cohorts or phases (e.g.,
#'               "01-02" or "09-10")
#' @param SUBPOP A character vector of the subpopulations to be plotted.  These
#'               must exactly match the subpopulations contained in the original
#'               codes_file.
#' @param save_output Set to TRUE to save the results. Note that the default
#'                    save location is the current working directory.
#' @param save_directory String giving the path of where to save the output table
#'                       and plot.
#' @param print_plot Boolean input. Default is FALSE. Set to TRUE to plot in the
#'                   current R session.
#'
#' @return exposureDF:  data.frame with rows representing the chemicals and columns
#'                      representing the different subpopulations
#'         Plot: The plot is printed for viewing, but you can choose where to save
#'               it or view it again using print(p)
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @export
#'
#'
plot_subpops <- function(path_to_lPs = NULL, cohort = NULL, SUBPOP = "Total", save_output = FALSE, save_directory = ".", print_plot = FALSE) {

  path_to_lPs = paste(save_directory, "/", cohort, sep = "")

  # Valid subpop?
  subpops <- c("Total", "Male", "Female", "0 - 5", "6 - 11 years", "12 - 19 years", "20 - 65 years", "66 years and older",
               "ReproAgeFemale", "BMI <= 30", "BMI > 30")
  try(if(!SUBPOP %in% subpops) stop(paste("please supply one of the following possible subpops:", subpops, sep = "")))

  # Valid cohort?
  cohorts <- c("99-00", "01-02", "03-04", "05-06", "07-08", "09-10", "11-12", "13-14", "15-16")
  try(if(!cohort %in% cohorts) stop(paste("please supply one of the following possible cohorts:", cohorts, sep = "")))

  print("Retrieving file names for the subpopulations:")
  print(unname(subpops))

  # Obtain file names for lPsamps
  if (is.null(path_to_lPs)){
    files <- list.files(path = ".", pattern = "lPsamps-gm_kg_day")
  } else {
    files <- list.files(path = path_to_lPs, pattern = "lPsamps-gm_kg_day")
  }

  # Read in files and build data.frame
  load(paste(path_to_lPs, files[grep(paste("_", SUBPOP, sep = ""), files)], sep = ""))
  samps <- as.matrix(lPsampsgm)
  result.med <- apply(samps, 2, mean)
  result.bar <- mean(apply(samps, 1, mean))
  # Get 95th Percentiles
  result.CI <- t(apply(samps, 2, function(x) quantile(x, probs = c(0.0275, 0.975))))
  dta <- data.frame(colnames(samps), exp(result.med))
  dta.low <- data.frame(colnames(samps), exp(result.CI[,1]))
  dta.up <- data.frame(colnames(samps), exp(result.CI[,2]))
  colnames(dta) <- c("CAS", SUBPOP)
  colnames(dta.low) <- c("CAS", paste(SUBPOP, "_low95", sep = ""))
  colnames(dta.up) <- c("CAS", paste(SUBPOP, "_up95", sep = ""))

  # Combine
  df <- data.frame("Chemical" = dta$CAS, "Median" = dta[,2], "low95" = dta.low[,2],
                   "up95" = dta.up[,2])

  # Sort from largest to smallest predicted exposure
  sorted <- sort(df[,2], decreasing = TRUE, index.return = TRUE)

  # Complete sort
  df$Chemical <- as.factor(df$Chemical)
  df$Chemical <- factor(df$Chemical, levels = levels(df$Chemical)[sorted$ix])

  # Generate plot
  p <- ggplot(df, aes(x=Chemical, y=Median)) +
    geom_point() +
    geom_errorbar(aes(ymin = low95, ymax = up95), width = 0.1) +
    ylab("Predicted Exposure (mg/kg BW/day)") +
    scale_y_log10() +
    theme_set(theme_gray(base_size = 14)) +
    theme(axis.text.x=element_text(angle=90, hjust=1))

  if (print_plot) {
    print(p)
  }

  result <- list(df = df,
                 Plot = p)

  if (save_output) {
    save(result, file=file.path(save_directory, paste("Exposures_", cohort, "_", SUBPOP, "_",
                                      format(Sys.time(), "%Y-%m-%d"), ".RData", sep = "")))
  }

  return(result)

}







