#' plot_subpops
#'
#' Plot the median exposure estimate for all subpopulations of each parent
#' chemical that have calculated exposure estimates.  Sorted from highest
#' to lowest predicted exposure.
#'
#' @param path_to_lPs String representing the location of the lPsamps outputs.
#'                    All subpopulations that are desired in the plot need to
#'                    be within the directory provided.
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
plot_subpops <- function(path_to_lPs = NULL, SUBPOP = "all", save_output = FALSE, save_directory = ".", print_plot = FALSE) {

  # Reference for population groups
  subpops <- c("Total", "Male", "Female", "0 - 5", "6 - 11 years", "12 - 19 years", "20 - 65 years", "66 years and older",
               "ReproAgeFemale", "BMI <= 30", "BMI > 30")
  names(subpops) <- c("Total", "Male", "Female", "0 - 5", "6 - 11 years", "12 - 19 years", "20 - 65 years", "66 years and older",
                      "ReproAgeFemale", "BMI_le_30", "BMI_gt_30")

  # Subpops the user wants plotted
  if (SUBPOP != "all"){
    ind <- SUBPOP %in% subpops
    subpops <- subpops[ind]
  }

  print("Retrieving file names for the subpopulations:")
  print(unname(subpops))

  # Obtain file names for lPsamps
  if (is.null(path_to_lPs)){
    files <- list.files(path = ".", pattern = "lPsamps-gm_kg_day")
  } else {
    files <- list.files(path = file.path(".", path_to_lPs), pattern = "lPsamps-gm_kg_day")
  }

  # Read in files and build data.frame
  for (i in 1:length(subpops)){
    load(files[grep(subpops[i], files)])
    samps <- as.matrix(lPsampsgm)
    ind <- grep("lP\\[", colnames(samps))
    result.med <- apply(samps[,ind], 2, mean)
    result.bar <- mean(apply(samps[,ind], 1, mean))
    dta <- data.frame(row.names(samps), exp(result.med))
    colnames(dta) <- c("CAS", subpops[i])
    if (i > 1){
      dta.1 <- merge(dta.1, dta, by = "CAS", all.x = TRUE)
    } else {
      dta.1 <- dta
    }
  }

  ####################################################
  #load("MCMCdata.RData")
  #dta.1 <- cbind("Name" = pred.data$Name, dta.1)
  #write.csv(dta.1, file = "./updatedPipeline/lPests_mgkgday_bySubpop.csv", row.names = FALSE)
  #rm(Bstart, Bstop, F, index, Measured, NBranches, Ndelta, Phi, Uses, mapping)
  ######################################################

  # Sort from largest to smallest predicted exposure
  ind <- apply(dta.1[2:length(subpops)], 1, function(x) max(x, na.rm = TRUE))
  sorted <- sort(ind, decreasing = TRUE, index.return = TRUE)

  exposureDF <- melt(dta.1, id.vars="CAS", value.name="loggm", variable.name="Demo")

  # Complete sort
  exposureDF$CAS <- factor(exposureDF$CAS, levels = levels(exposureDF$CAS)[sorted$ix])

  # Generate plot
  p <- ggplot(exposureDF, aes(x=CAS, y=loggm, group=Demo, color=Demo)) + geom_point() +
    ylab("Predicted Exposure (mg/kg BW/day)") +
    scale_y_log10(limits = c(1e-10, 10)) +
    theme_set(theme_gray(base_size = 18)) +
    theme(axis.text.x=element_text(angle=90, hjust=1))

  if (print_plot) {
    print(p)
  }

  result <- list(exposureDF = exposureDF,
                 Plot = p)

  if (save_output) {
    save(result, file=file.path(save_directory, paste("CompareSubpops_", SUBPOP, "_",
                                      format(Sys.time(), "%Y-%m-%d"), ".RData", sep = "")))
  }

  return(result)

}







