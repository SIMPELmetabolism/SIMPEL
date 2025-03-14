#' @title Plot barplots for MIDs comparison between the conditions at a selected time point
#'
#' @description This function takes an input time point and extract the MIDs for a barplot comparison plotting.
#'
#' @param data object (1) the "MIDs" or the "scaled_MIDs" object created using \link{get_table_objects_NA_corrected} function
#' @param time numeric (1), time point of interest
#' @param axisTitle string (1) y-axis text on the plot
#' @param plotTitle string (1) title to be used for barplots i.e. "Bin" or "Compound" column
#'
#' @seealso \link{get_table_objects_NA_corrected}, \link{MIDplot}
#'
#' @return The function returns a list that contains all figures generated. Plots will also be saved to the current directory.
#'
#' @export
#'
#' @examples
#' \donttest{try(MIDs_comp <- MIDs_comparison(test_13C15NGlutamine$MIDs_tableNAcorrected, time = 8, axisTitle = "MID", plotTitle = "Compound"))}

MIDs_comparison <- function(data, time, axisTitle = "MID", plotTitle = "Bin"){
  if(length(time) != 1){
    stop("Can only select a single time point")
  }
  data_cols <- colnames(data)[colnames(data) %like% paste0("X",time,"_")]
  if(length(data_cols) == 0){
    stop("Selected time point not available in the data")
  }
  if(plotTitle == "Compound" & !"Compound" %in% colnames(data)){
    stop("Compound column is not present in the data, choose Bin as the plot title")
  }
  if(plotTitle == "Compound"){
    data <- data %>%
      dplyr::select(-Bin) %>%
      dplyr::mutate(Bin = Compound)
  }
  data <- data %>%
    dplyr::select(all_of(c("mz", "rt", "polarity", "Isotopologue", "Bin", data_cols)))
  plot_data <- reshape2::melt(data[,c("Bin", "Isotopologue", data_cols)],
                              id.vars = c("Bin","Isotopologue"), variable.name = "category") %>%
    dplyr::mutate(category = sapply(strsplit(as.character(category), '[_]' ), `[` , 2)) %>%
    dplyr::group_by(Bin, Isotopologue, category) %>%
    summarise_at(vars(value), list(val = mean, se = ~sd(.)/sqrt(n())))
  plot_list <- list()
  for (i in unique(plot_data$Bin)) {
    current_bin <- plot_data %>%
      dplyr::filter(Bin == i) %>%
      dplyr::mutate(Isotopologue = factor(Isotopologue))
    # reorder the isotopologues
    if(grepl("M", current_bin$Isotopologue[1])){
      current_bin$Isotopologue = factor(current_bin$Isotopologue, 
                                        levels = levels(current_bin$Isotopologue)[order(as.numeric(gsub("M", "", levels(current_bin$Isotopologue))))])
    }else{
      calculate_total_isotopes <- function(factor_value) {
        sum(as.numeric(unlist(regmatches(factor_value, gregexpr("[0-9]+", factor_value)))))
      }
      current_bin$Isotopologue = factor(current_bin$Isotopologue, 
                                        levels = levels(current_bin$Isotopologue)[order(sapply(levels(current_bin$Isotopologue), calculate_total_isotopes))])
    }
    p <- ggplot(current_bin, aes(x = Isotopologue, y = val, fill = category)) +
      geom_col(position="dodge") +
      geom_errorbar(aes(ymin = pmax(0, val-se), ymax = pmin(val+se, 100)), width = 0.5, position = position_dodge(0.9), linewidth = 0.8) +
      labs(title = paste(i,"at time",time), x = "Isotopologue", y = axisTitle) +
      ylim(min(0, current_bin$val-current_bin$se, na.rm = TRUE), 1.05*max(current_bin$val+current_bin$se, na.rm = TRUE))
    plot_list <- append(plot_list, list(p))
  }
  allPlots <- marrangeGrob(plot_list, nrow = 2, ncol = 2)
  ggsave(paste0("MIDs_comparison_time",time,".pdf"), allPlots, width = 22, height = 11)
  return(plot_list)
}
