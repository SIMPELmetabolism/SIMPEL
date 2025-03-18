#' @title Plot barplots for comparison between the conditions at one or more time points
#'
#' @description This function takes input time point(s) and extract the label enrichment information for a side-by-side barplot comparison plotting.
#' When only a single time is supplied, the function also returns a table containing Bins that have significant different labeling between the two
#' categories at that time point.
#'
#' @param data object (1) The 'average_labeling' or the 'mol_equivalent_labeling' object created using \link{get_table_objects_NA_corrected} function
#' @param time numeric vector, time point(s) of interest. Default is the last time point
#' @param axisTitle string (1) y-axis text on the plot
#' @param plotTitle string (1) title to be used for barplots i.e. "Bin" or "Compound" column
#'
#' @seealso \link{get_table_objects_NA_corrected}, \link{label_enrichment_plot}
#'
#' @return The function returns a list that contains all figures generated and if a single time point is selected, a dataframe with compounds that are
#' labeled significantly different. Plots (and a table) will also be saved to the current directory.
#'
#' @export
#'
#' @examples
#' \donttest{try(barplots_comp <- barplots_comparison(test_13C15NGlutamine$average_labelingNAcorrected, time = "default", axisTitle = "% Labeling", plotTitle = "Compound"))}

barplots_comparison <- function(data, time = "default", axisTitle = "% Labeling", plotTitle = "Bin"){
  if(plotTitle == "Compound" & !"Compound" %in% colnames(data)){
    stop("Compound column is not present in the data, choose Bin as the plot title")
  }
  data_cols <- colnames(data)[colnames(data) %like% "_.+_"]
  times <- sapply(strsplit(data_cols, '[_]' ), `[` , 1) %>%
    sub(".", "", .) %>%
    as.numeric()
  if(identical(time, "default")){
    time <- max(times)
  }
  time_cols <- data_cols[times %in% time] # pull out columns with selected time points
  plot_list <- list()
  sig_diff_bins <- c()
  for (i in 1:nrow(data)){
    plot_data <- data.frame(value=unlist(data[i, time_cols])) %>%
      dplyr::mutate(time = as.numeric(sub(".","",sapply(strsplit(rownames(.), '[_]' ), `[` , 1))),
                    category = sapply(strsplit(rownames(.), '[_]' ), `[` , 2))
    # if a single time point was provided, check for significance
    if(length(time) == 1 &
       (t.test(value ~ category, data = plot_data)$p.value < 0.05)){
      sig_diff_bins <- append(sig_diff_bins, i)
    }
    plot_data <- plot_data %>%
      dplyr::group_by(time, category) %>%
      summarise_at(vars(value), list(val = ~mean(., na.rm = TRUE), se = ~sd(., na.rm = TRUE)/sqrt(length(na.omit(.))))) %>%
      ungroup()
    p <- ggplot(plot_data, aes(x=factor(time), y=val, fill=category)) +
      geom_col(position="dodge") +
      geom_errorbar(aes(ymin=pmax(0,val-se), ymax=pmin(val+se,100)), width=0.5,
                    position=position_dodge(0.9), linewidth=0.8) +
      labs(title=ifelse(plotTitle == "Bin", data$Bin[i], data$Compound[i]),
           x = "Time", y = axisTitle) +
      ylim(0, 1.1*max(plot_data$val+plot_data$se, na.rm = TRUE))
    plot_list <- append(plot_list, list(p))
  }
  allPlots <- marrangeGrob(plot_list, nrow = 3, ncol = 3)
  ggsave(paste0("labeling_barplots_time_", paste0(time, collapse = "&"),".pdf"), allPlots, width = 11, height = 11)
  if(length(time)==1){
    if(plotTitle == "Bin"){
      sig_bins_df <- data[sig_diff_bins, c("Bin", "mz", "rt")]
    }else{
      sig_bins_df <- data[sig_diff_bins, c("Bin", "Compound", "mz", "rt")]
    }
    write.csv(sig_bins_df, "labeling_barplots_sig_bins.csv", row.names = FALSE)
    return(list(allplots = plot_list, sig_bins = sig_bins_df))
  }else{
    return(list(allplots = plot_list))
  }
}
