#' @title Plot pool sizes along with label enrichment over time
#'
#' @description This function takes summation of intensities of each compound as a representation of the pool sizes
#' at each time, then visualize pool sizes with label incorporation at the same time.
#'
#' @param forBubblePlot object (1) data object for bubble plots, created using \link{get_table_objects_NA_corrected} function
#' @param labeling_table, object (1) label enrichment information, created using \link{get_table_objects_NA_corrected} function
#' @param plotTitle string (1) title to be used for bubble plots i.e. "Bin" or "Compound" column
#' @param labelLegend string (1) legend title to be used for labeling information. For average labeling, it could be "Enrichment (%)", for mol equivalents values, it could be "Enrichment (mol equivalents)"
#'
#' @seealso \link{get_table_objects_NA_corrected}, \link{label_enrichment_plot}
#'
#' @return The function returns a list that contains all figures generated and saves a pdf file to the current directory.
#'
#' @importFrom tidyr pivot_longer
#'
#' @export
#'
#' @examples
#' \donttest{try(bubble_plots <- get_bubble_plot(test_13C15NGlutamine$forBubblePlot,
#' test_13C15NGlutamine$average_labelingNAcorrected, plotTitle = "Compound", labelLegend = "Enrichment (%)"))}

get_bubble_plot <- function(forBubblePlot, labeling_table, plotTitle = "Bin", labelLegend = "Enrichment (%)"){
  if(plotTitle == "Compound"){
    if(! "Compound" %in% colnames(forBubblePlot)){
      stop("Compound column is not present in the data, choose Bin as the plot title")
    }
    forBubblePlot <- forBubblePlot %>%
      dplyr::mutate(Bin = Compound)
    labeling_table <- labeling_table %>%
      dplyr::mutate(Bin = Compound)
  }
  labeling_table <- labeling_table %>%
    pivot_longer(
      cols = matches("_.+_"),
      names_to = c("Time", "Category", "Replicate"),
      names_pattern = "X(\\d+\\.?\\d*)_(\\w+)_(\\d+)"
    ) %>%
    dplyr::mutate(Time = as.numeric(Time)) %>%
    dplyr::group_by(Bin, Time, Category) %>%
    dplyr::summarise(enrichment = mean(value, na.rm = TRUE))
  forBubblePlot <- forBubblePlot %>%
    dplyr::group_by(Bin) %>%
    dplyr::summarize(across(matches("_.+_"), ~ sum(., na.rm = TRUE))) %>%
    pivot_longer(
      cols = matches("_.+_"),
      names_to = c("Time", "Category", "Replicate"),
      names_pattern = "X(\\d+\\.?\\d*)_(\\w+)_(\\d+)"
    ) %>%
    dplyr::mutate(Time = as.numeric(Time)) %>%
    dplyr::group_by(Bin, Time, Category) %>%
    dplyr::summarise(poolSize = mean(value, na.rm = TRUE)) %>%
    dplyr::left_join(., labeling_table) %>%
    ungroup()

  min_enrichment <- min(forBubblePlot$enrichment, na.rm = TRUE)
  max_enrichment <- max(forBubblePlot$enrichment, na.rm = TRUE)
  min_poolSize <- min(forBubblePlot$poolSize, na.rm = TRUE)
  max_poolSize <- max(forBubblePlot$poolSize, na.rm = TRUE)

  round_to_two_digits <- function(x) {
    # Calculate the order of magnitude of the number
    magnitude <- 10^(floor(log10(abs(x))) - 1)
    # Round to the first two significant digits
    rounded <- round(x / magnitude) * magnitude
    return(rounded)
  }

  quantiles <- quantile(forBubblePlot$poolSize, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  pool_size_breaks <- sapply(quantiles, round_to_two_digits)
  pool_size_labels <- format(pool_size_breaks, scientific = TRUE)

  bubblePlots <- lapply(unique(forBubblePlot$Bin), function(bin){
    subdata = forBubblePlot %>%
      dplyr::filter(Bin == bin)
    # Creating the bubble plot
    ggplot(subdata, aes(x = factor(Time), y = Category, size = sqrt(poolSize), fill = enrichment)) +
      geom_point(shape = 21, color = "black", alpha = 0.7) +
      scale_size_continuous(range = c(1, 15),
                            limits = c(sqrt(min_poolSize), sqrt(max_poolSize)),
                            breaks = sqrt(pool_size_breaks),
                            labels = pool_size_labels) +
      scale_fill_gradient(low = "purple", high = "yellow",
                          limits = c(min_enrichment, max_enrichment)) +
      theme_minimal() +
      labs(title = paste0("Enrichment Dynamics - ", bin),
           x = "Time",
           y = "Category",
           size = "Pool Size",
           fill = labelLegend) +
      theme(legend.position = "bottom",
            legend.box = "vertical")
  })
  if(plotTitle == "Compound"){
    allPlots <- marrangeGrob(bubblePlots, nrow = 4, ncol = 2)
  }else{
    allPlots <- marrangeGrob(bubblePlots, nrow = 4, ncol = 3)
  }
  ggsave("BubblePlots.pdf", allPlots, width = 12, height = 14)
  return(bubblePlots)
}
