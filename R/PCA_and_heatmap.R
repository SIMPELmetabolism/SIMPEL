#' @title Perform global analysis of label enrichment and generate PCA plots and heatmaps
#'
#' @description This function performs Principle Component Analysis and plots the average_labeling or mol_equivalent_labeling data in a two dimensional space
#' to allow comparison of time course labeling enrichment data obtained from a non-stationary isotopic labeling experiment. In addition, this function
#' also plots the label enrichment for each of the compounds as a heat map to easily identify compounds that have significant labeling using global view
#'
#' @param mydata1 object (1) The 'average_labeling' or the 'mol_equivalent_labeling' object created using \link{get_table_objects_NA_corrected} function
#' @param heatMapCategories string (1) or (2) c("Category1") or c("Category1", "Category2") The category to be used for plotting a heatmap
#' @param PCMax numeric (1) maximum numbers of PC's to plot
#' @param labels string (1) label to be used for heat maps i.e. "Bin" or "Compound" column
#' @param outputName string (1) The name to be appended to the the output pdf
#'
#' @seealso \link{get_table_objects_NA_corrected}, \link{MIDplot}, \link{label_enrichment_plot}, \link{getClustersAndPlots}
#'
#' @return The function returns a list containing generated figures along with a pdf file with PCA plots comparing all PCs up to the user specified number and heat maps with label enrichment
#' for all compounds using the user specified Categories
#'
#' @importFrom gplots heatmap.2
#'
#' @export
#'
#' @examples
#' \donttest{try(PCA_Heatmap_Avg <- PCA_and_heatmap(test_13C15NGlutamine$average_labeling, PCMax = 3,
#' heatMapCategories = c("WT", "GAT"), labels="Bin", outputName = "Average"))}
#' \donttest{try(PCA_Heatmap_molEqui <- PCA_and_heatmap(test_13C15NGlutamine$mol_equivalent_labeling, PCMax = 3,
#' heatMapCategories = c("WT", "GAT"), labels="Bin", outputName = "mol_equivalents"))}

PCA_and_heatmap <- function(mydata1, PCMax=3, heatMapCategories, labels="Bin", outputName = "_")
{
  PCMax = PCMax
  allTotFigures = list()
  heatMapCategories = heatMapCategories
  mydata1 = as.data.frame(mydata1)
  outputName = outputName


  
  #remove bin column from table to do PCA on
  rownames(mydata1) =  pull(mydata1, labels)
  mydata1Backup = mydata1
  mydata1[labels] = NULL



  #make sure that there are now all zeroes columns or na containing columns
  #we're going to exclude the column labels whose rows do not include the numeric data
  vecToExclude = c("mz","polarity", "rt", "comp_result","Formula", "carbon", "nitrogen" , "total_isotopes", "Bin", "Compound", "CompoundPlaceholder", "Isotopologue" )

  #we'll only want to use the columns with numeric values
  columnsToUse = setdiff(colnames(mydata1), vecToExclude)

  #exclude the other columns
  mydata1 = mydata1[,colnames(mydata1) %in% columnsToUse]


  #exclude all of the nun-numeric columns
  mydata1 = mydata1[complete.cases(mydata1), ]
  mydata1 = mydata1[rowSums(mydata1) != 0 ,]


  #make a metadata table (DataFrameLabel) in order to do the PCA
  DataFrameLabel = colnames(mydata1)
  DataFrameLabel = as.data.frame(DataFrameLabel)
  colnames(DataFrameLabel) = c("Samples")

  #make sure unintended characters are removed with this function
  #and that we're including a column with the category + time as well as a separate
  #one with just the time
  DataFrameLabel$CategoryTime = paste(data_cleanII(as.character(DataFrameLabel$Samples)),data_clean(as.character(DataFrameLabel$Samples)), sep = "_")
  DataFrameLabel$Time = sub(".", "", data_clean(as.character(DataFrameLabel$Samples)))

  #just make sure that we're always having category come first
  #originally this was Category
  DataFrameLabel$Category = gsub("*.*_(.*)_.*", "\\1", DataFrameLabel$Samples)

  #get PCA object to determine maximum number of PC's
  prcompObject = prcomp(t(mydata1), center = TRUE, scale. = TRUE)

  if(PCMax > 0)
  {
    #now make the PCA autoplots
    for(x in 1:PCMax)
    {
      for(y in 1:PCMax)
      {
        if(x > y)
        {
          if( x <= ncol(prcompObject$x) & y <=  ncol(prcompObject$x))
          {
            #second plot, coloring by time and then Category as well
            autoplotList = list()
            autoplotList[[1]] = print(autoplot(prcomp(t(mydata1), center = TRUE, scale. = TRUE), 
                                               data =  DataFrameLabel, colour = 'Time', shape = 'Category', 
                                               frame.colour = 'CategoryTime', size = 3))
            allTotFigures = append(allTotFigures,autoplotList)
          }
        }
      }
    }
  }


  heatMapList = list()
  toPrintFile = paste(paste0(heatMapCategories, collapse = ""), outputName, "heatmap_and_PCA.pdf")
  for(i in 1:length(heatMapCategories))
  {
    catSubset = mydata1[,colnames(mydata1) %like% heatMapCategories[i]]
    heatmapName = paste("Heatmap for", heatMapCategories[i], sep = " ")
    heatmap.2(as.matrix(catSubset), Colv = F, scale = c("none"), dendrogram = "row", trace="none",
              col = cm.colors(16), keysize = 1.5, margins = c(8,12), cexCol = 0.6,
              main = heatmapName)
    heatMapList[[i]] = recordPlot()
  }

  returnList = c(heatMapList, allTotFigures)
  pdf(toPrintFile)
  print(lapply(returnList, function(x) x))
  dev.off()
  return(returnList)
}
