#' @title Plot time-course label enrichment both as average_labeling and mole_equivalents_labeling
#' in a non-stationary isotopic labeling experiment
#'
#' @description This functions allows plotting of the average_labeling and mol_equivalent_labeling objects for user interpretation of label enrichment over the course of the experiment
#'
#' @param mydata object (1) The 'average_labeling' or the 'mol_equivalent_labeling' object created using \link{get_table_objects_NA_corrected} function
#' @param Category string (1) name of the category by which you are binning, for the example data set - Lipidomics has two categories (i.e. tissue types),
#' "Cotyledon" and "EA" and the DualLabel dataset has two categories (i.e. genotypes), "WT" and "GAT"
#' @param yLim numeric (1) sets the limit of y-axis, default is set to maximum
#' @param xLim numeric (1) sets the limit of the x-axis, default is set to maximum
#' @param axisTitle string (1) The label to be used for the Y-Axis. Typically it is "\% labeling" for average_labeling and "mol equivalents of label" for mol_equivalent_labeling
#' @param plotTitle string (1) The title to be used for the plot - has to be one of the columns within the object.
#' Typically, it is "Compound", however, if the string is too long in the 9x9 output generated, "Bin" can be used as a substitute
#' @param plotTitle2 string (1) Additional title to be used for the plot, default is set to NULL
#' @param outputName string (1) The name to be appended to the the output pdf
#'
#' @seealso \link{get_table_objects_NA_corrected}, \link{MIDplot}, \link{PCA_and_heatmap}, \link{getClustersAndPlots}
#'
#' @return The function returns a list with all plots for label enrichment for all the compounds identified within the XCMS_data along with a pdf file with plots
#'
#' @export
#'
#' @examples
#' \donttest{try(Average_Labeling <- label_enrichment_plot(test_13C15NGlutamine$average_labeling,
#' Category = "WT", yLim=NULL, xLim=NULL, axisTitle=" % Labeling",
#' plotTitle= "Bin", outputName = "Average"))}
#' \donttest{try(Label_Enrichment <- label_enrichment_plot(test_13C15NGlutamine$mol_equivalent_labeling,
#' Category = "WT", yLim=NULL, xLim=NULL, axisTitle="mol equivalents of Label", plotTitle="Compound",
#' plotTitle2= "Bin", outputName = "mol Equi"))}

label_enrichment_plot <- function(mydata, Category, yLim=NULL,xLim=NULL, axisTitle="Labeling", plotTitle="Bin", plotTitle2=NULL, outputName = "_")
{
  #pull out the first field which correpsonds to timepoint of non-stationary labeling experiment
  data_clean <- function(x) sapply (strsplit(x , '[_]' ), `[` , 1)

  #read in the data table or use the object from previous function
  mydata = mydata
  mydata = na.omit(mydata)

  #set the axes limit
  yLim = yLim
  xLim = xLim

  #get the outputName
  outputName = outputName

  #set the axisTitle
  axisTitle = axisTitle
  plotTitle = plotTitle
  plotTitle2 = plotTitle2

  #workaround to the yLim resetting issue
  #we reset yLim to the max(avg + sd)*1.1
  #and then it is no longer null, so we need to know
  #if it was originally null
  YoriginallyNull = FALSE
  if(is.null(yLim) == TRUE)
  {
    YoriginallyNull = TRUE
  }

  XoriginallyNull = FALSE
  if(is.null(xLim) == TRUE)
  {
    XoriginallyNull = TRUE
  }

  #combine the multiple columns into a single column with unique names
  names = make.names(paste(mydata$Bin,mydata$rt,mydata$mz,sep = "_"), unique = TRUE)
  rownames(mydata) = names

  #what is the identifier for this group (i.e. the category, a condition or tissue group)
  Category = Category

  #pull out the condition/species specific information
  #by only those columns which match the Category
  mydata = mydata[,(!colnames(mydata) %like% "_") | (colnames(mydata) %like% Category)]

  #generate the pdf
  name = paste(Category,outputName,"label_enrichment.pdf", sep = "_")

  myListPlot = list()

  for(i in 1:nrow(mydata))
  {
    theVecOfInfo = mydata[i,]


    #if using id's for title
    #set title before we make it numeric
    title = theVecOfInfo[colnames(theVecOfInfo) == plotTitle][1,]


    if(!is.null(plotTitle2))
    {
      title = paste(theVecOfInfo[colnames(theVecOfInfo) == plotTitle][1,], theVecOfInfo[colnames(theVecOfInfo) == plotTitle2][1,], sep = "  ")
    }

    ##calculate each time point and each std dev

    #get the set of times from the column names
    myTimes = data_clean(as.character(colnames(mydata)))
    #pull out any non-numeric elements of the timepoints
    myTimes =  suppressWarnings(as.numeric(gsub("X","",as.character(myTimes))))
    #remove any non-valid timepoints
    myTimes = myTimes[is.na(myTimes) == FALSE]
    #get a unique vector of the timepoints
    myTimesUnique = unique(myTimes)

    #hold all the Avgs and standard deviations
    allMyInfoAvg = vector()
    allMyInfoSD = vector()

    ###get the average and the sd for each of the metabolites
    #for each time
    for(k in 1:length(myTimesUnique))
    {

      #get the time
      myTime = myTimesUnique[k]

      #pull out first field,i.e time
      myTimepointMatch = gsub("(*.*)_.*_.*", "\\1", names(theVecOfInfo))

      #pull out non-numeric characters
      myTimepointMatch = gsub("[^0-9.-]", "", myTimepointMatch)

      #rename the colnames with just the time information
      #I think that this is redundant, actually
      names(theVecOfInfo) = gsub("[^0-9.-]", "", names(theVecOfInfo))

      #make sure that it's numeric now
      myTimepointMatch = as.numeric(myTimepointMatch)

      #set the names to just the timepoint now, it's all we need
      names(theVecOfInfo) = as.numeric(myTimepointMatch)

      #calculate the mean for this timepoint
      myTimeMean = mean(as.numeric(theVecOfInfo[,as.numeric(colnames(theVecOfInfo)) %in% c(myTime)]))

      #calculate the SD for this timepoint
      myTimeSD = sd(as.numeric(theVecOfInfo[,as.numeric(colnames(theVecOfInfo)) %in% c(myTime)])) /
        sqrt(length(as.numeric(theVecOfInfo[,as.numeric(colnames(theVecOfInfo)) %in% c(myTime)])))

      allMyInfoAvg = c(allMyInfoAvg, myTimeMean)
      allMyInfoSD = c(allMyInfoSD, myTimeSD)
    }

    if(is.null(yLim) == TRUE)
    {
      # set yLim to 10% more than the max possible value
      yLim = (max(allMyInfoAvg) + max(c(0, allMyInfoSD), na.rm = TRUE)) * 1.1
      #print(yLim)
      #print("is yLim")
      #print(max(allMyInfoAvg+allMyInfoSD))
      #print("other potential max")
    }


    #if the user doesn't specify the x-limit, we're going
    #to use the max as the highest timepoint
    if(is.null(xLim) == TRUE)
    {
      xLim = max(myTimesUnique)
    }

    #set the y-limit to at least the value of the highest error bar
    if(yLim == "default")
    {
      yLim = max(allMyInfoAvg[!is.na(allMyInfoAvg)]) + allMyInfoAvg + allMyInfoSD

    }


    #declare the times which will be the x-axis of the plots
    x = myTimesUnique

    # put the data together for ggplot
    dataForPlot = data.frame(x, allMyInfoAvg, allMyInfoSD)

    #make the plots (time on the x-axis)
    #the label enrichment on the y-axis

    p = ggplot(dataForPlot, aes(x=x, y=allMyInfoAvg)) +
      geom_point(size = 2.5) +
      geom_errorbar(aes(ymin=allMyInfoAvg-allMyInfoSD, ymax=allMyInfoAvg+allMyInfoSD), width=0.4) +
      labs(title=title, x ="Time", y = axisTitle) +
      xlim(range(c(-0.4, xLim+0.4))) + ylim(range(c(min(0, allMyInfoAvg-allMyInfoSD), yLim)))

    myListPlot = append(myListPlot, list(p))

    #reset yLlim to null if needed
    if(YoriginallyNull == TRUE)
    {
      yLim = NULL
    }

    #reset xLlim to null if needed
    if(XoriginallyNull == TRUE)
    {
      xLim = NULL
    }

  }
  allPlots = marrangeGrob(myListPlot, nrow = 3, ncol = 3)
  ggsave(name, allPlots, width = 11, height = 11)

  return(myListPlot)
}
