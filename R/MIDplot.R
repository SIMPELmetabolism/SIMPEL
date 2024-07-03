#' @title Plot Mass Isotopologue Distribution, i.e. MIDs, in a non stationary isotopic
#' labeling experiment
#'
#' @description This function allows plotting of the MID and scaled_MID objects for user interpretation of isotopologue specific label enrichment over the course of the experiment
#'
#' @param myData object (1) The "MIDs" or the "scaled_MIDs" object created using \link{get_table_objects_NA_corrected} function
#' @param Category string (1) name of the category by which you are binning, for the example data set - Lipidomics has two categories (i.e. tissue types),
#' "Cotyledon" and "EA" and the DualLabel dataset has two categories (i.e. genotypes), "WT" and "GAT"
#' @param splitIsotopologue string (1) The isotopologue to be plotted separately in "split_MIDs". Generally if there is low amount of labeling in certain compounds,
#' splitting the MID plot by the unlabeled isotopologue (i.e. "C0N0" - set as default) will allow better visualization of label enrichment in the other isotopologues
#' @param axisTitle string (1) The label to be used for the Y-Axis. Typically it is "MIDs" for MIDs and "mol equivalents of label" for scaled_MIDs
#' @param plotTitle string (1) The title to be used for the plot - has to be one of the columns within the object. Typically, it is "Compound", however, if the string is too long in the 9x9 output generated, "Bin" can be used as a substitute
#' @param plotTitle2 string (1) Additional title to be used for the plot, default is set to NULL
#' @param yLimit numeric (1) sets the limit of y-axis, default is set to maximum
#' @param xLimit numeric (1) sets the limit of the x-axis, default is set to maximum
#' @param outputName string (1) The name to be appended to the the output pdf
#'
#' @seealso \link{get_table_objects_NA_corrected}, \link{label_enrichment_plot}, \link{PCA_and_heatmap}, \link{getClustersAndPlots}
#'
#' @return The function returns all plots generated along with two pdf files one with the MIDs of individual compounds as a single plot and the other with MIDs split by a user specified isotopologue for better visualization
#' of low labeled isotopologues
#'
#' @export
#'
#' @examples
#' \donttest{try(MIDs <- MIDplot(test_13C15NGlutamine$MIDs, Category = "WT",
#' splitIsotopologue = "C0N0", axisTitle = "MID", plotTitle = "Bin"))}
#' \donttest{try(scaled_MIDs <- MIDplot(test_DualLabel$scaled_MIDs, Category = "GAT",
#' axisTitle = "mol_equivalents", plotTitle = "Compound", plotTitle2 = "Bin", outputName = "scaled"))}

MIDplot <- function(myData, Category, splitIsotopologue = "default", axisTitle="MID", plotTitle="Bin", plotTitle2=NULL, yLimit="default", xLimit="default", outputName = "_")
{


  #pull out the first field which correpsonds to timepoint of non-stationary labeling experiment
  data_clean <- function(x) sapply (strsplit(x , '[_]' ), `[` , 1)


  #default options
  #splitIsotopologue = splitIsotopologue
  yLimit <- yLimit
  myData <- myData
  plotTitle <- plotTitle
  plotTitle2 <- plotTitle2
  axisTitle <- axisTitle
  Category <- Category
  outputName <- outputName
  if(splitIsotopologue == "default"){
    if(grepl('M', myData$Isotopologue[1])){ # if single-labeled, split on M0
      splitIsotopologue <- "M0"
    }else{
      whichLabelScheme <- which(sapply(c("C", "N", 'H'), function(x) grepl(x, myData$Isotopologue[1])))
      splitIsotopologue <- paste0(names(whichLabelScheme), 0, collapse = "")
    }
  }else{
    splitIsotopologue <- splitIsotopologue
  }

  # focus on the just the columns with the category/tissue of interest
  myData <- myData[,(!colnames(myData) %like% "_.+_") | (colnames(myData) %like% Category)]

  # make the rownames all the relevant identifying information (bin, retention time, m.z)
  # but make sure these are unique (just in case) using the make.names function
  names <- make.names(paste(myData$rt,myData$mz,myData$Bin, sep = "_"), unique = TRUE)
  rownames(myData) <- names
  # list of the plots
  p <- list()
  # plots with all of them together
  pAll <- list()
  # loop through each bin
  for(i in 1:length(unique(myData$Bin)))
  {

    theBin <- unique(myData$Bin)[i]
    df <- data.frame(Isotopologue=character(),
                     Time=numeric(),
                     Mean=numeric(),
                     stdDev=numeric())
    subsetBin <- subset(myData, Bin  == theBin)
    subsetBin$Isotopologue <- as.factor(make.unique(as.character(subsetBin$Isotopologue), sep = "."))
    # this is where we get the mass and retention time to format
    theMass <- subsetBin[subsetBin$Isotopologue == splitIsotopologue,]$mz
    retentTime <- subsetBin[subsetBin$Isotopologue == splitIsotopologue,]$rt

    # I've updated this on June 11th 2020 to make sure that the user can customize the column
    # they want to use in order to plot
    CompName <- subsetBin[subsetBin$Isotopologue == splitIsotopologue,][plotTitle]

    title <- paste(CompName)

    # if the user wants to include another column we'll include that in the title
    if(is.null(plotTitle2) == FALSE)
    {
      title <- paste(title, subsetBin[subsetBin$Isotopologue == splitIsotopologue,][plotTitle2], sep = " ")
    }

    # take all the data and put into a new table
    # we'll need to have this data formatted in order to be compatible with the ggplot2 function
    for(j in 1:nrow(subsetBin))
    {

      # get the row with all of the information
      myVecOfAllInfo <- subsetBin[j,]

      # determine how many timepoint

      # calculate each time point mean and each std dev
      myTimes <- data_clean(as.character(colnames(myData)))
      myTimes <- suppressWarnings(as.numeric(gsub("X","",as.character(myTimes))))
      myTimes <- myTimes[is.na(myTimes) == FALSE]
      myTimesUnique <- unique(myTimes)

      if(xLimit == "default")
      {
        xLim <- max(myTimesUnique)
      }

      isotopologueName <- as.character(myVecOfAllInfo$Isotopologue)

      # make sure that we get the information for each timepoint
      for(k in 1:length(myTimesUnique))
      {
        myTime <- myTimesUnique[k]

        # pull out first field
        myTimepointMatch <- gsub("(*.*)_.*_.*", "\\1", names(myVecOfAllInfo))

        # pull out non-numeric stuff
        myTimepointMatch <- gsub("[^0-9.-]", "", myTimepointMatch)
        names(myVecOfAllInfo) <- gsub("[^0-9.-]", "", names(myVecOfAllInfo))

        # make sure that it's numeric now
        myTimepointMatch <- as.numeric(myTimepointMatch)

        # set the names to just the timepoint now, it's all we need
        names(myVecOfAllInfo) <- as.numeric(myTimepointMatch)

        myTimeMean <- mean(as.numeric(myVecOfAllInfo[,as.numeric(colnames(myVecOfAllInfo)) %in% c(myTime)]))
        myTimeSD <- sd(as.numeric(myVecOfAllInfo[,as.numeric(colnames(myVecOfAllInfo)) %in% c(myTime)])) /
          sqrt(length(as.numeric(myVecOfAllInfo[,as.numeric(colnames(myVecOfAllInfo)) %in% c(myTime)])))
        df <- rbind(df, data.frame(Isotopologue=isotopologueName, Time=myTime, Mean=myTimeMean, stdDev=myTimeSD))
      }
    }

    calculate_total_isotopes <- function(factor_value) {
      # Extract numeric values from the factor string and sum them
      sum(as.numeric(unlist(regmatches(factor_value, gregexpr("[0-9]+", factor_value)))))
    }
    
    df$Isotopologue = as.factor(df$Isotopologue)

    # reorder Isotopologue to avoid plots showing M0,M1,M10,... (we want M0,M1,M2,...)
    if(grepl('M', df$Isotopologue[1])){
      df$Isotopologue = factor(df$Isotopologue, levels = levels(df$Isotopologue)[order(as.numeric(gsub("M", "", levels(df$Isotopologue))))])
    }else{
      df$Isotopologue = factor(df$Isotopologue, levels = levels(df$Isotopologue)[order(sapply(levels(df$Isotopologue), calculate_total_isotopes))])
    }

    # just M0
    smallDf <- subset(df, Isotopologue == splitIsotopologue)
    # everything except M0
    largerDF <- subset(df, Isotopologue != splitIsotopologue)

    # store all of these in the list of images
    # we will have to plot by default the entire range
    # as well as using the limits provided by the user if they want to use the option
    if(nrow(largerDF) > 0)
    {
      if(length(yLimit) > 1)
      {
        p[[i]] <- ggarrange(qplot(Time, Mean, data=smallDf, colour=Isotopologue, geom=c("line","point")) +
                             geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title) + ylim(yLimit) + ylab(axisTitle) ,qplot(Time, Mean, data=largerDF, colour=Isotopologue, geom=c("line","point")) +
                             geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev, width = .3)) + ggtitle(title) + ylim(yLimit) + ylab(axisTitle))
      }


      if(length(yLimit) == 1)
      {
        p[[i]] <- ggarrange(qplot(Time, Mean, data=smallDf, colour=Isotopologue, geom=c("line","point")) +
                             geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title) +  ylab(axisTitle),qplot(Time, Mean, data=largerDF, colour=Isotopologue, geom=c("line","point")) +
                             geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev, width = .3)) + ggtitle(title) + ylab(axisTitle))
      }
    }
    if(nrow(largerDF) == 0)
    {

      if(length(yLimit) == 1)
      {
        p[[i]] <- qplot(Time, Mean, data=df, colour=Isotopologue, geom=c("line","point")) + geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title) +  ylab(axisTitle)
      }

      if(length(yLimit) > 1)
      {
        p[[i]] <- qplot(Time, Mean, data=df, colour=Isotopologue, geom=c("line","point")) + geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title) + ylim(yLimit) + ylab(axisTitle)
      }

    }


    if(length(yLimit) == 1)
    {
      pAll[[i]] <- qplot(Time, Mean, data=df, colour=Isotopologue, geom=c("line","point")) + geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title) + ylab(axisTitle)

    }

    if(length(yLimit) > 1)
    {
      pAll[[i]] <- qplot(Time, Mean, data=df, colour=Isotopologue, geom=c("line","point")) + geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title) + ylab(axisTitle) + ylim(yLimit)
    }

  }


  # print out the side by side plots
  Split_MIDs <- marrangeGrob(p, nrow=2, ncol=1)

  splitMIDSname <- paste(Category, outputName, "split_MIDs.pdf",sep = "_")
  MIDSnonSplitname <- paste(Category, outputName, "MIDs.pdf",sep = "_")
  ggsave(splitMIDSname, Split_MIDs, width = 11, height = 11)
  # print out the all by all plots
  MIDs <- marrangeGrob(pAll, nrow=3, ncol=3)
  ggsave(MIDSnonSplitname, MIDs, width = 12, height = 12)

  return(list("Split_MID"=Split_MIDs, "MID"=MIDs, "p_split"=p, "p_all"=pAll,
              "split_filename"=splitMIDSname, "all_filename"=MIDSnonSplitname))
}
