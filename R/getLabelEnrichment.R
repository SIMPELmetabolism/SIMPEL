#' @title Calculate label enrichment
#'
#' @description  This function is going to calculate the label enrichment and is called in \link{get_table_objects_NA_corrected}
#'
#' @param MIDs_table The MIDs table
#' @param forMedianNormalization A dataframe used for calculation
#' @param columnsToSearch All columns that contain experiment data
#' @param subsetOfTableJustAnnotationData A subset of data with annotation information columns
#' @param compounds_data The compounds data
#'
#' @export

getLabelEnrichment = function(MIDs_table,forMedianNormalization,columnsToSearch,subsetOfTableJustAnnotationData,compounds_data)
{

  MIDs_table = MIDs_table
  forMedianNormalization = forMedianNormalization
  columnsToSearch = columnsToSearch
  subsetOfTableJustAnnotationData = subsetOfTableJustAnnotationData
  compounds_data = compounds_data

  #declare the label enrichment
  names <- c("Bin",columnsToSearch)
  label_enrichment <- data.frame(matrix(ncol = length(names), nrow = length(unique(forMedianNormalization$Bin))))
  x <-  names
  colnames(label_enrichment) <- x
  label_enrichment$Bin = unique(forMedianNormalization$Bin)

  label_enrichment$Compound = unique(forMedianNormalization$Compound)



  #remove the M0's from the MID tables
  MIDs_tableNoM0 = MIDs_table[MIDs_table$total_isotopes > 0,]


  #for each bin, sum up all of the MIDs
  for(column in 1:length(columnsToSearch))
  {
    columnSearch = columnsToSearch[column]

    for(noMIDBins in 1:length(label_enrichment$Bin))
    {
      myBinSearch = label_enrichment$Bin[noMIDBins]
      mySubsetToSum = MIDs_tableNoM0[MIDs_tableNoM0$Bin == myBinSearch, colnames(MIDs_tableNoM0) == columnSearch]

      MIDsSumToAdd = sum(mySubsetToSum)
      label_enrichment[label_enrichment$Bin == myBinSearch,colnames(label_enrichment) == columnSearch ] = MIDsSumToAdd
    }
  }


  #now, go back and add all of the annotation data to the average_labeling table
  #have a vector to store the formulas for each bin that we will add to the table
  formulaForColumn = vector()
  forColumnMz = vector()
  forColumnPolarity = vector()
  forColumnRT = vector()
  #mz, rt, polarity, Bin, Compound, Formula columns

  #look up the Formula associated with a bin
  for(lookUpBin in 1:length(unique(label_enrichment$Bin)))
  {
    binForFormula = unique(label_enrichment$Bin)[lookUpBin]
    correctFormula = unique(compounds_data[compounds_data$prefix == binForFormula,]$formula)
    #note that we want the first mass
    compMz = unique(subsetOfTableJustAnnotationData[subsetOfTableJustAnnotationData$Bin %in% binForFormula,]$mz)[1]
    compPolarity = as.character(unique(compounds_data[compounds_data$prefix == binForFormula,]$polarity))
    compRT = unique(compounds_data[compounds_data$prefix == binForFormula,]$rt)[1]


    if(is.na(correctFormula) == TRUE)
    {
      print(lookUpBin)
      print("gives an NA, unfortunately")
      print( binForFormula )
    }

    formulaForColumn = c(formulaForColumn, correctFormula)
    forColumnPolarity = c(forColumnPolarity, compPolarity)
    forColumnRT = c(forColumnRT, compRT)
    forColumnMz = c(forColumnMz, compMz )

  }

  label_enrichment$Formula = formulaForColumn
  label_enrichment$Compound = unique(forMedianNormalization$Compound)
  label_enrichment$polarity = forColumnPolarity
  label_enrichment$rt = forColumnRT
  label_enrichment$mz = forColumnMz

  return(label_enrichment)
}
