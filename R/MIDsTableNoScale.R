#' @title Calculate MIDs table and average labeling table
#'
#' @description This function is going to calculate the unscaled MIDs Table and is called in get_table_objects_NA_corrected()
#'
#' @param XCMS_data The dataframe after searching for isotopologues
#' @param label_scheme The labeling scheme of the experiments
#'
#' @return The function returns a list containing: an MIDs table, an average labeling table and a vector of
#' bins that do not have an M0
#'
#' @export

MIDsTableNoScale = function(XCMS_data, label_scheme)
{

  #get the xcms data as input
  MIDs_table = XCMS_data

  #get the bins
  allBins = MIDs_table$Bin
  allBins = unique(allBins)

  #have a vector to keep track of all
  #the compounds without an MID.
  #We don't use it here, but we may in the future
  vecOfNoMIDs = vector()

  #calculate the average labeling description
  average_labeling = XCMS_data

  #we're going to have a column by which we bin
  #and we don't need a lot of the information in here
  average_labeling = subset(average_labeling,  select = -c(carbon, nitrogen, hydrogen, total_isotopes, comp_result))

  #only one average per bin
  average_labeling = average_labeling[!duplicated(average_labeling$Bin),]

  #subset for each bin and creae the MIDs table
  for(i in 1:length(allBins))
  {

    myBin = allBins[i]

    #focus on just the xcms data associated with a bin
    #we will subsequently subset this by category (treatments/mutants) and then replicates
    AllMIDSubsBeforeCategories = subset(MIDs_table,  Bin == myBin)

    #determine the formula of whatever is in the bin
    myFormula = unique(AllMIDSubsBeforeCategories$Formula)

    #get the carbons and nitrogens in the formula
    carbonsFormula = get_element_count( myFormula )[['C']]
    nitrogenFormula = get_element_count( myFormula )[['N']]
    hydrogenFormula = get_element_count( myFormula )[['H']]

    #if we don't have any C or N's make sure to set the value to zero instead
    #of NULL, which will cause problems down the road
    if(is.null(carbonsFormula))
    {
      carbonsFormula = 0
    }


    if(is.null(nitrogenFormula))
    {
      nitrogenFormula = 0
    }


    if(is.null(hydrogenFormula))
    {
      hydrogenFormula = 0
    }

    #get all of the isotopologues
    isotopologueList = AllMIDSubsBeforeCategories$comp_result

    # To determine which entry is M0 based on current labeling scheme
    elementInLabel = c("C","N","H") %in% unlist(strsplit(label_scheme, ""))
    whichCols = c("carbon","nitrogen","hydrogen")[elementInLabel]
    whichCols = whichCols[whichCols!=0]
    # get the M0 for the bin
    M0 = subset(AllMIDSubsBeforeCategories, apply(AllMIDSubsBeforeCategories[,whichCols, drop=F], 1, sum)==0,
                select = -c(Bin, Formula, mz, rt, carbon, nitrogen, hydrogen, total_isotopes, comp_result, polarity))


    #get the total number of nitrogens
    nitrogens = AllMIDSubsBeforeCategories$nitrogen
    #get the total number of carbons
    carbons = AllMIDSubsBeforeCategories$carbon
    #get the total number of hydrogens
    hydrogens = AllMIDSubsBeforeCategories$hydrogen

    #now we've got everything that we need
    AllMIDSubsBeforeCategories = subset(MIDs_table,  Bin == myBin, select = -c(Bin, Formula, mz, rt, carbon, nitrogen, hydrogen, comp_result, polarity) )

    #if there's no M0 keep track of the bin so that we can remove it.  If there is no M0, the bin is likely noise so keep track of the bin so we can remove it
    if(nrow(M0) == 0)
    {
      vecOfNoMIDs = c(vecOfNoMIDs, myBin)
    }

    #if we have mutiple MIDs, we're gonna start binning them
    if(nrow(M0) > 0)
    {

      #get the replicates

      #see if we have replicates for multiple categories
      reps = unique(data_cleanIII(colnames(AllMIDSubsBeforeCategories)))

      #get these replicates
      reps = unique(reps)

      reps = reps[is.na(reps) == FALSE]

      #and get the categories as well, which may be treatment coniditions
      categories = unique(data_cleanII(colnames(AllMIDSubsBeforeCategories)[colnames(AllMIDSubsBeforeCategories) %like% "_.+_"]))
      categories = unique(categories)
      categories = categories[is.na(categories) == FALSE]

      #subset each row of MIDs by its replicate
      for(rep in 1:length(reps))
      {
        #in each category
        for(category in 1:length(categories))
        {


          justReps = colnames(AllMIDSubsBeforeCategories)[colnames(AllMIDSubsBeforeCategories) %like% paste0("_",reps[rep], sep = NULL)]

          #subset by the caetegorical variable (i.e. WT or MUT, treatment1 or treatment 2) and their associated reps
          repsAndCategories = justReps[justReps %like% categories[category]]


          whichRep = reps[rep]
          whichCategory = categories[category]
          AllMIDSubs = as.data.frame(subset(AllMIDSubsBeforeCategories, select = repsAndCategories))

          #this will take us across each timepoint for the given
          #category and replicate so that we can track how the MIDs change across time

          for(j in 1:ncol(AllMIDSubs))
          {
            myToGetPercent = AllMIDSubs[,j]
            justToSum = colSums(AllMIDSubs, na.rm = TRUE)[j]

            justPercent =  (myToGetPercent / justToSum) * 100

            #effectively "refill" the table
            #and also calculate the avgs information

            #store all the avgs information

            #will store all of the information at a single timepoint for all the MIDs
            vectorOfAvgsInfo = vector()
            allMIDsNumber = sum(AllMIDSubsBeforeCategories$total_isotopes)

            #this will take us across all the isotoplogues for a timepoint
            for(k in 1:length(justPercent))
            {
              nameOfColumn = colnames(AllMIDSubs)[j]
              toReplace = justPercent[k]

              myIsoBin = isotopologueList[k]

              #fill in all the elements in the jth row
              #by each of their columns
              MIDs_table[MIDs_table$comp_result == myIsoBin, colnames(MIDs_table) == nameOfColumn] = toReplace

              elementCount = list(carbons, nitrogens, hydrogens)
              numOfLabeled = sum(sapply(elementCount[elementInLabel], "[[", k))
              AvgsInfo = toReplace * numOfLabeled

              #sum up all of the isotopologues to create the mass distribution vector
              vectorOfAvgsInfo = c(vectorOfAvgsInfo, AvgsInfo)

            }

            #in order to calculate average labeling, first multiply each isotologue by the
            #number of labeled carbons and nitrogens in that isotopologue, then sum all the isotopologues
            #and divide the sum by total number of carbons and nitrogens in the compound

            CNHnumbers = c(carbonsFormula, nitrogenFormula, hydrogenFormula)
            CNHnumbers = sum(CNHnumbers[elementInLabel])

            #collapse all the isotologue info into one datapoint
            vectorOfAvgs = sum(vectorOfAvgsInfo) / CNHnumbers

            if(sum(CNHnumbers) == 0)
            {
              vectorOfAvgs = 0
            }

            #repopulate the average_labeling table with all of the updated data
            average_labeling[average_labeling$Bin == myBin, colnames(average_labeling) == nameOfColumn] = vectorOfAvgs

          }
        }
      }
    }
  }


  #initialize the isotopologue column
  MIDs_table$Isotopologue = MIDs_table$total_isotopes
  MIDs_table = as.data.frame(MIDs_table)
  #for each MID, go through the MID table and fill in the number of labeled carbon and nitrogen
  elementInLabel = c("C","N","H") %in% unlist(strsplit(label_scheme, ""))
  for(i in 1:nrow(MIDs_table))
  {
    carbonNum = MIDs_table[i,]$carbon
    nitrogenNum = MIDs_table[i,]$nitrogen
    hydrogenNum = MIDs_table[i,]$hydrogen
    elements = c("C", "N", "H")
    numbers = c(carbonNum, nitrogenNum, hydrogenNum)
    nameToAdd = paste0(paste0(elements[elementInLabel], numbers[elementInLabel]), collapse = "")
    MIDs_table[i, colnames(MIDs_table) == "Isotopologue"] = nameToAdd
  }


  listOfAllOutputs = list()
  listOfAllOutputs[[1]]= MIDs_table
  listOfAllOutputs[[2]] = average_labeling
  listOfAllOutputs[[3]] = vecOfNoMIDs
  # hold the MIDs_table before we adjust by proxy pool size
  return(listOfAllOutputs)
}
