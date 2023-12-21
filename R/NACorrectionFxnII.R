#' @title Perform natural abundance correction
#'
#' @description This function is going to do the NA correction and is called in get_table_objects_NA_corrected()
#'
#' @param MIDs_table The MIDs table calculated
#' @param label_scheme The labeling scheme of the experiments
#' @param output The output directory
#'
#' @return The function returns a dataframe of a corrected MIDs table
#'
#' @export

NACorrectionFxnII = function(MIDs_table, label_scheme, output="_")
{
  #we are going to make sure to place the isocorrector
  #data in the predetermined output file
  output = output

  MIDs_table = MIDs_table
  theMIDSFormat = MIDs_table

  #filepaths
  path.molecule = "MoleculeFile.csv"
  path.element = system.file("extdata","ElementFile_with_extra_elements.csv",package="SIMPEL")
  path.measurement = "MeasurementFile.csv"

  #new function
  get_element_count <- function(chemical_formula){
    chemical_formula %>%
      makeup() %>%
      as.list() %>%
      return()
  }


  #get the MIDs table that we are going to modify through the NA correction

  #remove any of the special characters that could cause an issue
  #get the MIDs table that we are going to modify through the NA correction

  #this is the MIDs template we'll be filling in
  theMIDSFormatFillin = theMIDSFormat

  #we've removed any non alphanumeric characters from the MIDs column
  theMIDSFormat$Compound = str_replace_all( theMIDSFormat$Compound, "[^[:alnum:]]", "-")

  #we'll have a column in the MIDs table we're going to fill in
  #with the cleaned column information
  theMIDSFormatFillin$CompoundPlaceholder = theMIDSFormat$Bin


  #we're going to need to make slight changes to the formulas
  #We won't be able to have chemicals of the formula CXNY
  #where C or Y  = 1 and the formula is instead provided as CN and not
  #C1N1
  #we'll store the adjusted formulas in this vector
  adjustedFormulas = vector()

  #reformat the chemical formulas to include the labeled atoms
  #get the formulas from the MIDs table that we are going to correct
  formulasToFix = theMIDSFormat$Formula

  #for each formula, we're going to recreate the formula
  #from the elements
  for(name in 1:length(formulasToFix))
  {
    #get a formula
    formulaUnique = formulasToFix[name]

    #get the name of each element, for example C, N
    formulaToFixElements = names(unlist(get_element_count(formulaUnique)))

    #get the number of each element
    elementNumber = unlist(get_element_count(formulaUnique))
    #we don't need the "names" attribute of the vector
    names(elementNumber) = NULL

    #Now, with elementNumber and elementName
    #we have everything we need in order to reassemble the formula correctly
    vecOfInfoToBeFormula = vector()
    elementAndCount = vector()

    #for each element in the formula:
    for(element in 1:length(formulaToFixElements))
    {
      #get the element name
      whichOne = formulaToFixElements[element]

      #get the element number
      howMany = elementNumber[element]

      #make a vector that contains the element and its quantity
      elementAndCount = c(whichOne,howMany)

      #add this to the vector that will contain all of the Element,Number pairs
      #for the formula. We're going to join this into a string
      #once we have this information for each element in the vector
      vecOfInfoToBeFormula = c(vecOfInfoToBeFormula, elementAndCount)
    }

    #all of the element, quantity information into a single string
    #to get the adjust formula
    adjustedFormula = paste(vecOfInfoToBeFormula, collapse = "")

    #store all of the adjusted formulas in a vector
    adjustedFormulas = c(adjustedFormulas, adjustedFormula)
  }

  #now, we have reformatted every formula to include information
  #about the number of elements lableed, add this column to the MIDs table
  #we want to work with for NA correction
  theMIDSFormat$Formula = adjustedFormulas

  ##In the case that a formula contained only on atom of on element
  #the formula will now contain a 1 following the element, as opposed
  #to no number indicating the number of elements in the formula


  ###Now that we have made sure all of the formulas are correct
  #for each MID, we will append the number of carbons and nitrogens that define the species
  #to be consistent with the formatting that IsocorrectoR requires


  #if there are X labeled carbons and Y labeled Nitrogens.  We'll have to do this for EVERY species
  #from C0N0 to CXNY.  Note, that there are some important issues that I found:
  #If only carbons are labeled, the M0 format has to be C0, not C0N0 and if no nitrogens are labled
  #M0 is going to have to have to be N0

  #in order to iterate through each labeled species, we'll need to get all of the MIDs together
  #so we add a placeholder column that we can sort the table on
  theMIDSFormat$CompoundPlaceholder = theMIDSFormat$Bin


  #we're going to make a reformatted MIDs table that will hold all the information
  #of the NA corrected values
  theMIDSReformat =  data.frame(matrix(nrow = 0, ncol = length(colnames(theMIDSFormat))))
  colnames(theMIDSReformat) = colnames(theMIDSFormat)
  newColumnFormulas = vector()

  #We're going to take our unreformmated MIDs table and sort it by our placeholder column (the metabolite name)
  theMIDSFormat = theMIDSFormat[order(theMIDSFormat$CompoundPlaceholder),]

  #then, we're going to go through all the MIDs for each metabolite
  #and put them in the proper format that IsoCorrectoR requires
  #making sure that we respect some of its formatting quirks which are described below
  MIDcompoundNameLast = ""

  #We will need to know if we are on the first line, because this will have to be the
  #M0 of a species.  This variable will keep track of that
  howManyTimesThrough = 1


  #for each row of an MID in the MID table
  for(MIDrow in 1:nrow(theMIDSFormat))
  {
    #extract the line with the MID's data
    myMIDLine = theMIDSFormat[MIDrow,]

    #will be of the form: name_CX.NY where X is the number of carbons and Y is the number
    #of nitrogens

    #get the N number from the line
    labeledNitrogensInMID = myMIDLine$nitrogen
    #get the C Number from the line
    labeledCarbonsInMID = myMIDLine$carbon
    #get the H number from the line
    labeledHydrogensInMID = myMIDLine$hydrogen

    #make sure that we have the parent compound name of the MID
    MIDcompoundName = myMIDLine$CompoundPlaceholder

    #get ALL of the information about the MID
    allMIDInfo = theMIDSFormat[theMIDSFormat$CompoundPlaceholder == MIDcompoundName,]

    myMIDTotCarbon = max(allMIDInfo$carbon)
    myMIDTotNitrogen = max(allMIDInfo$nitrogen)
    myMIDTotHydrogen = max(allMIDInfo$hydrogen)

    # first filter by labeling scheme
    elementInLabel = c("C", "N", "H") %in% unlist(strsplit(label_scheme, ""))
    elements = c("C", "N", "H")[elementInLabel]
    myMIDTot = c(myMIDTotCarbon, myMIDTotNitrogen, myMIDTotHydrogen)[elementInLabel]
    labeledElementsInMID = c(labeledCarbonsInMID, labeledNitrogensInMID, labeledHydrogensInMID)[elementInLabel]
    # then remove elements with 0 counts to create names
    if (sum(myMIDTot)==0){ # if no elements were labeled, we still need to assign some names
      reformattedName = paste0(MIDcompoundName, "_", elements[1], '0')
    }else{
      elements = elements[myMIDTot>0]
      #myMIDTot = myMIDTot[myMIDTot>0]
      labeledElementsInMID = labeledElementsInMID[myMIDTot>0]
      reformattedName = paste(MIDcompoundName, "_",
                              paste0(elements, labeledElementsInMID, collapse = "."), sep="")
    }


    #add to the MIDs table the reformatted name
    theMIDSFormat[MIDrow,]$Compound = reformattedName

    #we're going to keep track of the compound associated for both the current entry we're formatting and
    #the compound name associated with the last entry.  The reason for this is that if the current compound
    #name of the MID we're working on is different than the name of the last MID we worked on
    #then it means we just finished reformatting all of the MIDs of a metabolite
    #and we'll need to go back and check to see if we skipped any of the MIDs of that compound

    #the way we'll do this is create a dummy table with all 0's for every MID from C0N0 to CXNY
    #and see which one of those MIDs exist in the real data.  If any of them are skipped in the real
    #data, we'll take the row of zero values for that MID and included it in the reformatted table
    #since we can pass 0's to isocorrectoR, but there needs to be some measurement value

    #determine if we've moved on from one compound to another

    #if so, create the dummy MIDs - say that we get to MID10
    #create a set of dummy rows up to MID10, but discard those for which we already have MIDs
    #and any that we've skipped, add these zero values to the reformatted table

    #we've moved from one compound to another
    #determine if the compound we had just looked at
    #was missing any MIDs

    #determine if the metabolite we're on is the last one in the table
    locationMetab =


      if(((MIDcompoundNameLast != MIDcompoundName) & (MIDrow > 1)) || (howManyTimesThrough == nrow(theMIDSFormat)))
      {
        #create the subset of all the MIDs for the previous metabolite
        #we're going to determine if it is missing any MIDs and then replace them with zeroes
        mySubsetMIDs = theMIDSFormat[theMIDSFormat$CompoundPlaceholder == MIDcompoundNameLast,]


        #create the dummy table
        #create a dummy MID os zeroes for everything from M0 to MCXNY, where X is the maximum number of labeled
        #carbons an 'Y' is the maximum number of labeled nitrogens
        namesForDummyRows = vector()


        #include the baseline by default, we'll drop if if it's there
        #namesForDummyRows = c(namesForDummyRows,paste(MIDcompoundNameLast, "_", "C", 0, ".","N",0, sep = ""))


        #create a vector of 1:Y where Y is the number of labeled Nitrogens
        if(max(mySubsetMIDs$nitrogen) > 0)
        {
          allNitrogenLabels = c(1:max(mySubsetMIDs$nitrogen))
          allNitrogenLabels = c(allNitrogenLabels,0)
          allNitrogenLabels = unique(allNitrogenLabels)
        }

        if(max(mySubsetMIDs$nitrogen) == 0)
        {
          allNitrogenLabels = c(0)
        }

        #create a vector of 1:X where X is the number of labeled carbons
        if(max(mySubsetMIDs$carbon) > 0)
        {
          allCarbonLabels = c(1:max(mySubsetMIDs$carbon))
          allCarbonLabels = c(allCarbonLabels,0)
          allCarbonLabels = unique(allCarbonLabels)
        }

        if(max(mySubsetMIDs$carbon) == 0)
        {
          allCarbonLabels = c(0)
        }

        #create a vector of 1:Z where Z is the number of labeled Hydrogens
        if(max(mySubsetMIDs$hydrogen) > 0)
        {
          allHydrogenLabels = c(1:max(mySubsetMIDs$hydrogen))
          allHydrogenLabels = c(allHydrogenLabels,0)
          allHydrogenLabels = unique(allHydrogenLabels)
        }

        if(max(mySubsetMIDs$hydrogen) == 0)
        {
          allHydrogenLabels = c(0)
        }

        # create labels for the dummy table for any possible labeling schemes
        elementInLabel = c("C", "N", "H") %in% unlist(strsplit(label_scheme, ""))
        elements = c("C", "N", "H")[elementInLabel]
        numOfElements = list(allCarbonLabels, allNitrogenLabels, allHydrogenLabels)[elementInLabel]
        maxElementLabels = c(max(allCarbonLabels), max(allNitrogenLabels), max(allHydrogenLabels))[elementInLabel]
        # then remove elements with 0 counts to create names
        if (sum(maxElementLabels==0)){ # if no elements were labeled, we still need to assign some names
          namesForDummyRows = paste0(MIDcompoundNameLast, "_", elements[1], '0')
        }else{
          elements = elements[maxElementLabels>0]
          numOfElements = numOfElements[maxElementLabels>0]
          tempNames = lapply(c(1:length(elements)),
                             function(x) paste0(elements[x],unlist(numOfElements[x])))
          namesForDummyRows = paste0(MIDcompoundNameLast, "_",
                                     apply(expand.grid(tempNames) %>%
                                             arrange(Var1), 1, paste0, collapse="."))
        }


        namesForDummyRows = unique(namesForDummyRows)
        #Now that we have all of the MID species, create the dummy table of all zero values
        dummyTable = data.frame(matrix(nrow = length(namesForDummyRows), ncol = length(colnames(theMIDSFormat))))
        colnames(dummyTable) = colnames(theMIDSFormat)
        occursInFiles = 0

        dummyTable[is.na(dummyTable)] = 0
        dummyTable$Compound = namesForDummyRows

        #get the base formula for the compound whose MIDs we are adding
        formulaNameLast =  unique(mySubsetMIDs$Formula)
        #dummyTable$Formula = rep(formulaNameLast,length(dummyTable$Formula))
        dummyTable$Formula = rep(formulaNameLast,length(namesForDummyRows))


        #now that we have the dummy table we are going to want to keep only those rows that are not
        #present in the actual data, as these are the only MIDs that are skipped
        dummyTable = dummyTable[!(dummyTable$Compound %in% mySubsetMIDs$Compound),]

        #add to the dummy table which has data for the 'skipped' MIDs
        #the data for the present MIDs.  Now, the dummy table will contain all of the MIDs from
        #C0N0 to CXNY, even those that were 'skipped' in the data, as they have now been replaced by zeroes
        dummyTable = rbind(dummyTable,mySubsetMIDs)
        #now that we have made all the dummy rows, we can attach these dummy rose to the full table

        #take this table with all the MIDs from C0N0 to CXNY and bind it to the reformatted MIDs table
        theMIDSReformat = rbind(theMIDSReformat, dummyTable )

        #reset the modified MIDs vector
        newColumnFormulas = vector()
      }


    #update the name of the metabolite whose MID's we've just worked with
    MIDcompoundNameLast = MIDcompoundName
    howManyTimesThrough = howManyTimesThrough + 1
  }


  #we're going to call the reformatted MIDs table the "measurementFile" now
  #since that is the convention that isocorrectoR uses
  measureMentFile = theMIDSReformat


  #Now, we're just going to need to make sure that column names of the measurement table
  #are correct for isoCorrectoR
  newColumnCompounds = theMIDSReformat$Compound
  measureMentFile = cbind(newColumnCompounds,measureMentFile)
  colnames(measureMentFile)[1] = "Measurements/Samples"

  #we're going to have to take out the extra columns from the MIDs table for isocorrectoR to work
  #properly
  measureMentFile$CompoundPlaceholder = NULL
  measureMentFile$mz = NULL
  measureMentFile$rt = NULL
  measureMentFile$polarity = NULL
  measureMentFile$Isotopologue = NULL
  measureMentFile$Bin = NULL
  measureMentFile$total_isotopes = NULL
  measureMentFile$carbon = NULL
  measureMentFile$nitrogen = NULL
  measureMentFile$hydrogen = NULL
  measureMentFile$comp_result = NULL
  measureMentFile$Formula = NULL
  measureMentFile$Compound = NULL

  #remove the duplicates
  measureMentFile <-  measureMentFile[order(measureMentFile$"Measurements/Samples"), ]
  measureMentFile = measureMentFile[ !duplicated(measureMentFile$"Measurements/Samples"), ]

  #now that the extra columns are removed, we're going to write the reformatted MIDs
  #table as the 'measureMentFile'
  write.table(measureMentFile, file = "MeasurementFile.csv", sep = ",", row.names = FALSE)
  measureMentFilePath = "MeasurementFile.csv"

  #Now, the last thing to do is to create the moleculeFile for isocorrectoR
  #which will have a line for each metabolite describing the number of labeled C's and N's
  #This format is going to be, given a compound named test1 with formula CXNY, test1_LabCXLabNY


  #get the names of each metabolite
  uniqueCompoundFormulas =  unique(theMIDSFormat$CompoundPlaceholder)
  uniqueCompoundFormulas = uniqueCompoundFormulas[uniqueCompoundFormulas != 0]

  #store the labeled chemical formulas for the moleculeFile in a vector
  vecOfFormulas = vector()
  for(formulaIndex in 1:length(uniqueCompoundFormulas))
  {

    #get the metabolite name
    formulaUnique =  unique(theMIDSFormat[theMIDSFormat$CompoundPlaceholder == uniqueCompoundFormulas[formulaIndex],]$Formula)


    #get numberofLabeled carbons
    labeledC = max(theMIDSFormat[theMIDSFormat$CompoundPlaceholder == uniqueCompoundFormulas[formulaIndex],]$carbon)
    labeledN = max(theMIDSFormat[theMIDSFormat$CompoundPlaceholder == uniqueCompoundFormulas[formulaIndex],]$nitrogen)
    labeledH = max(theMIDSFormat[theMIDSFormat$CompoundPlaceholder == uniqueCompoundFormulas[formulaIndex],]$hydrogen)


    elementInLabel = c("C", "N", "H") %in% unlist(strsplit(label_scheme, ""))
    elements = c("LabC", "LabN", "LabH")[elementInLabel]
    labeledElements = c(labeledC, labeledN, labeledH)[elementInLabel]
    elements = elements[labeledElements > 0]
    labeledElements = labeledElements[labeledElements > 0]
    totalFormulaInfo = paste(formulaUnique,
                             paste0(elements, labeledElements, collapse = ""), sep="")


    #add the metabolite names which are now reformatted to be acceptable for isocorrectoR's
    #moleculeFile format into a vector
    vecOfFormulas = c(vecOfFormulas, totalFormulaInfo)
  }

  #create the molecule table now
  #by combining the column of metabolite names, the labeled formulas and the additional "" in the third column
  MoleculeTable = cbind(uniqueCompoundFormulas, vecOfFormulas, rep("", length(vecOfFormulas)))


  #make this a into a dataframe
  MoleculeTable = as.data.frame(MoleculeTable)

  #give it the proper columnm names
  colnames(MoleculeTable) = c("Molecule", "ms ion or ms/ms product ion", "ms/ms neutral loss")

  #write the molculeTable
  write.table(MoleculeTable, file = "MoleculeFile.csv", sep = ",", row.names = FALSE)


  # 2) run correction algorithm and save results in variable
  correctionResults <- IsoCorrectoR::IsoCorrection(MeasurementFile=path.measurement,
                                                   ElementFile=path.element,
                                                   MoleculeFile=path.molecule,UltraHighRes=TRUE, DirOut = output)

  #now, we're going to need to recreate the adjusted table by filling in the MIDs table
  NAcorrected = correctionResults$results$CorrectedFractions


  for(i in 1:nrow(NAcorrected))
  {
    myNAInfo = NAcorrected[i,]

    #get the compound
    myNAName = unlist(strsplit(rownames(myNAInfo), "_"))[1]
    myMIDNames =  unlist(strsplit(rownames(myNAInfo), "_"))[2]

    elementInLabel = c("C", "N", "H") %in% unlist(strsplit(label_scheme, ""))
    numOfElements = c(0, 0, 0) # initial value

    if(sum(elementInLabel)==1){
      MIDInfo = myMIDNames
    }else{
      MIDInfo = c(data_cleanDot(myMIDNames), data_cleanDotII(myMIDNames))
      MIDInfo[is.na(MIDInfo)] = 0
    }
    MIDInfo = gsub("[^0-9.-]", "", MIDInfo)
    numOfElements[which(elementInLabel)] = as.numeric(MIDInfo)


    for(j in 1:length(colnames(NAcorrected)))
    {

      sampleName = colnames(NAcorrected)[j]
      toReplace = myNAInfo[names(myNAInfo) == sampleName]

      if(is.na(toReplace)  == FALSE)
      {
        toReplace =toReplace * 100
      }


      if(is.na(toReplace) == TRUE)
      {
        toReplace =toReplace

      }

      theMIDSFormatFillin[theMIDSFormatFillin$CompoundPlaceholder == myNAName & theMIDSFormatFillin$nitrogen == numOfElements[2] & theMIDSFormatFillin$carbon == numOfElements[1] & theMIDSFormatFillin$hydrogen == numOfElements[3],colnames(theMIDSFormatFillin) == sampleName ] = toReplace
    }

    #get just the compound name
  }

  return(theMIDSFormatFillin)
}
