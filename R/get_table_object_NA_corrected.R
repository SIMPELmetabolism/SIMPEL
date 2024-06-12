#' @title Identifying single and dual labeled isotopologues from XCMS lists for chemical formulas provided.
#' Create natural abundance corrected post-processing data tables for pre-processed, stable isotope labeling based metabolomics datasets
#'
#' @description Function to process XCMS data for stable isotope enriched metabolomics experiments. This function first calculates all the
#' isotopologue m/z for a list of chemical formulae provided by the user as "compound_data" and identifies the matched m/z features within
#' the "xcms_data". In addition, natural abundance correction is achieved using the package IsoCorrectoR.
#' Four objects are created using this function, an MIDs table containing mass isotopologue distribution matrix for all the
#' compounds, a scaled_MIDs table where the MIDs are adjusted by pool sizes provided in the compounds data to account for differences in pool sizes of compounds,
#' an average_labeling table where the label enrichment for each compound is represented as a percentage and finally a mol_equivalent_labeling
#' table where the label enrichment is represented as mol equivalents of labeled compound at each time point. Four additional objects similar to the above
#' but corrected for natural abundance are also created.
#'
#' @param XCMS_data XCMS pre-processed data provided as user input in .csv format
#' @param ppm numeric (1) define the maximal tolerated deviation in m/z (as parts per million) from calculate labeled isotopologue within XCMS_data
#' @param rt_tolerance  numeric (1) the maximal tolerate retention time deviation in minutes (from specified value provided in the annotation file) when searching for labeled isotopologues within XCMS_data.
#' @param compounds_data User annotation file (in .txt format) containing a list of compounds for which isotopologues are to be identified.
#' @param output A title to append to the files if the user wishes to print out the analyses tables/objects.
#' @param label_scheme string (1) The labeling scheme of the experiments. It currently supports single and dual labeled settings that involve
#' 13C, 15N, and 2H. All available options are: "C", "N", "H", "CN", "CH", "NH".
#'
#' @seealso \link{MIDplot}, \link{label_enrichment_plot}, \link{PCA_and_heatmap}, \link{getClustersAndPlots}
#'
#' @return The function returns a class instance with eight objects 1. MIDs 2. scaled_MIDs 3. average_labeling 4. mol_equivalent_labeling as well as natural
#' abundance corrected versions of 1-4.
#' In addition, the function also exports the eight objects as .txt files containing the specific data tables for user verification into two
#' separate directories.
#'
#' @import stats
#' @import utils
#' @import data.table
#' @import purrr
#' @import stringr
#' @import ggfortify
#' @import ggplot2
#' @import ggpubr
#' @import gridExtra
#' @import grid
#' @import lattice
#' @import RColorBrewer
#' @import grDevices
#' @import tidyverse
#' @import dplyr
#' @import IsoCorrectoR
#'
#' @export
#'
#' @examples
#' olddir <- setwd(tempdir())
#' \donttest{try(test_13C15NGlutamine <- get_table_objects_NA_corrected(XCMS_data, compounds_data,
#' ppm = 5, rt_tolerance = 0.1, label_scheme = "CN", output = "13C15N_Glutamine"))}
#' setwd(olddir)

get_table_objects_NA_corrected <- function(XCMS_data, compounds_data, ppm=10, rt_tolerance=.1, output="sample_output", label_scheme){

  output = output

  #use the working directory that the user has specified via getwd()
  setWorkingDir = getwd()

  #create the output directory for the user which will hold all their data
  dir.create(file.path(setWorkingDir), showWarnings = FALSE)

  #we are now going to be working in the directory we have created to hold all
  #of the outputs
  dir.create(file.path(setWorkingDir, output), showWarnings = FALSE)

  setwd(file.path(setWorkingDir, output))

  #make the directories for the not NA_adjusted table
  dir.create(file.path(".","get_table_objects"), showWarnings = FALSE)

  #make the directories for the NA corrected table
  dir.create(file.path(".", "NA_corrected"), showWarnings = FALSE)


  #set the variables from the input, i.e. ppm
  ppm = ppm
  rt_tolerance = rt_tolerance
  XCMS_data = XCMS_data
  compounds_data = compounds_data
  output = output
  label_scheme = label_scheme

  #add the columns to the table that we want all the outputs to have
  XCMS_data <- XCMS_data %>%
    mutate(comp_result = NA)

  XCMS_data <- XCMS_data %>%
    mutate(carbon = NA)

  XCMS_data <- XCMS_data %>%
    mutate(nitrogen = NA)

  XCMS_data <- XCMS_data %>%
    mutate(hydrogen = NA)

  XCMS_data <- XCMS_data %>%
    mutate(total_isotopes = NA)

  XCMS_data <- XCMS_data %>%
    mutate(Bin = NA)

  XCMS_data <- XCMS_data %>%
    mutate(Compound = NA)

  XCMS_data <- XCMS_data %>%
    mutate(Formula = NA)



  #looping through the list of mass features in the XCMS data
  #to identfiy isotopologs and bins
  for(i in 1:nrow(compounds_data))
  {
    #for each compound we're gonna pull out formula and retention time
    #from the annotation file
    compounds_data[i, ] %>%
      unlist() %>%
      as.vector() %>%
      print()

    ##for each compound in the annotation file with rt and formula
    #determine its m/z
    comp_formula <- as.character(compounds_data[i, "formula"])

    #retention time and polarity
    r_time  <- as.numeric(compounds_data[i, "rt"])
    polarity = compounds_data$polarity[i]

    #send that m/z to look up all of the C and N isotopes
    #and get the upper and lower bounds on each m/z based on the specified ppm error

    #filter the compounds_data by polarity
    subset_comp_lookup_table = compounds_data[compounds_data$polarity == polarity,]
    comp_lookup_table <- get_comp_mz_lookup(subset_comp_lookup_table, comp_formula, r_time, ppm, polarity, label_scheme)

    #information for the bin as found in the prefix column of the annotation
    #file
    myBin = as.character(compounds_data[i, "prefix"])
    myFormula = as.character(compounds_data[i, "formula"])
    myCompound = as.character(compounds_data[i, "Compound"])

    #looping through the xcms data

    #for each feature in the XCMS data, id'd by RT and m/z
    #determine if it is one of the C or N or a combination of isotopes by whether or not
    #it falls within the error in ppm that are stored in comp_lookup_table
    for(j in 1:nrow(XCMS_data)){
      if(is.na(XCMS_data[j, "comp_result"]))
      {
        #only make sure to search the correct polarity
        if(XCMS_data[j, "polarity"] == polarity)
        {
          #print(j)
          x = as.numeric(XCMS_data[j, "mz"])
          y = as.numeric(XCMS_data[j, "rt"])

          #determine if this feature, id'd by m/z and RT is one of the C and N isotopes
          #whose ppm error is stored in comp_lookup_table
          #print("in the second loop")
          #add the number of N's and the number of C's
          #add the isotopologue columns
          val <- get_comp_stage(x, y, comp_lookup_table, r_time, rt_tolerance)

          if(is.null(val) == FALSE)
          {

            XCMS_data[j, "comp_result"] <- val$compound
            XCMS_data[j, "carbon"] <- val$carbon
            XCMS_data[j, "nitrogen"] <- val$nitrogen
            XCMS_data[j, "hydrogen"] <- val$hydrogen

            isotopeNumbers = val$carbon$carbon + val$nitrogen$nitrogen + val$hydrogen$hydrogen

            XCMS_data[j, "total_isotopes"] <- isotopeNumbers
            XCMS_data[j, "Bin"] <- myBin
            XCMS_data[j, "Compound"] <- myCompound
            XCMS_data[j, "Formula"] <- myFormula

          }

          if(is.null(val) == TRUE)
          {
            XCMS_data[j, "comp_result"] <- NA
            XCMS_data[j, "carbon"] <- 0
            XCMS_data[j, "nitrogen"] <- 0
            XCMS_data[j, "hydrogen"] <- 0
            XCMS_data[j, "total_isotopes"] <- 0
            XCMS_data[j, "Bin"] <- myBin
          }

        }
      }
    }
  }


  #create a subset table that only contains binned rows with
  #labeled data
  XCMS_data = XCMS_data[is.na(XCMS_data$comp_result) == FALSE,]

  # #get the xcms data as input
  MIDs_table = XCMS_data

  #run a master function to get the MIDs table, the averages table and also the list of
  #bins with incomplete sets of MIDs
  listOfAllOutputs = MIDsTableNoScale(MIDs_table, label_scheme)#, XCMS_data)

  #This is going to get us the 1st object of the list which will be the unscaled MIDs table
  MIDs_table = listOfAllOutputs[[1]]
  if(nchar(MIDs_table$Isotopologue[1]) < 4){ # if single-labeled, use M0,M1,...
    MIDs_table$Isotopologue = paste0('M', MIDs_table$total_isotopes)
  }
  #This is going to get us the second object of the list, which will be the averages table
  average_labeling = listOfAllOutputs[[2]]
  #This is going to get us the third object of the list, which will be any of the MIDs without
  #sufficient species present
  vecOfNoMIDs = listOfAllOutputs[[3]]


  #scale the MIDs table
  categories = unique(data_cleanII(colnames(MIDs_table)[colnames(MIDs_table) %like% "_.+_"]))
  MIDs_tableScaled = MIDs_table %>%
    group_by(Bin) %>%
    left_join(. , dplyr::select(compounds_data, c('prefix', contains('poolsize'))), by=join_by(Bin==prefix)) %>%
    mutate(across(contains(paste0(categories[1], "_")), ~.*eval(parse(text=paste0("poolsize_", categories[1]))))) %>%
    mutate(across(contains(paste0(categories[2], "_")), ~.*eval(parse(text=paste0("poolsize_", categories[2]))))) %>%
    dplyr::select(-contains('poolsize')) %>%
    as.data.frame()


  #Next we calculate label enrichment which is the nmole equivalent of labeled compound at a given timepoint
  #Label enrichment is the sum of labeled isotopologues of each bin within the pool size scaled MID table
  #note that this does not include M0 [unlabeled pool is removed]
  label_enrichment = MIDs_tableScaled %>%
    dplyr::filter(total_isotopes != 0) %>%
    select(colnames(.)[colnames(.) %like% "_.+_"], Bin) %>%
    group_by(Bin) %>%
    summarise(across(everything(), sum)) %>%
    left_join(. , MIDs_tableScaled %>%
                dplyr::filter(total_isotopes == 0) %>%
                dplyr::select(mz, rt, polarity, Bin, Compound, Formula),
              by=join_by(Bin)) %>%
    as.data.frame()


  #remove bins that have no M0s
  MIDs_table = MIDs_table[!MIDs_table$Bin %in% vecOfNoMIDs,]
  MIDs_tableScaled = MIDs_tableScaled[!MIDs_tableScaled$Bin %in% vecOfNoMIDs,]
  average_labeling = average_labeling[!average_labeling$Bin %in% vecOfNoMIDs,]
  label_enrichment = label_enrichment[!label_enrichment$Bin %in% vecOfNoMIDs,]


  # To get subsets of clean bins (bins with small se)
  data_cols = colnames(MIDs_table)[colnames(MIDs_table) %like% "_.+_"]
  clean_subset = MIDs_table %>%
    dplyr::filter(total_isotopes == 0) %>%
    dplyr::select(all_of(data_cols), Bin) %>%
    reshape2::melt(measure.vars = data_cols, variable.name = "sample") %>%
    dplyr::mutate(sample = as.character(sample), 
                  time = as.numeric(sub(".","",sapply(strsplit(sample, '[_]' ), `[` , 1))), 
                  category = sapply(strsplit(sample, '[_]' ), `[` , 2)) %>%
    dplyr::group_by(Bin, time, category) %>%
    summarise_at(vars(value), list(se = ~sd(.)/sqrt(length(.)), t_crit = function(x) qt(0.975,df=length(x)-1))) %>%
    dplyr::mutate(se = coalesce(se, 0)) %>%
    ungroup() %>%
    dplyr::group_by(Bin, category) %>%
    dplyr::summarise(clean = sum(se > 7.5) < 2)
  clean_bins = lapply(unique(clean_subset$category), function(category) clean_subset$Bin[clean_subset$category==category & clean_subset$clean])
  names(clean_bins) = unique(clean_subset$category)
  clean_objects = unlist(lapply(unique(clean_subset$category), function(x) 
                  return(list(MIDs_table[MIDs_table$Bin %in% clean_bins[[x]],],
                         MIDs_tableScaled[MIDs_tableScaled$Bin %in% clean_bins[[x]],],
                         average_labeling[average_labeling$Bin %in% clean_bins[[x]],],
                         label_enrichment[label_enrichment$Bin %in% clean_bins[[x]],]))), recursive = FALSE)
  combs = expand.grid(c("MIDs_table_clean_", "MIDs_tableScaled_clean_", "average_labeling_clean_", "label_enrichment_clean_"), unique(clean_subset$category))
  names(clean_objects) = paste0(combs$Var1, combs$Var2)


  if (!is.null(output))
  {
    write.csv(MIDs_table, file = file.path("./get_table_objects", paste(output, "MIDs.csv", sep = "_")), row.names = FALSE)
    write.csv(MIDs_tableScaled, file = file.path("./get_table_objects", paste(output, "scaled_MIDs.csv", sep = "_")), row.names = FALSE)
    write.csv(average_labeling, file = file.path("./get_table_objects", paste(output, "average_labeling.csv", sep = "_")), row.names = FALSE)
    write.csv(label_enrichment, file = file.path("./get_table_objects", paste(output, "mole_equivalents_of_label.csv", sep = "_")), row.names = FALSE)
    lapply(names(clean_objects), function(name) {
      write.csv(clean_objects[[name]], paste0("get_table_objects/", name, ".csv"), row.names = FALSE)
    })
  }



  #NA correct the MIDs now!
  # leave compounds that only have unlabeled isotopologues out of the correction process
  which_unlabeled = names(which(table(MIDs_table$Bin)==1))
  if(length(which_unlabeled) > 0){ # if there are such compounds
    MIDs_unlabeled = MIDs_table[MIDs_table$Bin %in% which_unlabeled,]
    MIDs_for_correction = MIDs_table[!MIDs_table$Bin %in% which_unlabeled,]
    # also extract same info for other tables
    MIDs_tableScaled_unlabeled = MIDs_tableScaled[MIDs_tableScaled$Bin %in% which_unlabeled,]
    average_labeling_unlabeled = average_labeling[average_labeling$Bin %in% which_unlabeled,]
  }else{
    MIDs_for_correction = MIDs_table
  }

  MIDS_NACorrected = NACorrectionFxnII(MIDs_for_correction, label_scheme, output)

  #get average and MIDs table
  listOfAllOutputsNAcorrected = MIDsTableNoScale(MIDS_NACorrected, label_scheme)


  MIDs_tableNAcorrected = listOfAllOutputsNAcorrected[[1]]
  if(nchar(MIDs_tableNAcorrected$Isotopologue[1]) < 4){ # if single-labeled, use M0,M1,...
    MIDs_tableNAcorrected$Isotopologue = paste0('M', MIDs_tableNAcorrected$total_isotopes)
  }
  average_labelingNAcorrected = subset(listOfAllOutputsNAcorrected[[2]], select=c(-Isotopologue))
  vecOfNoMIDsNAcorrected = listOfAllOutputsNAcorrected[[3]]


  #scale the MIDs table by pool sizes
  scaledMIDsTableNAcorrected = MIDs_tableNAcorrected %>%
    group_by(Bin) %>%
    left_join(. , dplyr::select(compounds_data, c('prefix',contains('poolsize'))), by=join_by(Bin==prefix)) %>%
    mutate(across(contains(paste0(categories[1], "_")), ~.*eval(parse(text=paste0("poolsize_", categories[1]))))) %>%
    mutate(across(contains(paste0(categories[2], "_")), ~.*eval(parse(text=paste0("poolsize_", categories[2]))))) %>%
    dplyr::select(-contains('poolsize')) %>%
    as.data.frame()


  #get the enrichment adjusted MIDs
  labelEnrichmentMIDsNAcorrected = scaledMIDsTableNAcorrected %>%
    dplyr::filter(total_isotopes != 0) %>%
    select(colnames(.)[colnames(.) %like% "_.+_"], Bin) %>%
    group_by(Bin) %>%
    summarise(across(everything(), sum)) %>%
    left_join(. , scaledMIDsTableNAcorrected %>%
                dplyr::filter(total_isotopes == 0) %>%
                dplyr::select(mz, rt, polarity, Bin, Compound, Formula),
              by=join_by(Bin)) %>%
    as.data.frame()


  # put all compounds that only had unlabeled isotopologues back to the tables
  if(length(which_unlabeled) > 0){ # if there are such compounds
    MIDs_unlabeled$CompoundPlaceholder = MIDs_unlabeled$Bin
    MIDs_tableNAcorrected = rbind(MIDs_tableNAcorrected, MIDs_unlabeled)
    MIDs_tableScaled_unlabeled$CompoundPlaceholder = MIDs_tableScaled_unlabeled$Bin
    scaledMIDsTableNAcorrected = rbind(scaledMIDsTableNAcorrected, MIDs_tableScaled_unlabeled)
    average_labeling_unlabeled$CompoundPlaceholder = average_labeling_unlabeled$Bin
    average_labelingNAcorrected = rbind(average_labelingNAcorrected, average_labeling_unlabeled)
  }


  #remove bins that have no M0s
  MIDs_tableNAcorrected = MIDs_tableNAcorrected[!MIDs_tableNAcorrected$Bin %in% vecOfNoMIDs,]
  scaledMIDsTableNAcorrected = scaledMIDsTableNAcorrected[!scaledMIDsTableNAcorrected$Bin %in% vecOfNoMIDs,]
  average_labelingNAcorrected = average_labelingNAcorrected[!average_labelingNAcorrected$Bin %in% vecOfNoMIDs,]
  labelEnrichmentMIDsNAcorrected = labelEnrichmentMIDsNAcorrected[!labelEnrichmentMIDsNAcorrected$Bin %in% vecOfNoMIDs,] 


  # subset by clean bins
  clean_objects_NAcorrected = unlist(lapply(unique(clean_subset$category), function(x) 
    return(list(MIDs_tableNAcorrected[MIDs_tableNAcorrected$Bin %in% clean_bins[[x]],],
                scaledMIDsTableNAcorrected[scaledMIDsTableNAcorrected$Bin %in% clean_bins[[x]],],
                average_labelingNAcorrected[average_labelingNAcorrected$Bin %in% clean_bins[[x]],],
                labelEnrichmentMIDsNAcorrected[labelEnrichmentMIDsNAcorrected$Bin %in% clean_bins[[x]],]))), recursive = FALSE)
  combs = expand.grid(c("MIDsTableNAcorrected_clean_", "MIDsTableScaledNAcorrected_clean_", "averageLabelingNAcorrected_clean_", "labelEnrichmentNAcorrected_clean_"), unique(clean_subset$category))
  names(clean_objects_NAcorrected) = paste0(combs$Var1, combs$Var2)


  if (!is.null(output))
  {
    write.csv(MIDs_tableNAcorrected, file =  file.path("NA_corrected", paste(output, "MIDsNAcorrected.csv", sep = "_")), row.names = FALSE)
    write.csv(scaledMIDsTableNAcorrected, file = file.path("NA_corrected", paste(output, "scaled_MIDsNAcorrected.csv", sep = "_")), row.names = FALSE)
    write.csv(average_labelingNAcorrected, file = file.path("NA_corrected", paste(output, "average_labelingNAcorrected.csv", sep = "_")), row.names = FALSE)
    write.csv(labelEnrichmentMIDsNAcorrected, file = file.path("NA_corrected", paste(output, "mole_equivalents_of_labelNAcorrected.csv", sep = "_")), row.names = FALSE)
    lapply(names(clean_objects_NAcorrected), function(name) {
      write.csv(clean_objects_NAcorrected[[name]], paste0("NA_corrected/", name, ".csv"), row.names = FALSE)
    })
  }

  #return all of outputs as a get_table_objects() object
  listReturn = list(MIDs = MIDs_table, scaled_MIDs = MIDs_tableScaled,
                    average_labeling = average_labeling, mol_equivalent_labeling = label_enrichment,
                    MIDS_Corrected = MIDs_tableNAcorrected, scaledMIDsCorrected = scaledMIDsTableNAcorrected,
                    averageLabeling_corrected = average_labelingNAcorrected, molEquivalent_corrected = labelEnrichmentMIDsNAcorrected)
  listReturn = c(listReturn, clean_objects, clean_objects_NAcorrected)
}
