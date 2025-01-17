#' @title Perform k-means clustering and create plots from non-stationary labeling experiments
#'
#' @description This function allows you to perform k-means clustering analyses using label enrichment data at both label_enrichment and MID levels
#' The function can be used with both average_labeling and mol_equivalent_labeling data. If clustering at MID level is performed, the objects
#' MIDs and scaled_MIDs should be used with average_labeling and mol_equivalent_labeling respectively. The k-means clustering returns two plots
#' one on a 2D plot in PCA space with clusters grouped by compounds and the other is a time scale representation of label enrichment pattern
#' for each of the clusters. The clusters with the steepest slope represent compounds that are closer to the labeled source and a lower slope
#' indicates compounds farther away from the labeled source. The slopes can be used to infer direction of label flow within a pathway context
#'
#' @param mydata1 object (1) The 'average_labeling' or the 'mol_equivalent_labeling' objects created using get_table_objects() function
#' @param mydata2 object (1) The 'MIDs' or 'scaled_MIDs' objects created using get_table_objects() function. Use 'MIDs' if mydata1 is 'average_labeling'
#' and use 'scaled_MIDs' if mydata1 is 'mol_equivalent_labeling'
#' @param Category string (1) name of the category by which you are binning, for the example data set - Lipidomics has two categories (i.e. tissue types),
#' "Cotyledon" and "EA" and the DualLabel dataset has two catergories (i.e. genotypes), "WT" and "GAT"
#' @param nClust user specified number of clusters to be used for k-means analysis. The default is set to 0 in which the number of clusters
#' is calculated using a formula sqrt(n/2), where n is the number of bins being clustered.
#' @param labels string (1) label to be used for datapoints in the 2D plot i.e. "Bin" or "Compound" column
#' @param doMIDs logical (1) whether or not clustering should be done at the MIDs level. Default is set to FALSE
#' @param outputName string (1) The name to be appended to the the output pdf
#'
#' @seealso \link{get_table_objects_NA_corrected}, \link{MIDplot}, \link{label_enrichment_plot}, \link{PCA_and_heatmap}
#'
#' @return The function returns a list that contains all figures generated and a dataframe with cluster information. It also
#' provides a pdf file containing clustering results saved into specified directory
#'
#' @importFrom drc drm
#' @importFrom motifcluster kmeanspp
#'
#' @export
#'
#' @examples
#' \donttest{try(clusters_moleEquivalents <- getClusterAndPlots(test_13C15NGlutamine$mol_equivalent_labeling,
#' test_13C15NGlutamine$scaled_MIDs, Category = "WT", nClust=0, labels="Bin",
#' doMIDs=FALSE, outputName = "mol Equi"))}
#' \donttest{try(clusters_AvgLabeling <- getClusterAndPlots(test_13C15NGlutamine$average_labeling,
#' test_13C15NGlutamine$MIDs, Category = "WT", nClust=0, labels="Bin",
#' doMIDs=FALSE, outputName = "Average"))}

getClustersAndPlots <- function(mydata1, mydata2, Category, nClust=0, labels="Bin", doMIDs=FALSE, fit=FALSE, outputName = "_")
{

  ##function that will create all the dataframes
  createDataframe = function(columnnames)
  {
    columnnames = columnnames
    lengthColNames = length(columnnames)
    df =  data.frame(matrix(nrow = 0, ncol = lengthColNames))
    colnames(df) = columnnames
    return(df)
  }

  listOfAllOutputs = list()
  #first things first, make sure that we've used the right rownames
  #for the input tables

  #how many iterations for the kmeans - 25 should be more than enough for convergence
  howmanyN = 25
  mydata1 = as.data.frame(mydata1)
  mydata2 = as.data.frame(mydata2)
  labels = labels
  outputName = outputName

  #labels = "Bin"
  rownames(mydata1) = pull(mydata1, labels)
  rownames(mydata2) = make.unique(pull(mydata2, labels))

  #get the category name, pare down these variables, note well
  Category = Category
  #how many clusters for the kmeans

  #if nClust is 0, we will implement the sqrt
  nClust = nClust

  #subset the table by just the columns matching our category of interest
  catSubset = mydata1[,colnames(mydata1) %like% Category]

  # move everything down parallelly so that initial labeling wouldn't be a large number to cause confusion
  initial_time <- min(as.numeric(gsub("X", "", data_clean(names(catSubset)))))
  last_time <- max(as.numeric(gsub("X", "", data_clean(names(catSubset))))) # for logistic fit prediction
  min_initial_time <- do.call(pmin, catSubset[, c(names(catSubset) %like% paste0("X",initial_time,"_"))]) # min of the initial time for each compound
  # catSubset <- catSubset - min_initial_time

  #there'll be some na's so set them to zero
  catSubset[is.na(catSubset)] = 0

  #if the user is using the default kmeans setting
  #we'll using the sqrt formula to determine the number of clusters
  if(nClust == 0)
  {
    nClust = min(8, as.integer(sqrt(nrow(catSubset) / 2)))
  }

  #go ahead and do the kmeans right away on the normalized table
  kmeans_object = kmeanspp(catSubset, k=nClust, nstart = howmanyN, iter.max = 100)

  # summarize each row for logistic function fitting
  times = data_clean(colnames(catSubset)) %>%
    sub(".", "", .) %>%
    as.numeric()
  every_compound = do.call(rbind,
                           lapply(1:nrow(catSubset), function(x){
                             current_comp = data.frame(time=times, value=unlist(catSubset[x,])) %>%
                               group_by(time) %>%
                               summarise_at(vars(value), list(val = mean)) %>%
                               mutate(bin = rownames(catSubset)[x], Cluster = as.factor(kmeans_object[["cluster"]][x]))
                           }))

  #bring the kmeans and the PCA in here
  #have one object to store all of the plots
  plotsList <- list()

  #this data frame is going to have the global information for a complete k-mean cluster
  #across all timepoints
  dfAllClusters = createDataframe(c("Time","mean","sd","Cluster"))
  dfAllClusters$Time = as.numeric(dfAllClusters$Time)
  dfAllClusters$mean = as.numeric(dfAllClusters$mean)
  dfAllClusters$sd = as.numeric(dfAllClusters$sd)
  dfAllClusters$Cluster = as.numeric(dfAllClusters$Cluster)


  #once we've done the PCA's
  #progress to analyzing the kmeans analysis

  #for each of the clusters, we will want to know all of the data
  #to plot

  #for each kmeans cluster in the normalized label enrichment table
  #we're going to subset by just the compounds in that cluster
  #then we're going to get all the MIDs for these compounds
  #and we will further cluster them

  #we will collect information for all the metabolites + MID's in each cluster
  #through plots that focus on all of the compounds in a single cluster
  #as well as all the clusters to reflect the metabolism of pooled, similar compounds
  for(i in 1:nClust)
  {
    #cluster specific information
    theCluster = i

    #we reuse i before making it up here, so have a placeholder for the cluster number
    #since we will need it after the resetting
    theClusterAvgNum = i

    #plot out everything associated with each cluster
    #this is the kmeans results on the label enrichment table
    toSubsetTable = names(kmeans_object$cluster[which(kmeans_object$cluster == i)])

    #for each cluster, pull out just those compounds which fall in that cluster
    for_plottingFurther = subset(catSubset, rownames(catSubset) %in% toSubsetTable)

    #get their MIDs
    #make sure that the naming is consistent with previous table


    #moved this from 201
    #this is softcoded to get all of the timepoints from the column names
    vecOfExpTimes = unique(data_clean(colnames(for_plottingFurther)))


    ####get all the timepoint specific information about each cluster
    for(timepoint in 1:length(vecOfExpTimes))
    {
      whichTimepoint = vecOfExpTimes[timepoint]
      #this is going to be the avg's and standard deviation for all of the compounds associated with a cluster
      myMeanToAdd = mean(rowMeans(for_plottingFurther[,sapply(strsplit(colnames(for_plottingFurther) , '[_]' ), `[` , 1) == whichTimepoint, drop = FALSE]))
      mySdToAdd = sd(rowMeans(for_plottingFurther[,sapply(strsplit(colnames(for_plottingFurther) , '[_]' ), `[` , 1) == whichTimepoint, drop = FALSE]))
      dfAllClusters = rbind(dfAllClusters, data.frame(Time=whichTimepoint, mean=myMeanToAdd, sd=mySdToAdd,Cluster=theClusterAvgNum))
    }


    #if we want to do the MIDs as well
    if(doMIDs == TRUE)
    {
      #we have to match up the label enrichment and the MIDs
      #only pull out those MIDs  of compounds within the cluster
      for_plottingMIDs <- mydata2[ which(mydata2[,colnames(mydata2) %in% c(labels)] %in% rownames(for_plottingFurther)), ]

      #get the bin info
      justBinInfo = for_plottingMIDs[labels]

      #make sure that it's unique
      justBinInfo = make.names(justBinInfo, unique = TRUE)

      for_plottingMIDs[labels] = NULL

      #get the number of clusters for the MIDs kmeans
      numberOfClusters = max(as.integer(sqrt(nrow(for_plottingMIDs) / 2)), 1)

      #make sure that our MIDs are just our Category of interest
      for_plottingMIDs = for_plottingMIDs[,colnames(for_plottingMIDs) %like% Category]

      #if we have any Na's set them to zero
      for_plottingMIDs[is.na(for_plottingMIDs)] = 0

      #now, do kmeans on the many MIDs associated
      #with the Avg concentrations that cluster together
      kmeans_MIDs_objects = kmeans(for_plottingMIDs, numberOfClusters, nstart = howmanyN, iter.max = 100)

      #as we iterate through clusters, this list will store our plots
      allAverageAndMIDSCluster = list()

      #store each MID cluster's information

      #for each of the clusters of MIDs
      #pool together all of the information
      #of the MIDs in that cluster
      dfAllClustersMIDs <- data.frame(
        Time=numeric(),
        mean=numeric(),
        sd=numeric(),
        Cluster=factor()
      )

      #within each cluster of MIDs
      #keep track of all of the abundance data
      dfMIDs <- data.frame(
        Time=numeric(),
        mean=numeric(),
        sd=numeric(),
        Compound=factor(),
        cluster=factor()
      )

      #build up a list of everything for MIDs clusters
      MidsPlotList = list()

      ##for each MIDs cluster
      #we're going to go through and collect the information to fill in the data frames
      for(j in 1:length(unique(kmeans_MIDs_objects$cluster)))
      {

        toSubsetTableII = names(kmeans_MIDs_objects$cluster[which(kmeans_MIDs_objects$cluster == j)])
        for_plottingSub = subset(mydata2, rownames(mydata2) %in% toSubsetTableII)
        for_plottingSub= for_plottingSub[,colnames(for_plottingSub) %like% Category]

        #theMIDs cluster is here
        theCluster = j

        ####get all the specific information about each compounds cluster
        for(timepoint in 1:length(vecOfExpTimes))
        {
          whichTimepoint = vecOfExpTimes[timepoint]
          myMeanToAdd =  mean(rowMeans(for_plottingSub[,sapply(strsplit(colnames(for_plottingSub) , '[_]' ), `[` , 1) == whichTimepoint]))
          mySdToAdd =  sd(as.vector(as.matrix(for_plottingSub[,sapply(strsplit(colnames(for_plottingSub) , '[_]' ), `[` , 1) == whichTimepoint])))

          #get all the MIDs associated with that
          dfAllClustersMIDs = rbind(dfAllClustersMIDs, data.frame(Time=whichTimepoint, mean=myMeanToAdd, sd=mySdToAdd,Cluster=theCluster))
        }


        #Get the compound specific information in each MIDs cluster
        for(k in 1:length(rownames(for_plottingSub)))
        {
          theCompound = rownames(for_plottingSub)[k]
          for(timepoint in 1:length(vecOfExpTimes))
          {
            whichTimepoint = vecOfExpTimes[timepoint]

            myMeanToAdd =  mean(as.numeric(for_plottingSub[rownames(for_plottingSub) == theCompound,sapply(strsplit(colnames(for_plottingSub) , '[_]' ), `[` , 1) == whichTimepoint]))
            mySdToAdd =  sd(as.numeric(for_plottingSub[rownames(for_plottingSub) == theCompound,sapply(strsplit(colnames(for_plottingSub) , '[_]' ), `[` , 1) == whichTimepoint]))
            #we're looking at all of the MIDs in a single cluster
            dfMIDs = rbind(dfMIDs, data.frame(Time=whichTimepoint, mean=myMeanToAdd, sd=mySdToAdd, Compound=theCompound, cluster = j))
          }

        }

        theTotalName = paste(Category, "compounds cluster", theClusterAvgNum, "MIDs cluster", j)
        whichMIDsCluster = j
        toSubsetPlotMIDs = dfMIDs[dfMIDs$cluster == j,]
        toSubsetPlotMIDs$Time = as.numeric(gsub("X","",as.character(toSubsetPlotMIDs$Time)))

        #we will re-order after sorting by the slope

        ###do the slope and sorting now
        #at the mids and averages level
        for(l in 1:length(unique(toSubsetPlotMIDs$Compound)))
        {


          #extract all of the MIDs
          myCompound = unique(toSubsetPlotMIDs$Compound)[l]

          #get the slope through the midpoint time
          midTime = unique(toSubsetPlotMIDs)$Time[length(unique(toSubsetPlotMIDs$Time)) / 2]
          startTime = unique(toSubsetPlotMIDs$Time)[1]
          p1 = subset(toSubsetPlotMIDs, subset=(Time=="0" & Compound== myCompound))$mean
          p2 = subset(toSubsetPlotMIDs, subset=(Time==midTime & Compound== myCompound))$mean


          theSlope = p2 - p1

          toSubsetPlotMIDs[toSubsetPlotMIDs$Compound ==  myCompound, "ClusterSlope"]  = rep(theSlope,length(toSubsetPlotMIDs$Compound[toSubsetPlotMIDs$Compound == myCompound]))

        }

        ##sort by the slope
        toSubsetPlotMIDs = toSubsetPlotMIDs[order(toSubsetPlotMIDs$ClusterSlope,decreasing=T),]

        labeledSortVec = vector()
        for(i in 1:length(unique(toSubsetPlotMIDs$Compound)))
        {
          labeledSortVec=  c(labeledSortVec,rep(i,length(unique(vecOfExpTimes))))
        }

        toSubsetPlotMIDs$slopeSorted = labeledSortVec

        #set the factor levels based off of the reordered by slope, now
        toSubsetPlotMIDs$SortNames = toSubsetPlotMIDs$Compound %>% factor(levels = unique(toSubsetPlotMIDs$Compound))

        #plot all of the MIDs from this cluster of MIDs associated with the Avgs Cluster
        MidsPlotList[[j]] =  ggplot(toSubsetPlotMIDs, aes(x=Time,y=mean,colour=SortNames,group=SortNames)) + geom_line() + ggtitle(theTotalName) + geom_ribbon(aes(ymax=mean + sd, ymin=mean - sd, linetype=NA, fill = factor(SortNames)), show.legend = F, alpha=.1)
      }


      dfAllClustersMIDs$Cluster = as.factor(dfAllClustersMIDs$Cluster)
      MIDsAllClustersPlotName = paste(" Global MIDs cluster for average cluster", theClusterAvgNum)

      #plot everything for the MIDs cluster, now
      dfAllClustersMIDs$Time = as.numeric(gsub("X","",as.character(dfAllClustersMIDs$Time)))


      #for each cluster of MIDs, calculate the collective slope
      for(j in 1:length(unique(dfAllClustersMIDs$Cluster)))
      {

        p1 = subset(dfAllClustersMIDs, subset=(Time==unique(dfAllClustersMIDs$Time)[1] & Cluster==j))$mean
        p2 = subset(dfAllClustersMIDs, subset=(Time==unique(dfAllClustersMIDs$Time)[(length(unique(dfAllClustersMIDs$Time)) / 2)] & Cluster==j))$mean

        theSlope = p2 - p1
        dfAllClustersMIDs[dfAllClustersMIDs$Cluster == j, "ClusterSlope"]  = rep(theSlope,length(dfAllClustersMIDs$Cluster[dfAllClustersMIDs$Cluster == j]))

      }


      dfAllClustersMIDs = dfAllClustersMIDs[order(dfAllClustersMIDs$ClusterSlope,decreasing=T),]


      labeledSortVec = vector()

      #each cluster will contains the information
      #of many compounds
      for(i in 1:length(unique(dfAllClustersMIDs$Cluster)))
      {
        labeledSortVec=  c(labeledSortVec,rep(i,  length(unique(vecOfExpTimes))))
      }

      dfAllClustersMIDs$slopeSorted = labeledSortVec

      #make sure that we have the legends sorted here, now
      dfAllClustersMIDs$SortNames = dfAllClustersMIDs$Cluster %>% factor(levels = unique(dfAllClustersMIDs$Cluster))

      listToAdd1 = list()
      listToAdd1[[1]] = ggplot(dfAllClustersMIDs, aes(x=Time,y=mean,colour=SortNames,group=SortNames)) + geom_line() + ggtitle(MIDsAllClustersPlotName) + geom_ribbon(aes(ymax=mean + sd, ymin=mean - sd, linetype=NA, fill = factor(SortNames)), show.legend = F, alpha=.1)

      #store all of the MIDs plots in a list
      allAverageAndMIDSCluster = append(listToAdd1, MidsPlotList)
      plotsList = append(plotsList, allAverageAndMIDSCluster)
    }
  }

  ######End the MIDs here !!!! #####



  #determine the cluster by slope
  #working with the averages
  dfAllClusters$ClusterSlope = dfAllClusters$Cluster
  dfAllClusters$Cluster = as.factor(dfAllClusters$Cluster)

  #determine which cluster has the best slope
  #subset by each cluster
  for(j in 1:length(unique(dfAllClusters$Cluster)))
  {

    p1 = 0
    p2 = 0

    p1 = subset(dfAllClusters, subset=(Time==unique(dfAllClusters$Time)[1] & Cluster==j))$mean
    p2 = subset(dfAllClusters, subset=(Time==unique(dfAllClusters$Time)[(length(unique(dfAllClusters$Time)) / 2)] & Cluster==j))$mean

    theSlope = 0
    theSlope = as.numeric(p2) - as.numeric(p1)

    dfAllClusters[dfAllClusters$Cluster == j, "ClusterSlope"]  =  theSlope
  }

  forAllAvgsCluster = paste(Category, ": K-means Clustering")


  ##for now, we want the ordering to be consistent with the kmeans biplot
  dfAllClusters = dfAllClusters[order(dfAllClusters$Cluster,decreasing=F),]


  #we've re-ordered after sorting by the slope
  allTimesVector = vector()
  for(time in 1:nClust)
  {
    allTimesVector  = c(allTimesVector,rep(time,  length(unique(vecOfExpTimes))))
  }
  dfAllClusters$slopeSorted = allTimesVector

  dfAllClusters$slopeSorted = as.factor(dfAllClusters$slopeSorted)
  dfAllClusters$Cluster = as.factor(dfAllClusters$Cluster)
  dfAllClusters$Time = as.numeric(gsub("X","",as.character(dfAllClusters$Time)))


  ##add the overall kmeans plot of the averages
  #allTotClusters = prepend(allTotClusters,allClusterForBeginning)
  listToAdd = list()
  listToAdd[[1]] = autoplot(kmeans_object, data = catSubset, label = TRUE, label.size = 3, label.hjust = 0.5, label.vjust = 1,
                            frame = TRUE, frame.type = 'norm', main = forAllAvgsCluster)
  listToAdd2 = list()
  listToAdd2[[1]] = ggplot(dfAllClusters, aes(x=Time,y=mean,colour=Cluster,group=Cluster)) + geom_line() + ggtitle(forAllAvgsCluster) + geom_ribbon(aes(ymax=mean + sd, ymin=mean - sd, linetype=NA, fill = factor(Cluster)), show.legend = F, alpha=.1)

  if(fit == TRUE){
    # add logistic curve fitting
    logistic_fit = do.call(rbind,
                           lapply(1:nClust, function(x){
                             subdata_cluster = every_compound %>%
                               dplyr::filter(Cluster == x) %>%
                               dplyr::group_by(time) %>%
                               dplyr::mutate(w = 1/abs(val - mean(val))) %>%
                               dplyr::mutate(w = (w-min(w)) / (max(w)-min(w))) %>%
                               ungroup()
                             # to avoid w being undefined in the case of only one of two compounds in a cluster
                             if(length(unique(subdata_cluster$bin)) == 1){ # if only one compound in a cluster
                               subdata_cluster$w = 1
                             }
                             if(length(unique(subdata_cluster$bin)) == 2){ # if only two compounds in a cluster
                               subdata_cluster$w = 0.5
                             }
                             # logistic fit
                             mL <- drm(val~time, weights = w, data = subdata_cluster, fct = L.5(), type = "continuous")
                             # predictions and confidence intervals
                             subdata_fit <- expand.grid(time = seq(initial_time, last_time, length = 100))
                             pm <- predict(mL, newdata = subdata_fit, interval = "confidence")
                             subdata_fit$p <- pm[,1]
                             subdata_fit$pmin <- pm[,2]
                             subdata_fit$pmax <- pm[,3]
                             subdata_fit$cluster <- as.factor(x)
                             return(subdata_fit)
                           })
    )
    listToAdd2[[2]] = ggplot(dfAllClusters, aes(x = Time, y = mean, colour = Cluster, group = Cluster)) +
      geom_ribbon(aes(ymax = mean+sd, ymin = mean-sd, linetype = NA, fill = factor(Cluster)), show.legend = F, alpha = .1) +
      geom_line(data = logistic_fit, aes(x = time, y = p, color = cluster, group = cluster)) +
      ggtitle(paste0(forAllAvgsCluster, "_Logistic Fit"))
  }

  autoplotAndGlobalClusters = append(listToAdd, listToAdd2)

  ##prepend our global kmeans plots to the full list of plot
  allTotClusters = append(autoplotAndGlobalClusters, plotsList)

  pdf(paste(Category, outputName,"kmeans_plots.pdf", sep = "_"))
  for(i in allTotClusters){print(i)}
  dev.off()

  clusterInfo = data.frame(cluster=allTotClusters[[1]][["data"]][["cluster"]],
                           compound=allTotClusters[[1]][["data"]][["rownames"]]) %>%
    merge(., mydata1[,c("mz","rt")], by.x="compound", by.y=0) %>%
    arrange(cluster)

  write.csv(clusterInfo, paste0(Category, "_", nClust, "clust_", outputName, ".csv"), row.names = F)

  allTotClusters = append(allTotClusters, list(clusterInfo))

  return(allTotClusters)
}
