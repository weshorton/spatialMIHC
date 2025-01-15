clusterNeighbors <- function(nn_lslsmat, sampleNames_v = "cohort", ks_lsv, seedClasses_v = "all") {
  #' Cluster K Neighbors
  #' @description Cluster neighborhood percent matrices using provided k values
  #' @param nn_lslsmat list of lists of neighbor matrices (count and pct) along with metadata
  #' @param sampleNames_v vector of names in nn_lslsmat to run. Usually just run for cohort.
  #' @param ks_lsv list of vectors of k's to test. One list element per sampleNames_v element, can have 1 or multiple k's to test per.
  #' @param seedClasses_v vector of cell classes to use as seeds. Default is all, can be any of the 'class' values
  #' @details It takes a long time to run the clustering so only choose a few ks.
  #' @return modified nn_lslsmat with columns for cluster results of each k in the 'pct' and 'meta' matrices
  #' @export
  
  ### Run for each sample
  for (i in 1:length(sampleNames_v)) {
    
    ### Get info
    currSampleName_v <- sampleNames_v[i]
    currKs_v <- ks_lsv[[currSampleName_v]]
    currMats_lsmat <- nn_lslsmat[[currSampleName_v]]
    currPct_mat <- currMats_lsmat$pct
    currClasses_v <- currMats_lsmat$meta[,"class"]
    seedName_v <- paste0(seedClasses_v, collapse = "_")
    
    ### Subset for seedClasses
    if (is.logical(all.equal(seedClasses_v, "all"))) seedClasses_v <- unique(currClasses_v)
    currCells_v <- names(currClasses_v[currClasses_v%in% seedClasses_v])
    currToCluster_mat <- currPct_mat[rownames(currPct_mat) %in% currCells_v,]
    
    ### Run kmeans on percent matrices
    clusters_lsv <- sapply(currKs_v, function(x) { 
      y <- kmeans(currToCluster_mat, centers = x, nstart = 10, iter.max = 50)
      return(y$cluster)}, simplify = F)
    names(clusters_lsv) <- paste0("k", currKs_v)
    
    ### Add them to the tables
    for (j in 1:length(clusters_lsv)) {
      
      ### Get info
      currClusterName_v <- names(clusters_lsv)[j]
      currClusters_v <- clusters_lsv[[currClusterName_v]]
      currColName_v <- paste0(currSampleName_v, "_", currClusterName_v, "clusters_onPct_", seedName_v, "_seedCells")
      
      ### Make column to add
      #currNewCol_v <- rep(NA, nrow(currPct_mat)); names(currNewCol_v) <- rownames(currPct_mat)
      currNewCol_v <- sapply(rownames(currPct_mat), function(x) {
        ifelse(x %in% names(currClusters_v), currClusters_v[x], NA)
      })
      
      ### Add to percent matrix
      currPct_mat <- cbind(currPct_mat, "newCol" = rep(NA, nrow(currPct_mat)))
      currPct_mat[rownames(currPct_mat %in% rownames(currToCluster_mat)), "newCol"] <- currClusters_v
      # currPct_mat <- cbind(currPct_mat, currClusters_v)
      colnames(currPct_mat)[colnames(currPct_mat) == "currClusters_v"] <- currColName_v
      ### Add to metadata
      #currMats_lsmat$meta <- cbind(currMats_lsmat$meta, currClusters_v)
      colnames(currMats_lsmat$meta)[colnames(currMats_lsmat$meta) == "currClusters_v"] <- currColName_v
    } # for j
    
    ### Add back to input
    nn_lslsmat[[currSampleName_v]] <- currMats_lsmat
    
  } # for i
  
  ### Return input
  return(nn_lslsmat)
  
} # clusterNeighbors