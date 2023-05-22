summarizeClusterNeighbors <- function(seurat_obj, clusters_v, clusterCol_v, idCol_v = "OBJECTID", nn_matrix) {
  #' Summarize Cluster Neighbors
  #' @description Summarize the neighbor composition of each cluster
  #' @param seurat_obj seurat object made by buildSeurat
  #' @param clusters_v clusters to iterate over. TO DO - get these from seurat object instead of having to supply them
  #' @param clusterCol_v name of column in seurat metadata to get clusters from
  #' @param idCol_v cell id column
  #' @param nn_matrix nn matrix made earlier
  #' @details to do
  #' @return summary table of neighbors within each cluster
  #' @export
  
  ### Sort clusters
  clusters_v <- sort(unique(clusters_v))
  
  ### Empty matrix
  summary_df <- NULL
  
  ### Iterate over each cluster and get summary
  for (j in 1:length(clusters_v)) {
    
    ### Get cluster
    currClust_v <- clusters_v[j]
    
    ### Get cells from that cluster
    currIndices_v <- which(seurat_obj[[clusterCol_v]] == currClust_v)
    currCells_v <- seurat_obj@meta.data[currIndices_v,idCol_v]
    
    ### Subset neighborhood matrix for just these cells
    currSub_mat <- nn_matrix[rownames(nn_matrix) %in% currCells_v,]
    
    ### Convert to percentage
    currSums_v <- colSums(currSub_mat)
    currSumPct_v <- round(currSums_v / sum(currSums_v) * 100, digits = 2)
    
    ### Add to output
    summary_df <- rbind(summary_df, currSumPct_v)
    
  } # for j
  
  ### Add row names
  summary_df <- as.data.frame(summary_df)
  rownames(summary_df) <- clusters_v
  
  ### Output
  return(summary_df)
  
}