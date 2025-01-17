addClusters <- function(seurat_ls, nn_lslsmat, sampleNames_v = "cohort") {
  #' Add Clusters
  #' @description Add neighborhood clusters to seurat objects
  #' @param seurat_ls list of seurat objects
  #' @param nn_lslsmat list of neighborhood lists with clusterNeighbors run
  #' @param sampleNames_v vector of sample names to run. Mainly just use 'cohort'. See details
  #' @details
  #' clusterNeighbors can be run on individual samples, but is most commonly run on the entire cohort of samples.
  #' The default of this function will take the cohort clusters from clusterNeighbors and add them to each individual seurat object in seurat_ls
  #' If an individual sample name from both seurat_ls and nn_lslsmat is also in sampleNames_v, then those individual clusters will also be added.
  #' @return seurat objects with new obj@meta.data columns.
  #' @export
  
  for (i in 1:length(sampleNames_v)) {
    
    ### Get info
    currSampleName_v <- sampleNames_v[i]
    currCols_v <- grep(paste0(currSampleName_v, "_k[0-9]+clusters_"), colnames(nn_lslsmat[[currSampleName_v]]$meta), value = T)
    currMeta_mat <- nn_lslsmat[[currSampleName_v]]$meta[,currCols_v,drop=F]
    
    ### Add to all if it's cohort
    if (currSampleName_v == "cohort") {
      
      for (j in 1:length(seurat_ls)) {
        ### Get info
        currSample_v <- names(seurat_ls)[j]
        currObj <- seurat_ls[[currSample_v]]
        ### Subset Meta
        currSubMeta_mat <- currMeta_mat[grep(currSample_v, rownames(currMeta_mat)),]
        rownames(currSubMeta_mat) <- gsub(paste0(currSample_v, "_"), "", rownames(currSubMeta_mat))
        ### Add, factorize, add back
        currObj <- AddMetaData(currObj, currSubMeta_mat, col.name = colnames(currSubMeta_mat))
        for (c_v in currCols_v) currObj@meta.data[[c_v]] <- factor(currObj@meta.data[[c_v]], levels = sort(unique(currObj@meta.data[[c_v]])))
        seurat_ls[[currSample_v]] <- currObj
      } # for j
      
      ranCohort_v <- T
      
      ### If it's an individual sample, add it to that sample only
    } else {
      
      ### If it's not in seurat_ls, can't add it
      if (!currSampleName_v %in% names(seurat_ls)) {
        warning(sprintf("Provided sample name <%s> not in names of seurat_ls. Please check spelling!.\n", currSampleName_v))
        if (!ranCohort_v) warning("Also haven't run 'cohort' (which should have been first, if present).\n")
        next
        
        ### Add to meta.data if it's in there
      } else {
        
        ### Meta will be null if wasn't run
        if (is.null(currMeta_mat)) {
          warning(sprintf("Neighborhood clusters not run for <%s>. Skipping.\n", currSampleName_v))
          next
        } # fi
        
        ### Get info
        currObj <- seurat_ls[[currSampleName_v]]
        ### Subset meta
        currSubMeta_mat <- currMeta_mat[grep(currSampleName_v, rownames(currMeta_mat)),]
        rownames(currSubMeta_mat) <- gsub(paste0(currSampleName_v, "_"), "", rownames(currSubMeta_mat))
        ### Add, factorize, add back
        currObj <- AddMetaData(currObj, currSubMeta_mat, col.name = colnames(currSubMeta_mat))
        for (c_v in currCols_v) currObj@meta.data[[c_v]] <- factor(currObj@meta.data[[c_v]], levels = sort(unique(currObj@meta.data[[c_v]])))
        seurat_ls[[currSampleName_v]] <- currObj
      } # fi in seurat_ls
      
    } # fi sample is cohort
    
  } # for i
  
  return(seurat_ls)
  
} # addClusters