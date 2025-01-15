wrangleNeighbors <- function(seurat_ls, neighbors_lsNN, sampleNames_v = NULL, overwrite_v = F, verbose_v = F,
                             cellIDCol_v = "ObjectNumber", classCol_v = "class", metaCols_v = c(cellIDCol_v, classCol_v), levels_v = NULL, 
                             neighborDir_v, suffix_v = "_neighbors.xlsx") {
  #' Wranlge Neighbors
  #' @description Either read in summarizeNeighbor results, or run it
  #' @param seurat_ls seurat object that has been wrangled for spatial mIHC (buildSeurat)
  #' @param neighbors_lsNN list of neighborhood objects output by findCellNeighbors
  #' @param sampleNames_v optional vector of names of samples in seurat_ls. If NULL (default), will use all samples
  #' @param overwrite_v logical indicating if existing neighbor data should be overwritten (TRUE). If FALSE (default), existing data will be read in, not re-run.
  #' @param verbose_v logical. If FALSE (default), no printing. If TRUE, will print for each sample whether it's being run or read in.
  #' @param cellIDCol_v character vector indicating which column
  #' @param classCol_v name of column that contains cell identity classes
  #' @param metaCols_v column names to extract for reference
  #' @param levels_v vector of factor levels to assign to classCol_v. If NULL (default), will assign levels automatically
  #' @param neighborDir_v directory to hold neighbor excel outputs. Will check here to read in results first.
  #' @param suffix_v file name suffix. Used to write and to search
  #' @return a list of either length(sampleNames_v) or length(seurat_ls). Each list element is a list of length 3, containing `mat`, `pct`, and `meta`
  #' 1. `mat` = one row per cell, one column per potential neighbor class, cell values are the number of the column's cell class that are neighbors of that row's cell
  #' 1. `pct` = same as `mat`, but values are percentages of total neighbors for that cell
  #' 1. `meta` = same number of rows as `mat` and `pct` (one per cell), and whatever columns are specified in metaCols_v.
  #' @export
  
  ### Wrangle Samples
  if (is.null(sampleNames_v)) sampleNames_v <- names(seurat_ls)
  
  ### Run for each
  out_lslsmat <- list()
  for (i in 1:length(sampleNames_v)) {
    
    ### Get file name
    currSampleName_v <- sampleNames_v[i]
    outFile_v <- file.path(neighborDir_v, paste0(currSampleName_v, "_neighbors.xlsx"))
    
    ### If it doesn't exist, have to run it.
    ### If it exists, read if overwrite
    if (!file.exists(outFile_v) | (file.exists(outFile_v) & overwrite_v)) {
      
      ### Update
      if (verbose_v) {
        if (!file.exists(outFile_v)) cat(sprintf("Neighbor file: %s does not exist. Running neighbor calculation.\n", outFile_v))
        if (file.exists(outFile_v) & overwrite_v) cat(sprintf("Neighbor file: %s exists, but overwrite == T. Running neighbor calculation.\n", outFile_v))
      } # fi verbose
      
      ### Make neighbor matrix
      nn_lsmat <- summarizeNeighbors(seurat_obj = seurat_ls[[currSampleName_v]], neighbors_nn = neighbors_lsNN[[currSampleName_v]], 
                                     levels_v = levels_v, metaCols_v = metaCols_v)
      
      ### Convert to data.table for export
      nn_lsdt <- sapply(names(nn_lsmat), function(name_v) {
        mat <- nn_lsmat[[name_v]]
        if (name_v != "meta") {
          out <- convertDFT(mat, newName_v = "ObjectNumber")
        } else {
          out <- as.data.table(mat)
        } # fi
        return(out) }, simplify = F, USE.NAMES = T)
      
      ### Write
      myOpenXWriteWkbk(nn_lsdt, file_v = outFile_v)
      
    } else {
      
      ### Update
      if (verbose_v) cat(sprintf("Neighbor file: %s exists. Loading file.\n", outFile_v))
      
      ## Load workbook
      nn_lsdt <- wrh.rUtils::readAllExcel(outFile_v)
      ## Convert to non-dt version
      nn_lsmat <- sapply(nn_lsdt, function(dt) {
        out <- as.matrix(convertDFT(dt, col_v = "ObjectNumber"))
        return(out)}, simplify = F, USE.NAMES = T)
    } # fi file exits
    
    ### Add to list
    out_lslsmat[[currSampleName_v]] <- nn_lsmat
    
  } # for i
  
  ### Return
  return(out_lslsmat)
  
} # wrangleNeighbors