summarizeNeighbors <- function(seurat_obj, neighbors_nn, classCol_v = "class", levels_v = NULL,
                               metaCols_v = c("ObjectNumber", "class") ) {
  #' Take the neighborhood output and determine which cell types they belong to
  #' @description Map the neighborhood$id values to cell types
  #' @param seurat_obj seurat object output by buildSeurat
  #' @param neighbors_nn neighborhood object output by findCellNeighbors
  #' @param classCol_v name of column that contains cell identity classes
  #' @param levels_v vector of factor levels to assign to classCol_v. If NULL (default) will assign levels automatically
  #' @param metaCols_v column names to extract for reference.
  #' @description To do
  #' @return list of data.frames. "mat" is counts while "pct" is percentages
  #' @export
  
  ### Double check objects
  if (!all.equal(colnames(seurat_obj), names(neighbors_nn$id))) {
    stop("Check agreement between seurat_obj columns and neighbors_nn$id names")
  }
  
  ### Classes must be factorized
  if (class(seurat_obj@meta.data[[classCol_v]]) != "factor") {
    cat("Cell type column 'class' is not a factor. Should have been converted at end of buildObject section.
        May want to check above to make sure ImageDimPlot has correct colors.")
    if (is.null(levels_v)) {
      seurat_obj@meta.data[[classCol_v]] <- factor(seurat_obj@meta.data[[classCol_v]])
    } else {
      seurat_obj@meta.data[[classCol_v]] <- factor(seurat_obj@meta.data[[classCol_v]], levels = levels_v)
    }
    Idents(seurat_obj) <- classCol_v
  }
  
  ### Map ids
  neighborClasses_lstab <- lapply(neighbors_nn$id, function(x) table(seurat_obj@meta.data[[classCol_v]][x]))
  nn_matrix <- do.call(rbind, neighborClasses_lstab)
  
  ### Calculate percentages
  nn_pct <- (nn_matrix / rowSums(nn_matrix)) * 100
  
  ### Conversion to percent causes cells with no neighbors to have NA
  ### instead of 0 across the board.
  nn_pct <- naTo0(nn_pct)
  #nn_pct <- nn_pct[rowSums(is.na(nn_pct)) != ncol(nn_pct),] # originally tried removing them
  
  ### Wrangle meta
  meta <- seurat_obj@meta.data[,metaCols_v]
  
  ### Output
  out_lsmat <- list("mat" = nn_matrix, "pct" = nn_pct, "meta" = meta)
  
}