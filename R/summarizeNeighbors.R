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
  
  ### Check levels
  ### Should probably change this in wrangleNeighbors, but the levels_v argument there
  ### is the named color vector for plotting, so I have to use the names for the factorization, not the values.
  ### To make this work with either, I need to check if they're colors and then grab the names, if so.
  checkHex_v <- grepl("\\#[0-9A-Z]{6}", levels_v)
  if (length(which(checkHex_v))) levels_v <- names(levels_v)
  
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
  } # fi not factor
  
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
  
  ### Add Neighbors
  tmp_pct <- nn_pct[,levels_v]
  colnames(tmp_pct) <- paste0("pctNeighbor_", colnames(tmp_pct))
  meta <- merge(meta, tmp_pct, by = 0, sort = F)
  rownames(meta) <- meta$Row.names; meta$Row.names <- NULL
  
  ### Output
  out_lsmat <- list("mat" = nn_matrix, "pct" = nn_pct, "meta" = meta)
  
}