checkNeighbors <- function(neighbors_nn, seurat_obj, cell_v, idCol_v = "OBJECTID", classCol_v = NULL, slideName_v, radius_v) {
  #' Check Neighbors
  #' @description Visual check of neighbor calculation
  #' @param neighbors_nn neighborhood object output by findCellNeighbors
  #' @param seurat_obj seurat object output by buildSeurat
  #' @param cell_v character vector name of the ID of a cell to check
  #' @param idCol_v name of column containing cell IDs
  #' @param classCol_v name of column that contains cell identity classes (not required)
  #' @param slideName_v name of the SpatialImage object
  #' @details placeholder
  #' @return Print a plot
  #' @export
  
  ### Select first cell if none provided
  if (is.null(cell_v)) {
    cell_v <- seurat_obj@meta.data[["OBJECTID"]][1]
  }
  
  ### Make sure that cell_v is character
  if (class(cell_v) != "character") {
    warning("cell_v argument is not a character. Converting to character. Make sure everything looks good.")
    cell_v <- as.character(cell_v)
  } # fi
  
  ### Make sure cell is in data
  if (!cell_v %in% seurat_obj@meta.data[[idCol_v]]) {
    warning(sprintf("cell_v %s is not in object. Choosing first cell instead\n", cell_v))
    cell_v <- seurat_obj@meta.data[["OBJECTID"]][1]
  }
  
  ### Split seurat object for easier manipulation
  spatial_df <- GetTissueCoordinates(seurat_obj, slideName_v)
  spatial_mat <- as.matrix(spatial_df[,c("x", "y")]); rownames(spatial_mat) <- spatial_df$cell
  meta_df <- seurat_obj@meta.data
  
  ### Get neighbor Ids and their coordinates
  ids_ls <- neighbors_nn$id
  idsToTest_v <- neighbors_nn$id[[as.character(cell_v)]]
  idsToTestCoords_df <- spatial_mat[idsToTest_v,,drop=F]
  
  ### Subset plotting data
  plot_df <- spatial_df[spatial_df$cell %in% rownames(idsToTestCoords_df),]
  
  ### Merge with colors
  if (!is.null(classCol_v)) {
    plot_df <- merge(plot_df, meta_df[,c(idCol_v, classCol_v)], by.x = "cell", by.y = idCol_v, sort = F)
  } # fi
  
  ### Make plot
  plot_gg <- ggplot(data = plot_df, aes(x = x, y = y)) +
    geom_point() +
    ggforce::geom_circle(aes(x0 = spatial_df[spatial_df$cell == cell_v, "x"], y0 = spatial_df[spatial_df$cell == cell_v, "y"], r = 50)) +
    geom_point(data = spatial_df[spatial_df$cell == cell_v,,drop=F], color = "red", size = 3) + coord_equal() + my_theme() +
    ggtitle(paste0("Cells within ", radius_v, "um of Cell ", cell_v))
  
  ### Add color
  if (!is.null(classCol_v)) {
    plot_gg <- plot_gg + geom_point(aes(color = !!sym(classCol_v)))
  } # fi
  
  ### Print plot
  print(plot_gg)
  
}
