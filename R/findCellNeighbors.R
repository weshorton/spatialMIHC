findCellNeighbors <- function(seurat_obj, slideName_v, radius_v = 50) {
  #' Find Cell Neighbors
  #' @description Use the spatial coordinates to find all neighbors within a certain radius
  #' @param seurat_obj seurat object that has an entry for a SpatialImage object (this will be true if made running buildSeurat)
  #' @param slideName_v name of the SpatialImage object
  #' @param radius_v size in um of radius within which to search for neighbors
  #' @details Something
  #' @export
  
  ### Create a matrix of the coordinates
  coords_df <- GetTissueCoordinates(seurat_obj, slideName_v)
  coords_mat <- as.matrix(coords_df[,c("x", "y")])
  rownames(coords_mat) <- coords_df$cell
  
  ### Find neighbors
  neighbors_nn <- frNN(x = coords_mat, eps = radius_v)
  
  ### Output
  return(neighbors_nn)
  
}