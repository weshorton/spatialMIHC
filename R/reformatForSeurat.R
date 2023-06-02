reformatForSeurat <- function(df, idCol_v = "ObjectNumber", 
                              coordCols_v = c("OBJECTID", "Location_Center_X", "Location_Center_Y"), 
                              exprCols_v, metaCols_v) {
  #' Reformat data for Seurat
  #' @description Reformat cleaned data.frames output from annotation scripts so they can be converted to Seurat objects
  #' @param df clean data.frame output by annotations script
  #' @param idCol_v column that contains the cell IDs. Used to make sure all data.frames are in same order.
  #' @param coordCols_v character vector of the 3 columns that will be extracted in order to form the spatial component. 
  #' These must represent the cell ID (idCol_v) and the X/Y coordinates for the cell's center IN THAT ORDER. 
  #' Default is c("ObjectNumber", "XCENTRE", "YCENTRE")
  #' @param exprCols_v character vector of the marker intensity columns. (TRY TO MAKE THIS AUTOMATIC)
  #' @param metaCols_v character vector of any columns to be added to the Seurat object metadata. 
  #' Needs to include at a minimum the annotated cell class column and the cell ID column. 
  #' Can include any other metadata that you may want to visualize downstream. (_func columns are common)
  #' @details To construct a Seurat object, we need:
  #' A data.frame of x/y coordinates. This is used for plotting and is analagous to the XY coordinates generated from a UMAP call
  #' A data.frame of marker intensities. This is analogous to gene expression data in a normal scRNAseq Seurat object
  #' A data.frame of cell classifications. This is added to the metadata of the seurat object. (Used for neighborhood analysis)
  #' @return a list of data.frames: "intensity", "spatial", "meta"
  #' @export
  
  ### Intensity
  intensity_df <- t(df[,exprCols_v])
  colnames(intensity_df) <- df[[idCol_v]]
  
  ### Spatial
  spatial_df <- df[,coordCols_v]
  colnames(spatial_df) <- c("cell", "x", "y")
  
  ### Meta
  meta_df <- df[,metaCols_v]
  rownames(meta_df) <- meta_df[[idCol_v]]
  
  ### Output
  out_lsdf <- list("intensity" = intensity_df, "spatial" = spatial_df, "meta" = meta_df)
  return(out_lsdf)
  
}
