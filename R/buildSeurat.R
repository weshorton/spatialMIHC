buildSeurat <- function(data_lsdf = NULL, intensity_df = NULL, spatial_df = NULL, meta_df = NULL,
                        assayName_v = "mIHC", slideName_v) {
  #' Build Seurat Object using spatial data
  #' @description Use the output of reformatForSeurat to build a spatial seurat object
  #' @param data_lsdf list of intensity, spatial, and meta dfs output by reformatForSeurat
  #' @param intensity_df data.frame of marker intensity values to add to seurat object (only used if data_lsdf is null)
  #' @param spatial_df data.frame of spatial values to add to seurat object (only used if data_lsdf is null)
  #' @param meta_df data.frame of metadata values to add to seurat object (only used if data_lsdf is null)
  #' @param assayName_v character vector used to name assay when creating object
  #' @param slideName_v character vector used to label the coordinate data
  #' @details Creates a seurat object with the given "count" data, metadata, and coordinates.
  #' @export
  
  ### Split list, if present
  if (!is.null(data_lsdf)) {
    intensity_df <- data_lsdf$intensity
    spatial_df <- data_lsdf$spatial
    meta_df <- data_lsdf$meta
  } else {
    if (is.null(intensity_df)) stop("Must set either data_lsdf or intensity_df")
    if (is.null(spatial_df)) stop("Must set either data_lsdf or spatial_df")
    if (is.null(meta_df)) stop("Must set either data_lsdf or meta_df")
  }
  
  ### Create object using the intensity data
  seurat_obj <- suppressWarnings(CreateSeuratObject(counts = intensity_df, assay = assayName_v))
  
  ### Add metadata
  seurat_obj <- AddMetaData(seurat_obj, meta_df)
  
  ### Create centroids object
  cents <- CreateCentroids(spatial_df)
  
  # ### List of segmentation info
  # segmentationsData_ls <- list("centroids" = cents,
  #                              "segmentation" = NULL)
  
  ### FOV coordinates
  coords <- CreateFOV(#coords = segmentationsData_ls$centroids,
                      coords = cents,
                      type = c("segmentation", "centroids"),
                      molecules = NULL,
                      assay = assayName_v)
  
  ### Add to object
  seurat_obj[[slideName_v]] <- coords
  
  ### Output
  return(seurat_obj)
}