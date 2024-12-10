plotROI <- function(seuratObj, fov_v, projName_v = "ap_dp20", colorBy_v = c("class", "RCN")) {
  #' Plot ROI
  #' @description make a dim plot of the ROI colored by population or neighborhood.
  #' @param seuratObj seurat object of mIHC data
  #' @param fov_v name of the seurat object's field of view. Should be the sample name, can be found in `seuratObj@images`
  #' @param projName_v name of project's color object from spatial package. Default to AP DP20
  #' @param colorBy_v vector of columns to color by. Will grep for "^colorBy_v[i]$" as well as "^colorBy_v[i][0-9]*_cluster"
  #' @details Use 'class' to plot the populations. Use 'RCN' to plot all RCNs, specify RCN09 to just do one.  
  #' If column to plot is neither class nor RCN, auto-generate colors.
  #' @return list of ggplots
  #' @export
  
  ### Get meta for easier use
  meta_dt <- as.data.table(seuratObj@meta.data)
  
  ### Wrangle colors to get columns to run.
  colorColumns_v <- NULL
  for (i in 1:length(colorBy_v)) {
    currCB_v <- colorBy_v[i]
    currGrep_v <- c(paste0("^", currCB_v, "$"),
                    paste0("^", currCB_v, "[0-9]*_cluster$"))
    currColumns_v <- c(grep(currGrep_v[1], colnames(meta_dt), value = T),
                       grep(currGrep_v[2], colnames(meta_dt), value = T))
    colorColumns_v <- c(colorColumns_v, currColumns_v)
  } # for i
  
  ### Make a plot for each column
  plot_lsgg <- list()
  for (i in 1:length(colorColumns_v)) {
    
    currColumn_v <- colorColumns_v[i]
    
    ## Get colors
    if (currColumn_v == "class") {
      objName_v <- paste0(projName_v, "Colors_v")
      if (objName_v %in% ls()) {
        currColors_v <- eval(as.name(objName_v))
      } else {
        warning(sprintf("Don't have colors for class.\nChecked: %s\n", objName_v))
        currColors_v <- NULL
      } # fi
    } else if (grepl("RCN", currColumn_v)) {
      currColors_v <- rcnColors_v[1:length(unique(meta_dt[[currColumn_v]]))]
    } else {
      currColors_v <- NULL
    } # fi
    
    ## Make plot
    currPlot_gg <- ImageDimPlot(object = seuratObj, fov = fov_v, axes = T,
                                group.by = currColumn_v, cols = currColors_v)
    
    ## Add title and theme
    currPlot_gg <- currPlot_gg + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0(fov_v, " Cells colored by\n", currColumn_v))
    
    ## Add to list
    plot_lsgg[[currColumn_v]] <- currPlot_gg
    
  } # for i
  
  ### Return all the plots
  return(plot_lsgg)
  
} # plotROI