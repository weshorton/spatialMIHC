plotROI <- function(seuratObj, fov_v, projName_v = "ap_dp20", colorBy_v = c("class", "RCN"), dark_v = F, size_v = 1) {
  #' Plot ROI
  #' @description make a dim plot of the ROI colored by population or neighborhood.
  #' @param seuratObj seurat object of mIHC data
  #' @param fov_v name of the seurat object's field of view. Should be the sample name, can be found in `seuratObj@images`
  #' @param projName_v name of project's color object from spatial package. Default to AP DP20
  #' @param colorBy_v vector of columns to color by. Will grep for "^colorBy_v[i]$" as well as "^colorBy_v[i][0-9]*_cluster"
  #' @param dark_v logical passed to dark.background of ImageDimPlot. FALSE: Plot cells on white background; TRUE: black background
  #' @param size_v size of points. ImageDimPlot default is 0.5, which is a bit small.
  #' @details Use 'class' to plot the populations. Use 'RCN' to plot all RCNs, specify RCN09 to just do one.  
  #' If column to plot is neither class nor RCN, auto-generate colors.
  #' @return list of ggplots
  #' @export
  
  ###
  ### Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Get meta for easier use
  meta_dt <- as.data.table(seuratObj@meta.data)
  
  ### Empty list for plots
  plot_lsgg <- list()
  
  ###
  ### Wrangle colors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
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
  
  ###
  ### Construct Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Make a plot for each column
  for (i in 1:length(colorColumns_v)) {
    
    currColumn_v <- colorColumns_v[i]
    
    ###
    ### Get colors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###
    
    ### Class
    if (currColumn_v == "class") {
      
      ### Construct name of color vector in package
      objName_v <- paste0(projName_v, "Colors_v")
      
      ### Grab if exists, or set to null for auto-generation
      if (objName_v %in% ls("package:spatialMIHC")) {
        currColors_v <- eval(as.name(objName_v))
      } else {
        warning(sprintf("Don't have colors for class.\nChecked: %s\n", objName_v))
        currColors_v <- NULL
      } # fi
      
    ## RCNs
    } else if (grepl("RCN", currColumn_v)) {
      
      currColors_v <- rcnColors_v[1:length(unique(meta_dt[[currColumn_v]]))]
      
    ## Auto-generate otherwise
    } else {
      currColors_v <- NULL
    } # fi currColumn_v
    
    ###
    ### Order legend
    ###
    
    if (grepl("RCN", currColumn_v)) {
      seuratObj@meta.data[[currColumn_v]] <- factor(seuratObj@meta.data[[currColumn_v]], levels = sort(unique(seuratObj@meta.data[[currColumn_v]])))
      #Idents(seuratObj) <- currColumn_v
    }
    
    ###
    ### Build Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###
    
    ### Make it
    currPlot_gg <- ImageDimPlot(object = seuratObj, 
                                fov = fov_v, 
                                axes = T,
                                group.by = currColumn_v, 
                                cols = currColors_v,
                                dark.background = dark_v,
                                size = size_v)
    
    ## Add title and theme
    currPlot_gg <- currPlot_gg + 
      ggtitle(paste0(fov_v, " Cells colored by\n", currColumn_v)) +
      theme(plot.title = element_text(hjust = 0.5, size = 32),
            axis.text = element_text(size = 16), legend.text = element_text(size = 16),
            legend.title = element_text(size = 20), legend.key.size = unit(1, 'cm'))
    
    ## Make white background
    if (!dark_v) {
      currPlot_gg <- currPlot_gg + theme(panel.background = element_rect(fill = "white", colour = NA), 
                                         panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    } # fi !dark_v
       
    
    ## Add to list
    plot_lsgg[[currColumn_v]] <- currPlot_gg
    
  } # for i
  
  ### Return all the plots
  return(plot_lsgg)
  
} # plotROI