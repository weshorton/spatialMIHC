fixImageDimPlot <- function(plot_gg) {
  #' Fix output of ImageDimPlot
  #' @description Edit ggplot object to avoid white circles around points
  #' @param plot_gg a ggplot/patchwork object output by ImageDimPlot()
  #' @details ImageDimPlot uses pch 21 as default point. This renders fine on Rstudio and markdown,
  #' but when printed to pdf, there is a white outline on each point, making the plot difficult to
  #' view. Have to update the colour parameter from "white" to match the fill
  #' @return object of classes "gtable", "gTree", "grob", and "gDesc". Note that print() will not display plot
  #' must use plot() on the output.
  #' @export
   
  
  ### Extract plot
  gg_build <- ggplot_build(plot_gg)
  
  ### Update color to be same as fill
  gg_build$data[[1]]$colour <- gg_build$data[[1]]$fill
  
  ### Reassemble
  gg_out <- ggplot_gtable(gg_build)
  
  ### Output
  return(gg_out)
  
} # fixImageDimPlot