plotClusterCompOfSamples <- function(meta_dt, xVar_v, classCol_v = "class", plotColGrep_v = NULL, plotCol_v = NULL, colors_lsv,
                                     annot_lsv = NULL, annotColors_lsv) {
  #' Plot Cluster Composition of Samples
  #' @description Stacked bar chart of each samples' composition of cluster assignments.
  #' @param meta_dt cohort metadata
  #' @param xVar_v column in meta_dt to use for x-axis. Usually corresponds to sample, treatment, timepoint, etc.
  #' @param classCol_v column name in meta_dt corresponding to cell classes
  #' @param plotColGrep_v regex to grep columns from meta_dt to plot
  #' @param plotCol_v column names in meta_dt to plot
  #' @param colors_v list of vectors of colors for bars. Each lest element name should correspond to a plot column supplied by plotCol_v or plotColGrep_v.
  #' @param annot_lsv list of annotations to add below the plot. e.g. annot_lsv = list("outName1" = "colName1", "outName2" = "colName2")
  #' @param annotColors_lsv list of colors for annotations. list element names are same as annot_lsv names; list elements are named color vectors, names are values of annot_lsv colNames in meta_dt
  #' @return list of stacked bar charts - combined, count, pct and data tables of counts and pcts
  #' @details can only use one of plotColGrep_v or plotCol_v
  #' @export
  
  # Check variable assignment
  if (is.null(plotColGrep_v) & is.null(plotCol_v)) stop("Must have one of plotColGrep_v or plotCol_v set.\n")
  if (!is.null(plotColGrep_v) & !is.null(plotCol_v)) {
    warning(sprintf("Both plotColGrep_v (%s) and plotCol_v (%s) are set. Using plotCol_v.\n", plotColGrep_v, paste(plotCol_v, collapse = "; ")))
    plotColGrep_v <- NULL
  } # fi
  
  # Convert class
  class_v <- class(meta_dt)
  if (is.logical(all.equal(class_v, c("matrix", "array"))) |
      is.logical(all.equal(class_v, "data.frame"))) {
    meta_dt <- convertDFT(meta_dt, newName_v = "cell")
  }
  
  # Get columns
  if (!is.null(plotColGrep_v)) {
    plotCol_v <- grep(plotColGrep_v, colnames(meta_dt), value = T)
  } # fi
  
  # Handle xVar_v
  if (length(xVar_v) > 1) {
    newCol_v <- paste(xVar_v, collapse = "_")
    meta_dt[ , (newCol_v) := do.call(paste, c(.SD, sep = "_")), .SDcols = xVar_v]
    xVar_v <- newCol_v
  } # fi xVar_v > 1
  
  # Get counts and pcts for each column
  counts_lsdt <- countPlots_lsgg <- pcts_lsdt <- pctPlots_lsgg <- comboPlots_lsgg <-  list()
  for (col_v in plotCol_v) {
    
    ## Subset and factorize meta, get colors
    currMeta_dt <- meta_dt[!is.na(get(col_v)),]
    currMeta_dt[[col_v]] <- factor(currMeta_dt[[col_v]], levels = sort(unique(currMeta_dt[[col_v]])))
    currColors_v <- colors_lsv[[col_v]]
    
    ## Get counts and pct for output
    currCounts_dt <- table(currMeta_dt[,mget(c(xVar_v, col_v))])
    counts_lsdt[[col_v]] <- convertDFT(as.data.frame.matrix((currCounts_dt)), newName_v = xVar_v)
    pcts_lsdt[[col_v]] <- convertDFT(as.data.frame.matrix(apply(currCounts_dt, 2, function(x) round(x / sum(x) * 100, digits = 2))), newName_v = xVar_v)
    
    if (!is.null(annot_lsv)) {
      
      ## Make count plot
      countPlots_lsgg[[col_v]] <- wrh.rUtils::annotatedBar(data_dt = currMeta_dt, x_v = xVar_v, position_v = "stack", fill_v = col_v,
                                                           fillColors_v = currColors_v, annot_lsv = annot_lsv, annotColors_lsv = annotColors_lsv,
                                                           title_v = paste0("Count of Cluster Assignments per ", xVar_v))$annotPlot
      
      ## Make pct plot
      pctPlots_lsgg[[col_v]] <- wrh.rUtils::annotatedBar(data_dt = currMeta_dt, x_v = xVar_v, position_v = "fill", fill_v = col_v,
                                                         fillColors_v = currColors_v, annot_lsv = annot_lsv, annotColors_lsv = annotColors_lsv,
                                                         title_v = paste0("Percentage of Cluster Assignments per ", xVar_v))$annotPlot
    } else {
      
      ## Make base plot
      base_gg <- ggplot(data = currMeta_dt, aes(x = !!sym(xVar_v), fill = !!sym(col_v))) +
        theme(plot.title = element_text(hjust = 0.5, size = 26), axis.text = element_text(size = 16), axis.title = element_text(size = 20)) +
        xlab(xVar_v) +
        scale_fill_manual(values = currColors_v, breaks = names(currColors_v))

      ## Make each plot
      countPlots_lsgg[[col_v]] <- base_gg + geom_bar(stat = "count") + ggtitle(paste0("Count of Cluster Assignments per ", xVar_v)) + ylab("Count")
      pctPlots_lsgg[[col_v]] <- base_gg + geom_bar(stat = "count", position = "fill") + ggtitle(paste0("Percentage of Cluster Assignments per ", xVar_v)) + ylab("Percent")
    } # fi
    
    ## Combine - stack them if sample
    if (grepl("[Ss]ample", xVar_v)) {
      combo_gg <- ggpubr::ggarrange(plotlist = list(countPlots_lsgg[[col_v]], pctPlots_lsgg[[col_v]]), nrow = 2, common.legend = T, legend = "bottom")
      comboPlots_lsgg[[col_v]] <- ggpubr::annotate_figure(p = combo_gg, top = ggpubr::text_grob(label = col_v, size = 30))
    } else {
      combo_gg <- ggpubr::ggarrange(plotlist = list(countPlots_lsgg[[col_v]], pctPlots_lsgg[[col_v]]), ncol = 2, common.legend = T, legend = "bottom")
      comboPlots_lsgg[[col_v]] <- ggpubr::annotate_figure(p = combo_gg, top = ggpubr::text_grob(label = col_v, size = 30))
    } # fi
    
  } # for col_v
  
  ### Return outputs
  return(list("comboPlot" = comboPlots_lsgg, "countPlot" = countPlots_lsgg, "pctPlot" = pctPlots_lsgg, "counts" = counts_lsdt, "pcts" = pcts_lsdt))
  
} # plotSeedCellCompOfClusters
