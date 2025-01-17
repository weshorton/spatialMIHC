plotSeedCellCompOfClusters <- function(meta_dt, classCol_v = "class", plotColGrep_v = NULL, plotCol_v = NULL, colors_v) {
  #' Plot Cell Composition of Clusters
  #' @description Stacked bar chart of each clusters' composition of seed cells.
  #' @param meta_dt cohort metadata
  #' @param classCol_v column name in meta_dt corresponding to cell classes
  #' @param plotColGrep_v regex to grep columns from meta_dt to plot
  #' @param plotCol_v column names in meta_dt to plot
  #' @param colors_v vector of colors for bars
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
  
  # Get counts and pcts for each column
  counts_lsdt <- countPlots_lsgg <- pcts_lsdt <- pctPlots_lsgg <- comboPlots_lsgg <-  list()
  for (col_v in plotCol_v) {
    
    ## Subset and factorize meta
    currMeta_dt <- meta_dt[!is.na(get(col_v)),]
    currMeta_dt[[col_v]] <- factor(currMeta_dt[[col_v]], levels = sort(unique(currMeta_dt[[col_v]])))
    
    ## Get counts and pct for output
    currCounts_dt <- table(currMeta_dt[,mget(c(classCol_v, col_v))])
    counts_lsdt[[col_v]] <- convertDFT(as.data.frame.matrix((currCounts_dt)), newName_v = classCol_v)
    pcts_lsdt[[col_v]] <- convertDFT(as.data.frame.matrix(apply(currCounts_dt, 2, function(x) round(x / sum(x) * 100, digits = 2))), newName_v = "Class")
    
    ## Make base plot
    base_gg <- ggplot(data = currMeta_dt, aes(x = !!sym(col_v), fill = !!sym(classCol_v))) +
      theme(plot.title = element_text(hjust = 0.5, size = 26), axis.text = element_text(size = 16), axis.title = element_text(size = 20)) +
      xlab("Cluster") +
      scale_fill_manual(values = colors_v, breaks = names(colors_v))
    
    ## Make each plot
    countPlots_lsgg[[col_v]] <- base_gg + geom_bar(stat = "count") + ggtitle("Count of Cells per Cluster") + ylab("Count")
    pctPlots_lsgg[[col_v]] <- base_gg + geom_bar(stat = "count", position = "fill") + ggtitle("Percentage of Cells per Cluster") + ylab("Percent")
    
    ## Combine
    combo_gg <- ggpubr::ggarrange(plotlist = list(countPlots_lsgg[[col_v]], pctPlots_lsgg[[col_v]]), ncol = 2, common.legend = T, legend = "bottom")
    comboPlots_lsgg[[col_v]] <- ggpubr::annotate_figure(p = combo_gg, top = ggpubr::text_grob(label = col_v, size = 30))
    
  } # for col_v
  
  ### Return outputs
  return(list("comboPlot" = comboPlots_lsgg, "countPlot" = countPlots_lsgg, "pctPlot" = pctPlots_lsgg, "counts" = counts_lsdt, "pcts" = pcts_lsdt))
  
} # plotSeedCellCompOfClusters