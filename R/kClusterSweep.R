kClusterSweep <- function(nn_pct, ks_v = 10, seed_v = 123, verbose = F) {
  #' Sweep K for kmeans
  #' @description Test multiple k values for kmeans clustering
  #' @param nn_pct neighbor matrix converted to percentages
  #' @param ks_v either a single numeric value indicating to test k = 2 to k = ks_v;
  #' or a range of values to test.
  #' @param seed_v seed to set prior to running kmeans() for reproducibility. Will use same seed for each run of k.
  #' @param verbose_v logical indicating whether or not to print message.
  #' @details Run kmeans clustering for specified values of k. Calculate mean silhouette width, mean RMSD (root mean-squared deviation),
  #' and within-cluster sum of squares. Plot silhouette and RMSD together and wcss as an elbow plot.
  #' @return prints 2 plots to console
  #' @import bluster
  #' @export
  
  ### Determine k
  if (length(ks_v) == 1) {
    if (verbose) cat(sprintf("Single k value provided. Going to test values from 2 to %s.\n", ks_v))
    ks_v <- 2:ks_v
  } else {
    cat(sprintf("Range of k values provided. Testing: %s\n", paste(ks_v, collapse = ' ')))
    if (verbose) if (1 %in% ks_v) warning("Testing a k of 1 will not provide a silhouette score.")
  } # fi
  
  ### Empty vector
  wcss <- vector("numeric", length(ks_v))
  meanSil <- vector("numeric", length(ks_v))
  meanRmsd <- vector("numeric", length(ks_v))
  tocs_v <- NULL
  
  ### k means for each
  for (i in 1:length(ks_v)) {
    k <- ks_v[i]
    if (nrow(nn_pct) < 20) next  # careful this is arbitrary
    set.seed(seed_v)
    tictoc::tic()
    kmeansMod <- kmeans(nn_pct, centers = k, nstart = 10, iter.max = 300)
    wcss[i] <- kmeansMod$tot.withinss
    tocs_v <- c(tocs_v, tictoc::toc(quiet = T)$callback_msg)
    if (k == 1) {
      meanSil[i] <- NA
      meanRmsd[i] <- NA
    } else {
      tictoc::tic()
      sil <- as.data.frame(bluster::approxSilhouette(x = nn_pct, clusters = unname(kmeansMod$cluster)))
      tocs_v <- c(tocs_v, tictoc::toc(quiet = T)$callback_msg)
      meanSil[i] <- mean(sil$width)
      tictoc::tic()
      rmsd <- clusterRMSD(x = nn_pct, clusters = unname(kmeansMod$cluster))
      tocs_v <- c(tocs_v, tictoc::toc(quiet = T)$callback_msg)
      meanRmsd[i] <- mean(rmsd)
    } # fi
    
  } # for k
  
  ### Build data.table for plot
  metric_dt <- data.table("K" = ks_v, "SilWidth" = meanSil, "RMSD" = meanRmsd)
  meltMetric_dt <- reshape2::melt(metric_dt, id.vars = "K")
  metric_gg <- ggplot(meltMetric_dt, aes(x = K, y = value, color = variable)) +
    geom_point() + geom_line() + theme_bw() + facet_wrap(~variable, ncol = 2, scales = "free_y") +
    scale_x_continuous(breaks = ks_v) +
    xlab("Number of Clusters") + ggtitle("Cluster QC Results") + theme(plot.title = element_text(hjust = 0.5))
    
  
  ### Plot
  print(plot(ks_v, wcss, type = "b", pch = 19))
  print(metric_gg)
  
}
