kClusterSweep <- function(nn_pct, ks_v = 10) {
  #' Test multiple k values for kmeans clustering
  #' @description to do
  #' @param nn_pct neighbor matrix converted to percentages
  #' @param ks_v either a single numeric value indicating to test k = 1 to k = ks_v;
  #' or a range of values to test.
  #' @details to do
  #' @return to do
  #' @import bluster
  #' @export
  
  ### Determine k
  if (length(ks_v) == 1) {
    cat(sprintf("Single k value provided. Going to test values from 1 to %s.\n", ks_v))
    ks_v <- 1:ks_v
  } else {
    cat(sprintf("Range of k values provided. Testing: %s\n", paste(ks_v, collapse = ' ')))
  } # fi
  
  ### Empty vector
  wcss <- vector("numeric", length(ks_v))
  #sil <- vector("numeric", length(ks_v))
  
  ### k means for each
  for (k in ks_v) {
    if (nrow(nn_pct) < 20) next  # careful this is arbitrary
    kmeansMod <- kmeans(nn_pct, centers = k, nstart = 10, iter.max = 300)
    wcss[k] <- kmeansMod$tot.withinss
    #sil[k] <- bluster::approxSilhouette(x = as.matrix(nn_pct), clusters = factor(kmeansMod$cluster))
  } # for k
  
  ### Plot
  print(plot(ks_v, wcss, type = "b", pch = 19))
  
}
