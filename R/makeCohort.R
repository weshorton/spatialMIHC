makeCohort <- function(nn_lslsmat) {
  #' Make Cohort
  #' @description make cohort-level objects
  #' @param nn_lslsmat list of list of neighborhood matrices
  #' @return list of cohort level objects
  #' @export
  
  ### Add sample prefix to rownames for all three objects
  for (i in 1:length(nn_lslsmat)) {
    currName_v <- names(nn_lslsmat)[i]
    rownames(nn_lslsmat[[currName_v]]$mat) <- paste(currName_v, rownames(nn_lslsmat[[currName_v]]$mat), sep = "_")
    rownames(nn_lslsmat[[currName_v]]$pct) <- paste(currName_v, rownames(nn_lslsmat[[currName_v]]$pct), sep = "_")
    nn_lslsmat[[currName_v]]$meta <- cbind(nn_lslsmat[[currName_v]]$meta, "sample" = currName_v)
    rownames(nn_lslsmat[[currName_v]]$meta) <- paste(currName_v, rownames(nn_lslsmat[[currName_v]]$meta), sep = "_")
  } # for i
  
  ### Combine each
  out_ls <- list("mat" = do.call(rbind, sapply(names(nn_lslsmat), function(x) nn_lslsmat[[x]]$mat)),
                 "pct" = do.call(rbind, sapply(names(nn_lslsmat), function(x) nn_lslsmat[[x]]$pct)),
                 "meta" = do.call(rbind, unname(sapply(names(nn_lslsmat), function(x) nn_lslsmat[[x]]$meta, simplify = F))))
  
  ### Output
  return(out_ls)
} # makeCohort