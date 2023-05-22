processcounts <- function(countspath){
  #' Process Counts
  #' @description
    #' Read in count and gate data
  #' @param countspath path to counts data
  #' @return data frame of counts
  #' @export
  
  ### Read
  countsdf <<- fread(countspath, header = TRUE)
  
  ### Get gates
  finalgates <<- countsdf$final_gate
  
  ### Get names
  countfilenames <<- as.list(tools::file_path_sans_ext(colnames(countsdf[,1:length(countsdf)]))) #column to match will = index + 1
  
  ### Make output
  count.tab <<- countsdf[,3:length(countsdf)]
  
  ### Return
  return(count.tab)
} # fi