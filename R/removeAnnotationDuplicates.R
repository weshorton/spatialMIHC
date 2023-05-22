removeAnnotationDuplicates <- function(df, idCol_v = "OBJECTID", classCol_v = "class", mainClass_v = "Myeloid cells") {
  #' Remove Annotation Duplicates
  #' @description It's possible for the intensity and combination of markers in certain cells to cause that cell to be
  #' assigned multiple labels. Must pick one. 
  #' @param df clean data.frame output by annotations script
  #' @param idCol_v column that contains the cell IDs
  #' @param classCol_v column that contains the cell type annotations
  #' @param mainClass_v so far it appears that all duplicate IDs belong to one main class and an alternate class. 
  #' Not sure about this part entirely. Will this always be the case, or can there be multiple two-class combinations?
  #' @details Take the cleaned, annotated data and search for duplicated cell IDs. Grab the cell annotations for these cells
  #' and split them. One annotation will be kept for subsequent analysis, the other will be removed and recorded.
  #' @return list of two data.frames: (1) "clean" contains the input data.frame with duplicates removed.
  #' (2) "dup" contains the duplicate IDs and their classes.
  #' @export
  
  ### Find Duplicated IDs
  dupIDs_df <- as.data.frame(table(df[[idCol_v]]))
  dupIDs_df <- dupIDs_df[dupIDs_df$Freq > 1,]
  dupIDs_v <- unique(dupIDs_df$Var1)
  
  ### Subset input data for these
  dupData_df <- df[df[[idCol_v]] %in% dupIDs_v,]
  
  ### Get identities of duplicated cells
  dupCalls_v <- dupData_df[[classCol_v]]
  print(table(dupCalls_v))
  
  ### Split the duplicated IDs into two equal data.frames. Merge together to compare
  ### Ideally one data.frame will have all one cell type and the other will be a mix
  dupData_df$which <- 1
  for (d_v in dupIDs_v) dupData_df[dupData_df[[idCol_v]] == d_v, "which"] <- c(1,2)
  dup1_df <- dupData_df[dupData_df$which == 1,]; uniq1_v <- unique(dup1_df[[classCol_v]])
  dup2_df <- dupData_df[dupData_df$which == 2,]; uniq2_v <- unique(dup2_df[[classCol_v]])
  dupCompareIDs_df <- merge(dup1_df[,c(idCol_v, classCol_v)], dup2_df[,c(idCol_v, classCol_v)], by = idCol_v, sort = F)
  
  ### Check dimensions
  if (nrow(dup1_df) != nrow(dup2_df)) stop("The two data.frames created by splitting the duplicate IDs are not the same size.\n")
  
  ### Determine which, if any, of the above has the unique class and also check if it's the specified main class.
  if (length(uniq1_v) == 1 & length(uniq2_v) == 1) {
    
    cat("Both sets of duplicate IDs have a single class (i.e. 1:1 class mapping)\n")
    if (uniq1_v == mainClass_v) {
      cat(sprintf("First data.frame matches main class %s.\n", mainClass_v))
      out_df <- dup1_df; other_df <- dup2_df
      outClass_v <- uniq1_v
    } else if (uniq2_v == mainClass_v) {
      cat(sprintf("Second data.frame matches main class %s.\n", mainClass_v))
      out_df <- dup2_df; other_df <- dup1_df
      outClass_v <- uniq2_v
    } else {
      warning(sprintf("Neither data.frame matches main class %s.\nUsing first data.frame class %s\n", mainClass_v, uniq1_v))
      outClass_v <- uniq1_v
    }
    
  } else if (length(uniq1_v) == 1 & length(uniq2_v) > 1) {
    
    cat("First data.frame has single class\n")
    if (uniq1_v == mainClass_v) {
      cat(sprintf("Selected class %s matches main class provided (%s).\n", uniq1_v, mainClass_v))
    } else {
      warning(sprintf("Selected class %s doesn't match main class provided (%s).\n", uniq1_v, mainClass_v))
    } # fi
    
    out_df <- dup1_df; other_df <- dup2_df
    outClass_v <- uniq1_v
    
  } else if (length(uniq1_v) > 1 & length(uniq2_v) == 1) {
    
    cat("Second data.frame has single class\n")
    if (uniq2_v == mainClass_v) {
      cat(sprintf("Selected class %s matches main class provided (%s).\n", uniq2_v, mainClass_v))
    } else {
      warning(sprintf("Selected class %s doesn't match main class provided (%s).\n", uniq2_v, mainClass_v))
    } # fi
    
    out_df <- dup2_df; other_df <- dup1_df
    outClass_v <- uniq2_v
    
  } else {
    
    warning("Both data.frames have >1 unique class. Keeping df1 and assigning df2 to 'other'. Be sure to check output\n")
    out_df <- dup1_df; other_df <- dup2_df
    outClass_v <- uniq1_v
    
  } # fi
  
  ### Remove the set we don't want
  df <- df[!(df[[idCol_v]] %in% dupData_df[[idCol_v]] & df[[classCol_v]] == outClass_v),]
  
  ### Outputs
  out_lsdf <- list("clean" = df, "dup" = dupCompareIDs_df)
  return(out_lsdf)
  cat("\n")
  
} # removeAnnotationDuplicates