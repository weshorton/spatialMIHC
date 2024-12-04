removeAnnotationDuplicates <- function(df, idCol_v = "ObjectNumber", classCol_v = "class", keep_v = 2, verbose_v = F) {
  #' Remove Annotation Duplicates
  #' @description It's possible for the intensity and combination of markers in certain cells to cause that cell to be
  #' assigned multiple labels. Must pick one. 
  #' @param df clean data.frame output by annotations script
  #' @param idCol_v column that contains the cell IDs
  #' @param classCol_v column that contains the cell type annotations
  #' @param keep_v value indicating which class to keep. Acceptable values are "first" for 1st class, "last" for last class, 
  #' or a numeric value (e.g. 2) indicating which to take.
  #' @param verbose_v logical indicating whether or not to print messages.
  #' @details Take the cleaned, annotated data and search for duplicated cell IDs. Grab the cell annotations for these cells
  #' and split them. One annotation will be kept for subsequent analysis, the other will be removed and recorded.
  #' @return list of three data.frames: (1) "clean" contains the input data.frame with duplicates removed.
  #' (2) "dup" contains the duplicate IDs and their classes. (3) is a summary of 2.
  #' @import data.table
  #' @export
  
  ### Check class
  class_v <- class(df)
  if (length(class_v) == 1) {
    if (class_v == "data.frame") {
      df <- as.data.table(df)
    } else {
      stop("Only single-length class allowed is data.frame\n")
    } # fi data.frame
  } else {
    if (!is.logical(all.equal(class_v, c("matrix", "array")))) {
      df <- as.data.table(df)
    } else if (!is.logical(all.equal(class_v, c("data.table", "data.frame")))) {
      stop("Must input data.frame, matrix, or data.table.\n")
    } # fi matrix/dt
  } # fi length
  
  ### Find duplicated IDs
  dupIDs_df <- df[, .N, by = idCol_v][N > 1]
  
  if (dupIDs_df[,.N] == 0 & verbose_v) {
    
    cat("No duplicates found.")
    
  } else {
    
    ###
    ### Wrangle Duplicate IDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###
  
    ### Subset input data for these
    dupData_df <- df[get(idCol_v) %in% dupIDs_df[[idCol_v]], ]
    
    ### New column to identify which occurrence of ID we have
    dupData_df[, which := 1:.N, by = idCol_v]
    
    ### Split into separate data.frames for each occurrence
    lsdf <- split(dupData_df, dupData_df$which)
    
    ### Merge
    dupCompareIDs_df <- suppressWarnings(mergeDTs(data_lsdt = lsdf, mergeCol_v = idCol_v, keepCol_v = classCol_v, sort = F))
    
    ### Summarize
    dupSummary_dt <- as.data.table(table(apply(dupCompareIDs_df[,!idCol_v,with=F], 1, function(x) paste(x, collapse = "-_-"))))
    
    ###
    ### Grab Correct Class
    ###
    
    if (keep_v == "first") {
      
      out_df <- lsdf[[1]]
      other_df <- rbindlist(lsdf[-1])
      
    } else if (keep_v == "last") {
      
      warning("Attempting to grab the last ID for all. This hasn't been thoroughly tested.")
      last_v <- length(lsdf)
      out_df <- lsdf[[last_v]]
      counter_v <- 1
      while (nrow(out_df) < nrow(lsdf[[1]])) {
        missing_v <- setdiff(lsdf[[(last_v-counter_v)]][[idCol_v]], out_df[[idCol_v]])
        out_df <- rbind(out_df, lsdf[[(last_v-counter_v)]][lsdf[[(last_v-counter_v)]][[idCol_v]] %in% missing_v,])
        lsdf[[(last_v-counter_v)]] <- lsdf[[(last_v-counter_v)]][!(lsdf[[(last_v-counter_v)]][[idCol_v]] %in% missing_v),]
        counter_v <- counter_v + 1
      } # while
      
      other_df <- rbindlist(lsdf[-last_v])
      
    } else {
      
      out_df <- lsdf[[keep_v]]
      other_df <- rbindlist(lsdf[-keep_v])
      
    } # fi
    
    ###
    ### Remove from main
    ###
    
    cols_v <- setdiff(colnames(other_df), "which")
    df <- df[!other_df, on=cols_v]
    
    ### Outputs
    out_lsdf <- list("clean" = df, "dup" = dupCompareIDs_df, "summary" = dupSummary_dt)
    return(out_lsdf)
    cat("\n")
  
  } # fi
  
} # removeAnnotationDuplicates


# Potential other way to do 'last' option
# last_v <- length(lsdf)
# out_df <- lsdf[[last_v]]
# other_lsdf <- lsdf[-last_v]
# missing_v <- setdiff(other_lsdf[[1]][[idCol_v]], out_df[[idCol_v]])
# for (i in seq_along(other_lsdf)) {
#   out_df <- rbind(out_df, other_lsdf[[i]][idCol_v %in% missing_v, ])
#   other_lsdf[[i]] <- other_lsdf[[i]][!(idCol_v %in% missing_v), ]
# }
# other_df <- rbindlist(other_lsdf)


### Original function:
# ### Find Duplicated IDs
# dupIDs_df <- as.data.frame(table(df[[idCol_v]]))
# dupIDs_df <- dupIDs_df[dupIDs_df$Freq > 1,]
# dupIDs_v <- unique(as.numeric(as.character(dupIDs_df$Var1)))
# 
# ### Subset input data for these
# dupData_df <- df[df[[idCol_v]] %in% dupIDs_v,]
# 
# ### Get identities of duplicated cells
# dupCalls_v <- dupData_df[[classCol_v]]

###
### Something ~~~~~~~~~~~~~
###

# ### New column to identify which occurrence of ID we have
# dupData_df$which <- 1
# for (i in 1:nrow(dupIDs_df)) dupData_df[dupData_df[[idCol_v]] == dupIDs_df$Var1[i], "which"] <- 1:dupIDs_df$Freq[i]
# 
# ### Dataframe for each ocurrence
# which_v <- unique(dupData_df$which)
# if (length(which_v) > 3) warning("At least some cells have >3 classes. Have only tested this for up to 3 classes.\n")
# lsdf <- uniq_lsv <- list()
# for (w_v in which_v) {
#   curr_df <- dupData_df[dupData_df$which == w_v,]
#   currUniq_v <- unique(curr_df[[classCol_v]])
#   lsdf[[paste0("V", w_v)]] <- curr_df
#   uniq_lsv[[paste0("V", w_v)]] <- currUniq_v
# }
# 
# ### Merge
# dupCompareIDs_df <- mergeDTs(data_lsdt = lsdf, mergeCol_v = idCol_v, keepCol_v = classCol_v, sort = F)
# 
# ### Summarize
# dupSummary_dt <- as.data.table(table(apply(dupCompareIDs_df[,!idCol_v,with=F], 1, function(x) paste(x, collapse = "-_-"))))
# 
# ###
# ### Grab Correct Class
# ###
# 
# if (keep_v == "first") {
#   
#   out_df <- lsdf[[1]]
#   other_df <- do.call(rbind, lsdf[-1])
#   
# } else if (keep_v == "last") {
#   
#   warning("Attempting to grab the last ID for all. This hasn't been thoroughly tested.")
#   last_v <- length(which_v)
#   out_df <- lsdf[[last_v]]
#   counter_v <- 1
#   while (nrow(out_df) < nrow(lsdf[[1]])) {
#     missing_v <- setdiff(lsdf[[(last_v-counter_v)]][[idCol_v]], out_df[[idCol_v]])
#     out_df <- rbind(out_df, lsdf[[(last_v-counter_v)]][lsdf[[(last_v-counter_v)]][[idCol_v]] %in% missing_v,])
#     lsdf[[(last_v-counter_v)]] <- lsdf[[(last_v-counter_v)]][!(lsdf[[(last_v-counter_v)]][[idCol_v]] %in% missing_v),]
#     counter_v <- counter_v + 1
#   } # while
#   
#   other_df <- do.call(rbind, lsdf[-last_v])
#   
# } else {
#   
#   out_df <- lsdf[[keep_v]]
#   other_df <- do.call(rbind, lsdf[-keep_v])
#   
# } # fi
# 
# ###
# ### Remove from main
# ###
# 
# cols_v <- setdiff(colnames(other_df), "which")
# df <- df[!other_df, on=cols_v]
# 
# ### Outputs
# out_lsdf <- list("clean" = df, "dup" = dupCompareIDs_df, "summary" = dupSummary_dt)
# return(out_lsdf)
# cat("\n")
