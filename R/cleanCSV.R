cleanCSV <- function(df, grep_v = "^X$|^V1$", cols_v = NULL) {
  #' Clean CSV data.frame
  #' @description Current processing method causes these files to contain two extra columns that we don't want. 
  #' Use this to remove them
  #' @param df clean data.frame output by annotations script
  #' @param grep_v character vector of regex pattern to grab columns that need removing. Use either this or cols_v.
  #' This is default. set to NULL if you want to use cols_v.
  #' @param cols_v character vector of columns to remove. Use either this or grep_v
  #' @details Current processing method causes the input files to contain two extra columns that we don't want. 
  #' Consider taking a detailed look at the processing and use data.table, write out tables without row.names, or
  #' something else to eliminate this issue. Using either grep_v or cols_v, remove specified columns from data.frame
  #' @return same data.frame as before, with specified columns removed
  #' @export
  
  ### Check arguments
  if (!is.null(grep_v) & !is.null(cols_v)) stop("Both grep_v and cols_v are set. Choose one or the other")
  if (is.null(grep_v) & is.null(cols_v)) stop("Neither grep_v nor cols_v are set. Must choose one.")
  
  ### Grab columns
  if (!is.null(grep_v)) cols_v <- grep(grep_v, colnames(df), value = T)
  
  ### Remove columns
  if (length(cols_v) > 0) {
    for (col_v in cols_v) {
      temp_v <- head(df[[col_v]])
      cat(sprintf("Removing column %s, which looks like: %s\n", col_v, paste(temp_v, collapse = " ")))
      df[[col_v]] <- NULL
    } # for col_v
  } else {
    cat(sprintf("Specified columns (%s) not found in df. Please check.\n", paste(cols_v, collapse = " ")))
  } # fi
  
  ### Return
  return(df)
}