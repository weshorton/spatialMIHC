get_gate_parents <- function(finalgates, Root){
  dfcol <<- data.frame(matrix(vector(),0,2,
                              dimnames=list(c(), c("parent", "current_gate"))),
                       stringsAsFactors=F)
  Areaindx <- as.numeric(match('Area', colnames(Root)))
  for(m in 1:length(finalgates)){
    gate <- str_split(finalgates[m], "_", simplify = TRUE)
    gate.m <- str_remove(gate, "[np]")
    markerdepth <- (length(gate))
    current_gate <- as.character(gate.m[markerdepth]) # marker is the last marker and current gate
    parentdf <- str_c(gate[-markerdepth], collapse = "_")
    if(current_gate == "Cells"){
      parentdf <- "Root"
      current_gate <- "Area"
      newg <- c(parentdf, current_gate)
      dfcol <<- rbind(dfcol,newg)
    } else {
      parentdf <- "Cellsp"
      newg <- c(parentdf, current_gate)
      dfcol <<- rbind(dfcol,newg)
    }
    #names(dfcol)<<- c("parent", "current_gate")
  }
  names(dfcol)<<- c("parent", "current_gate")
}