get_thresh_from_count<- function(parentgate, marker, gate_counts){
  parentdf <<- mget(ls(.GlobalEnv, pattern = parentgate), inherits = TRUE)
  marker_sort <<- rev(sort(as.vector(parentdf[[parentgate]][[marker]])))
  unqmarker <- sym(marker)
  new_gate_name <- paste0(parentgate,"_",marker)
  pos_th <- paste0("thrsh_",new_gate_name,"p") 
  if(!exists(pos_th)){
    threshold <<- as.numeric(marker_sort[gate_counts])
    if(length(threshold) == 0 | length(is.na(threshold)) == 0){
      threshold <<- 0
    }
    if(grepl(new_gate_name, "Root_Area", fixed=TRUE)){
      new_gate_name <- "Cells"
    }
    threshold_temp <- data.frame(parentgate, marker, threshold)
    assign(paste0("thrsh_",new_gate_name,"p"), threshold_temp, envir = .GlobalEnv)
  }
}