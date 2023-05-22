get_gates_from_count<- function(marker, parentgate, gate_counts){
  parentdf <<- mget(ls(.GlobalEnv, pattern = parentgate), inherits = TRUE)
  marker_sort <<- rev(sort(as.vector(parentdf[[parentgate]][[marker]])))
  unqmarker <- sym(marker)
  new_gate_name <- paste0(parentgate,"_",marker)
  
  negative <- paste0(new_gate_name,"n") 
  if(!exists(negative)){
    threshold <<- as.numeric(marker_sort[gate_counts])
    
    if(grepl(new_gate_name, "Root_Area", fixed=TRUE)){
      new_gate_name <- "Cells"
    }
    
    if(length(threshold) == 0 | length(is.na(threshold)) == 0){
      assign(paste0(new_gate_name,"n"), parentdf[[parentgate]], envir = .GlobalEnv)
      assign(paste0(new_gate_name,"p"), parentdf[[parentgate]][0,], envir = .GlobalEnv)
      threshold=0
    } else if (is.na(threshold)) {
      assign(paste0(new_gate_name,"n"), parentdf[[parentgate]], envir = .GlobalEnv)
      assign(paste0(new_gate_name,"p"), parentdf[[parentgate]][0,], envir = .GlobalEnv)
      threshold=0
    } else {
      assign(paste0(new_gate_name,"p"), filter(parentdf[[parentgate]], !!unqmarker >= threshold), envir = .GlobalEnv)
      assign(paste0(new_gate_name,"n"), filter(parentdf[[parentgate]], !!unqmarker < threshold), envir = .GlobalEnv)
    }
    threshold_temp <- data.frame(marker, parentgate, threshold)
    assign(paste0("thrsh_",new_gate_name,"p"), threshold_temp, envir = .GlobalEnv)
  }
}
