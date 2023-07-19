#!/usr/bin/env R


parse_mnp_FrozenFFPEstatus_table <- function(fn, prefix) {
  #fn = tmp$mnp_FrozenFFPEstatus_table[1]
  
  out <- read.table(fn) |> 
    dplyr::rename(predicted_array_type = V2) |> 
    dplyr::rename(predicted_sample_type = V3) |> 
    dplyr::rename_with( ~ paste0(prefix, .x)) |> 
    dplyr::rename(sentrix_id = 1)
  
  return(out)
}

