#!/usr/bin/env R

parse_mnpQCReport_csv <- function(fn, prefix) {
  a <- read.csv(fn,header=T, sep="\t") |> 
    dplyr::mutate(key=paste0(type, "_",name, "_",color, "_",check.if, "_",warning, "_",fail)) |> 
    dplyr::select(key, norm) |> 
    tibble::column_to_rownames('key') |> 
    `colnames<-`('value') |> 
    t() |> 
    as.data.frame() |> 
    dplyr::rename_with( ~ paste0(prefix, .x)) 
  
  return(a)
}

