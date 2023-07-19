#!/usr/bin/env R

parse_CNVPoncogeneScores_csv <- function(fn, prefix) {
  a <- read.csv(fn,header=T, sep="\t", stringsAsFactors = F) |> 
    dplyr::select(name, value) |> 
    tibble::column_to_rownames('name') |> 
    dplyr::rename(cnv_score = value) |> 
    t() |> 
    as.data.frame(stringsAsFactors = F) |> 
    dplyr::rename_with( ~ paste0(prefix, .x))
  
  return(a)
}

