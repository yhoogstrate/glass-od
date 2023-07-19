#!/usr/bin/env R


calc_ratio_detP <- function(sentrix_path) {
  files <- data.frame(files = sentrix_path)
  
  detP <- minfi::read.metharray(files) |> 
    minfi::detectionP(type = "m+u") |> 
    as.data.frame(stringsAsFactors=F) |> 
    dplyr::rename(pval = dplyr::everything()) |> # rename first column to fixed name
    dplyr::mutate(fail = pval > 0.01) |> # define probes that fail
    dplyr::group_by(fail) |> 
    dplyr::summarise(count = dplyr::n()) |> 
    dplyr::mutate(percentage.detP.signi = count / sum(count) * 100)
  
  out <- detP |> 
    dplyr::filter(fail == T) |> 
    dplyr::pull(percentage.detP.signi)
  
  gc() 
  
  return (out)
}

