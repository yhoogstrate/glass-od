#!/usr/bin/env R 


if(!exists('glass_od.metadata.idats')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


if(!exists('glass_nl.metadata.idats')) {
  source('scripts/load_GLASS-NL_metadata.R')
}


if(!exists('gsam.metadata.idats')) {
  source('scripts/load_G-SAM_metadata.R')
}




ratio_detP <- function(sentrix_path) {
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



# ratio_detP("data/GLASS_OD/Methylation_data/EPIC data Iris de Heer - data UMCU/EPIC data Iris de Heer/201496850071_R02C01")


#' do this for ALL possible samples - QC etc. is based on this
tmp <- glass_od.metadata.idats |> 
  dplyr::select(sentrix_id, channel_green) |> 
  dplyr::mutate(sentrix_path = gsub("_Grn.idat$","",channel_green)) |> 
  dplyr::mutate(channel_green = NULL) |> 
  dplyr::mutate(percentage.detP.signi = unlist(pbapply::pblapply(sentrix_path, ratio_detP)))



write.table(tmp |> 
          dplyr::mutate(channel_green = NULL) |> 
          dplyr::mutate(files = NULL)
        , file="output/tables/percentage_detP_probes.txt")


