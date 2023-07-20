#!/usr/bin/env R 

# load stuff ----


source('scripts/load_functions.R')


if(!exists('glass_od.metadata.idats')) {
  source('scripts/load_GLASS-OD_metadata.R')
}

if(!exists('glass_nl.metadata.idats')) {
  source('scripts/load_GLASS-NL_metadata.R')
}

if(!exists('gsam.metadata.idats')) {
  source('scripts/load_G-SAM_metadata.R')
}




#' test
#' ratio_detP("data/GLASS_OD/Methylation_data/EPIC data Iris de Heer - data UMCU/EPIC data Iris de Heer/201496850071_R02C01")


# calc detP and export ----

#' do this for ALL possible samples - QC etc. is based on this
tmp <- rbind(
  glass_od.metadata.idats |> dplyr::select(sentrix_id, channel_green),
  glass_nl.metadata.idats |> dplyr::select(sentrix_id, channel_green),
  gsam.metadata.idats |> dplyr::select(sentrix_id, channel_green)
  ) |> 
  dplyr::mutate(sentrix_path = gsub("_Grn.idat$","",channel_green)) |> 
  dplyr::mutate(channel_green = NULL) |> 
  dplyr::mutate(percentage.detP.signi = unlist(pbapply::pblapply(sentrix_path, calc_ratio_detP)))



write.table(tmp |> dplyr::mutate(sentrix_path = NULL), file="output/tables/percentage_detP_probes.txt")




