#!/usr/bin/env R

# load data ----


source('scripts/load_functions.R')


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}




# load libs ----



library(IlluminaHumanMethylationEPICmanifest)
data(IlluminaHumanMethylationEPICmanifest)
data('IlluminaHumanMethylationEPICmanifest')
library(minfi)
library(RFpurify)
data('IlluminaHumanMethylationEPICmanifest')
data('RFpurify_ABSOLUTE')
data('RFpurify_ESTIMATE')


# if(!require(minfiData)){
#   if(!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   BiocManager::install("minfiData")
# }

devtools::install_github('mwsill/RFpurify')


# analysis ----




fun <- function(idat_red, idat_grn, Basename) {
  
  df = data.frame(
    sample_red = idat_red,
    sample_grn = idat_grn,
    Basename = Basename
  )

  # df = data.frame(
  #   sample_red = sel$array_channel_red,
  #   sample_grn = sel$array_channel_green,
  #   Basename = sel$Basename
  # ) |> 
  #   head(n=1)
  
  #print(df)
  
  RGset <- minfi::read.metharray.exp(targets = df, force=T)
  MsetEx <- minfi::preprocessRaw(RGset)
  
  abs <- predict_purity(MsetEx, method="ABSOLUTE")
  #print(abs)
  
  est <- predict_purity(MsetEx, method="ESTIMATE")
  #print(est)
  
  
  out = data.frame(array_RFpurity_absolute=abs, array_RFpurity_estimate=est)
  
  
  return(out)
}



# 50 werkt
sel <- glass_od.metadata.array_samples |> 
  dplyr::filter(arraychip_version == "EPICv1") |> 
  #head(n=105) |> 
  dplyr::mutate(Basename = gsub("^(.+)_(Grn|Red).idat$","\\1", array_channel_red)) |> 
  dplyr::mutate(Slide = gsub("^.+/([^/_]+)_.+$","\\1", array_channel_red)) |> 
  dplyr::mutate(Array = gsub("^.+/[^_]+_([^_]+).+$","\\1", array_channel_red))

sel2 <- sel |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = fun(array_channel_red, array_channel_green, Basename)) |> 
  dplyr::ungroup() |> 
  dplyr::select(array_sentrix_id, tmp) |> 
  tidyr::unnest(tmp)



saveRDS(sel2, file="cache/analysis_RFpurity.Rds")


