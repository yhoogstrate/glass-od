#!/usr/bin/env R

# load data ----


source('scripts/load_functions.R')


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}




# load libs ----



#library(IlluminaHumanMethylationEPICmanifest)
#library(IlluminaHumanMethylation450kmanifest)
#library(tidyverse)
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




fun <- function(idat_red, idat_grn) {
  
  df = data.frame(
    sample_red = idat_red,
    sample_grn = idat_grn
  )
  
  RGset <- read.metharray.exp(targets = df)
  MsetEx <- preprocessRaw(RGset)
  
  abs <- predict_purity(MsetEx,method="ABSOLUTE")
  #print(abs)
  
  est <- predict_purity(MsetEx,method="ESTIMATE")
  #print(est)
  
  
  out <- rbind(out, slice %>% dplyr::mutate(absolute=abs,estimate=est))
}



sel <- glass_od.metadata.array_samples |> 
  head(n=3) |> 
  dplyr::mutate(Basename = gsub("^(.+)_(Grn|Red).idat$","\\1",fn)) |> 
  dplyr::mutate(Slide = gsub("^.+/([^/_]+)_.+$","\\1",fn)) |> 
  dplyr::mutate(Array = gsub("^.+/[^_]+_([^_]+).+$","\\1",fn))



df <- do.call(
  rbind,
  pbapply::pblapply(sel |> dplyr::pull(array_channel_red),
                    sel |> dplyr::pull(array_channel_green), FUN = fun)
) 
