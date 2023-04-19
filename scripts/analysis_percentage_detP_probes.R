#!/usr/bin/env R 

# data('IlluminaHumanMethylationEPICmanifest') #BUGresolver
#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")


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


tmp <- glass_od.metadata.idats |> 
  dplyr::select(sentrix_id, channel_green) |> 
  dplyr::mutate(files = gsub("_Grn.idat$","",channel_green)) |> 
  dplyr::mutate(percentage.detP.signi = unlist(pbapply::pblapply(files, ratio_detP)))



saveRDS(tmp |> 
          dplyr::mutate(channel_green = NULL) |> 
          dplyr::mutate(files = NULL)
        , file="output/tables/percentage_detP_probes.txt")


