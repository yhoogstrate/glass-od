#!/usr/bin/env R 


library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylation450kmanifest)
data(IlluminaHumanMethylation450kmanifest) # HACK!
#library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) # BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")



# load stuff ----


source('scripts/load_functions.R')


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}

if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}

if(!exists('gsam.metadata.array_samples')) {
  source('scripts/load_G-SAM_metadata.R')
}

if(!exists('tcga_laml.metadata.array_samples')) {
  source('scripts/load_TCGA-LAML_metadata.R')
}




#' test
#' ratio_detP("data/GLASS_OD/Methylation_data/EPIC data Iris de Heer - data UMCU/EPIC data Iris de Heer/201496850071_R02C01")


# calc detP and export ----

# skip running old files over and over
if(file.exists("output/tables/percentage_detP_probes.txt")) {
  old <- read.table("output/tables/percentage_detP_probes.txt") |> 
    dplyr::select(array_sentrix_id, array_percentage.detP.signi)
} else {
  old <- data.frame(array_sentrix_id = c(), array_percentage.detP.signi = c())
}



#' do this for ALL possible samples - QC etc. is based on this
tmp <- rbind(
  tcga_laml.metadata.array_samples |> dplyr::filter(array_type == "450k") |> dplyr::select(array_sentrix_id, array_channel_green),
  glass_od.metadata.array_samples |> dplyr::select(array_sentrix_id, array_channel_green),
  glass_nl.metadata.array_samples |> dplyr::select(array_sentrix_id, array_channel_green),
  gsam.metadata.array_samples |> dplyr::select(array_sentrix_id, array_channel_green)
  ) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_GLASS_OD_ALL_SAMPLES  +
                                        CONST_N_VALIDATION_ALL_SAMPLES +
                                        CONST_N_CATNON_ALL_SAMPLES +
                                        CONST_N_GLASS_NL_ALL_SAMPLES +
                                        CONST_N_GSAM_ALL_SAMPLES +
                                        194
    ))
    return(.)
  })() |> 
  dplyr::mutate(array_sentrix_path = gsub("_Grn.idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL) |> 
  
  dplyr::filter(array_sentrix_id %in% old$array_sentrix_id == F) |> 
  (function(.) {
    print(dim(.))
    return(.)
  })() |>
  dplyr::mutate(array_percentage.detP.signi = unlist(pbapply::pblapply(array_sentrix_path, calc_ratio_detP))) 



# write.table(tmp |> dplyr::mutate(array_sentrix_path = NULL), file="output/tables/percentage_detP_probes_TCGA_LAML.txt")
# write.table(tmp |> dplyr::mutate(array_sentrix_path = NULL), file="output/tables/percentage_detP_probes.txt")

# out <- rbind(
#   read.table("output/tables/percentage_detP_probes.txt") |> 
#     dplyr::rename_with( ~ paste0("array_", .x))
#   ,
#   tmp |> 
#     dplyr::select(array_sentrix_id, array_percentage.detP.signi)
# )


exp <- rbind(old, tmp |> dplyr::mutate(array_sentrix_path = NULL)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_GLASS_OD_ALL_SAMPLES  +
                                          CONST_N_VALIDATION_ALL_SAMPLES +
                                          CONST_N_CATNON_ALL_SAMPLES +
                                          CONST_N_GLASS_NL_ALL_SAMPLES +
                                          CONST_N_GSAM_ALL_SAMPLES +
                                          194
    ))
    return(.)
  })()

exp |> 
  write.table(file="output/tables/percentage_detP_probes.txt")


