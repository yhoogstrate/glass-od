#!/usr/bin/env R

devtools::install_github("hhhh5/ewastools")


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


# calc ewastools qc and export ----

# skip running old files over and over
if(file.exists("output/tables/ewastools.txt")) {
  old <- read.table("output/tables/ewastools.txt")
} else {
  old <- data.frame()
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
                                          CONST_N_OD_VALIDATION_ALL_SAMPLES +
                                          CONST_N_CATNON_ALL_SAMPLES +
                                          CONST_N_GLASS_NL_ALL_SAMPLES +
                                          CONST_N_GSAM_ALL_SAMPLES +
                                          194
    ))
    return(.)
  })() |> 
  dplyr::mutate(gsm = gsub("_Grn.idat$","", array_channel_green)) |>
  dplyr::select(array_sentrix_id, gsm)



todo <- tmp |> 
  dplyr::filter(array_sentrix_id %in% old$array_sentrix_id == F)
print(dim(todo))

ewastools_qc <- function(fn) {
  idat <- ewastools::read_idats(fn, quiet=T)
  
  metrics <- as.data.frame(ewastools::control_metrics(idat)) |> 
    dplyr::rename_with( ~ paste0("array_ewastools_qc_", .x))
  
  return(metrics)
}


data <- todo |> 
  #head(n=10) |> 
  dplyr::rowwise() |>  # before everything is loaded at once
  dplyr::mutate(tmp = ewastools_qc(gsm)) |> 
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(gsm = NULL)



stopifnot(colnames(old) == colnames(data))
exp <- rbind(old, data)


stopifnot(nrow(exp) > 0)


write.table(exp, "output/tables/ewastools.txt")


rm(tmp, todo, data, exp, ewastools_qc)


