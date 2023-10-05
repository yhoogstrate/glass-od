#!/usr/bin/env R
#' combined reading and normalization of GLASS-NL and GLASS-OD data

# load data ----


library(ggplot2)
library(minfi)

library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) # BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")



source('scripts/load_functions.R')


if(!exists('metadata.cg_probes.epic')) {
  source('scripts/load_probe_annotations.R')
}


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}

if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}

if(!exists('gsam.metadata.array_samples')) {
  source('scripts/load_G-SAM_metadata.R')
}



# metadata ----


#'@todo four replicates need to be erased after analysis
metadata.glass_od <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(215) |>  #@todo replicates
  dplyr::select(array_sentrix_id, array_channel_green, array_percentage.detP.signi)


metadata.glass_nl <- glass_nl.metadata.array_samples |> 
  filter_GLASS_NL_idats(218) |> 
  dplyr::select(array_sentrix_id, array_channel_green, array_percentage.detP.signi)


metadata.gsam <- gsam.metadata.array_samples |> 
  filter_GSAM_idats(77) |> 
  dplyr::select(array_sentrix_id, array_channel_green, array_percentage.detP.signi)





# integrate & make m-values ----


targets <- rbind(
  metadata.glass_od |> dplyr::mutate(dataset = "GLASS-OD"),
  metadata.glass_nl |> dplyr::mutate(dataset = "GLASS-NL"),
  metadata.gsam |> dplyr::mutate(dataset = "G-SAM")
) |> 
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL, array_channel_red = NULL) |> 
  dplyr::mutate(Sample_Name = paste0(dataset, "_", array_sentrix_id)) |>
  dplyr::mutate(Array = gsub("^.+_","", array_sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1", array_sentrix_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 215 + 218 + 77)
    return(.)
  })()


RGSet <- minfi::read.metharray.exp(targets = targets, force = T) #red/green channel together


detP <- minfi::detectionP(RGSet, type = "m+u")
proc <- minfi::preprocessNoob(RGSet, offset = 0, dyeCorr = TRUE, verbose = TRUE, dyeMethod="reference") 
rm(RGSet)
gc()


### mask ----


masked.value <- ifelse(detP > 0.01 , NA, 1) |> 
  data.table::as.data.table(keep.rownames = "probe_id") |> 
  dplyr::filter(probe_id %in% (
    metadata.cg_probes.epic |> 
      dplyr::filter(MASK_general == F) |> 
      dplyr::pull(probe_id)
  )) |> 
  tibble::column_to_rownames('probe_id')  |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 760405)
    assertthat::assert_that(ncol(.) == 510)
    return(.)
  })()

stopifnot(sum(is.na(masked.value)) > 0)




### m-values ----


mvalue <- minfi::ratioConvert(proc, what = "M") |> 
  assays() |> 
  purrr::pluck('listData') |> 
  purrr::pluck("M") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 865859)
    assertthat::assert_that(ncol(.) == 510)
    return(.)
  })() |> 
  data.table::as.data.table(keep.rownames = "probe_id") |> 
  dplyr::filter(probe_id %in% (
    metadata.cg_probes.epic |> 
      dplyr::filter(MASK_general == F) |> 
      dplyr::pull(probe_id)
  )) |> 
  tibble::column_to_rownames('probe_id') |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 760405)
    return(.)
  })()


stopifnot(sum(is.na(mvalue)) == 0)
stopifnot(targets$array_sentrix_id == colnames(mvalue))


saveRDS(mvalue, "cache/mvalues.HQ_samples.Rds")


rm(mvalue)
gc()





### beta.values ----


beta.value <- minfi::ratioConvert(proc, what = "beta") |> 
  assays() |> 
  purrr::pluck('listData') |> 
  purrr::pluck("Beta") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 865859)
    assertthat::assert_that(ncol(.) == 510)
    return(.)
  })() |> 
  data.table::as.data.table(keep.rownames = "probe_id") |> 
  dplyr::filter(probe_id %in% (
    metadata.cg_probes.epic |> 
      dplyr::filter(MASK_general == F) |> 
      dplyr::pull(probe_id)
  )) |> 
  tibble::column_to_rownames('probe_id') |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 760405)
    return(.)
  })()


stopifnot(sum(is.na(beta.value)) == 0)
stopifnot(targets$array_sentrix_id == colnames(beta.value))


saveRDS(beta.value, "cache/beta.values.HQ_samples.Rds")


rm(beta.value)
gc()




# export & cleanup ----




rm(targets, detP, proc)
gc()



saveRDS(masked.value, "cache/detP_masked_values.HQ_samples.Rds")



rm(masked.value)
gc()




