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


metadata.glass_od <- glass_od.metadata.array_samples |> 
  dplyr::filter(!is.na(array_sentrix_id)) |> 
  dplyr::filter(arraychip_version == "EPICv1") |> 
  dplyr::select(array_sentrix_id, array_channel_green) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (274 + CONST_N_OD_VALIDATION_INCLUDED_SAMPLES))
    return(.)
  })()


metadata.glass_nl <- glass_nl.metadata.array_samples |> 
  dplyr::filter(!is.na(array_sentrix_id)) |> 
  dplyr::select(array_sentrix_id, array_channel_green) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (235))
    return(.)
  })()


metadata.gsam <- gsam.metadata.array_samples |> 
  dplyr::filter(!is.na(array_sentrix_id)) |> 
  dplyr::select(array_sentrix_id, array_channel_green) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (79))
    return(.)
  })()






# integrate & make m-values ----


targets <- rbind(
  metadata.glass_od |> dplyr::mutate(dataset = "GLASS-OD"),
  metadata.glass_nl |> dplyr::mutate(dataset = "GLASS-NL"),
  metadata.gsam |> dplyr::mutate(dataset = "G-SAM")
  ) |> 
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL) |> 
  dplyr::mutate(Sample_Name = paste0(dataset, "_", array_sentrix_id)) |>
  dplyr::mutate(Array = gsub("^.+_", "", array_sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1", array_sentrix_id))


RGSet <- minfi::read.metharray.exp(targets = targets, force = T) #red/green channel together


#detP <- minfi::detectionP(RGSet, type = "m+u") # moved to external script
proc <- minfi::preprocessNoob(RGSet, offset = 0, dyeCorr = TRUE, verbose = TRUE, dyeMethod="single")  # dyeMethod="reference"

rm(RGSet)
gc()


### m-values for all samples ----


mvalue <- minfi::ratioConvert(proc, what = "M")

#stopifnot(rownames(mvalue) == rownames(detP))
#stopifnot(colnames(mvalue) == colnames(detP))


mvalue <-  mvalue |> 
  assays() |> 
  purrr::pluck('listData') |> 
  purrr::pluck("M")


#stopifnot(dim(mvalue) == dim(detP))


# mvalue.mask <- ifelse(detP > 0.01 , NA, 1) |> 
#   data.table::as.data.table(keep.rownames = "probe_id") |> 
#   dplyr::filter(probe_id %in% (
#     metadata.cg_probes.epic |> 
#       dplyr::filter(MASK_general == F) |> 
#       dplyr::pull(probe_id)
#   )) |> 
#   tibble::column_to_rownames('probe_id')


#dim(mvalue.mask)

mvalue <- mvalue |> 
  data.table::as.data.table(keep.rownames = "probe_id") |> 
  dplyr::rename_with(~ gsub("^GSM[0-9]+_([^_]+_[^_]+)$","\\1",.x)) |> # I Guess GEO hard re-codes internal identifiers
  dplyr::filter(probe_id %in% (
    metadata.cg_probes.epic |> 
      dplyr::filter(MASK_general == F) |> 
      dplyr::pull(probe_id)
  )) |> 
  tibble::column_to_rownames('probe_id')

dim(mvalue)






stopifnot(sum(is.na(mvalue)) == 0)
#stopifnot(sum(is.na(mvalue.mask)) > 0)


stopifnot(targets$array_sentrix_id == colnames(mvalue))


# barely an elbow - either ~0 or high
# detP.rowMeans <- rowMeans(detP)
# detP.rowMedians <- rowMedians(detP)
# plot(sort( detP.rowMeans), type="l")
# plot(sort( detP.rowMedians), type="l",ylim=c(0,0.000001))


# cleanup 

rm(proc) #detP,
gc()


saveRDS(mvalue, "cache/mvalues.all_samples.Rds")
#saveRDS(mvalue.mask, "cache/mvalues.all_samples.detP_mask.Rds")

rm(mvalue)# , mvalue.mask
gc()




