#!/usr/bin/env R
#' combined reading and normalization of GLASS-NL and GLASS-OD data

# load data ----


library(ggplot2)
library(minfi)

library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) # BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")



if(!exists('youri_gg_theme')) {
  source('scripts/youri_gg_theme.R')
}


if(!exists('glass_od.metadata.idats')) {
  source('scripts/load_GLASS-OD_metadata.R')
}

if(!exists('glass_nl.metadata.idats')) {
  source('scripts/load_GLASS-NL_metadata.R')
}

if(!exists('gsam.metadata.idats')) {
  source('scripts/load_G-SAM_metadata.R')
}


source('scripts/load_functions.R')


# metadata ----


metadata.glass_od <- glass_od.metadata.idats |> 
  dplyr::filter(!is.na(sentrix_id)) |> 
  dplyr::select(sentrix_id, channel_green, percentage.detP.signi, mnp_QC_predicted_sample_type) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (222))
    return(.)
  })()


metadata.glass_nl <- glass_nl.metadata.idats |> 
  dplyr::filter(!is.na(sentrix_id)) |> 
  dplyr::select(sentrix_id, channel_green, percentage.detP.signi, mnp_QC_predicted_sample_type) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (235))
    return(.)
  })()


metadata.gsam <- gsam.metadata.idats |> 
  dplyr::filter(!is.na(sentrix_id)) |> 
  dplyr::select(sentrix_id, channel_green, percentage.detP.signi, mnp_QC_predicted_sample_type) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (75))
    return(.)
  })()






# integrate & make m-values ----

targets <- rbind(
  metadata.glass_od |> dplyr::mutate(dataset = "GLASS-OD"),
  metadata.glass_nl |> dplyr::mutate(dataset = "GLASS-NL"),
  metadata.gsam |> dplyr::mutate(dataset = "G-SAM")
  ) |> 
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","",channel_green)) |> 
  dplyr::mutate(channel_green = NULL, channel_red = NULL) |> 
  dplyr::mutate(Sample_Name = paste0(dataset,"_",sentrix_id)) |>
  dplyr::mutate(Array = gsub("^.+_","",sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1",sentrix_id))


RGSet <- minfi::read.metharray.exp(targets = targets, force = T) #red/green channel together


detP <- minfi::detectionP(RGSet, type = "m+u")
proc <- minfi::preprocessNoob(RGSet, offset = 0, dyeCorr = TRUE, verbose = TRUE, dyeMethod="reference") 
rm(RGSet)
gc()


### m-values for all samples ----


mvalue <- minfi::ratioConvert(proc, what = "M")

stopifnot(rownames(mvalue) == rownames(detP))
stopifnot(colnames(mvalue) == colnames(detP))


mvalue <-  mvalue |> 
  assays() |> 
  purrr::pluck('listData') |> 
  purrr::pluck("M")


stopifnot(dim(mvalue) == dim(detP))


#mvalue.mask <- mvalue |> 
#  magrittr::multiply_by(ifelse(detP > 0.01 , NA, 1))

mvalue.mask <- ifelse(detP > 0.01 , NA, 1) |> 
  data.table::as.data.table(keep.rownames = "probeID") |> 
  dplyr::filter(probeID %in% (
    read.delim("~/mnt/neuro-genomic-1-ro/catnon/Methylation - EPIC arrays/EPIC.hg38.manifest.tsv.gz") |> 
      dplyr::filter(MASK_general == F) |> 
      dplyr::pull(probeID)
  )) |> 
  tibble::column_to_rownames('probeID')

dim(mvalue.mask)


mvalue <- mvalue |> 
  data.table::as.data.table(keep.rownames = "probeID") |> 
  dplyr::filter(probeID %in% (
    read.delim("~/mnt/neuro-genomic-1-ro/catnon/Methylation - EPIC arrays/EPIC.hg38.manifest.tsv.gz") |> 
      dplyr::filter(MASK_general == F) |> 
      dplyr::pull(probeID)
  )) |> 
  tibble::column_to_rownames('probeID')

dim(mvalue)







stopifnot(sum(is.na(mvalue)) == 0)
stopifnot(sum(is.na(mvalue.mask)) > 0)


stopifnot(targets$sentrix_id == colnames(mvalue))

# cleanup 

rm(detP, proc)
gc()


saveRDS(mvalue, "cache/mvalues.all_samples.Rds")
saveRDS(mvalue.mask, "cache/mvalues.all_samples.detP_mask.Rds")

rm(mvalue, mvalue.mask)
gc()




