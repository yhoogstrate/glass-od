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
  source('scripts/load_metadata.R')
}


source('scripts/load_functions.R')



# GLASS-NL ----


targets.glass_nl <- glass_nl.metadata.resections |> 
  dplyr::filter(!is.na(sentrix_id)) |> 
  dplyr::select(sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (214))
    return(.)
  })()


idats.glass_nl <-  list.files(path = "data/GLASS_NL/Methylation/Methylation Array Data/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/GLASS_NL/Methylation/Methylation Array Data/", filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (658))
    return(.)
  })() |>
  assertr::verify(file.exists(filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", filename)) |>
  dplyr::mutate(channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", filename)) |>
  tidyr::pivot_wider(id_cols = sentrix_id, names_from = channel, values_from = c(filename)) |>
  dplyr::rename(channel_green = Grn) |>
  dplyr::rename(channel_red = Red)


targets.glass_nl <- targets.glass_nl |> 
  dplyr::left_join(idats.glass_nl, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(channel_green)) |> 
  assertr::verify(!is.na(channel_red)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (214))
    return(.)
  })()



# GLASS-OD ----


targets.glass_od <- glass_od.metadata.idats |> 
  filter_GLASS_OD_idats(163) |> 
  dplyr::select(sentrix_id, channel_green, channel_red)


# combine ----

targets <- rbind(
  targets.glass_od |> dplyr::mutate(dataset = "GLASS-OD"),
  targets.glass_nl |> dplyr::mutate(dataset = "GLASS-NL")
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


### m-values ----


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


saveRDS(mvalue, "cache/mvalues.combi.Rds")
saveRDS(mvalue.mask, "cache/mvalues.combi.detP_mask.Rds")

rm(mvalue, mvalue.mask)
gc()



