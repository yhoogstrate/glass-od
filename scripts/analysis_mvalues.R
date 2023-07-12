#!/usr/bin/env R


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


# GLASS-OD ----

## load all samples, qc and corr w/ qc stats ----


targets <- glass_od.metadata.idats |> 
  dplyr::filter(is.na(reason_excluded_patient)) |> 
  dplyr::filter(is.na(reason_excluded_resection)) |> 
  dplyr::filter(is.na(reason_excluded_resection_isolation)) |> 
  #dplyr::filter(is.na(reason_excluded_array_sample)) # those are typically the replicates
  assertr::verify(!duplicated(resection_id)) |> 
  assertr::verify(sentrix_id != "204808700074_R04C01") |>  # 0003-R3-repA
  dplyr::rename(Sample_Name = resection_isolation_id) |>
  dplyr::mutate(Array = gsub("^.+_","",sentrix_id)) |> 
  dplyr::rename(Slide = methylation_array_chip_id) |> 
  dplyr::mutate(Basename = gsub("_Grn.idat$","", channel_green)) |> 
  dplyr::select(Sample_Name, sentrix_id, Array,Slide, Basename)

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


saveRDS(mvalue, "cache/mvalues.Rds")
saveRDS(mvalue.mask, "cache/mvalues_detP_mask.Rds")

rm(mvalue, mvalue.mask)
gc()



