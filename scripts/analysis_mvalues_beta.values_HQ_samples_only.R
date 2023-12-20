#!/usr/bin/env R
#' combined reading and normalization of GLASS-NL and GLASS-OD data

# load data ----


library(ggplot2)
library(minfi)

library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) # BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")



source('scripts/load_constants.R')
source('scripts/load_functions.R')
#source('scripts/load_palette.R')
#source('scripts/load_themes.R')



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
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |>
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(array_sentrix_id, array_channel_green, array_percentage.detP.signi)


metadata.catnon <- glass_od.metadata.array_samples |> 
  dplyr::filter(patient_study_name == "CATNON") |> 
  dplyr::select(array_sentrix_id, array_channel_green, array_percentage.detP.signi)


metadata.od_validation <- glass_od.metadata.array_samples |> 
  filter_OD_validation_idats(CONST_N_OD_VALIDATION_INCLUDED_SAMPLES) |>
  dplyr::select(array_sentrix_id, array_channel_green, array_percentage.detP.signi)


metadata.glass_nl <- glass_nl.metadata.array_samples |> 
  filter_GLASS_NL_idats(218) |> 
  dplyr::select(array_sentrix_id, array_channel_green, array_percentage.detP.signi)


metadata.gsam <- gsam.metadata.array_samples |> 
  filter_GSAM_idats(77) |> 
  dplyr::select(array_sentrix_id, array_channel_green, array_percentage.detP.signi)





# integrate & make m-values ----

all_targets <- rbind(
  metadata.glass_od |> dplyr::mutate(dataset = "GLASS-OD"),
  metadata.catnon |> dplyr::mutate(dataset = "CATNON"),
  metadata.od_validation |> dplyr::mutate(dataset = "OD-validation"),
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
    assertthat::assert_that(nrow(.) == CONST_N_SAMPLES)
    return(.)
  })() |> 
  dplyr::mutate(cache_mask =         paste0("cache/masks_hq/",array_sentrix_id,".Rds")) |> 
  dplyr::mutate(cache_intensities =    paste0("cache/intensities_hq/",array_sentrix_id,".Rds")) |> 
  dplyr::mutate(cache_intensities_m =  paste0("cache/intensities_m_hq/",array_sentrix_id,".Rds")) |> 
  dplyr::mutate(cache_intensities_um = paste0("cache/intensities_um_hq/",array_sentrix_id,".Rds")) |> 
  dplyr::mutate(cache_mvalues =       paste0("cache/mvalues/",array_sentrix_id,".Rds")) |> 
  dplyr::mutate(cache_betavalues =         paste0("cache/betavalues_hq/",array_sentrix_id,".Rds")) |> 
  dplyr::mutate(cached = file.exists(cache_mask) &
                  file.exists(cache_intensities) & file.exists(cache_intensities_m) & file.exists(cache_intensities_um) & 
                  file.exists(cache_mvalues) & file.exists(cache_betavalues))


message(paste0("excluded n=",sum(all_targets$cached),"/",nrow(all_targets)," samples that were already cached"))




targets <- all_targets |>  # only process new samples
  dplyr::filter(!cached)


#options(warn = 1) # nasty setting set by some other packages to 2, making an internal warning crash minfi
for(target_id in targets$array_sentrix_id) { # yep loop, minor overhead but safe that it does not paralellize
  print(target_id)
  
  target <- all_targets |> 
    dplyr::filter(array_sentrix_id %in% c( "201496850071_R02C01", "203293640061_R08C01", target_id)) |> 
    (function(.) {
      #print(dim(.))
      assertthat::assert_that(nrow(.) > 1) # we need  to include 2 or more samples. some idats do not have identical file size (force=T), and normalisation will only be conducted on the subset. force=t adds all of them
      # the actual problem was that running the normalisation with one sample differed in output as compared to two or more files as input
      return(.)
    })()
  
  
  # load
  
  RGSet <- minfi::read.metharray.exp(targets = as.data.frame(target), force = T) #red/green channel together
  
  
  # detP & mask
  detP <- minfi::detectionP(RGSet, type = "m+u") |> 
    data.table::as.data.table(keep.rownames = "probe_id") |> 
    tibble::column_to_rownames('probe_id') |>
    dplyr::rename_with(~ gsub("^GSM[0-9]+_([^_]+_[^_]+)$","\\1",.x)) |> # I Guess GEO hard re-codes internal identifiers
    dplyr::select(all_of(target_id)) |> 
    (function(.) {
      assertthat::assert_that(ncol(.) == 1) # we need  to include 2 or more samples. some idats do not have identical file size (force=T), and normalisation will only be conducted on the subset. force=t adds all of them
      return(.)
    })()
  
  masked.value <- ifelse(detP > 0.01 , NA, 1) |> 
    data.table::as.data.table(keep.rownames = "probe_id") |> 
    dplyr::filter(probe_id %in% (
      metadata.cg_probes.epic |> 
        dplyr::filter(MASK_general == F) |> 
        dplyr::pull(probe_id)
    )) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
      assertthat::assert_that(ncol(.) == 2)
      return(.)
    })()
  
  stopifnot(sum(is.na(masked.value)) > 0)
  
  saveRDS(
    masked.value |> tibble::tibble(),
    file = target |> dplyr::filter(array_sentrix_id == target_id) |> dplyr::pull(cache_mask)
  )
  
  rm(detP, masked.value)
  
  # normalise
  
  proc <- minfi::preprocessNoob(RGSet, offset = 0, dyeCorr = T, verbose = TRUE, dyeMethod="single")  #dyeMethod="reference"
  
  
  # intensities
  
  
  intensities <- log2(minfi::getMeth(proc) + minfi::getUnmeth(proc)) |> 
    data.table::as.data.table(keep.rownames = "probe_id") |> 
    dplyr::rename_with(~ gsub("^GSM[0-9]+_([^_]+_[^_]+)$","\\1",.x)) |> # I Guess GEO hard re-codes internal identifiers
    dplyr::select(probe_id, target_id) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES)
      assertthat::assert_that(ncol(.) == 2)
      return(.)
    })() |> 
    dplyr::filter(probe_id %in% (
      metadata.cg_probes.epic |> 
        dplyr::filter(MASK_general == F) |> 
        dplyr::pull(probe_id)
    )) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
      return(.)
    })()
  
  
  stopifnot(sum(is.na(intensities)) == 0)
  stopifnot(colnames(intensities) == c("probe_id", target_id))
  
  saveRDS(
    intensities  |> tibble::tibble(),
    file = target |> dplyr::filter(array_sentrix_id == target_id) |> dplyr::pull(cache_intensities)
  )
  
  rm(intensities)
  
  
  # intensities M
  
  intensities_m <- log2(minfi::getMeth(proc)) |> 
    data.table::as.data.table(keep.rownames = "probe_id") |> 
    dplyr::rename_with(~ gsub("^GSM[0-9]+_([^_]+_[^_]+)$","\\1",.x)) |> # I Guess GEO hard re-codes internal identifiers
    dplyr::select(probe_id, target_id) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES)
      assertthat::assert_that(ncol(.) == 2)
      return(.)
    })() |> 
    dplyr::filter(probe_id %in% (
      metadata.cg_probes.epic |> 
        dplyr::filter(MASK_general == F) |> 
        dplyr::pull(probe_id)
    )) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
      return(.)
    })()
  
  
  stopifnot(sum(is.na(intensities_m)) == 0)
  stopifnot(colnames(intensities_m) == c("probe_id", target_id))
  
  saveRDS(
    intensities_m  |> tibble::tibble(),
    file = target |> dplyr::filter(array_sentrix_id == target_id) |> dplyr::pull(cache_intensities_m)
  )
  
  rm(intensities_m)
  
  
  # intensities UM
  
  intensities_um <- log2(minfi::getUnmeth(proc)) |> 
    data.table::as.data.table(keep.rownames = "probe_id") |> 
    dplyr::rename_with(~ gsub("^GSM[0-9]+_([^_]+_[^_]+)$","\\1",.x)) |> # I Guess GEO hard re-codes internal identifiers
    dplyr::select(probe_id, target_id) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES)
      assertthat::assert_that(ncol(.) == 2)
      return(.)
    })() |> 
    dplyr::filter(probe_id %in% (
      metadata.cg_probes.epic |> 
        dplyr::filter(MASK_general == F) |> 
        dplyr::pull(probe_id)
    )) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
      return(.)
    })()
  
  
  stopifnot(sum(is.na(intensities_um)) == 0)
  stopifnot(colnames(intensities_um) == c("probe_id", target_id))
  
  saveRDS(
    intensities_um  |> tibble::tibble(),
    file = target |> dplyr::filter(array_sentrix_id == target_id) |> dplyr::pull(cache_intensities_um)
  )
  
  rm(intensities_um)

  
  # M-value
  
  mvalue <- minfi::ratioConvert(proc, what = "M") |> 
    assays() |> 
    purrr::pluck('listData') |> 
    purrr::pluck("M") |> 
    data.table::as.data.table(keep.rownames = "probe_id") |> 
    dplyr::rename_with(~ gsub("^GSM[0-9]+_([^_]+_[^_]+)$","\\1",.x)) |> # I Guess GEO hard re-codes internal identifiers
    dplyr::select(probe_id, target_id) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES)
      assertthat::assert_that(ncol(.) == 2)
      return(.)
    })() |> 
    dplyr::filter(probe_id %in% (
      metadata.cg_probes.epic |> 
        dplyr::filter(MASK_general == F) |> 
        dplyr::pull(probe_id)
    )) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
      return(.)
    })()
  
  
  stopifnot(sum(is.na(mvalue)) == 0)
  stopifnot(colnames(mvalue) == c("probe_id", target_id))

  
  saveRDS(
    mvalue  |> tibble::tibble(),
    file = target |> dplyr::filter(array_sentrix_id == target_id) |> dplyr::pull(cache_mvalues)
  )
  
  
  rm(mvalue)
  
  
  # beta value
  
  beta.value <- minfi::ratioConvert(proc, what = "beta") |> 
    assays() |> 
    purrr::pluck('listData') |> 
    purrr::pluck("Beta") |> 
    data.table::as.data.table(keep.rownames = "probe_id") |> 
    dplyr::rename_with(~ gsub("^GSM[0-9]+_([^_]+_[^_]+)$","\\1",.x)) |> # I Guess GEO hard re-codes internal identifiers
    dplyr::select(probe_id, target_id) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES)
      assertthat::assert_that(ncol(.) == 2)
      return(.)
    })() |> 
    dplyr::filter(probe_id %in% (
      metadata.cg_probes.epic |> 
        dplyr::filter(MASK_general == F) |> 
        dplyr::pull(probe_id)
    )) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
      return(.)
    })()
  
  
  stopifnot(sum(is.na(beta.value)) == 0)
  stopifnot(colnames(beta.value) == c("probe_id", target_id))
  
  
  saveRDS(
    beta.value |> tibble::tibble(),
    file = target |> dplyr::filter(array_sentrix_id == target_id) |> dplyr::pull(cache_betavalues),
    compress=T
  )
  
  
  rm(beta.value)
}


# export ----
## mask ----


exp <- all_targets |>
  dplyr::mutate(mask_cached = file.exists(cache_mask)) |> 
  assertr::verify(mask_cached) |> 
  dplyr::pull(cache_mask) |> 
  pbapply::pblapply(readRDS) |> 
  purrr::reduce(function(x, y) dplyr::left_join(x, y, by = 'probe_id')) |>  # considerable faster than reshape::merge_all
  tibble::column_to_rownames('probe_id')


saveRDS(exp, "cache/masks_hq/masks_hq.Rds")


rm(exp)



## intensities ----


exp <- all_targets |>
  dplyr::mutate(intensities_cached = file.exists(cache_intensities)) |> 
  assertr::verify(intensities_cached) |> 
  dplyr::pull(cache_intensities) |> 
  pbapply::pblapply(readRDS) |> 
  purrr::reduce(function(x, y) dplyr::left_join(x, y, by = 'probe_id')) |> 
  tibble::column_to_rownames('probe_id')


saveRDS(exp, "cache/intensities_hq/intensities_hq.Rds")


rm(exp)



## intensities_m ----


exp <- all_targets |>
  dplyr::mutate(intensities_m_cached = file.exists(cache_intensities_m)) |> 
  assertr::verify(intensities_m_cached) |> 
  dplyr::pull(cache_intensities_m) |> 
  pbapply::pblapply(readRDS) |> 
  purrr::reduce(function(x, y) dplyr::left_join(x, y, by = 'probe_id')) |> 
  tibble::column_to_rownames('probe_id')


saveRDS(exp, "cache/intensities_m_hq/intensities_m_hq.Rds")


rm(exp)




## intensities_um ----


exp <- all_targets |>
  dplyr::mutate(intensities_um_cached = file.exists(cache_intensities_um)) |> 
  assertr::verify(intensities_um_cached) |> 
  dplyr::pull(cache_intensities_um) |> 
  pbapply::pblapply(readRDS) |> 
  purrr::reduce(function(x, y) dplyr::left_join(x, y, by = 'probe_id')) |> 
  tibble::column_to_rownames('probe_id')


saveRDS(exp, "cache/intensities_um_hq/intensities_um_hq.Rds")


rm(exp)



## m-values ----


exp <- all_targets |>
  dplyr::mutate(mvalues_cached = file.exists(cache_mvalues)) |> 
  assertr::verify(mvalues_cached) |> 
  dplyr::pull(cache_mvalues) |> 
  pbapply::pblapply(readRDS) |> 
  purrr::reduce(function(x, y) dplyr::left_join(x, y, by = 'probe_id')) |> 
  tibble::column_to_rownames('probe_id')


saveRDS(exp, "cache/mvalues/mvalues_hq.Rds")


rm(exp)



## beta-values ----


exp <- all_targets |>
  dplyr::mutate(betavalues_cached = file.exists(cache_betavalues)) |> 
  assertr::verify(betavalues_cached) |> 
  dplyr::pull(cache_betavalues) |> 
  pbapply::pblapply(readRDS) |> 
  purrr::reduce(function(x, y) dplyr::left_join(x, y, by = 'probe_id')) |> 
  tibble::column_to_rownames('probe_id')


saveRDS(exp, "cache/betavalues_hq/betavalues.hq")


rm(exp)



# 
# RGSet <- minfi::read.metharray.exp(targets = targets, force = T) #red/green channel together
# 
# 
# detP <- minfi::detectionP(RGSet, type = "m+u")
# proc <- minfi::preprocessNoob(RGSet, offset = 0, dyeCorr = T, verbose = TRUE, dyeMethod="single")  #dyeMethod="reference"
# #proc <- minfi::preprocessNoob(RGSet, offset = 0, dyeCorr = T, verbose = TRUE, dyeMethod="reference") 
# rm(RGSet)
# gc()
# 
# 
# 
# 
# 
# ### mask 
# 
# 
# masked.value <- ifelse(detP > 0.01 , NA, 1) |> 
#   data.table::as.data.table(keep.rownames = "probe_id") |> 
#   dplyr::filter(probe_id %in% (
#     metadata.cg_probes.epic |> 
#       dplyr::filter(MASK_general == F) |> 
#       dplyr::pull(probe_id)
#   )) |> 
#   tibble::column_to_rownames('probe_id') |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
#     assertthat::assert_that(ncol(.) == CONST_N_SAMPLES)
#     return(.)
#   })()
# 
# stopifnot(sum(is.na(masked.value)) > 0)
# 
# 
# 
# 
# ## intensities ---
# 
# 
# intensities <- log2(minfi::getMeth(proc) + minfi::getUnmeth(proc)) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == CONST_N_PROBES)
#     assertthat::assert_that(ncol(.) == CONST_N_SAMPLES)
#     return(.)
#   })() |> 
#   data.table::as.data.table(keep.rownames = "probe_id") |> 
#   dplyr::filter(probe_id %in% (
#     metadata.cg_probes.epic |> 
#       dplyr::filter(MASK_general == F) |> 
#       dplyr::pull(probe_id)
#   )) |> 
#   tibble::column_to_rownames('probe_id') |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
#     return(.)
#   })()
# 
# 
# stopifnot(sum(is.na(intensities)) == 0)
# stopifnot(targets$array_sentrix_id == colnames(intensities))
# 
# 
# saveRDS(intensities, "cache/intensities.HQ_samples.Rds")
# 
# 
# rm(intensities)
# gc()
# 
# 
# 
# ## Meth intensities ---
# 
# 
# meth_intensities <- log2(minfi::getMeth(proc)) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == CONST_N_PROBES)
#     assertthat::assert_that(ncol(.) == CONST_N_SAMPLES)
#     return(.)
#   })() |> 
#   data.table::as.data.table(keep.rownames = "probe_id") |> 
#   dplyr::filter(probe_id %in% (
#     metadata.cg_probes.epic |> 
#       dplyr::filter(MASK_general == F) |> 
#       dplyr::pull(probe_id)
#   )) |> 
#   tibble::column_to_rownames('probe_id') |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
#     return(.)
#   })()
# 
# 
# stopifnot(sum(is.na(meth_intensities)) == 0)
# stopifnot(targets$array_sentrix_id == colnames(meth_intensities))
# 
# 
# saveRDS(meth_intensities, "cache/meth_intensities.HQ_samples.Rds")
# 
# 
# rm(meth_intensities)
# gc()
# 
# 
# 
# ## Unmeth intensities ---
# 
# 
# unmeth_intensities <- log2(minfi::getUnmeth(proc)) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == CONST_N_PROBES)
#     assertthat::assert_that(ncol(.) == CONST_N_SAMPLES)
#     return(.)
#   })() |> 
#   data.table::as.data.table(keep.rownames = "probe_id") |> 
#   dplyr::filter(probe_id %in% (
#     metadata.cg_probes.epic |> 
#       dplyr::filter(MASK_general == F) |> 
#       dplyr::pull(probe_id)
#   )) |> 
#   tibble::column_to_rownames('probe_id') |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
#     return(.)
#   })()
# 
# 
# stopifnot(sum(is.na(unmeth_intensities)) == 0)
# stopifnot(targets$array_sentrix_id == colnames(unmeth_intensities))
# 
# 
# saveRDS(unmeth_intensities, "cache/unmeth_intensities.HQ_samples.Rds")
# 
# 
# rm(unmeth_intensities)
# gc()
# 
# 
# 
# ### m-values ---
# 
# 
# mvalue <- minfi::ratioConvert(proc, what = "M") |> 
#   assays() |> 
#   purrr::pluck('listData') |> 
#   purrr::pluck("M") |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == CONST_N_PROBES)
#     assertthat::assert_that(ncol(.) == CONST_N_SAMPLES)
#     return(.)
#   })() |> 
#   data.table::as.data.table(keep.rownames = "probe_id") |> 
#   dplyr::filter(probe_id %in% (
#     metadata.cg_probes.epic |> 
#       dplyr::filter(MASK_general == F) |> 
#       dplyr::pull(probe_id)
#   )) |> 
#   tibble::column_to_rownames('probe_id') |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
#     return(.)
#   })()
# 
# 
# stopifnot(sum(is.na(mvalue)) == 0)
# stopifnot(targets$array_sentrix_id == colnames(mvalue))
# 
# 
# saveRDS(mvalue, "cache/mvalues.HQ_samples.Rds")
# 
# 
# rm(mvalue)
# gc()
# 
# 
# 
# 
# 
# ### beta.values ---
# 
# 
# beta.value <- minfi::ratioConvert(proc, what = "beta") |> 
#   assays() |> 
#   purrr::pluck('listData') |> 
#   purrr::pluck("Beta") |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == CONST_N_PROBES)
#     assertthat::assert_that(ncol(.) == CONST_N_SAMPLES)
#     return(.)
#   })() |> 
#   data.table::as.data.table(keep.rownames = "probe_id") |> 
#   dplyr::filter(probe_id %in% (
#     metadata.cg_probes.epic |> 
#       dplyr::filter(MASK_general == F) |> 
#       dplyr::pull(probe_id)
#   )) |> 
#   tibble::column_to_rownames('probe_id') |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
#     return(.)
#   })()
# 
# 
# stopifnot(sum(is.na(beta.value)) == 0)
# stopifnot(targets$array_sentrix_id == colnames(beta.value))
# 
# 
# saveRDS(beta.value, "cache/beta.values.HQ_samples.Rds")
# 
# 
# rm(beta.value)
# gc()
# 
# 
# 
# # export & cleanup ---
# 
# 
# 
# 
# rm(targets, detP, proc)
# gc()
# 
# 
# 
# saveRDS(masked.value, "cache/detP_masked_values.HQ_samples.Rds")
# 
# 
# 
# rm(masked.value)
# gc()
# 
# 


