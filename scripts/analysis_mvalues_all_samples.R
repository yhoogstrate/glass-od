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
    assertthat::assert_that(nrow(.) == (CONST_N_GLASS_OD_ALL_SAMPLES +
                                        CONST_N_CATNON_ALL_SAMPLES + 
                                        CONST_N_OD_VALIDATION_INCLUDED_SAMPLES))
    return(.)
  })()



metadata.glass_nl <- glass_nl.metadata.array_samples |> 
  dplyr::filter(!is.na(array_sentrix_id)) |> 
  dplyr::select(array_sentrix_id, array_channel_green) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_NL_ALL_SAMPLES)
    return(.)
  })()


metadata.gsam <- gsam.metadata.array_samples |> 
  dplyr::filter(!is.na(array_sentrix_id)) |> 
  dplyr::select(array_sentrix_id, array_channel_green) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GSAM_ALL_SAMPLES)
    return(.)
  })()






# integrate & make m-values ----


all_targets <- rbind(
  metadata.glass_od |> dplyr::mutate(dataset = "GLASS-OD"), # includes validation + catnon
  metadata.glass_nl |> dplyr::mutate(dataset = "GLASS-NL"),
  metadata.gsam |> dplyr::mutate(dataset = "G-SAM")
  ) |> 
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL) |> 
  dplyr::mutate(Sample_Name = paste0(dataset, "_", array_sentrix_id)) |>
  dplyr::mutate(Array = gsub("^.+_", "", array_sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1", array_sentrix_id)) |> 
  dplyr::mutate(cache_mvalues = paste0("cache/mvalues/",array_sentrix_id,".Rds")) |> 
  dplyr::mutate(cached = file.exists(cache_mvalues))


message(paste0("excluded n=",sum(all_targets$cached),"/",nrow(all_targets)," samples that were already cached"))


targets <- all_targets |> # don't re-do stuff
  dplyr::filter(!cached)



# run normalisation per sample, saves lots of mem, requires 'single' dye correction
for(target_id in targets$array_sentrix_id) { # for loop, to ensure things don't start running in parallel and spoil the memory
  print(target_id)
  
  target <- all_targets |> 
    dplyr::filter(array_sentrix_id %in% c( "201496850071_R02C01", "203293640061_R08C01", target_id)) |> 
    (function(.) {
      #print(dim(.))
      assertthat::assert_that(nrow(.) > 1) # we need  to include 2 or more samples. some idats do not have identical file size (force=T), and normalisation will only be conducted on the subset. force=t adds all of them
                                           # the actual problem was that running the normalisation with one sample differed in output as compared to two or more files as input
      return(.)
    })()
  
  
  
  
  RGSet <- minfi::read.metharray.exp(targets = target, force = T, verbose=F) #red/green channel together
  proc <- minfi::preprocessNoob(RGSet, offset = 0, dyeCorr = TRUE, verbose = TRUE, dyeMethod="single")  # dyeMethod="reference"
  
  rm(RGSet)
  gc()
  
  
  ## m-values for all samples
  
  mvalue <- minfi::ratioConvert(proc, what = "M") |> 
    assays() |> 
    purrr::pluck('listData') |> 
    purrr::pluck("M") |> 
    data.table::as.data.table(keep.rownames = "probe_id") |> 
    dplyr::rename_with(~ gsub("^GSM[0-9]+_([^_]+_[^_]+)$","\\1",.x)) |> # I Guess GEO hard re-codes internal identifiers
    dplyr::select(probe_id, target_id) |> 
    (function(.) {
      assertthat::assert_that(ncol(.) == 2) # we need  to include 2 or more samples. some idats do not have identical file size (force=T), and normalisation will only be conducted on the subset. force=t adds all of them
      return(.)
    })() |> 
    dplyr::filter(probe_id %in% (
      metadata.cg_probes.epic |> 
        dplyr::filter(MASK_general == F) |> 
        dplyr::pull(probe_id)
    ))
  
  stopifnot(sum(is.na(mvalue)) == 0)
  stopifnot(colnames(mvalue) == c("probe_id", target_id))

  # cleanup 
  
  rm(proc) #detP,

  saveRDS(mvalue, file=paste0("cache/mvalues/",target_id,".Rds"))
  
  rm(mvalue, target, target_id)
}





# export ----


exp <- all_targets |>
  assertr::verify(file.exists(cache_mvalues)) |> 
  dplyr::pull(cache_mvalues) |> 
  pbapply::pblapply(readRDS) |> 
  purrr::reduce(function(x, y) dplyr::full_join(x, y, by = 'probe_id')) |> 
  tibble::column_to_rownames('probe_id') |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
    assertthat::assert_that(ncol(.) == nrow(all_targets))
    return(.)
  })()



saveRDS(exp, "cache/mvalues/mvalues_all.Rds")


rm(exp)


