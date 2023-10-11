#!/usr/bin/env R 

# load ----


source('scripts/load_functions.R')


if(!exists('metadata.cg_probes.epic')) {
  source('scripts/load_probe_annotations.R')
}



# EPIC: all hq ----


data.mvalues.hq_samples <- readRDS("cache/mvalues.HQ_samples.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == (510)) # 4 replicates still need to be removed still
    assertthat::assert_that(nrow(.) == (760405))
    return(.)
  })()

data.mvalues.mask.hq_samples <- readRDS("cache/detP_masked_values.HQ_samples.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == (510)) # 4 replicates still need to be removed still
    assertthat::assert_that(nrow(.) == (760405))
    return(.)
  })()

stopifnot(colnames(data.mvalues.hq_samples) == colnames(data.mvalues.mask.hq_samples))
stopifnot(rownames(data.mvalues.hq_samples) == rownames(data.mvalues.mask.hq_samples))




## probe table - mask / detP counts ----


data.mvalues.probes <- data.mvalues.mask.hq_samples |> 
  is.na() |> 
  rowSums() |> 
  as.data.frame() |> 
  dplyr::rename(n_na = 1) |> 
  dplyr::mutate(detP_good_probe = n_na == 0) |> 
  tibble::rownames_to_column('probe_id')



## probe table - probe annotation ----


data.mvalues.probes <- data.mvalues.probes |> 
  dplyr::left_join(metadata.cg_probes.epic, by=c('probe_id'='probe_id'), suffix=c('','')) |>  # + annotation of ALL probes
  assertr::verify(MASK_general == F)


## filter for good probes ----


data.mvalues.good_probes <- data.mvalues.probes |> 
  dplyr::filter(detP_good_probe) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017))
    return(.)
  })() |> 
  assertr::verify(MASK_general == F) |> 
  dplyr::pull(probe_id)




# 450K ----
## AD: from Beta-value exported files ----



if(!exists('ad_bmc_clin_epi.metadata.array_samples')) {
  source('scripts/load_AD_BMC_Clin_Epi.R')
}



tmp.1 <- read.delim("data/GLASS_OD/DNA Methylation - 450K arrays - Alzheimer/SampleMethFinalReport_ctrl-bkg.txt", skip=8, header=T) |> 
  tibble::column_to_rownames("TargetID") |> 
  dplyr::select(contains("_Beta")) |> 
  dplyr::mutate_all(function(x){return (  log2( `x` / (1 - `x`))  )}) |>  # Beta to M-values
  tibble::rownames_to_column("probe_id")



tmp.2 <- read.delim("data/GLASS_OD/DNA Methylation - 450K arrays - Alzheimer/SampleMethFinalReport_nonorm.txt", skip=8, header=T) |> 
  tibble::column_to_rownames("TargetID") |> 
  dplyr::select(contains("_Beta")) |> 
  dplyr::mutate_all(function(x){return (  log2( `x` / (1 - `x`))  )}) |>  # Beta to M-values
  tibble::rownames_to_column("probe_id")


stopifnot(tmp.1$probe_id == tmp.2$probe_id)
stopifnot(nrow(tmp.1) == nrow(tmp.2))


data.mvalues.alzheimer.dirty <- tmp.1.m |> 
  dplyr::left_join(tmp.2.m, by=c('probe_id'='probe_id'), suffix=c('','')) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::rename_with(~ gsub("\\.AVG_Beta$","",.x)) |> 
  dplyr::rename_with(~ gsub("^X","",.x)) |> 
  dplyr::rename_with(~ gsub("\\.","_",.x)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == (26 + 12))
    assertthat::assert_that(all(colnames(.) %in% tmp$DNAm_id))
    
    return(.)
  })() |> 
  dplyr::select(ad_bmc_clin_epi.metadata.array_samples |> dplyr::filter(is.na(reason_excluded)) |>  dplyr::pull(DNAm_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == (26 + 12 - 1))
    
    return(.)
  })()


rm(tmp.1, tmp.2)





