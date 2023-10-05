#!/usr/bin/env R 

# load ----


source('scripts/load_functions.R')


if(!exists('metadata.cg_probes.epic')) {
  source('scripts/load_probe_annotations.R')
}



# all hq ----


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




