#!/usr/bin/env R 

# load ----


source('scripts/load_functions.R')


if(!exists('metadata.cg_probes.epic')) {
  source('scripts/load_probe_annotations.R')
}



# all hq ----


data.beta.values.hq_samples <- readRDS("cache/beta.values.HQ_samples.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == (163+218+73))
    assertthat::assert_that(nrow(.) == (760405))
    return(.)
  })()

data.beta.values.mask.hq_samples <- readRDS("cache/detP_masked_values.HQ_samples.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == (163+218+73))
    assertthat::assert_that(nrow(.) == (760405))
    return(.)
  })()

stopifnot(colnames(data.beta.values.hq_samples) == colnames(data.beta.values.mask.hq_samples))
stopifnot(rownames(data.beta.values.hq_samples) == rownames(data.beta.values.mask.hq_samples))


## probe table - mask / detP counts ----


data.beta.values.probes <- data.beta.values.mask.hq_samples |> 
  is.na() |> 
  rowSums() |> 
  as.data.frame() |> 
  dplyr::rename(n_na = 1) |> 
  dplyr::mutate(detP_good_probe = n_na == 0) |> 
  tibble::rownames_to_column('probe_id')


## probe table - probe annotation ----


data.beta.values.probes <- data.beta.values.probes |> 
  dplyr::left_join(metadata.cg_probes.epic, by=c('probe_id'='probe_id'), suffix=c('','')) |>  # + annotation of ALL probes
  assertr::verify(MASK_general == F)


## filter for good probes ----


data.beta.values.good_probes <- data.beta.values.probes |> 
  dplyr::filter(detP_good_probe) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (694299))
    return(.)
  })() |> 
  assertr::verify(MASK_general == F) |> 
  dplyr::pull(probe_id)



