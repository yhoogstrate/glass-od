#!/usr/bin/env R 

# load ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')
#source('scripts/load_palette.R')
#source('scripts/load_themes.R')


if(!exists('metadata.cg_probes.epic')) {
  source('scripts/load_probe_annotations.R')
}



# all hq ----


data.beta.values.hq_samples <- readRDS("cache/betavalues_hq/betavalues_hq.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_SAMPLES)
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
    return(.)
  })()

data.beta.values.mask.hq_samples <- readRDS("cache/masks_hq/masks_hq.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_SAMPLES)
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
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
  dplyr::filter(grepl("^cg", probe_id)) |> # non CG probes https://knowledge.illumina.com/microarray/general/microarray-general-troubleshooting-list/000005501
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
    return(.)
  })() |> 
  assertr::verify(MASK_general == F) |> 
  dplyr::pull(probe_id)




