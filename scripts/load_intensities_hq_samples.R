#!/usr/bin/env R 

# load ----


source('scripts/load_functions.R')


if(!exists('metadata.cg_probes.epic')) {
  print("First loading probe annotations")
  source('scripts/load_probe_annotations.R')
  print("Done")
}



# EPIC: all hq ----


data.intensities.hq_samples <- readRDS("cache/intensities.HQ_samples.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == (505)) # 4 replicates still need to be removed still
    assertthat::assert_that(nrow(.) == (760405))
    return(.)
  })()

data.intensities.mask.hq_samples <- readRDS("cache/detP_masked_values.HQ_samples.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == (505)) # 4 replicates still need to be removed still
    assertthat::assert_that(nrow(.) == (760405))
    return(.)
  })()

stopifnot(colnames(data.intensities.hq_samples) == colnames(data.intensities.mask.hq_samples))
stopifnot(rownames(data.intensities.hq_samples) == rownames(data.intensities.mask.hq_samples))




## probe table - mask / detP counts ----


data.intensities.probes <- data.intensities.mask.hq_samples |> 
  is.na() |> 
  rowSums() |> 
  as.data.frame() |> 
  dplyr::rename(n_na = 1) |> 
  dplyr::mutate(detP_good_probe = n_na == 0) |> 
  tibble::rownames_to_column('probe_id')



## probe table - probe annotation ----


data.intensities.probes <- data.intensities.probes |> 
  dplyr::left_join(metadata.cg_probes.epic, by=c('probe_id'='probe_id'), suffix=c('','')) |>  # + annotation of ALL probes
  assertr::verify(MASK_general == F)



## filter for good probes ----


data.intensities.good_probes <- data.intensities.probes |> 
  dplyr::filter(detP_good_probe) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693060))
    return(.)
  })() |> 
  assertr::verify(MASK_general == F) |> 
  dplyr::pull(probe_id)



## DMP outcomes ----

### prim rec ----


# fn <- "cache/analysis_differential__primary_recurrence__partial_paired_nc__stats.Rds"
# if(file.exists(fn)) {
#   data.intensities.probes <- data.intensities.probes |> 
#     dplyr::left_join(
#       readRDS(fn) |> 
#         dplyr::rename_with(~paste0("DMP__primary_recurrence__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)),
#       by=c('probe_id'='probe_id'), suffix=c('','') )
# } else {
#   warning("DMP result primary - recurrence is missing")
# }
# 
# rm(fn)


### g3 g2 ----


# fn <- "cache/analysis_differential__g2_g3__partial_paired_nc__stats.Rds"
# if(file.exists(fn)) {
#   data.intensities.probes <- data.intensities.probes |> 
#     dplyr::left_join(
#       readRDS(fn) |> 
#         dplyr::rename_with(~paste0("DMP__g2_g3__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)),
#       by=c('probe_id'='probe_id'), suffix=c('','') )
# } else {
#   warning("DMP result g2 - g3 is missing")
# }
# 
# rm(fn)

