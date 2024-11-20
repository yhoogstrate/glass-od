#!/usr/bin/env R 

# load ----



source('scripts/load_constants.R')
source('scripts/load_functions.R')
#source('scripts/load_palette.R')
#source('scripts/load_themes.R')



if(!exists('metadata.cg_probes.epic')) {
  print("First loading probe annotations")
  source('scripts/load_probe_annotations.R')
  print("Done")
}



# EPIC: all hq ----





#data.intensities.combined.hq_samples <- readRDS("cache/intensities_hq/intensities_hq.Rds") |> 
data.intensities.combined.hq_samples <- readRDS("cache/intensities_raw_hq/intensities_raw_hq.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_SAMPLES) # 4 replicates still need to be removed still
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
    return(.)
  })()


#data.intensities.methylated.hq_samples <- readRDS("cache/intensities_m_hq/intensities_m_hq.Rds") |> 
  data.intensities.methylated.hq_samples <- readRDS("cache/intensities_raw_m_hq/intensities_raw_m_hq.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_SAMPLES) # 4 replicates still need to be removed still
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
    return(.)
  })()


#data.intensities.unmethylated.hq_samples <- readRDS("cache/intensities_um_hq/intensities_um_hq.Rds") |> 
data.intensities.unmethylated.hq_samples <- readRDS("cache/intensities_raw_um_hq/intensities_raw_um_hq.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_SAMPLES) # 4 replicates still need to be removed still
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
    return(.)
  })()


data.intensities.mask.hq_samples <- readRDS("cache/masks_hq/masks_hq.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_SAMPLES) # 4 replicates still need to be removed still
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
    return(.)
  })()


stopifnot(colnames(data.intensities.combined.hq_samples) == colnames(data.intensities.mask.hq_samples))
stopifnot(rownames(data.intensities.combined.hq_samples) == rownames(data.intensities.mask.hq_samples))




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
  dplyr::filter(grepl("^cg", probe_id)) |> # non CG probes https://knowledge.illumina.com/microarray/general/microarray-general-troubleshooting-list/000005501
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
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


# glass_od.metadata.array_samples |> 
#   dplyr::filter(grepl("Valor",patient_center_name)) |> 
#   dplyr::select(array_sentrix_id, resection_id, patient_id, resection_number, patient_study_name, isolation_pathology_number , contains("note"), contains("reaso"), patient_suspected_noncodel) |> 
#   View()


