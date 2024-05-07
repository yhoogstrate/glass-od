#!/usr/bin/env R

# load data ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')
#source('scripts/load_palette.R')
#source('scripts/load_themes.R')


if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}


data <- data.mvalues.hq_samples |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_SAMPLES)
    return(.)
  })()
data.mask <- data.mvalues.mask.hq_samples


stopifnot(nrow(data) == nrow(data.mask))
stopifnot(ncol(data) == ncol(data.mask))


rm(data.mvalues.hq_samples)
rm(data.mvalues.mask.hq_samples)


gc()




data <- data |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
    return(.)
  })()

rm(data.mask)
gc()


# median meth all probes ----


median.overall.meth <- data |> 
  t() |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  as.data.frame() |>
  dplyr::rename(median.overall.methylation = 1) |> 
  assertr::verify(!is.na(median.overall.methylation)) |> 
  tibble::rownames_to_column('sentrix_id')



# median meth GLASS-NL probes ----


median.glass_nl.meth <- data |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% (read.csv('data/GLASS_NL/Metadata/(Epi)genetic_data/ProbeSelection_IvR_FDR 1e-9 Delta 1_06072022.csv') |> dplyr::pull(Probe_ID))) |> 
  tibble::column_to_rownames('probe_id') |> 
  t() |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  as.data.frame() |>
  dplyr::rename(median.glass_nl_supervised.methylation = 1) |> 
  assertr::verify(!is.na(median.glass_nl_supervised.methylation)) |> 
  tibble::rownames_to_column('sentrix_id')




# merge & export ----


export <- median.overall.meth |> 
  dplyr::left_join(median.glass_nl.meth, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_SAMPLES)
    return(.)
  })()


saveRDS(export, file="cache/analysis_median_methylation.Rds")



rm(data, export)
gc()



