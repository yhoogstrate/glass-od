#!/usr/bin/env R

# load data ----


if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}


data <- data.mvalues.hq_samples |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == (456-2))
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
    assertthat::assert_that(nrow(.) == (694299))
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




# clean-up ----



# merge & export ----


export <- median.overall.meth |> 
  dplyr::left_join(median.glass_nl.meth, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))


saveRDS(export, file="cache/analysis_median_methylation.Rds")



rm(data, export, targets)
gc()


# sandbox test plot ----
# 
# if(!exists('glass_od.metadata.idats')) {
#   source('scripts/load_GLASS-OD_metadata.R')
# }
# 
# if(!exists('glass_nl.metadata.idats')) {
#   source('scripts/load_GLASS-NL_metadata.R')
# }
# 
# if(!exists('gsam.metadata.idats')) {
#   source('scripts/load_G-SAM_metadata.R')
# }
# 
# 
# 
# 
# plt <- readRDS('cache/analysis_median_methylation.Rds') |>
#   dplyr::mutate(dataset = dplyr::case_when(
#     sentrix_id %in% glass_od.metadata.idats$sentrix_id ~ "GLASS-OD",
#     sentrix_id %in% glass_nl.metadata.idats$sentrix_id ~ "GLASS-NL",
#     sentrix_id %in% gsam.metadata.idats$sentrix_id ~ "G-SAM",
#     T ~ "error"
#   )) 
# 
# 
# 
# ggplot(plt, aes(x=median.overall.methylation, y=median.glass_nl_supervised.methylation, col=dataset)) +
#   geom_point()
# 
# 
# 
# 
# plt.glnl <- readRDS('cache/analysis_median_methylation.Rds') |>
#   dplyr::mutate(dataset = dplyr::case_when(
#     sentrix_id %in% glass_od.metadata.idats$sentrix_id ~ "GLASS-OD",
#     sentrix_id %in% glass_nl.metadata.idats$sentrix_id ~ "GLASS-NL",
#     sentrix_id %in% gsam.metadata.idats$sentrix_id ~ "G-SAM",
#     T ~ "error"
#   )) |> 
#   #dplyr::filter(dataset == "GLASS-NL") |> 
#   dplyr::left_join(
#     glass_nl.metadata.idats |>
#       dplyr::select(sentrix_id, Sample_Type, mnp_predictBrain_v12.8_cal_class, WHO_Classification2021),
#     by=c('sentrix_id'='sentrix_id'),suffix=c('','')
#   )
# 
# ggplot(plt.glnl, aes(x=median.overall.methylation, y=median.glass_nl_supervised.methylation, col=dataset)) +
#   geom_point()
# 
# 
