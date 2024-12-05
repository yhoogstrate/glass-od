#!/usr/bin/env R

# load libs, config & db ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')



# normal cortex samples ----
## idat(s) ----


normal_cortex.metadata.array_samples <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - Normal Brain/ctrl_brain_PMC6355837_GSE111165/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |> 
  data.frame(array_filename = _) |>
  dplyr::mutate(array_filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - Normal Brain/ctrl_brain_PMC6355837_GSE111165/", array_filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 28)
    return(.)
  })() |> 
  assertr::verify(file.exists(array_filename)) |>
  dplyr::mutate(array_sentrix_id = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", array_filename)) |>
  dplyr::mutate(array_sentrix_id = gsub("GSM[0-9]+_", "", array_sentrix_id)) |>
  dplyr::mutate(array_channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", array_filename)) |>
  tidyr::pivot_wider(id_cols = array_sentrix_id, names_from = array_channel, values_from = c(array_filename)) |>
  dplyr::rename(array_channel_green = Grn) |>
  dplyr::rename(array_channel_red = Red) |>
  dplyr::mutate(array_channel_green_filesize = file.info(array_channel_green)$size) |> 
  dplyr::mutate(array_channel_red_filesize   = file.info(array_channel_red)$size) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 14)
    return(.)
  })() |>
  assertr::verify(!is.na(array_channel_green)) |>
  assertr::verify(!is.na(array_channel_red)) |> 
  dplyr::mutate(array_control_id = dplyr::case_when(array_sentrix_id == '200392810022_R03C01' ~ 's107',
                                                    array_sentrix_id == '200392810022_R04C01' ~ 's108',
                                                    array_sentrix_id == '200607090086_R07C01' ~ 's110',
                                                    array_sentrix_id == '200392810089_R01C01' ~ 's112',
                                                 T ~ as.character(NA)))



## mnp probabilities ----



tmp.ls <- c(list.files(
  path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/non-malignant-cortex-and-spike-ins/", pattern = "*_scores_cal.csv", recursive = TRUE)
)


tmp <- tmp.ls |>
  data.frame(tmp_filename = _) |>
  dplyr::mutate(array_mnp_predictBrain_filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/non-malignant-cortex-and-spike-ins/", tmp_filename)) |>
  dplyr::mutate(array_mnp_predictBrain_version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", tmp_filename)) |>
  dplyr::mutate(array_sentrix_id = gsub("^.+[^0-9]([0-9]{10,12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_filename)) |>
  dplyr::mutate(tmp_uniq = paste0(array_sentrix_id, array_mnp_predictBrain_version)) |> 
  assertr::verify(!duplicated(tmp_uniq)) |>
  dplyr::mutate(tmp_uniq = NULL) |> 
  dplyr::mutate(tmp_filename = NULL) |> 
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |>  # only one version per sentrix_id 
  dplyr::filter(array_sentrix_id %in% normal_cortex.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 14)
    return(.)
  })() |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_reportBrain_csv(array_mnp_predictBrain_filename, paste0("array_mnp_predictBrain_", array_mnp_predictBrain_version, "_"))) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_class)) |>  # mnp version is hardcoded here
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_A_IDH_LG)) |>
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_A_IDH_HG)) |>
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_OLIGOSARC_IDH))



normal_cortex.metadata.array_samples <- normal_cortex.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_class)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 14)
    return(.)
  })()



rm(tmp, tmp.ls)



## mnp cnv bins ----


tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/non-malignant-cortex-and-spike-ins/", pattern = "*.bins.igv", recursive = TRUE) |> 
  data.frame(array_heidelberg_cnvp_bins = _) |> 
  dplyr::mutate(array_mnp_CNVP_v12.8_v5.2_CNVP_bins = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/non-malignant-cortex-and-spike-ins/", array_heidelberg_cnvp_bins)) |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+[^0-9]([0-9]{10,12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", array_heidelberg_cnvp_bins)) |> 
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  dplyr::filter(array_sentrix_id %in% normal_cortex.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 14)
    return(.)
  })()


normal_cortex.metadata.array_samples <- normal_cortex.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_heidelberg_cnvp_bins)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 14)
    return(.)
  })()


rm(tmp)



## Heidelberg 12.8 bin-based tumor purity calls ----


tmp <- readRDS("cache/analysis_tumor_purity_EPIC_bin-based__normal_cortex.Rds") |> 
  assertr::verify(!is.na(array_methylation_bins_1p19q_median.lfc)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_median.lfc)) |> 
  assertr::verify(!is.na(array_methylation_bins_1p19q_purity)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_purity)) |> 
  assertr::verify(!is.na(array_methylation_bins_1p19q_sd)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_sd)) |> 
  dplyr::filter(array_sentrix_id %in% normal_cortex.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 14)
    return(.)
  })()


normal_cortex.metadata.array_samples <- normal_cortex.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)


## sandbox stats ----


normal_cortex.metadata.array_samples |>
  dplyr::filter(array_control_id == "s107") |>
  dplyr::select(contains("_cal_")) |>
  dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = NULL) |> 
  t() |> 
  as.data.frame() |> 
  dplyr::mutate(V1 = round(V1,3)) |> 
  dplyr::arrange(-V1) |> 
  head(n=3)



normal_cortex.metadata.array_samples |>
  dplyr::filter(array_control_id == "s108") |>
  dplyr::select(contains("_cal_")) |>
  dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = NULL) |> 
  t() |> 
  as.data.frame() |> 
  dplyr::mutate(V1 = round(V1,3)) |> 
  dplyr::arrange(-V1) |> 
  head(n=3)



normal_cortex.metadata.array_samples |>
  dplyr::filter(array_control_id == "s110") |>
  dplyr::select(contains("_cal_")) |>
  dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = NULL) |> 
  t() |> 
  as.data.frame() |> 
  dplyr::mutate(V1 = round(V1,3)) |> 
  dplyr::arrange(-V1) |> 
  head(n=3)



normal_cortex.metadata.array_samples |>
  dplyr::filter(array_control_id == "s112") |>
  dplyr::select(contains("_cal_")) |>
  dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = NULL) |> 
  t() |> 
  as.data.frame() |> 
  dplyr::mutate(V1 = round(V1,3)) |> 
  dplyr::arrange(-V1) |> 
  head(n=3)





# dilution series ----
## mnp probabilities ----


tmp.ls <- c(list.files(
  path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/non-malignant-cortex-and-spike-ins/", pattern = "*_scores_cal.csv", recursive = TRUE)
)


cortex_dilution.metadata.array_samples <- tmp.ls |>
  data.frame(tmp_filename = _) |>
  dplyr::mutate(array_mnp_predictBrain_filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/non-malignant-cortex-and-spike-ins/", tmp_filename)) |>
  dplyr::mutate(array_mnp_predictBrain_version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", tmp_filename)) |>
  dplyr::mutate(array_sentrix_id = gsub("^.+[^0-9]([0-9]{10,12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_filename)) |>
  dplyr::mutate(tmp_uniq = paste0(array_sentrix_id, array_mnp_predictBrain_version)) |> 
  assertr::verify(!duplicated(tmp_uniq)) |>
  dplyr::mutate(tmp_uniq = NULL) |> 
  dplyr::mutate(tmp_filename = NULL) |> 
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |>  # only one version per sentrix_id 
  dplyr::filter(array_sentrix_id %in% normal_cortex.metadata.array_samples$array_sentrix_id == F) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 160)
    return(.)
  })() |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_reportBrain_csv(array_mnp_predictBrain_filename, paste0("array_mnp_predictBrain_", array_mnp_predictBrain_version, "_"))) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_class)) |>  # mnp version is hardcoded here
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_A_IDH_LG)) |>
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_A_IDH_HG)) |>
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_OLIGOSARC_IDH)) |> 
  dplyr::mutate(array_tumor_id = gsub("^.+/28([0-9]{3})[0-9]+_R0([0-9]+)C[0-9]+_scores_cal.csv$", "0\\1-R\\2", array_mnp_predictBrain_filename)) |> 
  dplyr::mutate(array_fraction_normal_brain = as.numeric(gsub("^.+/28[0-9]{4}([0-9]{2})[0-9]+_R[0-9]+C[0-9]+_scores_cal.csv$", "\\1", array_mnp_predictBrain_filename))) |> 
  dplyr::mutate(array_control_id = gsub(".+/28[0-9]{7}([0-9]+)_R0[0-9]+C[0-9]+_scores_cal.csv$", "s\\1", array_mnp_predictBrain_filename))
  #dplyr::select(array_sentrix_id, array_tumor_id, array_fraction_normal_brain, array_control_id)



rm(tmp.ls)



## mnp cnv bins ----


tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/non-malignant-cortex-and-spike-ins/", pattern = "*.bins.igv", recursive = TRUE) |> 
  data.frame(array_heidelberg_cnvp_bins = _) |> 
  dplyr::mutate(array_mnp_CNVP_v12.8_v5.2_CNVP_bins = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/non-malignant-cortex-and-spike-ins/", array_heidelberg_cnvp_bins)) |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+[^0-9]([0-9]{10,12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", array_heidelberg_cnvp_bins)) |> 
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  dplyr::filter(array_sentrix_id %in% cortex_dilution.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 160)
    return(.)
  })()


cortex_dilution.metadata.array_samples <- cortex_dilution.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_heidelberg_cnvp_bins)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 160)
    return(.)
  })()


rm(tmp)



## Heidelberg 12.8 bin-based tumor purity calls ----


tmp <- readRDS("cache/analysis_tumor_purity_EPIC_bin-based__normal_cortex.Rds") |> 
  assertr::verify(!is.na(array_methylation_bins_1p19q_median.lfc)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_median.lfc)) |> 
  assertr::verify(!is.na(array_methylation_bins_1p19q_purity)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_purity)) |> 
  assertr::verify(!is.na(array_methylation_bins_1p19q_sd)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_sd)) |> 
  dplyr::filter(array_sentrix_id %in% cortex_dilution.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 160)
    return(.)
  })()


cortex_dilution.metadata.array_samples <- cortex_dilution.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)




## actual used spike-in percentages ----


tmp <- readLines("scripts/analysis_mix_idats.sh") |> 
  as.data.frame() |> 
  dplyr::rename(V1 = 1) |> 
  dplyr::filter(!grepl("^#", V1)) |> 
  dplyr::filter(grepl("idat-tools[ ]mix", V1)) |> 
  dplyr::mutate(array_fraction_normal_brain_cli = as.numeric(gsub("^.+ -r ([0-9\\.]+) .+$","\\1",V1)) * 100.0) |> 
  dplyr::mutate(sentrix_id_tumor = gsub('^.+/([0-9]+_[RC0-9]+)_.+.+/([A-Z0-9]+_[0-9]+_[RC0-9]+)_.+$','\\1',V1)) |> 
  dplyr::mutate(sentrix_id_control = gsub('^.+/([0-9]+_[RC0-9]+)_.+.+/([A-Z0-9]+_[0-9]+_[RC0-9]+)_.+$','\\2',V1)) |> 
  dplyr::mutate(array_sentrix_id = gsub('^.+/([0-9]+_[RC0-9]+)_.+.+/([A-Z0-9]+_[0-9]+_[RC0-9]+)_.+/([0-9]+_[RC0-9]+)_.+$','\\3',V1)) |> 
  dplyr::mutate(V1 = NULL) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 160 * 2)
    return(.)
  })() |> 
  dplyr::distinct() |> # collapse green and red channels
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 160)
    return(.)
  })() |> 
  assertr::verify(array_sentrix_id %in% cortex_dilution.metadata.array_samples$array_sentrix_id) |> 
  dplyr::mutate(sentrix_id_tumor = NULL) |> 
  dplyr::mutate(sentrix_id_control = NULL)


cortex_dilution.metadata.array_samples <- cortex_dilution.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)



## relevel ----


cortex_dilution.metadata.array_samples <- cortex_dilution.metadata.array_samples |> 
  dplyr::mutate(array_tumor_id = factor(array_tumor_id, levels=c("0017-R3", "0008-R2", "0121-R3", "0054-R3"))) |> 
  dplyr::mutate(array_control_id = factor(array_control_id, levels=c("s107", "s108", "s110", "s112"))) |> 
  dplyr::arrange(array_tumor_id, array_control_id, array_fraction_normal_brain) 




## sandbox ----


tmp <- normal_cortex.metadata.array_samples |> 
  dplyr::filter(!is.na(array_control_id)) |> 
  tibble::column_to_rownames('array_control_id') |> 
  dplyr::select(contains("array_mnp_predictBrain_v12.8_cal_")) |> 
  dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = NULL) |> 
  round(3) |> 
  t()


tmp |>
  as.data.frame() |> 
  dplyr::select(s107) |> 
  dplyr::arrange(-s107) |> 
  head(n=5)


tmp |>
  as.data.frame() |> 
  dplyr::select(s108) |> 
  dplyr::arrange(-s108) |> 
  head(n=5)


tmp |>
  as.data.frame() |> 
  dplyr::select(s110) |> 
  dplyr::arrange(-s110) |> 
  head(n=5)


tmp |>
  as.data.frame() |> 
  dplyr::select(s112) |> 
  dplyr::arrange(-s112) |> 
  head(n=5)



cortex_dilution.metadata.array_samples |>
  dplyr::select(array_tumor_id, array_control_id, array_fraction_normal_brain, array_fraction_normal_brain_cli) |> 
  View()


