#!/usr/bin/env R

# GLASS-NL ----


source('scripts/load_functions.R')


## idats ----


glass_nl.metadata.idats <-  list.files(path = "data/GLASS_NL/Methylation/Methylation Array Data/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/GLASS_NL/Methylation/Methylation Array Data/", filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (658)) # no idea how there can be 658 idats in here?
    return(.)
  })() |>
  assertr::verify(file.exists(filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", filename)) |>
  dplyr::mutate(channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", filename)) |>
  tidyr::pivot_wider(id_cols = sentrix_id, names_from = channel, values_from = c(filename)) |>
  dplyr::rename(channel_green = Grn) |>
  dplyr::rename(channel_red = Red) |> 
  assertr::verify(!is.na(sentrix_id)) |> 
  assertr::verify(!is.na(channel_green)) |> 
  assertr::verify(!is.na(channel_red))



## study identifier: Sample_Name & resection metadata ----


tmp <- read.csv("data/GLASS_NL/Metadata/Samples/Master Datasheet_ALL METHODS_27012023.csv") |>
  dplyr::select( Surgery_ID, GLASS_ID,  Sample_Name, Sample_Sex, Sample_Type, Resectie, Sample_ID, Recurrent_Select.Meth, Matched_Pair.Meth) |> 
  dplyr::rename(sentrix_id = Sample_ID) |> 
  dplyr::rename(resection_number = Resectie) |> 
  dplyr::filter(!is.na(sentrix_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (235)) 
    return(.)
  })() |> 
  dplyr::mutate(patient_id = as.factor(gsub("^.+_[0]*([0-9]+)$", "\\1", GLASS_ID)))


glass_nl.metadata.idats <- glass_nl.metadata.idats |>
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |>
  dplyr::filter(!is.na(Sample_Name)) |>  # exclude hundreds present in directory but absent in metadata
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (235)) 
    return(.)
  })()


rm(tmp)




## WHO classes ----


tmp <- read.csv("data/GLASS_NL/Metadata/Samples/WHOclassification_03052022.csv")  |> 
  dplyr::select(Surgery_ID, WHO_Classification2021) |> # join on Surgery_ID
  dplyr::mutate(WHO_Classification2021 = ifelse(WHO_Classification2021 == "Therapy Effects", as.character(NA), WHO_Classification2021)) |> 
  assertr::verify(!is.na(Surgery_ID)) |> 
  dplyr::filter(Surgery_ID %in% glass_nl.metadata.idats$Surgery_ID) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (235)) 
    return(.)
  })()



glass_nl.metadata.idats <- glass_nl.metadata.idats |>
  dplyr::left_join(tmp, by=c('Surgery_ID'='Surgery_ID'), suffix=c('','')) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (235)) 
    return(.)
  })()


rm(tmp)


## Percentage detP probes ----
#' from: scripts/analysis_percentage_detP_probes.R

tmp <- read.table("output/tables/percentage_detP_probes.txt") |> 
  assertr::verify(!is.na(percentage.detP.signi) & is.numeric(percentage.detP.signi)) |> 
  assertr::verify(glass_nl.metadata.idats$sentrix_id %in% sentrix_id)

glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(percentage.detP.signi)) 

rm(tmp)



#' from: scripts/analysis_unsupervised_qc.R

tmp <- readRDS('cache/unsupervised_qc_outliers_all.Rds') |> 
  assertr::verify(!is.na(qc.pca.comp1)) |> 
  assertr::verify(!is.na(qc.pca.detP.outlier)) |> 
  assertr::verify(glass_nl.metadata.idats$sentrix_id %in% sentrix_id)


glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(qc.pca.comp1)) |> 
  assertr::verify(!is.na(qc.pca.detP.outlier))

rm(tmp)




## Heidelberg 12.8 reportBrain files ----


tmp <- c(
  list.files(path = "data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/",
             pattern = "*_scores_cal.csv", recursive = TRUE)) |>
  data.frame(filename = _) |>
  dplyr::mutate(mnp_predictBrain_filename = paste0("data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/", filename)) |>
  dplyr::mutate(mnp_predictBrain_version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", filename)) |>
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |>  # only one version per sentrix_id desired
  assertr::verify(sentrix_id %in% glass_nl.metadata.idats$sentrix_id) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_reportBrain_csv(mnp_predictBrain_filename, paste0("mnp_predictBrain_", mnp_predictBrain_version, "_"))) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 235)
    return(.)
  })() |> 
  assertr::verify(!is.na(mnp_predictBrain_v12.8_cal_class)) |>  # version is hardcoded here
  assertr::verify(!is.na(mnp_predictBrain_v12.8_cal_A_IDH_LG)) |>
  assertr::verify(!is.na(mnp_predictBrain_v12.8_cal_A_IDH_HG)) |>
  dplyr::mutate(A_IDH_HG__A_IDH_LG_lr = log(mnp_predictBrain_v12.8_cal_A_IDH_HG / mnp_predictBrain_v12.8_cal_A_IDH_LG)) |> 
  dplyr::mutate(A_IDH_HG__A_IDH_LG_lr_neat = log(
    (mnp_predictBrain_v12.8_cal_A_IDH_HG / (1-mnp_predictBrain_v12.8_cal_A_IDH_HG))
    / 
      (mnp_predictBrain_v12.8_cal_A_IDH_LG / (1-mnp_predictBrain_v12.8_cal_A_IDH_LG))
  ))



glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(mnp_predictBrain_v12.8_cal_class)) |> 
  assertr::verify(!is.na(A_IDH_HG__A_IDH_LG_lr)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 235)
    return(.)
  })()

rm(tmp)



## Heidelberg 12.8 Frozen ~ FFPE status ----


tmp <- list.files(path = "data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/",  pattern = "_ffpe_frozen.txt", recursive = TRUE) |> 
  data.frame(mnp_QC_FrozenFFPEstatus_table = _) |> 
  dplyr::mutate(mnp_QC_FrozenFFPEstatus_table = paste0("data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/", mnp_QC_FrozenFFPEstatus_table)) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_FrozenFFPEstatus_table(mnp_QC_FrozenFFPEstatus_table, "mnp_QC_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% glass_nl.metadata.idats$sentrix_id) |> 
  dplyr::mutate(mnp_QC_FrozenFFPEstatus_table = NULL) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 235)
    return(.)
  })()



glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(mnp_QC_predicted_array_type)) |> 
  assertr::verify(!is.na(mnp_QC_predicted_sample_type)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 235)
    return(.)
  })() 

rm(tmp)



## Heidelberg 12.8 CNVP segment files ----

tmp <- query_mnp_12.8_CNVP_segment_csv(
  "data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/",
  235, glass_nl.metadata.idats$sentrix_id)


glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(heidelberg_cnvp_segments))  |> 
  assertr::verify(!is.na(heidelberg_cnvp_version))
rm(tmp)




## QC PCA outlier ----


tmp <- readRDS("cache/unsupervised_qc_outliers_all.Rds")

glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(qc.pca.comp1)) |> 
  assertr::verify(!is.na(qc.pca.detP.outlier))

rm(tmp)



## RF purity calls ----


tmp <- read.table("data/GLASS_NL/Methylation/Analysis/RFpurity/purities_RFpurity.txt") |> 
  dplyr::mutate(sentrix_id = gsub("^.+/([0-9]{10,}_R[0-9]+C[0-9]+).+\\.idat$","\\1",fn)) |> 
  dplyr::select(sentrix_id, absolute,  estimate) |> 
  assertr::verify(!is.na(sentrix_id)) |> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(glass_nl.metadata.idats$sentrix_id %in% sentrix_id) |> 
  dplyr::rename(RFpurity.absolute = absolute) |> 
  dplyr::rename(RFpurity.estimate = estimate)


glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(RFpurity.absolute)) |> 
  assertr::verify(!is.na(RFpurity.estimate))

rm(tmp)




## QC Purity PCA outlier ----
#' PCA was performed on only 218 HQ arrays
#' these outliers in PC3 (val >300) were CONTR_ and sometimes difficultly classifiable

# PC3_low_purity_samples <- c("203175700013_R08C01", "203189480016_R03C01", "203189480016_R05C01", "203986510092_R04C01", "203986510125_R04C01", "203989100024_R08C01",
#                             "203989100035_R02C01", "203989100096_R03C01", "203989100096_R05C01", "203989100142_R02C01", "203991400003_R01C01", "203991400003_R06C01",
#                             "204073520032_R07C01", "204073570005_R02C01", "203519500055_R03C01", "203519500055_R05C01")
# 
# glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
#   dplyr::mutate(qc.pca.pc3purity.outlier = sentrix_id %in% PC3_low_purity_samples) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == 235)
#     return(.)
#   })()
# 
# rm(PC3_low_purity_samples)
# 
# 

## ++ below: re-build because mvalue normalisation ++ ----

## Median methylation levels ----


tmp <- readRDS("cache/analysis_median_methylation.Rds")

glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(
    (qc.pca.detP.outlier == F & !is.na(median.overall.methylation)) |
    (qc.pca.detP.outlier == T & is.na(median.overall.methylation))
  ) |>
  assertr::verify(
    (qc.pca.detP.outlier == F & !is.na(median.glass_nl_supervised.methylation)) |
    (qc.pca.detP.outlier == T & is.na(median.glass_nl_supervised.methylation))
  )

rm(tmp)


## A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV ----


tmp <- readRDS(file="cache/analysis_A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV.Rds") |> 
  assertr::verify(!is.na(A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV)) |> 
  dplyr::filter(sentrix_id %in% glass_nl.metadata.idats$sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 218)
    return(.)
  })()


glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))


## unsupervised PCA ----


tmp <- readRDS(file="cache/analysis_unsupervised_PCA_GLASS-NL_x.Rds") |> 
  assertr::verify(!is.na(PC1)) |> 
  assertr::verify(!is.na(PC2)) |> 
  assertr::verify(!is.na(PC3)) |> 
  assertr::verify(!is.na(PC4)) |> 
  assertr::verify(!is.na(PC5)) |> 
  assertr::verify(!is.na(PC218)) |> 
  dplyr::filter(sentrix_id %in% glass_nl.metadata.idats$sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 218)
    return(.)
  })()


glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))




## unsupervised PCA [GLASS-OD + GLASS-NL combi] ----


# tmp <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_GLASS-NL_combined.Rds") |> 
#   dplyr::rename_with(~ gsub("^PC","PC.GLASS_OD_NL_combined.",.x), .cols = matches("^PC[0-9]", perl = T)) |> 
#   dplyr::filter(sentrix_id %in% glass_nl.metadata.idats$sentrix_id) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == 218)
#     return(.)
#   })()
# 
# 
# glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
#   dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))
# rm(tmp)

tmp <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_GLASS-NL_combined_no_1P19Q.Rds") |>
  dplyr::rename_with(~ gsub("^PC","PC.GLASS_OD_NL_combined_excl_1P19Q.",.x), .cols = matches("^PC[0-9]", perl = T)) |>
  dplyr::filter(sentrix_id %in% glass_nl.metadata.idats$sentrix_id) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 218)
    return(.)
  })()


glass_nl.metadata.idats <- glass_nl.metadata.idats |>
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))
rm(tmp)




