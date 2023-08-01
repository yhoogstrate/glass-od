#!/usr/bin/env R

# G-SAM ----


source('scripts/load_functions.R')


## get idats ----

gsam.metadata.idats <-  list.files(path = "data/G-SAM/DNA Methylation - EPIC arrays/", pattern = "_(Grn|Red).idat$", recursive = T) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/G-SAM/DNA Methylation - EPIC arrays/", filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 182)# just freeze the number to terminate on unexpected behavior
    return(.)
  })() |> 
  assertr::verify(file.exists(filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", filename)) |>
  dplyr::mutate(channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", filename)) |>
  tidyr::pivot_wider(id_cols = sentrix_id, names_from = channel, values_from = c(filename)) |>
  dplyr::rename(channel_green = Grn) |>
  dplyr::rename(channel_red = Red) |> 
  dplyr::filter(!grepl("/MET2017-126-014/", channel_green)) |> # stored there for historical reasons - IDH-mutant loss study
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 75)
    return(.)
  })() 



## link sample names ----


tmp <- read.csv("data/G-SAM/DNA Methylation - EPIC arrays/MET2022-350-014/MET2022-350-014_IdH.csv", skip=8) |> 
  dplyr::filter(!is.na(Sentrix_ID)) |> 
  assertr::verify(grepl("^[0-9]{12}_[A-Z][0-9]{2}[A-Z][0-9]{2}$", Column2)) |> 
  dplyr::rename(sentrix_id = Column2) |> 
  dplyr::mutate(study = gsub("^(....).+$","\\1",Sample_Name)) |> 
  dplyr::filter(study %in% c("MINT","GLSO") == F) |> 
  assertr::verify(study == "GSAM") |> 
  dplyr::select(sentrix_id, Sample_Name) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == nrow(gsam.metadata.idats)) # == 75
    return(.)
  })() |> 
  assertr::verify(sentrix_id %in% gsam.metadata.idats$sentrix_id) |> 
  dplyr::mutate(resection = paste0("R",gsub("GSAM_...(.)_.+$","\\1",Sample_Name))) |> 
  dplyr::mutate(patient = gsub("GSAM_(...)._.+$","\\1", Sample_Name)) |> 
  dplyr::mutate(IDH = patient %in% c("BAW","CAV","CBG","CDF","DAB","EAF","EBD","ECB","FAD","FAL","JAB","JAD","JAF","KAC")) # IDH mut according to "data/gsam/output/tables/dna/idh_mutations.txt" - EAF also by MNP brain Classifier





gsam.metadata.idats <- gsam.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(Sample_Name))
rm(tmp)





## Add survival data ----


tmp <- read.csv('data/G-SAM/administratie/GSAM_combined_clinical_molecular.csv',stringsAsFactors=F) |> 
  dplyr::rename(patient = studyID) |> 
  dplyr::arrange(patient) |> 
  dplyr::mutate(initialMGMT = NULL) |> 
  dplyr::mutate(gender = ifelse(patient %in% c('AAT', 'AAM', 'AZH', 'HAI', 'FAG'),"Male",gender)) |>  # there's a number of samples of which the gender does not fit with the omics data - omics data determined genders are the corrected ones
  dplyr::mutate(gender = as.factor(gender)) |> 
  dplyr::mutate(survival.events = dplyr::case_when(
    status == "Deceased" ~ 1,
    status == "Censored" ~ 0,
    T ~ as.numeric(NA))) |> 
  dplyr::rename(os.event = survival.events) |> 
  dplyr::mutate(survival.months = survivalDays / 365.0 * 12.0) |> 
  dplyr::select(patient, gender, 
                status, os.event,
                survivalDays, 
                progressionFreeDays,
                survivalFromSecondSurgeryDays
                ) |> 
  dplyr::filter(patient %in% gsam.metadata.idats$patient) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == length(unique(gsam.metadata.idats$patient)))
    return(.)
  })()


gsam.metadata.idats <- gsam.metadata.idats |> 
  dplyr::left_join(tmp, by=c('patient'='patient'), suffix=c('',''))
rm(tmp)





## Percentage detP probes ----
#' from: scripts/analysis_percentage_detP_probes.R & scripts/analysis_unsupervised_qc.R


tmp <- readRDS('cache/unsupervised_qc_outliers_all.Rds') |> 
  assertr::verify(!is.na(qc.pca.comp1)) |> 
  assertr::verify(!is.na(qc.pca.detP.outlier)) |> 
  assertr::verify(gsam.metadata.idats$sentrix_id %in% sentrix_id)


gsam.metadata.idats <- gsam.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(qc.pca.comp1)) |> 
  assertr::verify(!is.na(qc.pca.detP.outlier))

rm(tmp)



## Heidelberg 12.8 reportBrain files ----


tmp <- c(
  list.files(path = "data/G-SAM/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/",
             pattern = "*_scores_cal.csv", recursive = TRUE)) |>
  data.frame(filename = _) |>
  dplyr::mutate(mnp_predictBrain_filename = paste0("data/G-SAM/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", filename)) |>
  dplyr::mutate(mnp_predictBrain_version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", filename)) |>
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |>  # only one version per sentrix_id desired
  assertr::verify(sentrix_id %in% gsam.metadata.idats$sentrix_id) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_reportBrain_csv(mnp_predictBrain_filename, paste0("mnp_predictBrain_", mnp_predictBrain_version, "_"))) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 75)
    return(.)
  })() |> 
  assertr::verify(!is.na(mnp_predictBrain_v12.8_cal_class)) |>  # version is hardcoded here
  assertr::verify(!is.na(mnp_predictBrain_v12.8_cal_A_IDH_LG)) |>
  assertr::verify(!is.na(mnp_predictBrain_v12.8_cal_A_IDH_HG)) |>
  dplyr::mutate(A_IDH_HG__A_IDH_LG_lr = log(mnp_predictBrain_v12.8_cal_A_IDH_HG / mnp_predictBrain_v12.8_cal_A_IDH_LG)) 



gsam.metadata.idats <- gsam.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(mnp_predictBrain_v12.8_cal_class)) |> 
  assertr::verify(!is.na(A_IDH_HG__A_IDH_LG_lr)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 75)
    return(.)
  })() 

rm(tmp)



## Heidelberg 12.8 Frozen ~ FFPE status ----


tmp <- list.files(path = "data/G-SAM/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/",  pattern = "_ffpe_frozen.txt", recursive = TRUE) |> 
  data.frame(mnp_QC_FrozenFFPEstatus_table = _) |> 
  dplyr::mutate(mnp_QC_FrozenFFPEstatus_table = paste0("data/G-SAM/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", mnp_QC_FrozenFFPEstatus_table)) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_FrozenFFPEstatus_table(mnp_QC_FrozenFFPEstatus_table, "mnp_QC_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% gsam.metadata.idats$sentrix_id) |> 
  dplyr::mutate(mnp_QC_FrozenFFPEstatus_table = NULL) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 75)
    return(.)
  })()


gsam.metadata.idats <- gsam.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(mnp_QC_predicted_array_type)) |> 
  assertr::verify(!is.na(mnp_QC_predicted_sample_type)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 75)
    return(.)
  })()

rm(tmp)





## Heidelberg 12.8 CNVP segment files ----

tmp <- query_mnp_12.8_CNVP_segment_csv(
  "data/G-SAM/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/",
  75, gsam.metadata.idats$sentrix_id)



gsam.metadata.idats <- gsam.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(heidelberg_cnvp_segments)) |> 
  assertr::verify(!is.na(heidelberg_cnvp_version)) 
rm(tmp)



## QC PCA outlier ----


tmp <- readRDS("cache/unsupervised_qc_outliers_all.Rds")

gsam.metadata.idats <- gsam.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(qc.pca.comp1)) |> 
  assertr::verify(!is.na(qc.pca.detP.outlier))

rm(tmp)


## ++ below: re-build because mvalue normalisation ++ ----


## Median methylation levels ----


tmp <- readRDS("cache/analysis_median_methylation.Rds") |> 
  assertr::verify(!is.na(median.overall.methylation)) |> 
  assertr::verify(!is.na(median.glass_nl_supervised.methylation)) |> 
  dplyr::filter(sentrix_id %in% gsam.metadata.idats$sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 73) # only HQ samples
    return(.)
  })()



gsam.metadata.idats <- gsam.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))


rm(tmp)



## A_IDH_HG__A_IDH_LG_lr__lasso_fit ----


tmp <- readRDS(file="cache/analysis_A_IDH_HG__A_IDH_LG_lr__lasso_fit.Rds") |> 
  dplyr::filter(!is.na(A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  dplyr::filter(sentrix_id %in% gsam.metadata.idats$sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 73)
    return(.)
  })()


gsam.metadata.idats <- gsam.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))






