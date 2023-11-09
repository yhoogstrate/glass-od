#!/usr/bin/env R

# G-SAM ----


source('scripts/load_functions.R')


## get idats ----


gsam.metadata.array_samples <-  list.files(path = "data/G-SAM/DNA Methylation - EPIC arrays/", pattern = "_(Grn|Red).idat$", recursive = T) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/G-SAM/DNA Methylation - EPIC arrays/", filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 190)# just freeze the number to terminate on unexpected behavior
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
    assertthat::assert_that(nrow(.) == 79)
    return(.)
  })()  |> 
  dplyr::rename_with( ~ paste0("array_", .x)) 




## link sample names ----


tmp <- rbind(
  read.csv("data/G-SAM/DNA Methylation - EPIC arrays/EPIC2023-387-plate1/STS/EPIC2023-387-plate1.csv", skip=8) |> 
    dplyr::filter(!is.na(Sentrix_ID)) |> 
    dplyr::filter(grepl("^GSAM", Sample_Name)) |> 
    dplyr::mutate(array_sentrix_id = paste0(Sentrix_ID, "_", Sentrix_Position )) |> 
    dplyr::mutate(resection = paste0("R",gsub("GSAM_...(.)_.+$","\\1",Sample_Name))) |> 
    dplyr::mutate(patient_id = gsub("GSAM_(...)._.+$","\\1", Sample_Name)) |> 
    dplyr::select(array_sentrix_id, resection, patient_id) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == 4)
      return(.)
    })()
,
  read.csv("data/G-SAM/DNA Methylation - EPIC arrays/MET2022-350-014/MET2022-350-014_IdH.csv", skip=8) |> 
    dplyr::filter(!is.na(Sentrix_ID)) |> 
    dplyr::filter(grepl("^GSAM", Sample_Name)) |> 
    dplyr::mutate(array_sentrix_id = paste0(Sentrix_ID, "_", Sentrix_Position )) |> 
    dplyr::mutate(resection = paste0("R",gsub("GSAM_...(.)_.+$","\\1",Sample_Name))) |> 
    dplyr::mutate(patient_id = gsub("GSAM_(...)._.+$","\\1", Sample_Name)) |> 
    dplyr::select(array_sentrix_id, resection, patient_id) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == 75) # == 75
      return(.)
    })()
  ) |> 
  dplyr::mutate(IDH = patient_id %in% c("BAW","CAV","CBG","CDF","DAB","EAF","EBD","ECB","FAD","FAL","JAB","JAD","JAF","KAC")) |> # IDH mut according to "data/gsam/output/tables/dna/idh_mutations.txt" - EAF also by MNP brain Classifier
  dplyr::mutate(resection_id = paste0(patient_id, gsub("^R","",resection))) |> 
  assertr::verify(!duplicated(resection_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == nrow(gsam.metadata.array_samples)) # == 75
    return(.)
  })()



gsam.metadata.array_samples <- gsam.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(resection_id))
rm(tmp)





## Add survival data ----


tmp <- read.csv('data/G-SAM/administratie/GSAM_combined_clinical_molecular.csv',stringsAsFactors=F) |> 
  dplyr::rename(patient_id = studyID) |> 
  dplyr::arrange(patient_id) |> 
  dplyr::mutate(initialMGMT = NULL) |> 
  dplyr::mutate(gender = ifelse(patient_id %in% c('AAT', 'AAM', 'AZH', 'HAI', 'FAG'),"Male",gender)) |>  # there's a number of samples of which the gender does not fit with the omics data - omics data determined genders are the corrected ones
  dplyr::mutate(gender = as.factor(gender)) |> 
  dplyr::mutate(survival.events = dplyr::case_when(
    status == "Deceased" ~ 1,
    status == "Censored" ~ 0,
    T ~ as.numeric(NA))) |> 
  dplyr::rename(os.event = survival.events) |> 
  dplyr::mutate(survival.months = survivalDays / 365.0 * 12.0) |> 
  dplyr::select(patient_id, gender, 
                status, os.event,
                survivalDays, 
                progressionFreeDays,
                survivalFromSecondSurgeryDays
                ) |> 
  dplyr::filter(patient_id %in% gsam.metadata.array_samples$patient_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == length(unique(gsam.metadata.array_samples$patient_id)))
    return(.)
  })() |> 
  dplyr::rename_with( ~ paste0("patient_", .x), .cols=!matches("^patient_id$",perl = T))


gsam.metadata.array_samples <- gsam.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('patient_id'='patient_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(patient_gender))
rm(tmp)





## Percentage detP probes ----
#' from: scripts/analysis_percentage_detP_probes.R


tmp <- read.table("output/tables/percentage_detP_probes.txt") |> 
  assertr::verify(!is.na(array_percentage.detP.signi) & is.numeric(array_percentage.detP.signi)) |> 
  assertr::verify(gsam.metadata.array_samples$array_sentrix_id %in% array_sentrix_id)

gsam.metadata.array_samples <- gsam.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_percentage.detP.signi)) 

rm(tmp)



## Unsupervised PCA QC all samples ----


tmp <- readRDS("cache/unsupervised_qc_outliers_all.Rds") |> 
  assertr::verify(!is.na(array_qc.pca.comp1)) |> 
  assertr::verify(is.numeric(array_qc.pca.comp1)) |> 
  assertr::verify(!is.na(array_qc.pca.detP.outlier)) |> 
  assertr::verify(is.logical(array_qc.pca.detP.outlier)) |> 
  assertr::verify(gsam.metadata.array_samples$array_sentrix_id %in% array_sentrix_id) |> 
  dplyr::filter(array_sentrix_id %in% gsam.metadata.array_samples$array_sentrix_id)


gsam.metadata.array_samples <- gsam.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_qc.pca.comp1)) |> 
  assertr::verify(!is.na(array_qc.pca.detP.outlier))


rm(tmp)



## Heidelberg 12.8 reportBrain files ----


tmp <- c(
  list.files(path = "data/G-SAM/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/",
             pattern = "*_scores_cal.csv", recursive = TRUE)) |>
  data.frame(tmp_filename = _) |>
  dplyr::mutate(array_mnp_predictBrain_filename = paste0("data/G-SAM/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", tmp_filename)) |>
  dplyr::mutate(array_mnp_predictBrain_version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", tmp_filename)) |>
  dplyr::mutate(array_sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_filename)) |>
  dplyr::mutate(tmp_filename = NULL) |> 
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |>  # only one version per sentrix_id desired
  assertr::verify(array_sentrix_id %in% gsam.metadata.array_samples$array_sentrix_id) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_reportBrain_csv(array_mnp_predictBrain_filename, paste0("array_mnp_predictBrain_", array_mnp_predictBrain_version, "_"))) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 79)
    return(.)
  })() |> 
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_class)) |>  # version is hardcoded here
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_A_IDH_LG)) |>
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_A_IDH_HG)) |>
  
  dplyr::mutate(array_A_IDH_HG__A_IDH_LG_lr_v12.8 = log(array_mnp_predictBrain_v12.8_cal_A_IDH_HG / array_mnp_predictBrain_v12.8_cal_A_IDH_LG)) |>
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_LG_lr_v12.8))


gsam.metadata.array_samples <- gsam.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_class)) |> 
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_LG_lr_v12.8)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 79)
    return(.)
  })() 

rm(tmp)



## Heidelberg 12.8 Frozen ~ FFPE status ----


tmp <- list.files(path = "data/G-SAM/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/",  pattern = "_ffpe_frozen.txt", recursive = TRUE) |> 
  data.frame(array_mnp_QC_v12.8_FrozenFFPEstatus_table = _) |> 
  dplyr::mutate(array_mnp_QC_v12.8_FrozenFFPEstatus_table = paste0("data/G-SAM/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", array_mnp_QC_v12.8_FrozenFFPEstatus_table)) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_FrozenFFPEstatus_table(array_mnp_QC_v12.8_FrozenFFPEstatus_table, "array_mnp_QC_v12.8_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(array_mnp_QC_v12.8_FrozenFFPEstatus_table = NULL) |> 
  
  assertr::verify(!is.na(array_sentrix_id)) |> 
  assertr::verify(!duplicated(array_sentrix_id)) |>
  
  assertr::verify(array_sentrix_id %in% gsam.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 79)
    return(.)
  })()


gsam.metadata.array_samples <- gsam.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_QC_v12.8_predicted_array_type)) |> 
  assertr::verify(!is.na(array_mnp_QC_v12.8_predicted_sample_type)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 79)
    return(.)
  })()

rm(tmp)




## Heidelberg 12.8 CNVP segment files ----


tmp <- query_mnp_12.8_CNVP_segment_csv(
  "data/G-SAM/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/",
  79,
  gsam.metadata.array_samples$array_sentrix_id,
  "array_mnp_CNVP_v12.8_v5.2_"
  )



gsam.metadata.array_samples <- gsam.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_CNVP_v12.8_v5.2_CNVP_segments))  |> 
  assertr::verify(!is.na(array_mnp_CNVP_v12.8_v5.2_CNVP_version)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 79)
    return(.)
  })()
rm(tmp)




## ++ below: re-build because mvalue normalisation ++ ----


## Median methylation levels ----


tmp <- readRDS("cache/analysis_median_methylation.Rds") |> 
  dplyr::rename_with( ~ paste0("array_", .x)) |> 
  assertr::verify(!is.na(array_median.overall.methylation)) |> 
  assertr::verify(!is.na(array_median.glass_nl_supervised.methylation)) |> 
  dplyr::filter(array_sentrix_id %in% gsam.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 77) # only HQ samples
    return(.)
  })()


gsam.metadata.array_samples <- gsam.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))


rm(tmp)


## A_IDH_HG__A_IDH_LG_lr__lasso_fit ----


tmp <- readRDS(file="cache/analysis_A_IDH_HG__A_IDH_LG_lr__lasso_fit.Rds") |>
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  dplyr::filter(array_sentrix_id %in% gsam.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 77)
    return(.)
  })()


gsam.metadata.array_samples <- gsam.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(sum(is.na(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) == 2) # 2 of too poor quality, never used because of odd probes etc.



## epiTOC2 ----


tmp <- readRDS("cache/analysis_EPITOC2.Rds") |> 
  dplyr::filter(array_sentrix_id %in% gsam.metadata.array_samples$array_sentrix_id) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 77)
    return(.)
  })()


gsam.metadata.array_samples <- gsam.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))
rm(tmp)






