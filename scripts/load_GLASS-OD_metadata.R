#!/usr/bin/env R

# load libs, config & db ----


source('scripts/load_functions.R')


metadata.db.con <- DBI::dbConnect(RSQLite::SQLite(), "../glass-od-clinical-database/glass-od-clinical-database.db")


# 1. patient level ----
# pat-47 - obviously a codel

# pat-27 - could be with additional 19q deletion, could be not - 206119350033_R04C01
# pat-56 - 1P & 19Q stable
# pat-76 - 1P stable, 19Q partial deletion like in P-98.
# pat-93 - 1P & 19Q stable
# pat-96 - partial 1Q gain, no 1P and no 19Q deletion visible
# pat-97 - 1P & 19Q stable - no CODEL
# pat-98 - (only one resection available) - 1P and 19Q arms seem partially deleted, but not at centromeres so no classical centromere fusion. This one may be codel-definition dependent, FISH could reveal co-localization of 1Q-19P.


## temporary file - incomplete as of yet

glass_od.metadata.patients <- DBI::dbReadTable(metadata.db.con, 'view_patients') |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 102)
    return(.)
  })()

patients_without_array_samples <- DBI::dbReadTable(metadata.db.con, 'view_check_patients_without_array_samples')
stopifnot(sort(patients_without_array_samples$patient_id) == c(1)) # x-checked, patients currently missing samples

glass_od.metadata.patients <- glass_od.metadata.patients |> 
  dplyr::filter(patient_id %in% patients_without_array_samples$patient_id == F) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 101)
    return(.)
  })() |> 
  dplyr::filter(is.na(reason_excluded)) |> # 7 non(-canonical) codels
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 91)
    return(.)
  })()



rm(patients_without_array_samples)



# 2. resection level ----


glass_od.metadata.resections <- DBI::dbReadTable(metadata.db.con, 'view_resections') |> 
  dplyr::filter(patient_id %in% glass_od.metadata.patients$patient_id) |> 
  dplyr::filter(is.na(reason_excluded_resection)) |> 
  assertr::verify(patient_id %in% c('26','63','27', '56', '76', '93', '96', '97', '98', '85') == F) |> # hard coded non-codels
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 200)
    return(.)
  })() |> 
  assertr::verify(is.numeric(resection_number)) |> 
  assertr::verify(is.na(resection_tumor_grade) | resection_tumor_grade %in% c(2,3))



# 3. idat level ----
## a. load all idat files ----


glass_od.metadata.idats <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays/", filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (414 + 28 + 2))
    return(.)
  })() |> # equals in-pipe stopifnot(nrow(.) == 404)
  assertr::verify(file.exists(filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", filename)) |>
  dplyr::mutate(channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", filename)) |>
  tidyr::pivot_wider(id_cols = sentrix_id, names_from = channel, values_from = c(filename)) |>
  dplyr::rename(channel_green = Grn) |>
  dplyr::rename(channel_red = Red) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == ((414 + 28 + 2) / 2))
    return(.)
  })() |>
  assertr::verify(!is.na(channel_green)) |>
  assertr::verify(!is.na(channel_red)) |>
  assertr::verify(sentrix_id != "206119350032_R01C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(sentrix_id != "206119350032_R02C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(sentrix_id != "206119350032_R03C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(sentrix_id != "206119350032_R04C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(sentrix_id != "206119350032_R05C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(sentrix_id != "204808700073_R07C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(sentrix_id != "204808700074_R06C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(sentrix_id != "204808700074_R07C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(sentrix_id != "204808700074_R08C01") # sample not part in GLASS-OD - removal confirmed by Iris



## b. load all idat metadata from access ----


glass_od.metadata.array_samples <- DBI::dbReadTable(metadata.db.con, 'view_array_samples')


# setdiff(glass_od.metadata.array_samples$sentrix_id, glass_od.metadata.idats$sentrix_id)
# setdiff(glass_od.metadata.idats$sentrix_id, glass_od.metadata.array_samples$sentrix_id)

stopifnot(length(setdiff(glass_od.metadata.array_samples$sentrix_id, glass_od.metadata.idats$sentrix_id)) == 0)
stopifnot(length(setdiff(glass_od.metadata.idats$sentrix_id, glass_od.metadata.array_samples$sentrix_id)) == 0)



## link patient identifier ----


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })() |> 
  dplyr::left_join(glass_od.metadata.array_samples, by=c('sentrix_id' = 'sentrix_id'), suffix=c('','')) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })() |> 
  assertr::verify(!is.na(resection_id))



rm(glass_od.metadata.array_samples)



## Percentage detP probes ----

# from: scripts/analysis_percentage_detP_probes.R


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(read.table("output/tables/percentage_detP_probes.txt"), by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(percentage.detP.signi))


## Heidelberg 11b4[+12.5] QC full ----

# quite some of these files contain odd N/A's

tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", pattern = "_qc_full.txt", recursive = TRUE) |> 
  data.frame(heidelberg_qc_report_full = _) |> 
  dplyr::mutate(heidelberg_qc_report_full = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", heidelberg_qc_report_full)) |>
  dplyr::mutate(heidelberg_qc_report_version = gsub("^.+qc_(v[^_\\/]+)[_/].+$","\\1", heidelberg_qc_report_full)) |>
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", heidelberg_qc_report_full)) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_QCReport_csv(heidelberg_qc_report_full, "qc_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id)


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(heidelberg_qc_report_version)) |> 
  assertr::verify(!is.na(heidelberg_qc_report_full)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })() |> 
  dplyr::mutate(heidelberg_qc_report_full = NULL)

rm(tmp)





## Heidelberg 11b4[+12.5] reportBrain files ----
#' needed to correlate LGC


tmp <- c(list.files(
  path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", pattern = "*_scores_cal.csv", recursive = TRUE)
) |>
  data.frame(filename = _) |>
  dplyr::mutate(mnp_predictBrain_filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", filename)) |>
  dplyr::mutate(mnp_predictBrain_version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", filename)) |>
  dplyr::mutate(filename = NULL) |> 
  tidyr::pivot_wider(id_cols = sentrix_id,
                     names_from = mnp_predictBrain_version, 
                     values_from = c(mnp_predictBrain_filename), 
                     names_prefix = "mnp_predictBrain_") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })() |> 
  dplyr::rename(`mnp_predictBrain_v12.5_filename` = `mnp_predictBrain_v12.5`) |> 
  dplyr::rename(`mnp_predictBrain_v2.0.1_filename` = `mnp_predictBrain_v2.0.1`) |> 
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |>
  assertr::verify(file.exists(`mnp_predictBrain_v12.5_filename`)) |>
  assertr::verify(file.exists(`mnp_predictBrain_v2.0.1_filename`)) |>
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(`tmp_v2.0.1` = parse_mnp_reportBrain_csv(`mnp_predictBrain_v2.0.1_filename`, paste0("mnp_predictBrain_v2.0.1_"))) |>
  dplyr::mutate(tmp_v12.5 = parse_mnp_reportBrain_csv(`mnp_predictBrain_v12.5_filename`, paste0("mnp_predictBrain_v12.5_"))) |>
  dplyr::ungroup() |> 
  tidyr::unnest(`tmp_v2.0.1`) |> 
  tidyr::unnest(tmp_v12.5) |> 
  assertr::verify(!is.na(mnp_predictBrain_v2.0.1_cal_O_IDH)) |> 
  assertr::verify(!is.na(mnp_predictBrain_v2.0.1_cal_A_IDH)) |> 
  assertr::verify(!is.na(mnp_predictBrain_v2.0.1_cal_A_IDH_HG)) |> 
  assertr::verify(!is.na(mnp_predictBrain_v12.5_cal_O_IDH))  |> 
  assertr::verify(!is.na(mnp_predictBrain_v12.5_cal_A_IDH_LG))  |> 
  assertr::verify(!is.na(mnp_predictBrain_v12.5_cal_A_IDH_HG)) |>  
  assertr::verify(!is.na(mnp_predictBrain_v12.5_cal_OLIGOSARC_IDH)) |> 
  
  dplyr::mutate(A_IDH_HG__A_IDH_lr = log(mnp_predictBrain_v2.0.1_cal_A_IDH_HG / mnp_predictBrain_v2.0.1_cal_A_IDH)) |> 
  assertr::verify(!is.na(A_IDH_HG__A_IDH_lr)) |> 
  
  dplyr::mutate(A_IDH_HG__A_IDH_lr_neat = log(
    
    (mnp_predictBrain_v2.0.1_cal_A_IDH_HG / (1-mnp_predictBrain_v2.0.1_cal_A_IDH_HG))
    / 
      (mnp_predictBrain_v2.0.1_cal_A_IDH / (1-mnp_predictBrain_v2.0.1_cal_A_IDH))
    
  )) 



glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(A_IDH_HG__A_IDH_lr)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })()

rm(tmp)


## Heidelberg 11b4[+12.5] rs_gender ----

tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", pattern = "*.mix_gender.csv", recursive = TRUE) |> 
  data.frame(heidelberg_rs_gender_report = _) |> 
  dplyr::mutate(heidelberg_rs_gender_report = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", heidelberg_rs_gender_report)) |>
  assertr::verify(file.exists(heidelberg_rs_gender_report)) |> 
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", heidelberg_rs_gender_report)) |> 
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_RsGender_csv(heidelberg_rs_gender_report, "mnp_rsGender_11b4_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(heidelberg_rs_gender_report = NULL) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })()



glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(mnp_rsGender_11b4_predicted))
rm(tmp)



## Heidelberg 12.8 reportBrain files ----

tmp <- c(list.files(
  path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "*_scores_cal.csv", recursive = TRUE)
) |>
  data.frame(filename = _) |>
  dplyr::mutate(mnp_predictBrain_filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", filename)) |>
  dplyr::mutate(mnp_predictBrain_version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", filename)) |>
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |>  # only one version per sentrix_id 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_reportBrain_csv(mnp_predictBrain_filename, paste0("mnp_predictBrain_", mnp_predictBrain_version, "_"))) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })() |> 
  assertr::verify(!is.na(mnp_predictBrain_v12.8_cal_class)) |>  # version is hardcoded here
  assertr::verify(!is.na(mnp_predictBrain_v12.8_cal_A_IDH_LG)) |>
  assertr::verify(!is.na(mnp_predictBrain_v12.8_cal_A_IDH_HG)) |>
  dplyr::mutate(A_IDH_HG__A_IDH_LG_lr = log(mnp_predictBrain_v12.8_cal_A_IDH_HG / mnp_predictBrain_v12.8_cal_A_IDH_LG))  |>
  dplyr::mutate(A_IDH_HG__A_IDH_LG_lr_neat = log(
    (mnp_predictBrain_v12.8_cal_A_IDH_HG / (1-mnp_predictBrain_v12.8_cal_A_IDH_HG))
    / 
    (mnp_predictBrain_v12.8_cal_A_IDH_LG / (1-mnp_predictBrain_v12.8_cal_A_IDH_LG))
    )) 



glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(mnp_predictBrain_v12.8_cal_class)) |> 
  assertr::verify(!is.na(A_IDH_HG__A_IDH_LG_lr)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })() 

rm(tmp)



## Heidelberg 12.8 Frozen ~ FFPE status ----

tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/",  pattern = "_ffpe_frozen.txt", recursive = TRUE) |> 
  data.frame(mnp_QC_FrozenFFPEstatus_table = _) |> 
  dplyr::mutate(mnp_QC_FrozenFFPEstatus_table = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", mnp_QC_FrozenFFPEstatus_table)) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_FrozenFFPEstatus_table(mnp_QC_FrozenFFPEstatus_table, "mnp_QC_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id) |> 
  dplyr::mutate(mnp_QC_FrozenFFPEstatus_table = NULL) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })()


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(mnp_QC_predicted_array_type)) |> 
  assertr::verify(!is.na(mnp_QC_predicted_sample_type)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })()

rm(tmp)



## Heidelberg 12.8 CNVP segment files ----

tmp <- query_mnp_12.8_CNVP_segment_csv("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", 222, glass_od.metadata.idats$sentrix_id)

glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(heidelberg_cnvp_segments)) |> 
  assertr::verify(!is.na(heidelberg_cnvp_version)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })()
rm(tmp)




## Heidelberg 12.8 CNVP bins files ----


tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "*.bins.igv", recursive = TRUE) |> 
  data.frame(heidelberg_cnvp_bins = _) |> 
  dplyr::mutate(heidelberg_cnvp_bins = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", heidelberg_cnvp_bins)) |> 
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", heidelberg_cnvp_bins)) |> 
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id)


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(heidelberg_cnvp_bins)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })()
rm(tmp)


## Heidelberg 12.8 CNVP ongene scores ----


tmp <- c(
  list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "*.detail.txt", recursive = TRUE)
) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", filename)) |>
  assertr::verify(file.exists(filename)) |> 
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", filename)) |>
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_CNVPoncogeneScores_csv (filename, "mnp_CNVP_12.8_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(filename = NULL)


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(`mnp_CNVP_12.8_CDKN2A/B`)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })()

rm(tmp)






## Heidelberg 12.8 predictMGMT ----

tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "_mgmt.csv", recursive = TRUE) |> 
  data.frame(heidelberg_mgmt_report = _) |> 
  dplyr::mutate(heidelberg_mgmt_report = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", heidelberg_mgmt_report)) |>
  assertr::verify(file.exists(heidelberg_mgmt_report)) |> 
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", heidelberg_mgmt_report)) |>
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_predictMGMT_csv(heidelberg_mgmt_report, "mnp_MGMT_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(heidelberg_mgmt_report = NULL) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })()


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(mnp_MGMT_Estimated))
rm(tmp)



## Heidelberg 12.8 rs_gender ----

tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "*.mix_gender.csv", recursive = TRUE) |> 
  data.frame(heidelberg_rs_gender_report = _) |> 
  dplyr::mutate(heidelberg_rs_gender_report = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", heidelberg_rs_gender_report)) |>
  assertr::verify(file.exists(heidelberg_rs_gender_report)) |> 
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", heidelberg_rs_gender_report)) |> 
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_RsGender_csv(heidelberg_rs_gender_report, "mnp_rsGender_12.8_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(heidelberg_rs_gender_report = NULL) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })()
  


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(`mnp_rsGender_12.8_chrYintensity`))
rm(tmp)



## QC PCA outlier ----


tmp <- readRDS("cache/unsupervised_qc_outliers_all.Rds")

glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(qc.pca.comp1)) |> 
  assertr::verify(!is.na(qc.pca.detP.outlier))

rm(tmp)


## Median methylation levels ----


tmp <- readRDS("cache/analysis_median_methylation.Rds")

glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))

glass_od.metadata.idats |> 
  filter_GLASS_OD_idats(163) |> 
  assertr::verify(!is.na(median.overall.methylation)) |> 
  assertr::verify(!is.na(median.glass_nl_supervised.methylation))


rm(tmp)



#plot(glass_od.metadata.idats$A_IDH_HG__A_IDH_LG_lr, glass_od.metadata.idats$median.glass_nl_supervised.methylation, xlim=c(-8,16)) 
#plot(glass_od.metadata.idats$A_IDH_HG__A_IDH_LG_lr_neat, glass_od.metadata.idats$median.glass_nl_supervised.methylation, xlim=c(-8,16))

#plot(glass_od.metadata.idats$A_IDH_HG__A_IDH_LG_lr, glass_od.metadata.idats$median.overall.methylation, xlim=c(-8,16))
#plot(glass_od.metadata.idats$A_IDH_HG__A_IDH_LG_lr_neat, glass_od.metadata.idats$median.overall.methylation, xlim=c(-8,16))

# plot(glass_od.metadata.idats$A_IDH_HG__A_IDH_LG_lr, glass_od.metadata.idats$median.glass_nl_supervised.methylation)
# plot(glass_od.metadata.idats$A_IDH_HG__A_IDH_LG_lr, glass_od.metadata.idats$median.overall.methylation)

# plot(glass_od.metadata.idats$A_IDH_HG__A_IDH_lr, glass_od.metadata.idats$median.glass_nl_supervised.methylation)
# plot(glass_od.metadata.idats$A_IDH_HG__A_IDH_lr, glass_od.metadata.idats$median.overall.methylation)



# plot(glass_od.metadata.idats$A_IDH_HG__A_IDH_LG_lr, glass_od.metadata.idats$median.glass_nl_supervised.methylation)
# plot(glass_od.metadata.idats$A_IDH_HG__A_IDH_LG_lr, glass_od.metadata.idats$median.overall.methylation)

# plot(glass_od.metadata.idats$A_IDH_HG__A_IDH_lr, glass_od.metadata.idats$median.glass_nl_supervised.methylation)
# plot(glass_od.metadata.idats$A_IDH_HG__A_IDH_lr, glass_od.metadata.idats$median.overall.methylation)



## bin-based tumor purity calls ----


tmp <- readRDS("cache/analysis_tumor_purity_EPIC_bin-based.Rds")
glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))

rm(tmp)


## A_IDH_HG__A_IDH_LG_lr__lasso_fit ----


tmp <- readRDS(file="cache/GLASS-OD__A_IDH_HG__A_IDH_LG_lr__lasso_fit.Rds")


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'),suffix=c('',''))


stopifnot(glass_od.metadata.idats |>
  dplyr::filter(!is.na(A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  nrow() == 163)


# cleanup db connection ----

DBI::dbDisconnect(metadata.db.con)
rm(metadata.db.con)



