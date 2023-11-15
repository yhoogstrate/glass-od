#!/usr/bin/env R

# load libs, config & db ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')
#source('scripts/load_palette.R')
#source('scripts/load_themes.R')


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
    assertthat::assert_that(nrow(.) == 127)
    return(.)
  })() |> 
  assertr::verify(is.na(patient_diagnosis_date) | grepl("^[0-9]{1,2} (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) [0-9]{4}$", patient_diagnosis_date)) |> 
  assertr::verify(is.na(patient_IDH_mutation) | (patient_IDH_mutation %in% c(
    "IDH1 R132H",
    "IDH1 R132S",
    
    "IDH2 R172K",
    "IDH2 R172M"
  ))) |> 
  dplyr::mutate(patient_suspected_noncodel = dplyr::case_when(
    patient_suspected_noncodel == "true" ~ TRUE,
    patient_suspected_noncodel == "false" ~ FALSE,
    TRUE ~ as.logical(NA)
  ))



patients_without_array_samples <- DBI::dbReadTable(metadata.db.con, 'view_check_patients_without_array_samples')
stopifnot(sort(patients_without_array_samples$patient_id) == c("0001")) # x-checked, patients currently missing samples



glass_od.metadata.patients <- glass_od.metadata.patients |> 
  dplyr::filter(patient_id %in% patients_without_array_samples$patient_id == F) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) ==  126) # + 10x astro
    return(.)
  })() |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> # 7 non(-canonical) codels
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 126) # + 10x astro
    return(.)
  })()



rm(patients_without_array_samples)



# 2. resection level ----


glass_od.metadata.resections <- DBI::dbReadTable(metadata.db.con, 'view_resections') |> 
  dplyr::filter(patient_id %in% glass_od.metadata.patients$patient_id) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::mutate(patient_suspected_noncodel = dplyr::case_when(
    patient_suspected_noncodel == "true" ~ TRUE,
    patient_suspected_noncodel == "false" ~ FALSE,
    TRUE ~ as.logical(NA)
  )) |> 
  assertr::verify((patient_id %in% c('0026','0063','0027', '0056', '0076', '0093', '0096', '0097', '0098', '0085') & patient_suspected_noncodel == T) | (patient_suspected_noncodel == F)) |> 
  dplyr::mutate(patient_id = as.factor(patient_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 277) # + 22 x astro
    return(.)
  })() |> 
  assertr::verify(is.numeric(resection_number)) |> 
  assertr::verify(is.na(resection_tumor_grade) | resection_tumor_grade %in% c(2, 3)) |> 
  assertr::verify(is.na(resection_date) | grepl("^[0-9]{1,2} (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) [0-9]{4}$", resection_date))




# 3. idat level ----
## a. load all idat files ----


glass_od.metadata.array_samples <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(array_filename = _) |>
  dplyr::mutate(array_filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays/", array_filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (550))
    return(.)
  })() |> # equals in-pipe stopifnot(nrow(.) == 404)
  assertr::verify(file.exists(array_filename)) |>
  dplyr::mutate(array_sentrix_id = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", array_filename)) |>
  dplyr::mutate(array_channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", array_filename)) |>
  tidyr::pivot_wider(id_cols = array_sentrix_id, names_from = array_channel, values_from = c(array_filename)) |>
  dplyr::rename(array_channel_green = Grn) |>
  dplyr::rename(array_channel_red = Red) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == ((550) / 2))
    return(.)
  })() |>
  assertr::verify(!is.na(array_channel_green)) |>
  assertr::verify(!is.na(array_channel_red)) |>
  assertr::verify(array_sentrix_id != "206119350032_R01C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(array_sentrix_id != "206119350032_R02C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(array_sentrix_id != "206119350032_R03C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(array_sentrix_id != "206119350032_R04C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(array_sentrix_id != "206119350032_R05C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(array_sentrix_id != "204808700073_R07C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(array_sentrix_id != "204808700074_R06C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(array_sentrix_id != "204808700074_R07C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(array_sentrix_id != "204808700074_R08C01") |> # sample not part in GLASS-OD - removal confirmed by Iris

  assertr::verify(array_sentrix_id != "207331540058_R01C01") |> # non-GLASS_OD files included in a run shared with GLASS_OD samples (MINT, G-SAM, CATNON/HOX & IDH-inhibitor)
  assertr::verify(array_sentrix_id != "207331540058_R02C01") |> # non-GLASS_OD files included in a run shared with GLASS_OD samples (MINT, G-SAM, CATNON/HOX & IDH-inhibitor)
  assertr::verify(array_sentrix_id != "207331540058_R03C01") |> # non-GLASS_OD files included in a run shared with GLASS_OD samples (MINT, G-SAM, CATNON/HOX & IDH-inhibitor)
  assertr::verify(array_sentrix_id != "207331540058_R04C01") |> # non-GLASS_OD files included in a run shared with GLASS_OD samples (MINT, G-SAM, CATNON/HOX & IDH-inhibitor)
  assertr::verify(array_sentrix_id != "207513900021_R04C01") |> # non-GLASS_OD files included in a run shared with GLASS_OD samples (MINT, G-SAM, CATNON/HOX & IDH-inhibitor)
  assertr::verify(array_sentrix_id != "207513900021_R02C01") |> # non-GLASS_OD files included in a run shared with GLASS_OD samples (MINT, G-SAM, CATNON/HOX & IDH-inhibitor)
  assertr::verify(array_sentrix_id != "207513900021_R03C01") |> # non-GLASS_OD files included in a run shared with GLASS_OD samples (MINT, G-SAM, CATNON/HOX & IDH-inhibitor)
  assertr::verify(array_sentrix_id != "207331540058_R05C01") |> # non-GLASS_OD files included in a run shared with GLASS_OD samples (MINT, G-SAM, CATNON/HOX & IDH-inhibitor)
  assertr::verify(array_sentrix_id != "207331540058_R06C01") |> # non-GLASS_OD files included in a run shared with GLASS_OD samples (MINT, G-SAM, CATNON/HOX & IDH-inhibitor)
  assertr::verify(array_sentrix_id != "207331540058_R07C01") |> # non-GLASS_OD files included in a run shared with GLASS_OD samples (MINT, G-SAM, CATNON/HOX & IDH-inhibitor)
  assertr::verify(array_sentrix_id != "207331540058_R08C01") # non-GLASS_OD files included in a run shared with GLASS_OD samples (MINT, G-SAM, CATNON/HOX & IDH-inhibitor)




## b. load all idat metadata from access ----


tmp <- DBI::dbReadTable(metadata.db.con, 'view_array_samples') |> 
  dplyr::mutate(patient_suspected_noncodel = dplyr::case_when(
    patient_suspected_noncodel == "true" ~ TRUE,
    patient_suspected_noncodel == "false" ~ FALSE,
    TRUE ~ as.logical(NA)
  )) |> 
  assertr::verify(arrayplate_id != "GLASS-NL") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(sum(is.na(.$arraychip_date)) <= 6)
    return(.)
  })() |> 
  
  dplyr::mutate(patient_diagnosis_date = as.Date(patient_diagnosis_date, format = "%d %b %Y")) |> 
  dplyr::mutate(patient_birth_date = as.Date(patient_birth_date, format = "%d %b %Y")) |> 
  dplyr::mutate(resection_date = as.Date(resection_date, format = "%d %b %Y")) |> 
  dplyr::mutate(arraychip_date = as.Date(arraychip_date, format = "%d %b %Y")) |> 
  
  assertr::verify(is.na(patient_diagnosis_date) | is.na(resection_date) | patient_diagnosis_date <= resection_date) |> 
  assertr::verify(is.na(arraychip_date) | is.na(resection_date) | resection_date < arraychip_date) |> 
  assertr::verify(is.na(patient_birth_date) | is.na(resection_date) | patient_birth_date < resection_date) |> 
  
  dplyr::mutate(time_between_resection_and_array = arraychip_date - resection_date) |> 
  dplyr::mutate(time_between_birth_and_resection = resection_date - patient_birth_date)



stopifnot(length(setdiff(glass_od.metadata.array_samples$array_sentrix_id, tmp$array_sentrix_id)) == 0)
stopifnot(length(setdiff(tmp$array_sentrix_id, glass_od.metadata.array_samples$array_sentrix_id)) == 0)


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (275))
    return(.)
  })() |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id' = 'array_sentrix_id'), suffix=c('','')) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })() |> 
  assertr::verify(!is.na(resection_id))
  



rm(tmp)



## HE coupes ----


tmp <- DBI::dbReadTable(metadata.db.con, 'HE_coupes') |> 
  dplyr::rename_with( ~ paste0("HE_coupe_", .x), .cols=!matches("^resection_id$",perl = T)) |> 
  dplyr::mutate(file_exists = file.exists(paste0("data/GLASS_OD/Stainings/H&E Slides/",HE_coupe_filename))) |> 
  assertr::verify(file_exists) |> 
  dplyr::mutate(file_exists = NULL) |> 
  dplyr::mutate(join_prefix = gsub("\\.ndpi$", "", HE_coupe_filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 223)
    
    print(length(unique(.$resection_id)))
    assertthat::assert_that(length(unique(.$resection_id)) == 184)
    return(.)
  })()


tmp.thumbnails_full <- list.files(path = "output/figures/HE_thumbnails/full/", pattern = "*.png$", recursive = TRUE) |>
  data.frame(HE_coupe_thumbnail_full = _) |>
  dplyr::mutate(HE_coupe_thumbnail_full = paste0("output/figures/HE_thumbnails/full/", HE_coupe_thumbnail_full)) |>
  dplyr::mutate(join_prefix = gsub("^.+/(.+)_[0-9]+x[0-9]+px\\.png$", "\\1", HE_coupe_thumbnail_full)) |> 
  assertr::verify(!duplicated(join_prefix)) |> 
  assertr::verify(join_prefix %in% tmp$join_prefix) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 223)
    return(.)
  })()


tmp.thumbnails_roi <- list.files(path = "output/figures/HE_thumbnails/roi/", pattern = "*.jpg$", recursive = TRUE) |>
  data.frame(HE_coupe_thumbnail_roi = _) |>
  dplyr::mutate(HE_coupe_thumbnail_roi = paste0("output/figures/HE_thumbnails/roi/", HE_coupe_thumbnail_roi)) |>
  dplyr::mutate(join_prefix = gsub("^.+/(.+)\\.ndpi_[0-9]+\\.jpg$", "\\1", HE_coupe_thumbnail_roi)) |> 
  assertr::verify(!duplicated(join_prefix)) |> 
  assertr::verify(join_prefix %in% tmp$join_prefix) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 200)
    return(.)
  })()



tmp <- tmp |> 
  dplyr::left_join(tmp.thumbnails_full, by=c('join_prefix'='join_prefix'), suffix=c('','')) |> 
  dplyr::left_join(tmp.thumbnails_roi, by=c('join_prefix'='join_prefix'), suffix=c('','')) |> 
  dplyr::mutate(join_prefix = NULL) |> 
  dplyr::filter(HE_coupe_primary_coupe == "yes") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 184)
    return(.)
  })() |> 
  assertr::verify(!is.na(HE_coupe_thumbnail_full)) |> 
  assertr::verify(!is.na(HE_coupe_thumbnail_roi))


rm(tmp.thumbnails_full)






glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('resection_id'='resection_id'), suffix=c('',''))

rm(tmp)





## Percentage detP probesC ----
#' from: scripts/analysis_percentage_detP_probes.R


tmp <- read.table("output/tables/percentage_detP_probes.txt") |> 
  assertr::verify(!is.na(array_percentage.detP.signi) & is.numeric(array_percentage.detP.signi)) |> 
  dplyr::rename(array_sentrix_id = array_sentrix_id) |> 
  assertr::verify(glass_od.metadata.array_samples$array_sentrix_id %in% array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (783))
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_percentage.detP.signi))


rm(tmp)



## QC PCA QC all samples outlier ----


tmp <- readRDS('cache/unsupervised_qc_outliers_all.Rds') |> 
  assertr::verify(!is.na(array_qc.pca.comp1)) |> 
  assertr::verify(is.numeric(array_qc.pca.comp1)) |> 
  assertr::verify(!is.na(array_qc.pca.detP.outlier)) |> 
  assertr::verify(is.logical(array_qc.pca.detP.outlier)) |> 
  assertr::verify(glass_od.metadata.array_samples$array_sentrix_id %in% array_sentrix_id) |> 
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (275))
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_qc.pca.comp1)) |> 
  assertr::verify(!is.na(array_qc.pca.detP.outlier))

rm(tmp)





## Heidelberg 11b4[+12.5] QC full ----

# quite some of these files contain odd N/A's

tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", pattern = "_qc_full.txt", recursive = TRUE) |> 
  data.frame(array_mnp_qc_report_full = _) |> 
  dplyr::mutate(array_mnp_qc_report_full = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", array_mnp_qc_report_full)) |>
  dplyr::mutate(array_mnp_qc_report_version = gsub("^.+qc_(v[^_\\/]+)[_/].+$","\\1", array_mnp_qc_report_full)) |>
  dplyr::mutate(array_sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", array_mnp_qc_report_full)) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_QCReport_csv(array_mnp_qc_report_full, "array_qc_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_qc_report_version)) |> 
  assertr::verify(!is.na(array_mnp_qc_report_full)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })() |> 
  dplyr::mutate(array_mnp_qc_report_full = NULL)

rm(tmp)





## Heidelberg 11b4[+12.5] reportBrain files ----
#' needed to correlate LGC


tmp <- c(list.files(
  path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", pattern = "*_scores_cal.csv", recursive = TRUE)
) |>
  data.frame(tmp_filename = _) |>
  dplyr::mutate(array_mnp_predictBrain_filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", tmp_filename)) |>
  dplyr::mutate(array_mnp_predictBrain_version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", tmp_filename)) |>
  dplyr::mutate(array_sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_filename)) |>
  dplyr::mutate(tmp_filename = NULL) |> 
  tidyr::pivot_wider(id_cols = array_sentrix_id,
                     names_from = array_mnp_predictBrain_version, 
                     values_from = c(array_mnp_predictBrain_filename), 
                     names_prefix = "array_mnp_predictBrain_") |> 
  dplyr::rename(`array_mnp_predictBrain_v12.5_filename` = `array_mnp_predictBrain_v12.5`) |> 
  dplyr::rename(`array_mnp_predictBrain_v2.0.1_filename` = `array_mnp_predictBrain_v2.0.1`) |> 
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |>
  assertr::verify(file.exists(`array_mnp_predictBrain_v12.5_filename`)) |>
  assertr::verify(file.exists(`array_mnp_predictBrain_v2.0.1_filename`)) |>
  assertr::verify(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(`tmp_v2.0.1` = parse_mnp_reportBrain_csv(`array_mnp_predictBrain_v2.0.1_filename`, paste0("array_mnp_predictBrain_v2.0.1_"))) |>
  dplyr::mutate(tmp_v12.5 = parse_mnp_reportBrain_csv(`array_mnp_predictBrain_v12.5_filename`, paste0("array_mnp_predictBrain_v12.5_"))) |>
  dplyr::ungroup() |> 
  tidyr::unnest(`tmp_v2.0.1`) |> 
  tidyr::unnest(tmp_v12.5) |> 
  assertr::verify(!is.na(array_mnp_predictBrain_v2.0.1_cal_O_IDH)) |> 
  assertr::verify(!is.na(array_mnp_predictBrain_v2.0.1_cal_A_IDH)) |> 
  assertr::verify(!is.na(array_mnp_predictBrain_v2.0.1_cal_A_IDH_HG)) |> 
  assertr::verify(!is.na(array_mnp_predictBrain_v12.5_cal_O_IDH))  |> 
  assertr::verify(!is.na(array_mnp_predictBrain_v12.5_cal_A_IDH_LG))  |> 
  assertr::verify(!is.na(array_mnp_predictBrain_v12.5_cal_A_IDH_HG)) |>  
  assertr::verify(!is.na(array_mnp_predictBrain_v12.5_cal_OLIGOSARC_IDH)) |> 
  dplyr::mutate(array_A_IDH_HG__A_IDH_lr_v2.0.1 = log(array_mnp_predictBrain_v2.0.1_cal_A_IDH_HG / array_mnp_predictBrain_v2.0.1_cal_A_IDH)) |> 
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_lr_v2.0.1)) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()
  

glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_lr_v2.0.1)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()


rm(tmp)



## Heidelberg 11b4[+12.5] rs_gender ----


tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", pattern = "*.mix_gender.csv", recursive = TRUE) |> 
  data.frame(tmp_heidelberg_rs_gender_report = _) |> 
  dplyr::mutate(tmp_heidelberg_rs_gender_report = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", tmp_heidelberg_rs_gender_report)) |>
  assertr::verify(file.exists(tmp_heidelberg_rs_gender_report)) |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_heidelberg_rs_gender_report)) |> 
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_RsGender_csv(tmp_heidelberg_rs_gender_report, "array_mnp_rsGender_11b4_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(tmp_heidelberg_rs_gender_report = NULL) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_rsGender_11b4_predicted))

rm(tmp)



## Heidelberg 12.8 reportBrain files ----


tmp <- c(list.files(
  path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "*_scores_cal.csv", recursive = TRUE)
) |>
  data.frame(tmp_filename = _) |>
  dplyr::mutate(array_mnp_predictBrain_filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", tmp_filename)) |>
  dplyr::mutate(array_mnp_predictBrain_version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", tmp_filename)) |>
  dplyr::mutate(array_sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_filename)) |>
  dplyr::mutate(tmp_filename = NULL) |> 
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |>  # only one version per sentrix_id 
  assertr::verify(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_reportBrain_csv(array_mnp_predictBrain_filename, paste0("array_mnp_predictBrain_", array_mnp_predictBrain_version, "_"))) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 

  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_class)) |>  # version is hardcoded here
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_A_IDH_LG)) |>
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_A_IDH_HG)) |>
  
  dplyr::mutate(array_A_IDH_HG__A_IDH_LG_lr_v12.8 = log(array_mnp_predictBrain_v12.8_cal_A_IDH_HG / array_mnp_predictBrain_v12.8_cal_A_IDH_LG))  |>
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_LG_lr_v12.8)) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()



  #dplyr::mutate(array_A_IDH_HG__O_IDH_lr = log(array_mnp_predictBrain_v12.8_cal_A_IDH_HG / array_mnp_predictBrain_v12.8_cal_O_IDH))  |>
  #assertr::verify(!is.na(array_A_IDH_HG__O_IDH_lr)) |> 
  #dplyr::mutate(array_A_IDH_HGoligsarc__O_IDH_lr = log((array_mnp_predictBrain_v12.8_cal_A_IDH_HG + array_mnp_predictBrain_v12.8_cal_OLIGOSARC_IDH) / array_mnp_predictBrain_v12.8_cal_O_IDH))  |>
  #assertr::verify(!is.na(array_A_IDH_HGoligsarc__O_IDH_lr))
  # dplyr::mutate(A_IDH_HG__A_IDH_LG_lr_neat = log((mnp_predictBrain_v12.8_cal_A_IDH_HG / (1-mnp_predictBrain_v12.8_cal_A_IDH_HG)) /  (mnp_predictBrain_v12.8_cal_A_IDH_LG / (1-mnp_predictBrain_v12.8_cal_A_IDH_LG)) ))  |> 
  # assertr::verify(!is.na(A_IDH_HG__A_IDH_LG_lr_neat))



glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_class)) |> 
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_LG_lr_v12.8)) |> 
  (function(.) {
    print(paste0(sum(!is.na(.$array_mnp_predictBrain_v12.8_cal_class)), "/", nrow(.)))
    assertthat::assert_that(sum(!is.na(.$array_mnp_predictBrain_v12.8_cal_class)) == 275) # 0118-R4 still failing
    return(.)
  })() 

rm(tmp)




## Heidelberg 12.8 Frozen ~ FFPE status ----


tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/",  pattern = "_ffpe_frozen.txt", recursive = TRUE) |> 
  data.frame(array_mnp_QC_v12.8_FrozenFFPEstatus_table = _) |> 
  dplyr::mutate(array_mnp_QC_v12.8_FrozenFFPEstatus_table = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", array_mnp_QC_v12.8_FrozenFFPEstatus_table)) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_FrozenFFPEstatus_table(array_mnp_QC_v12.8_FrozenFFPEstatus_table, "array_mnp_QC_v12.8_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(array_mnp_QC_v12.8_FrozenFFPEstatus_table = NULL) |> 
  
  assertr::verify(!is.na(array_sentrix_id)) |> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_QC_v12.8_predicted_array_type)) |> 
  assertr::verify(!is.na(array_mnp_QC_v12.8_predicted_sample_type) ) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()

rm(tmp)



## Heidelberg 12.8 CNVP segment files ----


tmp <- query_mnp_12.8_CNVP_segment_csv(
  "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", 
  275, 
  glass_od.metadata.array_samples$array_sentrix_id,
  "array_mnp_CNVP_v12.8_v5.2_"
  ) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()



glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_CNVP_v12.8_v5.2_CNVP_segments) )  |> 
  assertr::verify(!is.na(array_mnp_CNVP_v12.8_v5.2_CNVP_version) ) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()
rm(tmp)




## Heidelberg 12.8 CNVP bins files ----


tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "*.bins.igv", recursive = TRUE) |> 
  data.frame(array_heidelberg_cnvp_bins = _) |> 
  dplyr::mutate(array_mnp_CNVP_v12.8_v5.2_CNVP_bins = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", array_heidelberg_cnvp_bins)) |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", array_heidelberg_cnvp_bins)) |> 
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_heidelberg_cnvp_bins) ) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()


rm(tmp)



## Heidelberg 12.8 bin-based tumor purity calls ----


tmp <- readRDS("cache/analysis_tumor_purity_EPIC_bin-based.Rds") |> 
  assertr::verify(!is.na(array_methylation_bins_1p19q_median.lfc)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_median.lfc)) |> 
  assertr::verify(!is.na(array_methylation_bins_1p19q_purity)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_purity)) |> 
  assertr::verify(!is.na(array_methylation_bins_1p19q_sd)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_sd)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275) # heidi still running
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)



## Heidelberg 12.8 CNVP ongene scores ----


tmp <- c(
  list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "*.detail.txt", recursive = TRUE)
) |>
  data.frame(tmp_filename = _) |>
  dplyr::mutate(tmp_filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", tmp_filename)) |>
  assertr::verify(file.exists(tmp_filename)) |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_filename)) |>
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_CNVPoncogeneScores_csv (tmp_filename, "array_mnp_CNVP_12.8_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(tmp_filename = NULL) |> 
  assertr::verify(!is.na(`array_mnp_CNVP_12.8_CDKN2A/B`)) |> 
  assertr::verify(is.numeric(`array_mnp_CNVP_12.8_CDKN2A/B`)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(`array_mnp_CNVP_12.8_CDKN2A/B`)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()

rm(tmp)



## Heidelberg 12.8 predictMGMT ----


tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "_mgmt.csv", recursive = TRUE) |> 
  data.frame(tmp_heidelberg_mgmt_report = _) |> 
  dplyr::mutate(tmp_heidelberg_mgmt_report = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", tmp_heidelberg_mgmt_report)) |>
  assertr::verify(file.exists(tmp_heidelberg_mgmt_report)) |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_heidelberg_mgmt_report)) |>
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_predictMGMT_csv(tmp_heidelberg_mgmt_report, "array_mnp_MGMT_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(tmp_heidelberg_mgmt_report = NULL) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_MGMT_Estimated)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()

rm(tmp)



## Heidelberg 12.8 rs_gender ----


tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "*.mix_gender.csv", recursive = TRUE) |> 
  data.frame(tmp_heidelberg_rs_gender_report = _) |> 
  dplyr::mutate(tmp_heidelberg_rs_gender_report = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", tmp_heidelberg_rs_gender_report)) |>
  assertr::verify(file.exists(tmp_heidelberg_rs_gender_report)) |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_heidelberg_rs_gender_report)) |> 
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_RsGender_csv(tmp_heidelberg_rs_gender_report, "array_mnp_rsGender_12.8_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(tmp_heidelberg_rs_gender_report = NULL)  |> 
  assertr::verify(!is.na(array_mnp_rsGender_12.8_chrYintensity) & is.numeric(array_mnp_rsGender_12.8_chrYintensity)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 275)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)




## ++ below: re-build because mvalue normalisation ++ ----

## Median methylation levels ----


tmp <- readRDS("cache/analysis_median_methylation.Rds") |> 
  dplyr::rename_with( ~ paste0("array_", .x)) |> 
  assertr::verify(!is.na(array_median.overall.methylation)) |> 
  assertr::verify(!is.na(array_median.glass_nl_supervised.methylation)) |> 
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 215) # only HQ samples - should become 210?
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))


rm(tmp)




## A_IDH_HG__A_IDH_LG_lr__lasso_fit ----



tmp <- readRDS(file="cache/analysis_A_IDH_HG__A_IDH_LG_lr__lasso_fit.Rds") |> 
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 215) # only HQ samples - should become 210?
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)



## unsupervised PCA ----


tmp <- readRDS(file="cache/analysis_unsupervised_PCA_GLASS-OD_x.Rds") |> 
  #dplyr::rename_with( ~ paste0("array_", .x)) |> 
  assertr::verify(!is.na(array_PC1)) |> 
  assertr::verify(!is.na(array_PC2)) |> 
  assertr::verify(!is.na(array_PC3)) |> 
  assertr::verify(!is.na(array_PC163)) |> 
  assertr::verify(!is.na(array_PC210)) |> 
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 210)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)



## unsupervised PCA [GLASS-OD + GLASS-NL combi] ----

# tmp <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_GLASS-NL_combined.Rds") |> 
#   dplyr::rename_with(~ gsub("^PC","PC.GLASS_OD_NL_combined.",.x), .cols = matches("^PC[0-9]", perl = T)) |> 
#   dplyr::filter(sentrix_id %in% glass_od.metadata.array_samples$sentrix_id) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == 163)
#     return(.)
#   })()
# 
# 
# glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
#   dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))
# rm(tmp)

tmp <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_GLASS-NL_combined_no_1P19Q.Rds") |>
  dplyr::rename_with( ~ paste0("array_", .x)) |> 
  dplyr::rename_with(~ gsub("^array_PC","array_PC.GLASS_OD_NL_combined_excl_1P19Q.",.x), .cols = matches("^array_PC[0-9]", perl = T)) |>
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 163)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)



## epiTOC2 ----


#' array_epiTOC2_tnsc: the estimated cumulative number of stem-cell divisions per stem-cell per year and per sample using the full epiTOC2 model
#' array_epiTOC2_tnsc2: the estimated cumulative number of stem-cell divisions per stem-cell per year and per sample using an approximation of epiTOC2 which assumes all epiTOC2 CpGs have beta-values exactly 0 in the fetal stage
#' array_epiTOC2_hypoSC: the HypoClock score over the 678 solo-WCGWs - QC associated?
#' array_epiTOC2_pcgtAge: this is the mitotic-score obtained using our previous epiTOC model


tmp <- readRDS("cache/analysis_EPITOC2.Rds") |> 
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 215)# only HQ samples - should become 210?
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)




## dnaMethyAge ----


tmp <- readRDS("cache/analysis_dnaMethyAge.Rds") |>
  dplyr::rename_with( ~ paste0("array_", .x)) |> 
  dplyr::mutate(array_dnaMethyAge__epiTOC2 = NULL) |> # superseded by `array_epiTOC2_tnsc`
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 215)# only HQ samples - should become 210?
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)


#plot(glass_od.metadata.array_samples$array_epiTOC2_hypoSC , glass_od.metadata.array_samples$dnaMethyAge__epiTOC2)
#plot(glass_od.metadata.array_samples$array_epiTOC2_pcgtAge , glass_od.metadata.array_samples$dnaMethyAge__epiTOC2)
#plot(glass_od.metadata.array_samples$array_epiTOC2_tnsc , glass_od.metadata.array_samples$dnaMethyAge__epiTOC2) ~= 1.0
#plot(glass_od.metadata.array_samples$array_epiTOC2_tnsc2 , glass_od.metadata.array_samples$dnaMethyAge__epiTOC2)


## RepliTali ----
#' scripts/analysis_EPITOC2.R


tmp <- readRDS("cache/analysis_RepliTali.Rds") |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::rename(array_RepliTali = RepliTali) |> 
  assertr::verify(is.numeric(array_RepliTali)) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 510)# only HQ samples - should become 210?
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)


## GLASS-NL median methylation 1300 probes signature ----


tmp <- readRDS(file="cache/analysis_progression_signatures__GLASS-NL_prim_rec_signature.Rds") |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 510)# only HQ samples - should become 210?
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)



## GLASS-OP prog. signatures ----
# 
# tmp <- readRDS(file="cache/analysis_progression_signatures__GLASS-OD__g2_g3.Rds")
# 
# glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
#   dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))
# rm(tmp)
# 
# 
# tmp <- readRDS(file="cache/analysis_progression_signatures__GLASS-OD__primary_recurrence.Rds")
# 
# 
# glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
#   dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))
# rm(tmp)
# 


# cleanup db connection ----

DBI::dbDisconnect(metadata.db.con)
rm(metadata.db.con)



