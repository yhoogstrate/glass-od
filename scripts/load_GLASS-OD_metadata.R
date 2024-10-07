#!/usr/bin/env R

# load libs, config & db ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')
#source('scripts/load_palette.R')
#source('scripts/load_themes.R')


metadata.db.con <- DBI::dbConnect(RSQLite::SQLite(), "../glass-od-clinical-database/glass-od-clinical-database.db")

# check whether views exist
tmp <- DBI::dbListTables(metadata.db.con) |> 
  as.data.frame() |> 
  dplyr::rename(table = 1) |> 
  dplyr::filter(grepl("^view_", table)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) >= 5)
    return(.)
  })()

rm(tmp)



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
    assertthat::assert_that(nrow(.) == 127 + 57 + 7 + 4 + 20 + 21 + 7)
    return(.)
  })()  |> 
  assertr::verify(is.na(patient_diagnosis_date) | grepl("^[0-9]{1,2} (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) [0-9]{4}$", patient_diagnosis_date)) |> 
  assertr::verify(is.na(patient_IDH_mutation) | (patient_IDH_mutation %in% c(
    "IDH1 R132H",
    "IDH1 R132G",
    "IDH1 R132S",
    "IDH1 R132L",
    
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
    assertthat::assert_that(nrow(.) ==  126 + 57 + 7 + 4 + 20 + 21 + 7) # + 10x astro
    return(.)
  })() |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> # 7 non(-canonical) codels
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 126 + 57 + 7 + 4 + 20 + 21 + 7) # + 10x astro
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
  assertr::verify((patient_id %in% c(
    "0026", "0027", "0056", "0063", "0076", "0085", "0093", "0096", "0097",
    "0098", "0131", "0132", "0133", "GLSS-LX-0170"
    
  ) & patient_suspected_noncodel == T) | (patient_suspected_noncodel == F)) |> 
  dplyr::mutate(patient_id = as.factor(patient_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_OD_VALIDATION_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + 6) # + 57 oligosarcoma-paper + 1x vali + 1x catnon - 3 resections without arrays (to complete db)
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
    assertthat::assert_that(nrow(.) == (2 * (CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_OD_VALIDATION_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES)))
    return(.)
  })() |> # equals in-pipe stopifnot(nrow(.) == 404)
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
    assertthat::assert_that(nrow(.) == (CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_OD_VALIDATION_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES))
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
    assertthat::assert_that(sum(is.na(.$arraychip_date)) <= (6 + 57 + 17 + 4 + 25 + 21))
    return(.)
  })() |> 
  
  assertr::verify(is.na(patient_diagnosis_date) | grepl("^[0-9]{1,2} (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) [0-9]{4}$", patient_diagnosis_date)) |> 
  assertr::verify(is.na(patient_birth_date) | grepl("^[0-9]{1,2} (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) [0-9]{4}$", patient_birth_date)) |> 
  assertr::verify(is.na(patient_last_follow_up_date) | grepl("^[0-9]{1,2} (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) [0-9]{4}$", patient_last_follow_up_date)) |> 
  assertr::verify(is.na(resection_date) | grepl("^[0-9]{1,2} (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) [0-9]{4}$", resection_date)) |> 
  assertr::verify(is.na(arraychip_date) | grepl("^[0-9]{1,2} (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) [0-9]{4}$", arraychip_date)) |> 
  
  dplyr::mutate(patient_diagnosis_date = as.Date(patient_diagnosis_date, format = "%d %b %Y")) |> 
  dplyr::mutate(patient_birth_date = as.Date(patient_birth_date, format = "%d %b %Y")) |> 
  dplyr::mutate(patient_last_follow_up_date = as.Date(patient_last_follow_up_date, format = "%d %b %Y")) |> 
  dplyr::mutate(resection_date = as.Date(resection_date, format = "%d %b %Y")) |> 
  dplyr::mutate(arraychip_date = as.Date(arraychip_date, format = "%d %b %Y")) |> 
  
  assertr::verify(is.na(patient_last_follow_up_event) | grepl("^yes|no$", patient_last_follow_up_event)) |> 
  assertr::verify(is.na(resection_treatment_status_chemo) | grepl("^yes|no$", resection_treatment_status_chemo)) |> 
  assertr::verify(is.na(resection_treatment_status_radio) | grepl("^yes|no$", resection_treatment_status_radio)) |> 
  
  dplyr::mutate(patient_last_follow_up_event = patient_last_follow_up_event == "yes") |> 
  dplyr::mutate(resection_treatment_status_chemo = resection_treatment_status_chemo == "yes") |> 
  dplyr::mutate(resection_treatment_status_radio = resection_treatment_status_radio == "yes") |> 
  
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
    assertthat::assert_that(nrow(.) == (CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_OD_VALIDATION_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES))
    return(.)
  })() |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id' = 'array_sentrix_id'), suffix=c('','')) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_OD_VALIDATION_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES))
    return(.)
  })() |> 
  assertr::verify(!is.na(resection_id)) |> 
  assertr::verify(
    (array_channel_green_filesize < 10000000 & arraychip_version == "450k") |
    (array_channel_green_filesize > 10000000 & arraychip_version == "EPICv1") |
    (arraychip_version %in% c("450k", "EPICv1") == F)
  )
  


rm(tmp)


## c. calc survival ----


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::mutate(time_between_resection_and_last_event = as.Date(patient_last_follow_up_date) - as.Date(resection_date)) |> 
  dplyr::mutate(patient_last_follow_up_event = 1 * patient_last_follow_up_event)  # from boolean to 1 = T , 0 = F

#glass_od.metadata.array_samples |> 
#  dplyr::select(array_sentrix_id, resection_id, time_between_resection_and_last_event, resection_date, patient_last_follow_up_date, patient_last_follow_up_event) |> 
#  View()




## HE coupes ----


tmp <- DBI::dbReadTable(metadata.db.con, 'HE_coupes') |> 
  dplyr::rename_with( ~ paste0("HE_coupe_", .x), .cols=!matches("^resection_id$",perl = T)) |> 
  dplyr::mutate(file_exists = file.exists(paste0("data/GLASS_OD/Stainings/H&E Slides/",HE_coupe_filename))) |> 
  assertr::verify(file_exists) |> 
  dplyr::mutate(file_exists = NULL) |> 
  dplyr::mutate(join_prefix = gsub("\\.ndpi$", "", HE_coupe_filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 233)
    
    print(length(unique(.$resection_id)))
    assertthat::assert_that(length(unique(.$resection_id)) == 193)
    return(.)
  })() |> 
  dplyr::full_join( # scan for files on disc but missing in the DB
    list.files(path = "data/GLASS_OD/Stainings/H&E Slides/", pattern = "*.ndpi$", recursive = TRUE) |> 
      data.frame(HE_coupe_filename = _) |> 
      (function(.) {
        print(dim(.))
        assertthat::assert_that(nrow(.) == 233)
        return(.)
      })(),
    by=c('HE_coupe_filename'='HE_coupe_filename')
  ) |> 
  assertr::verify(!is.na(HE_coupe_md5sum))




tmp.thumbnails_full <- list.files(path = "output/figures/HE_thumbnails/full/", pattern = "*.png$", recursive = TRUE) |>
  data.frame(HE_coupe_thumbnail_full = _) |>
  dplyr::mutate(HE_coupe_thumbnail_full = paste0("output/figures/HE_thumbnails/full/", HE_coupe_thumbnail_full)) |>
  dplyr::mutate(join_prefix = gsub("^.+/(.+)_[0-9]+x[0-9]+px\\.png$", "\\1", HE_coupe_thumbnail_full)) |> 
  assertr::verify(!duplicated(join_prefix)) |> 
  assertr::verify(join_prefix %in% tmp$join_prefix) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 233)
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
    assertthat::assert_that(nrow(.) == 210)
    return(.)
  })()



tmp <- tmp |> 
  dplyr::left_join(tmp.thumbnails_full, by=c('join_prefix'='join_prefix'), suffix=c('','')) |> 
  dplyr::left_join(tmp.thumbnails_roi, by=c('join_prefix'='join_prefix'), suffix=c('','')) |> 
  dplyr::mutate(join_prefix = NULL) |> 
  dplyr::filter(HE_coupe_primary_coupe == "yes") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 193)
    return(.)
  })() |> 
  assertr::verify(!is.na(HE_coupe_thumbnail_full)) |> 
  assertr::verify(!is.na(HE_coupe_thumbnail_roi))


rm(tmp.thumbnails_full, tmp.thumbnails_roi)






glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('resection_id'='resection_id'), suffix=c('',''))


rm(tmp)



## Ki-67 stainings ----



tmp <- DBI::dbReadTable(metadata.db.con, 'stainings_KI67') |> 
  dplyr::arrange(resection_id, pathology_number) |> 
  dplyr::rename_with( ~ paste0("staining_KI67_", .x), .cols=!matches("^resection_id$",perl = T)) |> 
  dplyr::mutate(file_exists = file.exists(paste0("data/GLASS_OD/Stainings/KI67/",staining_KI67_filename))) |> 
  assertr::verify(file_exists) |> 
  assertr::verify(!is.na(staining_KI67_md5sum)) |> 
  assertr::verify(!duplicated(staining_KI67_md5sum)) |> 
  dplyr::mutate(file_exists = NULL) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 151)
    return(.)
  })()


tmp.verify <- list.files(path = "data/GLASS_OD/Stainings/KI67/", pattern = "*.ndpi$", recursive = TRUE) |> 
  data.frame(disc.filename = _) |> 
  assertr::verify(disc.filename %in% tmp$staining_KI67_filename)


tmp <- tmp |>
  dplyr::filter(staining_KI67_primary_coupe == "yes" & is.na(staining_KI67_reason_excluded)) |> 
  assertr::verify(!duplicated(resection_id)) |> 
  dplyr::mutate(file_exists = NULL) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (135))
    return(.)
  })() |> 
  assertr::verify(sum(resection_id %in% glass_od.metadata.array_samples$resection_id == F) == 1) # 0003-R1


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('resection_id'='resection_id'), suffix=c('',''))


rm(tmp, tmp.verify)




## Percentage detP probesC ----
#' from: scripts/analysis_percentage_detP_probes.R


tmp <- read.table("output/tables/percentage_detP_probes.txt") |> 
  assertr::verify(!is.na(array_percentage.detP.signi) & is.numeric(array_percentage.detP.signi)) |> 
  dplyr::rename(array_sentrix_id = array_sentrix_id) |> 
  assertr::verify(glass_od.metadata.array_samples$array_sentrix_id %in% array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (  (CONST_N_GLASS_OD_ALL_SAMPLES)  +
                                          CONST_N_CATNON_ALL_SAMPLES +
                                          (CONST_N_OD_VALIDATION_ALL_SAMPLES) +
                                          CONST_N_GLASS_NL_ALL_SAMPLES +
                                          (CONST_N_GSAM_ALL_SAMPLES) +
                                          194))
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
  assertr::verify(glass_od.metadata.array_samples |> dplyr::filter(arraychip_version == "EPICv1") |> dplyr::pull(`array_sentrix_id`) %in% array_sentrix_id) |> 
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 373)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify((arraychip_version == "EPICv1" &!is.na(array_qc.pca.comp1)) | arraychip_version != "EPICv1" ) |> 
  assertr::verify((arraychip_version == "EPICv1" &!is.na(array_qc.pca.detP.outlier)) | arraychip_version != "EPICv1" ) |> 
  assertr::verify(ifelse(resection_id == "0002-R1", (!is.na(array_qc.pca.detP.outlier)), T))



rm(tmp)



## minfi dyeCorrection stats ----


tmp <- readRDS("cache/analysis_dyeCorrection_rate.Rds") |> 
  dplyr::rename(array_sentrix_id = sample) |> 
  dplyr::rename_with(~paste0("array_minfi_dyeCorrection_", .x), .cols=!matches("^array_sentrix_id$",perl = T)) |> 
  assertr::verify(glass_od.metadata.array_samples$arraychip_version != "EPICv1" | glass_od.metadata.array_samples$array_sentrix_id %in% array_sentrix_id)


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))


rm(tmp)



## ewastools ----


tmp <- read.table("output/tables/ewastools.txt") |> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(glass_od.metadata.array_samples$array_sentrix_id %in% array_sentrix_id)


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))


rm(tmp)



## Heidelberg 11b4[+12.5] QC full ----

# quite some of these files contain odd N/A's



tmp.ls <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", pattern = "_qc_full.txt", recursive = TRUE)


tmp <- tmp.ls |> 
  data.frame(array_mnp_qc_report_full = _) |> 
  dplyr::mutate(array_mnp_qc_report_full = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", array_mnp_qc_report_full)) |>
  dplyr::mutate(array_mnp_qc_report_version = gsub("^.+qc_(v[^_\\/]+)[_/].+$","\\1", array_mnp_qc_report_full)) |>
  dplyr::mutate(array_sentrix_id = gsub("^.+[^0-9]([0-9]{10,12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", array_mnp_qc_report_full)) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_QCReport_csv(array_mnp_qc_report_full, "array_qc_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + 1 + 5 +
                              89 # validation files completed so far
                              ) # only one of the validation set so far & 1 CATNON
    return(.)
  })()




# glass_od.metadata.array_samples |> 
#   dplyr::filter(	
#     arraychip_version != "450k" ) |> 
#   dplyr::select(isolation_id, array_sentrix_id, 
#                 patient_study_name,
#                 arraychip_version,
#                 array_mnp_predictBrain_v2.0.1_cal_class,
#                 array_mnp_predictBrain_v12.8_cal_class) |> 
#   View()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_qc_report_version) | arraychip_version != "EPICv1") |> 
  assertr::verify(!is.na(array_mnp_qc_report_full) | arraychip_version != "EPICv1") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + CONST_N_OD_VALIDATION_ALL_SAMPLES) # only one of the validation set so far & 1 CATNON
    return(.)
  })() |> 
  dplyr::mutate(array_mnp_qc_report_full = NULL)


rm(tmp)



## Heidelberg 11b4[+12.5] reportBrain files ----
#' needed to correlate LGC


tmp.ls <- c(list.files(
  path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", pattern = "*_scores_cal.csv", recursive = TRUE)
)


tmp <- tmp.ls |> 
  data.frame(tmp_filename = _) |>
  dplyr::mutate(array_mnp_predictBrain_filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", tmp_filename)) |>
  dplyr::mutate(array_mnp_predictBrain_version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", tmp_filename)) |>
  dplyr::mutate(array_sentrix_id = gsub("^.+[^0-9]([0-9]{10,12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_filename)) |>
  dplyr::mutate(tmp_filename = NULL) |> 
  dplyr::mutate(tmp_uniq = paste0(array_sentrix_id, array_mnp_predictBrain_version)) |> 
  assertr::verify(!duplicated(tmp_uniq)) |>
  dplyr::mutate(tmp_uniq = NULL) |> 
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
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + 1 + 5+
                              89 # validaton files completed so far
                            )
    return(.)
  })()



glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_lr_v2.0.1) | arraychip_version != "EPICv1") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + CONST_N_OD_VALIDATION_ALL_SAMPLES) # only one of the validation set so far & 1 CATNON
    return(.)
  })()


rm(tmp, tmp.ls)




## Heidelberg 11b4[+12.5] rs_gender ----


tmp.ls <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", pattern = "*.mix_gender.csv", recursive = TRUE)

tmp <- tmp.ls |> 
  data.frame(tmp_heidelberg_rs_gender_report = _) |> 
  dplyr::mutate(tmp_heidelberg_rs_gender_report = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v11b4_sample_report__v3.3__125/", tmp_heidelberg_rs_gender_report)) |>
  assertr::verify(file.exists(tmp_heidelberg_rs_gender_report)) |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+[^0-9]([0-9]{10,12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_heidelberg_rs_gender_report)) |> 
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
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + 1 + 5 +
                              89 # validation files completed so far
                            )
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_rsGender_11b4_predicted) | arraychip_version != "EPICv1")


rm(tmp, tmp.ls)



## Heidelberg 12.8 reportBrain files ----


tmp.ls <- c(list.files(
  path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "*_scores_cal.csv", recursive = TRUE)
)

tmp <- tmp.ls |>
  data.frame(tmp_filename = _) |>
  dplyr::mutate(array_mnp_predictBrain_filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", tmp_filename)) |>
  dplyr::mutate(array_mnp_predictBrain_version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", tmp_filename)) |>
  dplyr::mutate(array_sentrix_id = gsub("^.+[^0-9]([0-9]{10,12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_filename)) |>
  dplyr::mutate(tmp_uniq = paste0(array_sentrix_id, array_mnp_predictBrain_version)) |> 
  assertr::verify(!duplicated(tmp_uniq)) |>
  dplyr::mutate(tmp_uniq = NULL) |> 
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
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + 1 + 5 +
                              89 # validation files not yet completed
                            )
    return(.)
  })()



glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_class) | arraychip_version != "EPICv1") |> 
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_LG_lr_v12.8) | arraychip_version != "EPICv1") |> 
  (function(.) {
    print(paste0(sum(!is.na(.$array_mnp_predictBrain_v12.8_cal_class)), "/", nrow(.)))
    assertthat::assert_that(sum(!is.na(.$array_mnp_predictBrain_v12.8_cal_class)) == (CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + 1  + 5 +
             89 # validation files not yet completed
                                                                      )) # 0118-R4 still failing
    return(.)
  })() 




rm(tmp, tmp.ls)




## Heidelberg 12.8 Frozen ~ FFPE status ----


tmp.ls <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/",  pattern = "_ffpe_frozen.txt", recursive = TRUE)

tmp <- tmp.ls |> 
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
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + 1 + 5 +
                              89) # validation files not yet completed
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_QC_v12.8_predicted_array_type) | arraychip_version != "EPICv1") |> 
  assertr::verify(!is.na(array_mnp_QC_v12.8_predicted_sample_type) | arraychip_version != "EPICv1" ) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + CONST_N_OD_VALIDATION_ALL_SAMPLES  )
    return(.)
  })()


rm(tmp.ls, tmp)



## Heidelberg 12.8 CNVP segment files ----


tmp <- query_mnp_12.8_CNVP_segment_csv(
  "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", 
  CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + 1 + 89 + 5, 
  glass_od.metadata.array_samples$array_sentrix_id,
  "array_mnp_CNVP_v12.8_v5.2_"
  )


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_CNVP_v12.8_v5.2_CNVP_segments)  | arraychip_version != "EPICv1")  |> 
  assertr::verify(!is.na(array_mnp_CNVP_v12.8_v5.2_CNVP_version) | arraychip_version != "EPICv1" ) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + CONST_N_OD_VALIDATION_ALL_SAMPLES)
    return(.)
  })()


rm(tmp)




## Heidelberg 12.8 CNVP bins files ----


tmp <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "*.bins.igv", recursive = TRUE) |> 
  data.frame(array_heidelberg_cnvp_bins = _) |> 
  dplyr::mutate(array_mnp_CNVP_v12.8_v5.2_CNVP_bins = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", array_heidelberg_cnvp_bins)) |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+[^0-9]([0-9]{10,12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", array_heidelberg_cnvp_bins)) |> 
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + 1 + 89 + 5)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_heidelberg_cnvp_bins) | arraychip_version != "EPICv1") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + CONST_N_OD_VALIDATION_ALL_SAMPLES)
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
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + 1 + 89 + 5)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)



## Heidelberg 12.8 CNVP ongene scores ----


tmp.ls <- c(
  list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "*.detail.txt", recursive = TRUE)
)

tmp <- tmp.ls |> 
  data.frame(tmp_filename = _) |>
  dplyr::mutate(tmp_filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", tmp_filename)) |>
  assertr::verify(file.exists(tmp_filename)) |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+[^0-9]([0-9]{10,12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_filename)) |>
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
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + 1 + 89 + 5)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(`array_mnp_CNVP_12.8_CDKN2A/B`) | arraychip_version != "EPICv1") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + CONST_N_OD_VALIDATION_ALL_SAMPLES)
    return(.)
  })()


rm(tmp, tmp.ls)




## Heidelberg 12.8 predictMGMT ----


tmp.ls <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "_mgmt.csv", recursive = TRUE)


tmp <- tmp.ls |> 
  data.frame(tmp_heidelberg_mgmt_report = _) |> 
  dplyr::mutate(tmp_heidelberg_mgmt_report = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", tmp_heidelberg_mgmt_report)) |>
  assertr::verify(file.exists(tmp_heidelberg_mgmt_report)) |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+[^0-9]([0-9]{10,12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_heidelberg_mgmt_report)) |>
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
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + 1 + 89 + 5)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_MGMT_Estimated)| arraychip_version != "EPICv1") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + CONST_N_OD_VALIDATION_ALL_SAMPLES)
    return(.)
  })()


rm(tmp, tmp.ls)




## Heidelberg 12.8 rs_gender ----


tmp.ls <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "*.mix_gender.csv", recursive = TRUE)


tmp <- tmp.ls |> 
  data.frame(tmp_heidelberg_rs_gender_report = _) |> 
  dplyr::mutate(tmp_heidelberg_rs_gender_report = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", tmp_heidelberg_rs_gender_report)) |>
  assertr::verify(file.exists(tmp_heidelberg_rs_gender_report)) |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+[^0-9]([0-9]{10,12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_heidelberg_rs_gender_report)) |> 
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
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + 1 + 89 + 5)
    return(.)
  })()



glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_rsGender_12.8_dist)| arraychip_version != "EPICv1") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + CONST_N_OD_VALIDATION_ALL_SAMPLES)
    return(.)
  })()



rm(tmp, tmp.ls)


## MethylScape Bethesda Classifier v2 ----


tmp <- readRDS( "cache/MethylScape_Bethesda_Classifier_v2.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_INCLUDED_SAMPLES)
    return(.)
  })() |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+_([0-9]+_[0-9A-Z]+)$","\\1",array_sentrix_id)) |> 
  assertr::verify(!is.na(array_sentrix_id)) |> 
  assertr::verify(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  assertr::verify(!is.na(array_methylscape_bethesda_class)) |> 
  assertr::verify(is.numeric(array_methylscape_bethesda_class_score))
  #dplyr::mutate(array_methylscape_bethesda_class = ifelse(array_methylscape_bethesda_class == "O_SARC_IDH", "OLIGOSARC_IDH" , array_methylscape_bethesda_class))




# if err, then:
# 
# tmp.ls <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays - MethylScape Bethesda classifier/", pattern = "*.html", recursive = TRUE)
# 
# tmp <- parse_MethylScape_Bv2(paste0("data/GLASS_OD/DNA Methylation - EPIC arrays - MethylScape Bethesda classifier/",tmp.ls)) |>
#  (function(.) {
#    print(dim(.))
#    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_INCLUDED_SAMPLES)
#    return(.)
#  })()
# 
# saveRDS(tmp, "cache/MethylScape_Bethesda_Classifier_v2.Rds")
# rm(tmp.ls)



glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_CATNON_ALL_SAMPLES + CONST_N_OD_VALIDATION_ALL_SAMPLES)
    return(.)
  })()



rm(tmp)



## RFpurity ----


tmp <- readRDS("cache/analysis_RFpurity.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 373) # CONST_N_SAMPLES
    return(.)
  })() |> 
  assertr::verify(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id)


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) 

rm(tmp)




## ++ below: re-build because mvalue normalisation ++ ----

## Median methylation levels ----


tmp <- readRDS("cache/analysis_median_methylation.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_SAMPLES) # CONST_N_SAMPLES
    return(.)
  })() |> 
  dplyr::rename_with( ~ paste0("array_", .x)) |> 
  assertr::verify(!is.na(array_median.overall.methylation)) |> 
  assertr::verify(!is.na(array_median.glass_nl_supervised.methylation)) |> 
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_INCLUDED_SAMPLES + CONST_N_OD_VALIDATION_INCLUDED_SAMPLES + CONST_N_CATNON_ALL_SAMPLES)  # CONST_N_GLASS_OD_INCLUDED_SAMPLES
    return(.)
  })()

glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))  |> 
  assertr::verify(ifelse(resection_id == "0002-R1", (!is.na(array_median.glass_nl_supervised.methylation)), T))



rm(tmp)




## A_IDH_HG__A_IDH_LG_lr__lasso_fit ----


tmp <- readRDS(file="cache/analysis_A_IDH_HG__A_IDH_LG_lr__lasso_fit.Rds") |> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_INCLUDED_SAMPLES + CONST_N_OD_VALIDATION_INCLUDED_SAMPLES + CONST_N_GSAM_INCLUDED_SAMPLES + CONST_N_CATNON_ALL_SAMPLES)
    return(.)
  })() |>
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_INCLUDED_SAMPLES + CONST_N_OD_VALIDATION_INCLUDED_SAMPLES + CONST_N_CATNON_ALL_SAMPLES)
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(ifelse(resection_id == "0002-R1", (!is.na(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)), T))



rm(tmp)



## unsupervised PCA ----


tmp <- readRDS(file="cache/analysis_unsupervised_PCA_GLASS-OD_x.Rds") |> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_INCLUDED_SAMPLES)
    return(.)
  })() |>
  
  assertr::verify(!is.na(array_PC1)) |> 
  assertr::verify(!is.na(array_PC2)) |> 
  assertr::verify(!is.na(array_PC3)) |> 
  assertr::verify(!is.na(array_PC163)) |> 
  assertr::verify(!is.na(array_PC212)) |> 
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_INCLUDED_SAMPLES) # CONST_N_GLASS_OD_INCLUDED_SAMPLES
    return(.)
  })()


tmp$array_sentrix_id[tmp$array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id == F]


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(ifelse(resection_id == "0002-R1", (!is.na(array_PC2)), T))


rm(tmp)


## unsupervised PCA [GL-OD + OD-vali] ----


tmp <- readRDS(file="cache/analysis_unsupervised_PCA_GLASS-OD_AND_VALIDATION_x.Rds") |> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_GLASS_OD_INCLUDED_SAMPLES + CONST_N_OD_VALIDATION_INCLUDED_SAMPLES))
    return(.)
  })() |>
  assertr::verify(!is.na(array_GLASS_OD_VALIDATION_PC1)) |> 
  assertr::verify(!is.na(array_GLASS_OD_VALIDATION_PC2)) |> 
  assertr::verify(!is.na(array_GLASS_OD_VALIDATION_PC298)) |> 
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_GLASS_OD_INCLUDED_SAMPLES + CONST_N_OD_VALIDATION_INCLUDED_SAMPLES))
    return(.)
  })()


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(ifelse(resection_id == "0002-R1", (!is.na(array_PC2)), T))


rm(tmp)





## unsupervised PCA [GLASS-OD + GLASS-NL combi] ----


# tmp <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_GLASS-NL_combined_no_1P19Q.Rds") |>
#   dplyr::rename_with( ~ paste0("array_", .x)) |> 
#   dplyr::rename_with(~ gsub("^array_PC","array_PC.GLASS_OD_NL_combined_excl_1P19Q.",.x), .cols = matches("^array_PC[0-9]", perl = T)) |>
#   dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |>
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == 163)
#     return(.)
#   })()
# 
# 
# glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
#   dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))
# 
# rm(tmp)



## epiTOC2 ----


#' array_epiTOC2_tnsc: the estimated cumulative number of stem-cell divisions per stem-cell per year and per sample using the full epiTOC2 model
#' array_epiTOC2_tnsc2: the estimated cumulative number of stem-cell divisions per stem-cell per year and per sample using an approximation of epiTOC2 which assumes all epiTOC2 CpGs have beta-values exactly 0 in the fetal stage
#' array_epiTOC2_hypoSC: the HypoClock score over the 678 solo-WCGWs - QC associated?
#' array_epiTOC2_pcgtAge: this is the mitotic-score obtained using our previous epiTOC model


tmp <- readRDS("cache/analysis_EPITOC2.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_SAMPLES)
    return(.)
  })() |> 
  
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id)


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(ifelse(resection_id == "0002-R1", (!is.na(array_epiTOC2_tnsc)), T))

rm(tmp)




## dnaMethyAge ----


tmp <- readRDS("cache/analysis_dnaMethyAge.Rds") |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_SAMPLES)
    return(.)
  })() |> 
  dplyr::rename_with( ~ paste0("array_", .x)) |> 
  dplyr::mutate(array_dnaMethyAge__epiTOC2 = NULL) |> # superseded by `array_epiTOC2_tnsc`
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id)



glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(ifelse(resection_id == "0002-R1", (!is.na(array_dnaMethyAge__HannumG2013)), T))

rm(tmp)


#plot(glass_od.metadata.array_samples$array_epiTOC2_hypoSC , glass_od.metadata.array_samples$dnaMethyAge__epiTOC2)
#plot(glass_od.metadata.array_samples$array_epiTOC2_pcgtAge , glass_od.metadata.array_samples$dnaMethyAge__epiTOC2)
#plot(glass_od.metadata.array_samples$array_epiTOC2_tnsc , glass_od.metadata.array_samples$dnaMethyAge__epiTOC2) ~= 1.0
#plot(glass_od.metadata.array_samples$array_epiTOC2_tnsc2 , glass_od.metadata.array_samples$dnaMethyAge__epiTOC2)


## RepliTali ----
#' scripts/analysis_EPITOC2.R


tmp <- readRDS("cache/analysis_RepliTali.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_SAMPLES)
    return(.)
  })() |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |>
  dplyr::rename(array_RepliTali = RepliTali) |> 
  assertr::verify(is.numeric(array_RepliTali))



glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(ifelse(resection_id == "0002-R1", (!is.na(array_RepliTali)), T))

rm(tmp)


## GLASS-NL median methylation 1300 probes signature ----
#' generated by: scripts/analysis_progression_signatures.R

tmp <- readRDS(file="cache/analysis_progression_signatures__GLASS-NL_prim_rec_signature.Rds") |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_SAMPLES)
    return(.)
  })() |> 
  dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id)


glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(ifelse(resection_id == "0002-R1", (!is.na(array_GLASS_NL_g2_g3_sig)), T))

rm(tmp)



## GLASS-OP prog. signatures ----
#' generated by: scripts/analysis_progression_signatures.R


# tmp <- readRDS(file="cache/analysis_progression_signatures__GLASS-OD__g2_g3.Rds") |>
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == 505) # CONST_N_SAMPLES
#     return(.)
#   })() |> 
#   dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |>
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == 210) # CONST_N_GLASS_OD_INCLUDED_SAMPLES
#     return(.)
#   })()
# 
# glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
#   dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))
# 
# rm(tmp)
# 
# 
# 
# tmp <- readRDS(file="cache/analysis_progression_signatures__GLASS-OD__primary_recurrence.Rds") |>
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == 505) # CONST_N_SAMPLES
#     return(.)
#   })() |> 
#   dplyr::filter(array_sentrix_id %in% glass_od.metadata.array_samples$array_sentrix_id) |>
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == 210 ) # CONST_N_GLASS_OD_INCLUDED_SAMPLES
#     return(.)
#   })()
# 
# 
# glass_od.metadata.array_samples <- glass_od.metadata.array_samples |>
#   dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))
# 
# rm(tmp)




# 4. proteomics level ----


glass_od.metadata.proteomics <- DBI::dbReadTable(metadata.db.con, 'view_proteomics') |> 
  assertr::verify(resection_id %in% glass_od.metadata.array_samples$resection_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_PROTEOMICS)
    return(.)
  })()



tmp <-
  readxl::read_xlsx('data/GLASS_OD/Protein - Tobias Weiss/Broad overview samples sent to Tobias Weiss.xlsx') |>
  dplyr::rename(proteomics_id = `Sample ID`) |> 
  dplyr::filter(!is.na(proteomics_id)) |> # two empty rows
  dplyr::rename(proteomics_box_location__box_row_column = `Box location\r\n(Box, Row, Column)`) |> 
  dplyr::rename(proteomics_notes_zurich = `Note`) |> 
  dplyr::rename(proteomics_material_source = `Type of starting material`) |> 
  dplyr::select(proteomics_id, proteomics_box_location__box_row_column, proteomics_material_source, proteomics_notes_zurich) |> 
  assertr::verify(proteomics_id %in% glass_od.metadata.proteomics$proteomics_id) |> 
  assertr::verify(glass_od.metadata.proteomics$proteomics_id %in% proteomics_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_ALL_PROTEOMICS)
    return(.)
  })()



glass_od.metadata.proteomics <- glass_od.metadata.proteomics |> 
  dplyr::left_join(tmp, by=c('proteomics_id'='proteomics_id'), suffix=c('',''))

rm(tmp)


# not needed
#tmp <- readxl::read_xlsx('data/GLASS_OD/Protein - Tobias Weiss/Oligodendroglioma_protein matrices.xlsx')



# cleanup db connection ----


DBI::dbDisconnect(metadata.db.con)
rm(metadata.db.con)


