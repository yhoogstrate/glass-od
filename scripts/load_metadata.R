#!/usr/bin/env R

# load libs, config & db ----


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
    assertthat::assert_that(nrow(.) == 101)
    return(.)
  })()

patients_without_array_samples <- DBI::dbReadTable(metadata.db.con, 'view_check_patients_without_array_samples')
stopifnot(sort(patients_without_array_samples$patient_id) == c(1)) # x-checked, patients currently missing samples

glass_od.metadata.patients <- glass_od.metadata.patients |> 
  dplyr::filter(patient_id %in% patients_without_array_samples$patient_id == F) |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == 100)
    return(.)
  })() |> 
  dplyr::filter(is.na(reason_excluded)) |> # 7 non(-canonical) codels
  (function(.) {
    assertthat::assert_that(nrow(.) == 93)
    return(.)
  })()



rm(patients_without_array_samples)



# 2. resection level ----


glass_od.metadata.resections <- DBI::dbReadTable(metadata.db.con, 'view_resections') |> 
  dplyr::filter(patient_id %in% glass_od.metadata.patients$patient_id) |> 
  dplyr::filter(is.na(reason_excluded_resection)) |> 
  assertr::verify(patient_id %in% c(27, 56, 76, 93, 96, 97, 98) == F) |> # hard coded non-codels
  (function(.) {
    assertthat::assert_that(nrow(.) == 205)
    return(.)
  })() |> 
  assertr::verify(is.numeric(resection_number))



# 3. idat level ----
## a. load all idat files ----


glass_od.metadata.idats <- list.files(path = "data/GLASS_OD/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/GLASS_OD/", filename)) |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == (414 + 28))
    return(.)
  })() |> # equals in-pipe stopifnot(nrow(.) == 404)
  assertr::verify(file.exists(filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", filename)) |>
  dplyr::mutate(channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", filename)) |>
  tidyr::pivot_wider(id_cols = sentrix_id, names_from = channel, values_from = c(filename)) |>
  dplyr::rename(channel_green = Grn) |>
  dplyr::rename(channel_red = Red) |>
  (function(.) {
    assertthat::assert_that(nrow(.) == ((414 + 28) / 2))
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
    assertthat::assert_that(nrow(.) == 221)
    return(.)
  })() |> 
  dplyr::left_join(glass_od.metadata.array_samples, by=c('sentrix_id' = 'sentrix_id'), suffix=c('','')) |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == 221)
    return(.)
  })() |> 
  assertr::verify(!is.na(resection_id))



rm(glass_od.metadata.array_samples)



## add percentage detP probes ----

# from: scripts/analysis_percentage_detP_probes.R

tmp <- read.table("output/tables/percentage_detP_probes.txt") |> 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id)

glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(percentage.detP.signi))




## heidelberg reportBrain files ----


w <- function(fn, prefix) {
  
  a = read.csv(fn) |> 
    tibble::column_to_rownames('X') |> 
    `colnames<-`('pval') |> 
    dplyr::arrange(-pval) |> 
    dplyr::mutate(pval = round(pval * 100,1)) 
  
  top <- a |> 
    tibble::rownames_to_column('class') |> 
    dplyr::slice_head(n=1) |> 
    dplyr::pull(class)
  
  a<- a |> 
    t() |> 
    as.data.frame() |> 
    dplyr::mutate(class = top) |> 
    dplyr::rename_with( ~ paste0(prefix,"cal_", .x)) 
  
  return(a)
}




tmp <- c(
  list.files(path = "data/GLASS_OD/Heidelberg_classifier_output/", pattern = "*_scores_cal.csv", recursive = TRUE)
) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/GLASS_OD/Heidelberg_classifier_output/", filename)) |>
  dplyr::mutate(basename = gsub("^.+/([^/]+)$", "\\1", filename)) |>
  dplyr::mutate(version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", filename)) |>
  dplyr::select(-basename) |>
  tidyr::pivot_wider(id_cols = sentrix_id, names_from = version, values_from = c(filename), names_prefix = "heidelberg_reportBrain_") |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = w(heidelberg_reportBrain_v12.5, "predictBrain_12.5_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp)



stopifnot(tmp$sentrix_id %in% glass_od.metadata.idats$sentrix_id)
stopifnot(duplicated(tmp$sentrix_id) == F)


stopifnot(nrow(glass_od.metadata.idats) == 221)


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))
rm(tmp, w)


stopifnot(nrow(glass_od.metadata.idats) == 221)


stopifnot(!is.na(glass_od.metadata.idats$heidelberg_reportBrain_v12.5))
stopifnot(!is.na(glass_od.metadata.idats$heidelberg_reportBrain_v2.0.1))





## heidelberg CNV files ----

### segment file ----


tmp <- list.files(path = "data/GLASS_OD/Heidelberg_classifier_output/", pattern = "*.seg", recursive = TRUE) |> 
  data.frame(heidelberg_cnvp_segments = _) |> 
  dplyr::mutate(heidelberg_cnvp_segments = paste0("data/GLASS_OD/Heidelberg_classifier_output/", heidelberg_cnvp_segments)) |> 
  dplyr::mutate(heidelberg_cnvp_version = gsub("^.+/cnvp_([^/]+)/.+$","\\1", heidelberg_cnvp_segments)) |> 
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", heidelberg_cnvp_segments)) |> 
  assertr::verify(!is.na(duplicated(sentrix_id)))


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))
rm(tmp)


stopifnot(!is.na(glass_od.metadata.idats$heidelberg_cnvp_segments))
stopifnot(!is.na(glass_od.metadata.idats$heidelberg_cnvp_version))



### bins file ----


tmp <- list.files(path = "data/GLASS_OD/Heidelberg_classifier_output/", pattern = "*.bins.igv", recursive = TRUE) |> 
  data.frame(heidelberg_cnvp_bins = _) |> 
  dplyr::mutate(heidelberg_cnvp_bins = paste0("data/GLASS_OD/Heidelberg_classifier_output/", heidelberg_cnvp_bins)) |> 
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", heidelberg_cnvp_bins)) |> 
  assertr::verify(!is.na(duplicated(sentrix_id)))


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))
rm(tmp)


stopifnot(!is.na(glass_od.metadata.idats$heidelberg_cnvp_bins))



## heidelberg qc full ----


x <- function(fn, prefix) {
  a <- read.csv(fn,header=T, sep="\t") |> 
    dplyr::mutate(key=paste0(type, "_",name, "_",color, "_",check.if, "_",warning, "_",fail)) |> 
    dplyr::select(key, norm) |> 
    tibble::column_to_rownames('key') |> 
    `colnames<-`('value') |> 
    t() |> 
    as.data.frame() |> 
    dplyr::rename_with( ~ paste0(prefix, .x)) 
  
  return(a)
}


tmp <- list.files(path = "data/GLASS_OD/Heidelberg_classifier_output/", pattern = "_qc_full.txt", recursive = TRUE) |> 
  data.frame(heidelberg_qc_report_full = _) |> 
  dplyr::mutate(heidelberg_qc_report_full = paste0("data/GLASS_OD/Heidelberg_classifier_output/", heidelberg_qc_report_full)) |>
  dplyr::mutate(basename = gsub("^.+/([^/]+)$", "\\1", heidelberg_qc_report_full)) |>
  dplyr::mutate(heidelberg_qc_report_version = gsub("^.+qc_(v[^_\\/]+)[_/].+$","\\1", heidelberg_qc_report_full)) |>
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", heidelberg_qc_report_full)) |>
  dplyr::select(-basename) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = x(heidelberg_qc_report_full, "qc_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp)


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))
rm(tmp, x)

stopifnot(!is.na(glass_od.metadata.idats$heidelberg_qc_report_full))
glass_od.metadata.idats$heidelberg_qc_report_full <- NULL # already parsed



## predictMGMT ----


y <- function(fn, prefix) {

  a <- read.csv(fn,header=T) |> 
    dplyr::mutate(X=NULL, Status = NULL) |>  # also returns status if it is inconfident 
    dplyr::mutate(status = dplyr::case_when(
      Estimated < Cutoff & CI_Lower < Cutoff & CI_Upper < Cutoff ~ "unmethylated",
      Estimated > Cutoff & CI_Lower > Cutoff & CI_Upper > Cutoff ~ "methylated",
      T ~ as.character(NA)
    )) |> 
    dplyr::rename_with( ~ paste0(prefix, .x)) 
  
  return(a)
}



tmp <- list.files(path = "data/GLASS_OD/Heidelberg_classifier_output/", pattern = "_mgmt.csv", recursive = TRUE) |> 
  data.frame(heidelberg_mgmt_report = _) |> 
  dplyr::mutate(heidelberg_mgmt_report = paste0("data/GLASS_OD/Heidelberg_classifier_output/", heidelberg_mgmt_report)) |>
  dplyr::mutate(basename = gsub("^.+/([^/]+)$", "\\1", heidelberg_mgmt_report)) |>
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", heidelberg_mgmt_report)) |>
  dplyr::select(-basename) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = y(heidelberg_mgmt_report, "mgmt_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp)


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))
rm(tmp, y)

stopifnot(!is.na(glass_od.metadata.idats$heidelberg_mgmt_report))
glass_od.metadata.idats$heidelberg_mgmt_report <- NULL # already parsed



## rs_gender ----


z <- function(fn, prefix) {

  a <- read.csv(fn,header=T) |> 
    dplyr::mutate(idat = NULL, array = NULL) |>  # also returns status if it is inconfident 
    dplyr::mutate(predicted = ifelse(predicted == F, "F", predicted)) |> 
    dplyr::rename_with( ~ paste0(prefix, .x)) 
  
  return(a)
}


tmp <- list.files(path = "data/GLASS_OD/Heidelberg_classifier_output/", pattern = "*.mix_gender.csv", recursive = TRUE) |> 
  data.frame(heidelberg_rs_gender_report = _) |> 
  dplyr::mutate(heidelberg_rs_gender_report = paste0("data/GLASS_OD/Heidelberg_classifier_output/", heidelberg_rs_gender_report)) |>
  dplyr::mutate(basename = gsub("^.+/([^/]+)$", "\\1", heidelberg_rs_gender_report)) |>
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", heidelberg_rs_gender_report)) |>
  dplyr::select(-basename) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = z(heidelberg_rs_gender_report, "rs_gender_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp)


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))
rm(tmp, z)

stopifnot(!is.na(glass_od.metadata.idats$heidelberg_rs_gender_report))
glass_od.metadata.idats$heidelberg_rs_gender_report <- NULL # already parsed


## QC PCA ----




tmp <- readRDS("cache/unsupervised_qc_qc.outliers.Rds")
glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))

rm(tmp)



# cleanup db connection ----

DBI::dbDisconnect(metadata.db.con)
rm(metadata.db.con)


# add tumor purity calls? ----



