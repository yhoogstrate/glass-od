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


glass_od.metadata.idats <- list.files(path = "data/GLASS_OD/Methylation data - EPIC arrays/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/GLASS_OD/Methylation data - EPIC arrays/", filename)) |> 
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



## add percentage detP probes ----

# from: scripts/analysis_percentage_detP_probes.R

tmp <- read.table("output/tables/percentage_detP_probes.txt") |> 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id)

glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(percentage.detP.signi))




## heidelberg reportBrain files ----


w <- function(fn, prefix) {
  
  a <- read.csv(fn) |> 
    tibble::column_to_rownames('X') |> 
    dplyr::rename(pval = 1) |> 
    dplyr::arrange(-pval) |> 
    dplyr::mutate(pval = round(pval * 100,1)) 
  
  top <- a |> 
    tibble::rownames_to_column('class') |> 
    dplyr::slice_head(n=1) |> 
    dplyr::pull(class)
  
  a <- a |> 
    t() |> 
    as.data.frame() |> 
    dplyr::mutate(class = top) |> 
    dplyr::rename_with( ~ paste0(prefix,"cal_", .x)) 
  
  return(a)
}




tmp <- c(
  list.files(path = "data/GLASS_OD/Methylation data - EPIC arrays - brain classifier/", pattern = "*_scores_cal.csv", recursive = TRUE)
) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/GLASS_OD/Methylation data - EPIC arrays - brain classifier/", filename)) |>
  dplyr::mutate(basename = gsub("^.+/([^/]+)$", "\\1", filename)) |>
  dplyr::mutate(version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", filename)) |>
  dplyr::select(-basename) |>
  tidyr::pivot_wider(id_cols = sentrix_id, names_from = version, values_from = c(filename), names_prefix = "heidelberg_reportBrain_") |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = w(heidelberg_reportBrain_v12.5, "predictBrain_12.5_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id)




glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(heidelberg_reportBrain_v12.5)) |> 
  assertr::verify(!is.na(heidelberg_reportBrain_v2.0.1)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })()
  
rm(tmp, w)







## QC PCA ----



tmp <- readRDS("cache/unsupervised_qc_qc.outliers.Rds")
glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))

rm(tmp)





# add tumor purity calls? ----


tmp <- readRDS("cache/analysis_tumor_purity_EPIC_bin-based.Rds")
glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))

rm(tmp)


# cleanup db connection ----

DBI::dbDisconnect(metadata.db.con)
rm(metadata.db.con)



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


# tmp <- read.csv('data/GLASS_NL/Methylation/Metadata/Datasheet4.csv') # former file

tmp <- read.csv("data/GLASS_NL/Metadata/Samples/Master Datasheet_ALL METHODS_27012023.csv") |>
  dplyr::select( Surgery_ID, GLASS_ID,  Sample_Name, Sample_Sex, Sample_Type, Resectie, Sample_ID, Recurrent_Select.Meth, Matched_Pair.Meth) |> 
  dplyr::rename(sentrix_id = Sample_ID) |> 
  dplyr::filter(!is.na(sentrix_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (235)) # no idea how there can be 658 idats in here?
    return(.)
  })()


glass_nl.metadata.idats <- glass_nl.metadata.idats |>
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |>
  dplyr::filter(!is.na(Sample_Name)) |>  # exclude hundreds present in directory but absent in metadata
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (235)) # no idea how there can be 658 idats in here?
    return(.)
  })()


rm(tmp)




## Type IDH-mut & CN stata (CDK etc.) ----
#tmp <- read.csv("data/GLASS_NL/Metadata/(Epi)genetic_data/(Epi)genetic data_GLASS-NL_01092021.csv") |> 
#  dplyr::mutate(X=NULL)



## WHO classes ----


tmp <- read.csv("data/GLASS_NL/Metadata/Samples/WHOclassification_03052022.csv")  |> 
  dplyr::select(Surgery_ID, WHO_Classification2021) |> # join on Surgery_ID
  dplyr::mutate(WHO_Classification2021 = ifelse(WHO_Classification2021 == "Therapy Effects", as.character(NA), WHO_Classification2021)) |> 
  assertr::verify(!is.na(Surgery_ID)) |> 
  dplyr::filter(Surgery_ID %in% glass_nl.metadata.idats$Surgery_ID) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (235)) # no idea how there can be 658 idats in here?
    return(.)
  })()



glass_nl.metadata.idats <- glass_nl.metadata.idats |>
  dplyr::left_join(tmp, by=c('Surgery_ID'='Surgery_ID'), suffix=c('','')) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (235)) # no idea how there can be 658 idats in here?
    return(.)
  })()


rm(tmp)



## heidelberg 12.8 reportBrain files ----


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
  dplyr::mutate(tmp = parse_reportBrain_csv(mnp_predictBrain_filename, paste0("predictBrain_", mnp_predictBrain_version, "_"))) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 214)
    return(.)
  })() |> 
  assertr::verify(!is.na(predictBrain_v12.8_cal_class)) |>  # version is hardcoded here
  assertr::verify(!is.na(predictBrain_v12.8_cal_A_IDH_LG)) |>
  assertr::verify(!is.na(predictBrain_v12.8_cal_A_IDH_HG)) |>
  dplyr::mutate(A_IDH_HG__A_IDH_LG_lr = log(predictBrain_v12.8_cal_A_IDH_HG / predictBrain_v12.8_cal_A_IDH_LG)) 



glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  #assertr::verify(!is.na(predictBrain_v12.8_cal_class)) |> 
  #assertr::verify(!is.na(A_IDH_HG__A_IDH_LG_lr)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 235)
    return(.)
  })() |> 
  (function(.) {
    assertthat::assert_that(sum(!is.na(.$A_IDH_HG__A_IDH_LG_lr)) == 214)
    return(.)
  })() 

rm(tmp)



## Heidelberg Frozen ~ FFPE status ----


tmp <- list.files(path = "data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/",  pattern = "_ffpe_frozen.txt", recursive = TRUE) |> 
  data.frame(mnpQC_FrozenFFPEstatus_table = _) |> 
  dplyr::mutate(mnpQC_FrozenFFPEstatus_table = paste0("data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/", mnpQC_FrozenFFPEstatus_table)) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnpFrozenFFPEstatus_table(mnpQC_FrozenFFPEstatus_table, "mnpQC_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% glass_nl.metadata.idats$sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 214)
    return(.)
  })()



glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  #assertr::verify(!is.na(mnpQC_predicted_array_type)) |> 
  #assertr::verify(!is.na(mnpQC_predicted_sample_type)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 235)
    return(.)
  })() |> 
  (function(.) {
    assertthat::assert_that(sum(!is.na(.$mnpQC_predicted_array_type)) == 214)
    return(.)
  })() 

rm(tmp)


