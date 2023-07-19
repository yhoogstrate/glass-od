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
  dplyr::mutate(tmp = parse_reportBrain_csv(mnp_predictBrain_filename, paste0("mnp_predictBrain_", mnp_predictBrain_version, "_"))) |>
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



## Heidelberg 12.8 Frozen ~ FFPE status ----


tmp <- list.files(path = "data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/",  pattern = "_ffpe_frozen.txt", recursive = TRUE) |> 
  data.frame(mnpQC_FrozenFFPEstatus_table = _) |> 
  dplyr::mutate(mnpQC_FrozenFFPEstatus_table = paste0("data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/", mnpQC_FrozenFFPEstatus_table)) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_FrozenFFPEstatus_table(mnpQC_FrozenFFPEstatus_table, "mnpQC_")) |>
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



## Heidelberg 12.8 CNVP segment files ----

tmp <- query_Heidelberg_12_8_CNVP_segment_files(
  "data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/",
  214, # 235 to become
  glass_nl.metadata.idats$sentrix_id)



glass_nl.metadata.idats <- glass_nl.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))
  #assertr::verify(!is.na(heidelberg_cnvp_segments)) 
  #assertr::verify(!is.na(heidelberg_cnvp_version)) 
rm(tmp)



