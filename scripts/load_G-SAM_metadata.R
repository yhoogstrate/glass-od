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


## link sample names


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
  assertr::verify(sentrix_id %in% gsam.metadata.idats$sentrix_id)



gsam.metadata.idats <- gsam.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(Sample_Name))
rm(tmp)



## heidelberg 12.8 reportBrain files ----


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
  dplyr::mutate(tmp = parse_reportBrain_csv(mnp_predictBrain_filename, paste0("predictBrain_", mnp_predictBrain_version, "_"))) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 75)
    return(.)
  })() |> 
  assertr::verify(!is.na(predictBrain_v12.8_cal_class)) |>  # version is hardcoded here
  assertr::verify(!is.na(predictBrain_v12.8_cal_A_IDH_LG)) |>
  assertr::verify(!is.na(predictBrain_v12.8_cal_A_IDH_HG)) |>
  dplyr::mutate(A_IDH_HG__A_IDH_LG_lr = log(predictBrain_v12.8_cal_A_IDH_HG / predictBrain_v12.8_cal_A_IDH_LG)) 



gsam.metadata.idats <- gsam.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(predictBrain_v12.8_cal_class)) |> 
  assertr::verify(!is.na(A_IDH_HG__A_IDH_LG_lr)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 75)
    return(.)
  })() 

rm(tmp)



## CNVP segment files ----


tmp <- list.files(path = "data/G-SAM/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/",
                  pattern = "*.seg", recursive = TRUE) |> 
  data.frame(heidelberg_cnvp_segments = _) |> 
  dplyr::mutate(heidelberg_cnvp_segments = paste0("data/G-SAM/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", heidelberg_cnvp_segments)) |> 
  assertr::verify(file.exists(heidelberg_cnvp_segments)) |> 
  dplyr::mutate(heidelberg_cnvp_version = gsub("^.+/cnvp_([^/]+)/.+$","\\1", heidelberg_cnvp_segments)) |> 
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", heidelberg_cnvp_segments)) |> 
  assertr::verify(!is.na(sentrix_id)) |> 
  assertr::verify(!duplicated(sentrix_id)) |> # only one version per sample needed
  assertr::verify(sentrix_id %in% gsam.metadata.idats$sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 75)
    return(.)
  })()

gsam.metadata.idats <- gsam.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))
rm(tmp)





