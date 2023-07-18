#!/usr/bin/env R

# G-SAM ----

## get idats ----

gsam.metadata <-  list.files(path = "data/G-SAM/", pattern = "_(Grn|Red).idat$", recursive = T) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/G-SAM/", filename)) |> 
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


tmp <- read.csv("data/G-SAM/MET2022-350-014/MET2022-350-014_IdH.csv", skip=8) |> 
  dplyr::filter(!is.na(Sentrix_ID)) |> 
  assertr::verify(grepl("^[0-9]{12}_[A-Z][0-9]{2}[A-Z][0-9]{2}$", Column2)) |> 
  dplyr::rename(sentrix_id = Column2) |> 
  dplyr::mutate(study = gsub("^(....).+$","\\1",Sample_Name)) |> 
  dplyr::filter(study %in% c("MINT","GLSO") == F) |> 
  assertr::verify(study == "GSAM") |> 
  dplyr::select(sentrix_id, Sample_Name) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == nrow(gsam.metadata)) # == 75
    return(.)
  })() |> 
  assertr::verify(sentrix_id %in% gsam.metadata$sentrix_id)



gsam.metadata <- gsam.metadata |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))


## heidelberg reportBrain ----


tmp <- list.files(path = "data/G-SAM/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "*.seg", recursive = TRUE) |> 
  data.frame(heidelberg_cnvp_segments = _) |> 
  dplyr::mutate(heidelberg_cnvp_segments = paste0("data/G-SAM/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", heidelberg_cnvp_segments)) |> 
  assertr::verify(file.exists(heidelberg_cnvp_segments)) |> 
  dplyr::mutate(heidelberg_cnvp_version = gsub("^.+/cnvp_([^/]+)/.+$","\\1", heidelberg_cnvp_segments)) |> 
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", heidelberg_cnvp_segments)) |> 
  assertr::verify(!is.na(sentrix_id)) |> 
  assertr::verify(!duplicated(sentrix_id)) |> # only one version per sample needed
  assertr::verify(sentrix_id %in% gsam.metadata$sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 75)
    return(.)
  })()





