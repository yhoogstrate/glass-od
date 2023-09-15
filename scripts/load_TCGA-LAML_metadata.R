#!/usr/bin/env


# load libs ----


source('scripts/load_functions.R')


# load idats ----


tcga_laml.metadata.array_samples <- list.files(path = "data/TCGA-LAML/DNA Methylation - EPIC arrays/GDCdownload/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(array_filename_GDC = _) |> 
  dplyr::mutate(array_filename_GDC = paste0("data/TCGA-LAML/DNA Methylation - EPIC arrays/GDCdownload/", array_filename_GDC)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (776))
    return(.)
  })() |>
  assertr::verify(file.exists(array_filename_GDC)) |>
  dplyr::mutate(array_sentrix_hash = gsub("^.+/([a-f0-9\\-]+)_noid.+$","\\1",array_filename_GDC)) |> 
  dplyr::mutate(array_sample_hash = gsub("^.+/([a-f0-9\\-]+)/.+$","\\1",array_filename_GDC)) |> 
  dplyr::arrange(array_sentrix_hash) |> 
  dplyr::group_by(array_sentrix_hash) |> 
  dplyr::add_count(array_sentrix_hash) |> 
  assertr::verify(n == 2) |> 
  dplyr::mutate(n = NULL) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(array_channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", array_filename_GDC)) |> 
  dplyr::mutate(array_channel = dplyr::recode(array_channel, `Grn` = 'green', `Red`='red')) |> 
  tidyr::pivot_wider(id_cols = c(array_sentrix_hash), names_from = array_channel, values_from = c(array_filename_GDC, array_sample_hash)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (776 / 2))
    return(.)
  })() |>
  dplyr::mutate(array_sentrix_hash = NULL) |> 
  dplyr::rename(array_channel_green_GDC = array_filename_GDC_green) |> 
  dplyr::rename(array_channel_red_GDC = array_filename_GDC_red) |> 
  assertr::verify(!is.na(array_channel_green_GDC)) |>
  assertr::verify(!is.na(array_channel_red_GDC)) |> 
  dplyr::mutate(array_sentrix_id = unlist(pbapply::pblapply(array_channel_green_GDC, parse_sentrix_id))) |> # recover the sentrix_id's
  assertr::verify(!is.na(array_sentrix_id)) |> 
  assertr::verify(!duplicated(array_sentrix_id))



# link metadata ----


tmp <- read.table("data/TCGA-LAML/Administration/query.txt") |> 
  dplyr::rename(array_sample_hash = id) |> 
  dplyr::mutate(array_type = as.factor(dplyr::case_when(
    file_size > 700000 & file_size < 750000 ~ "27k",
    file_size > 8000000 & file_size < 8500000 ~ "450k",
    T ~ as.character(NA)
  ))) |> 
  assertr::verify(!is.na(array_type)) |> 
  dplyr::select(array_sample_hash, cases, cases.submitter_id, sample.submitter_id, array_type)


tcga_laml.metadata.array_samples <- tcga_laml.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sample_hash_green'='array_sample_hash'), suffix=c('','')) |> 
  assertr::verify(!is.na(cases.submitter_id)) |> 
  assertr::verify(grepl("^TCGA-AB-[0-9]{4}$",cases.submitter_id)) |> 
  assertr::verify(!is.na(sample.submitter_id)) |> 
  assertr::verify(grepl("^TCGA-AB-[0-9]{4}-03A$", sample.submitter_id)) |> 
  assertr::verify(!is.na(array_type))



# restructure idats w/ symlinks ----


tcga_laml.metadata.array_samples <- tcga_laml.metadata.array_samples |>
  dplyr::mutate(array_channel_green = paste0("data/TCGA-LAML/DNA Methylation - EPIC arrays/restructured/",array_sentrix_id,"_Grn.idat")) |>
  dplyr::mutate(array_channel_red = paste0("data/TCGA-LAML/DNA Methylation - EPIC arrays/restructured/",array_sentrix_id,"_Red.idat"))


# 
# if(!file.exists('data/TCGA-LAML/DNA Methylation - EPIC arrays/restructured')) {
#   dir.create("data/TCGA-LAML/DNA Methylation - EPIC arrays/restructured")
# }
# 
# 
# tmp <- tcga_laml.metadata.array_samples |> 
#   dplyr::select(array_filename_Grn, array_filename_Red) |> 
#   tidyr::pivot_longer(c(array_filename_Grn, array_filename_Red))
# unlink(tmp$value)
#  
# file.symlink(gsub("data/TCGA-LAML/DNA Methylation - EPIC arrays/GDCdownload","../GDCdownload",tcga_laml.metadata.array_samples$array_channel_green_GDC), tcga_laml.metadata.array_samples$array_channel_green)
# file.symlink(gsub("data/TCGA-LAML/DNA Methylation - EPIC arrays/GDCdownload","../GDCdownload",tcga_laml.metadata.array_samples$array_channel_red_GDC), tcga_laml.metadata.array_samples$array_channel_red)


# add mutations ----

tmp <- read.table("data/TCGA-LAML/Administration/mutations.txt",header=T) |> 
  dplyr::mutate(STUDY_ID = NULL) |> 
  dplyr::mutate(cases.submitter_id = gsub("^(TCGA-AB-[0-9]{4}).*$","\\1",SAMPLE_ID)) |> 
  dplyr::mutate(SAMPLE_ID = NULL) |> 
  dplyr::mutate(IDH = ifelse(!is.na(IDH1) | !is.na(IDH2), T, F))


tcga_laml.metadata.array_samples <- tcga_laml.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('cases.submitter_id'='cases.submitter_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(IDH))


# add GDC / cBioportal Clinical xml ----



#fn = "data/TCGA-LAML/Clinical/GDCdownload/01cba3f5-2fdb-4290-b3b6-b11f54a358aa/genome.wustl.edu_clinical.TCGA-AB-2969.xml"
ff <- function(fn) {
  dat <- XML::xmlToDataFrame(fn)
  
  split_1 <- dat[,1:9] # file related?
  split_2 <- dat[,10:ncol(dat)] # patient related?
  
  stopifnot(is.na(split_1[2,]))
  stopifnot(is.na(split_2[1,]))
  
  split_1 <- split_1[1,]
  split_2 <- split_2[2,]
  
  merged <- cbind(split_1, split_2)
  
  return(merged)
}



tmp <- list.files(path = "data/TCGA-LAML/Clinical/GDCdownload/", pattern = "*.xml$", recursive = TRUE) |>
  data.frame(patient_filename_clinical = _) |> 
  dplyr::mutate(patient_filename_clinical = paste0("data/TCGA-LAML/Clinical/GDCdownload/", patient_filename_clinical)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (200))
    return(.)
  })() |> 
  dplyr::mutate(tmp = pbapply::pblapply(patient_filename_clinical, ff)) |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(tmp = NULL) |> 
  dplyr::select(bcr_patient_barcode, vital_status, days_to_death, days_to_last_followup) |> 
  dplyr::mutate(days_to_death = ifelse(days_to_death == "",NA, as.numeric(days_to_death))) |> 
  dplyr::mutate(days_to_last_followup = ifelse(days_to_last_followup %in% c("", 0),NA, as.numeric(days_to_last_followup))) |> 
  dplyr::mutate(day_to_last_event = ifelse(!is.na(days_to_death), days_to_death, days_to_last_followup)) |> 
  dplyr::mutate(last_event_status = dplyr::case_when(
    !is.na(days_to_death) & is.na(days_to_last_followup) ~ 1,
    is.na(days_to_death) & !is.na(days_to_last_followup) ~ 0,
    T ~ NA
  ))
rm(ff)


stopifnot(tcga_laml.metadata.array_samples$cases.submitter_id %in% tmp$bcr_patient_barcode)



tcga_laml.metadata.array_samples <- tcga_laml.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('cases.submitter_id'='bcr_patient_barcode'), suffix=c('',''))

rm(tmp)


# load data test ----

# 
# library(ggplot2)
# library(minfi)
# 
# library(IlluminaHumanMethylationEPICmanifest)
# library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) # BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
# 
# 
# targets <- tcga_laml.metadata.array_samples |> 
#   dplyr::filter(array_type == "450k") |> 
#   dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","",array_filename_Grn))
# 
# 
# 
# RGSet <- minfi::read.metharray.exp(targets = targets, force = T) #red/green channel together

