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






## heidelberg CNV files ----

### segment file ----


tmp <- list.files(path = "data/GLASS_OD/Methylation data - EPIC arrays - brain classifier/", pattern = "*.seg", recursive = TRUE) |> 
  data.frame(heidelberg_cnvp_segments = _) |> 
  dplyr::mutate(heidelberg_cnvp_segments = paste0("data/GLASS_OD/Methylation data - EPIC arrays - brain classifier/", heidelberg_cnvp_segments)) |> 
  dplyr::mutate(heidelberg_cnvp_version = gsub("^.+/cnvp_([^/]+)/.+$","\\1", heidelberg_cnvp_segments)) |> 
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", heidelberg_cnvp_segments)) |> 
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id)


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




### bins file ----


tmp <- list.files(path = "data/GLASS_OD/Methylation data - EPIC arrays - brain classifier/", pattern = "*.bins.igv", recursive = TRUE) |> 
  data.frame(heidelberg_cnvp_bins = _) |> 
  dplyr::mutate(heidelberg_cnvp_bins = paste0("data/GLASS_OD/Methylation data - EPIC arrays - brain classifier/", heidelberg_cnvp_bins)) |> 
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


### ongene scores ----


v <- function(fn, prefix) {
  a <- read.csv(fn,header=T, sep="\t", stringsAsFactors = F) |> 
    dplyr::select(name, value) |> 
    tibble::column_to_rownames('name') |> 
    dplyr::rename(cnv_score = value) |> 
    t() |> 
    as.data.frame(stringsAsFactors = F) |> 
    dplyr::rename_with( ~ paste0(prefix, .x))
  
  return(a)
}


tmp <- c(
  list.files(path = "data/GLASS_OD/Methylation data - EPIC arrays - brain classifier/", pattern = "*.detail.txt", recursive = TRUE)
) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/GLASS_OD/Methylation data - EPIC arrays - brain classifier/", filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", filename)) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = v(filename, "heidelberg_cnvp_")) |>
  tidyr::unnest(tmp) |> 
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id) |> 
  dplyr::mutate(filename = NULL)


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(`heidelberg_cnvp_CDKN2A/B`)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })()

rm(tmp, v)



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


tmp <- list.files(path = "data/GLASS_OD/Methylation data - EPIC arrays - brain classifier/", pattern = "_qc_full.txt", recursive = TRUE) |> 
  data.frame(heidelberg_qc_report_full = _) |> 
  dplyr::mutate(heidelberg_qc_report_full = paste0("data/GLASS_OD/Methylation data - EPIC arrays - brain classifier/", heidelberg_qc_report_full)) |>
  dplyr::mutate(basename = gsub("^.+/([^/]+)$", "\\1", heidelberg_qc_report_full)) |>
  dplyr::mutate(heidelberg_qc_report_version = gsub("^.+qc_(v[^_\\/]+)[_/].+$","\\1", heidelberg_qc_report_full)) |>
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", heidelberg_qc_report_full)) |>
  dplyr::select(-basename) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = x(heidelberg_qc_report_full, "qc_")) |>
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
  
rm(tmp, x)





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



tmp <- list.files(path = "data/GLASS_OD/Methylation data - EPIC arrays - brain classifier/", pattern = "_mgmt.csv", recursive = TRUE) |> 
  data.frame(heidelberg_mgmt_report = _) |> 
  dplyr::mutate(heidelberg_mgmt_report = paste0("data/GLASS_OD/Methylation data - EPIC arrays - brain classifier/", heidelberg_mgmt_report)) |>
  dplyr::mutate(basename = gsub("^.+/([^/]+)$", "\\1", heidelberg_mgmt_report)) |>
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", heidelberg_mgmt_report)) |>
  dplyr::select(-basename) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = y(heidelberg_mgmt_report, "mgmt_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id)


glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(heidelberg_mgmt_report)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 222)
    return(.)
  })() |> 
  dplyr::mutate(heidelberg_mgmt_report = NULL)
rm(tmp, y)





## rs_gender ----


z <- function(fn, prefix) {

  a <- read.csv(fn,header=T) |> 
    dplyr::mutate(idat = NULL, array = NULL) |>  # also returns status if it is inconfident 
    dplyr::mutate(predicted = ifelse(predicted == F, "F", predicted)) |> 
    dplyr::rename_with( ~ paste0(prefix, .x)) 
  
  return(a)
}


tmp <- list.files(path = "data/GLASS_OD/Methylation data - EPIC arrays - brain classifier/", pattern = "*.mix_gender.csv", recursive = TRUE) |> 
  data.frame(heidelberg_rs_gender_report = _) |> 
  dplyr::mutate(heidelberg_rs_gender_report = paste0("data/GLASS_OD/Methylation data - EPIC arrays - brain classifier/", heidelberg_rs_gender_report)) |>
  dplyr::mutate(basename = gsub("^.+/([^/]+)$", "\\1", heidelberg_rs_gender_report)) |>
  dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", heidelberg_rs_gender_report)) |>
  dplyr::select(-basename) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = z(heidelberg_rs_gender_report, "rs_gender_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  assertr::verify(!is.na(sentrix_id))|> 
  assertr::verify(!duplicated(sentrix_id)) |> 
  assertr::verify(sentrix_id %in% glass_od.metadata.idats$sentrix_id)


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





# add tumor purity calls? ----


tmp <- readRDS("cache/analysis_tumor_purity_EPIC_bin-based.Rds")
glass_od.metadata.idats <- glass_od.metadata.idats |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))

rm(tmp)


# cleanup db connection ----

DBI::dbDisconnect(metadata.db.con)
rm(metadata.db.con)



# GLASS-NL ----


glass_nl.metadata.resections <- read.csv("data/GLASS_NL/Metadata/Cleaned_clinical/metadata_2022/Samplesheet_GLASS_RNAseq__ALL.csv") |> 
  dplyr::mutate(institute = gsub("^.+_(.+)_.+$","\\1",GLASS_ID)) |> 
  dplyr::rename(genomescan.sid = GS_ID) |> 
  dplyr::mutate(rid = paste0(gsub("^(.+_)[^_]+$","\\1",GLASS_ID),Sample_Name)) |> 
  dplyr::rename(Exclude.by.Wies.on.complete.pair = Exclude) |> 
  dplyr::mutate(Sample_Type = dplyr::case_when(Sample_Type == "I" ~ "initial",
                                        Sample_Type == "R" ~ "recurrent",
                                        T ~ as.character(NA))) |> 
  dplyr::mutate(Sample_Type = factor(Sample_Type, levels=c('initial','recurrent','X'))) |> 
  dplyr::mutate(Customer_ID = NULL) # horribly confusing set of identifiers, some of which match Sample_Name but on different samples






tmp.1 <- read.csv('data/GLASS_NL/Methylation/Metadata/Datasheet4.csv') |> 
  dplyr::mutate(X=NULL) |> 
  dplyr::mutate(GLASS_ID = NULL,
                Surgery_ID = NULL,
                Sample_Sex = NULL,
                Recurrent_Type = NULL,
                Sample_Type = NULL,
                Sample_Resection = NULL,
                
                Sample_Plate = NULL,
                Basename = NULL,
                Array = NULL,
                Slide = NULL
                )


tmp.2 <- data.frame(Heidelberg.segment.file = Sys.glob("data/GLASS_NL/Methylation/Heidelberg/Heidelberg_unzip/*/cnvp_v3.0/*.segments.seg")) |> 
  dplyr::mutate(Sample_ID = gsub("^.+_unzip/([^/]+)_Run.+$","\\1",Heidelberg.segment.file))


tmp.3 <- read.csv("data/GLASS_NL/Metadata/(Epi)genetic_data/(Epi)genetic data_GLASS-NL_01092021.csv") |> 
  dplyr::mutate(X=NULL)


parse_predictbrain_csv <- function(file, suffix) {
  # file = 'data/GLASS_NL/Methylation/Heidelberg/Heidelberg_unzip/203989100107_R01C01_Run_78486/predictBrain_v2.1/203989100107_R01C01_scores.csv'
  tmp <- read.csv(file)
  colnames(tmp)[2] <- 'score'
  
  out <- list(
    A_IDH =    tmp |>  dplyr::filter(grepl('^A_IDH$',X))    |>  dplyr::pull(score) |>  as.numeric(),
    A_IDH_HG = tmp |>  dplyr::filter(grepl('^A_IDH_HG$',X)) |> dplyr::pull(score) |>  as.numeric(),
    O_IDH =    tmp |>  dplyr::filter(grepl('^O_IDH$',X))    |>  dplyr::pull(score) |>  as.numeric(),
    GBM_MES =  tmp  |>  dplyr::filter(grepl('^GBM_MES$',X))  |>  dplyr::pull(score) |>  as.numeric()
  )
  
  names(out) <- paste0(names(out), suffix)
  
  return(as.data.frame(out))
}


tmp.4 <- data.frame(predictBrain.scores.file = Sys.glob("data/GLASS_NL/Methylation/Heidelberg/Heidelberg_unzip/*/predictBrain_v2.1/*_scores.csv")) |> 
  dplyr::mutate(Sample_ID = gsub("^.+_v2.1/([^/]+)_scores.csv$","\\1", predictBrain.scores.file)) |> 
  dplyr::mutate(stats = pbapply::pblapply(predictBrain.scores.file, parse_predictbrain_csv, suffix='')) |> 
  tidyr::unnest(stats) |>  
  dplyr::mutate(predictBrain.scores.file = NULL) |> 
  tibble::column_to_rownames('Sample_ID') |> 
  dplyr::mutate_all(as.numeric) |> 
  tibble::rownames_to_column('Sample_ID')


tmp.5 <- data.frame(predictBrain.scores.file = Sys.glob("data/GLASS_NL/Methylation/Heidelberg/Heidelberg_unzip/*/predictBrain_v2.1/*_scores_cal.csv")) |> 
  dplyr::mutate(Sample_ID = gsub("^.+_v2.1/([^/]+)_scores_cal.csv$","\\1", predictBrain.scores.file)) |> 
  dplyr::mutate(stats = pbapply::pblapply(predictBrain.scores.file, parse_predictbrain_csv, suffix='_cal')) |> 
  tidyr::unnest(stats) |> 
  dplyr::mutate(predictBrain.scores.file = NULL) |> 
  tibble::column_to_rownames('Sample_ID') |> 
  dplyr::mutate_all(as.numeric) |> 
  tibble::rownames_to_column('Sample_ID') |> 
  dplyr::mutate(IDH_HG_IDH_ratio = log(A_IDH_HG_cal/A_IDH_cal))



stopifnot(duplicated(tmp.1$Sample_ID) == FALSE)
stopifnot(duplicated(tmp.2$Sample_ID) == FALSE)
stopifnot(duplicated(tmp.3$Sample_ID) == FALSE)
stopifnot(tmp.1$Sample_ID %in% tmp.3$Sample_ID)
stopifnot(tmp.3$Sample_ID %in% tmp.1$Sample_ID)


tmp <- tmp.1 |> 
  dplyr::left_join(tmp.2, by=c('Sample_ID'='Sample_ID')) |> 
  dplyr::left_join(tmp.3, by=c('Sample_ID'='Sample_ID')) |>  
  dplyr::left_join(tmp.4, by=c('Sample_ID'='Sample_ID')) |>  
  dplyr::left_join(tmp.5, by=c('Sample_ID'='Sample_ID'))  |> 
  dplyr::rename(sentrix_id = Sample_ID) |> 
  dplyr::rename(methylation.sub.diagnosis = sub.diagnosis)


stopifnot(sum(is.na(tmp$Heidelberg.segment.file)) <= 1) # one file missing so far, that's known


rm(tmp.1,tmp.2,tmp.3,tmp.4,tmp.5, parse_predictbrain_csv)



glass_nl.metadata.resections <- glass_nl.metadata.resections |> 
  dplyr::left_join(tmp, by=c('Sample_Name'='Sample_Name')) 

rm(tmp)




## idats ----


tmp <-  list.files(path = "data/GLASS_NL/Methylation/Methylation Array Data/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/GLASS_NL/Methylation/Methylation Array Data/", filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (658))
    return(.)
  })() |>
  assertr::verify(file.exists(filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", filename)) |>
  dplyr::mutate(channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", filename)) |>
  tidyr::pivot_wider(id_cols = sentrix_id, names_from = channel, values_from = c(filename)) |>
  dplyr::rename(channel_green = Grn) |>
  dplyr::rename(channel_red = Red) |> 
  assertr::verify(is.na(sentrix_id) == is.na(channel_green)) |>
  assertr::verify(is.na(sentrix_id) == is.na(channel_red))


glass_nl.metadata.resections <- glass_nl.metadata.resections |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (217))
    return(.)
  })() |>
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  dplyr::filter(!is.na(sentrix_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (214))
    return(.)
  })()




# G-SAM ----

## idats ----


gsam.metadata <- read.csv("data/G-SAM/MET2022-350-014/MET2022-350-014_IdH.csv", skip=8) |> 
  dplyr::filter(!is.na(Sentrix_ID)) |> 
  assertr::verify(grepl("^[0-9]{12}_[A-Z][0-9]{2}[A-Z][0-9]{2}$", Column2)) |> 
  dplyr::rename(sentrix_id = Column2) |> 
  dplyr::mutate(study = gsub("^(....).+$","\\1",Sample_Name)) |> 
  dplyr::filter(study %in% c("MINT","GLSO") == F) |> 
  assertr::verify(study == "GSAM") |> 
  dplyr::select(sentrix_id, Sample_Name)


tmp <-  list.files(path = "data/G-SAM/", pattern = "_(Grn|Red).idat$", recursive = T) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/G-SAM/", filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (240))
    return(.)
  })() |> 
  assertr::verify(file.exists(filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", filename)) |>
  dplyr::mutate(channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", filename)) |>
  tidyr::pivot_wider(id_cols = sentrix_id, names_from = channel, values_from = c(filename)) |>
  dplyr::rename(channel_green = Grn) |>
  dplyr::rename(channel_red = Red) |> 
  dplyr::filter(!grepl("/MET2017-126-014/", channel_green)) |> # stored there for historical reasons - IDH-mutant loss study
  dplyr::filter(!grepl("/GLSO/", channel_green)) |> # stored there for historical reasons
  dplyr::filter(!grepl("/MINT/", channel_green)) |> # stored there for historical reasons
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 75)
    return(.)
  })()


stopifnot(nrow(gsam.metadata) == nrow(tmp))
stopifnot(gsam.metadata$sentrix_id %in% tmp$sentrix_id)


gsam.metadata <- gsam.metadata |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))


# heidelberg reportBrain ----




