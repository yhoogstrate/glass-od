#!/usr/bin/env R



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


glass_od.metadata.patients <- read.table("data/GLASS_OD/Clinical_data/20230208_SentrixID_patient_link.txt", sep = "\t", header = T) |>
  dplyr::rename(GLASS_OD_patient = GLASS_OD_nr) |>
  dplyr::mutate(GLASS_OD_patient = paste0(GLASS_OD_patient, " [tmp-id]")) |>
  dplyr::mutate(Sentrix_ID = ifelse(Sentrix_ID == "", NA, Sentrix_ID)) |>
  dplyr::mutate(GLASS_OD_patient = ifelse(GLASS_OD_patient == "", NA, GLASS_OD_patient)) |> 
  dplyr::select(GLASS_OD_patient) |> 
  dplyr::distinct(GLASS_OD_patient) |> 
  dplyr::filter(grepl("NO_GLASS_OD", GLASS_OD_patient) == F)



# 2. resection level ----

# no overall operations / resection table available, yet?
# some resections have 1 or more arrays



## per probe level? ----




# 3. idat level ----


glass_od.metadata.idat <- list.files(path = "data/GLASS_OD/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/GLASS_OD/", filename)) |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == 404)
    return(.)
  })() |> # equals in-pipe stopifnot(nrow(.) == 404)
  assertr::verify(file.exists(filename)) |>
  dplyr::mutate(Sentrix_ID = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", filename)) |>
  dplyr::mutate(channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", filename)) |>
  tidyr::pivot_wider(id_cols = Sentrix_ID, names_from = channel, values_from = c(filename)) |>
  dplyr::rename(channel_green = Grn) |>
  dplyr::rename(channel_red = Red) |>
  (function(.) {
    assertthat::assert_that(nrow(.) == (404 / 2))
    return(.)
  })() |>
  assertr::verify(!is.na(channel_green)) |>
  assertr::verify(!is.na(channel_red)) |>
  assertr::verify(Sentrix_ID != "206119350032_R01C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(Sentrix_ID != "206119350032_R02C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(Sentrix_ID != "206119350032_R03C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(Sentrix_ID != "206119350032_R04C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(Sentrix_ID != "206119350032_R05C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(Sentrix_ID != "204808700073_R07C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(Sentrix_ID != "204808700074_R06C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(Sentrix_ID != "204808700074_R07C01") |> # sample not part in GLASS-OD - removal confirmed by Iris
  assertr::verify(Sentrix_ID != "204808700074_R08C01") # sample not part in GLASS-OD - removal confirmed by Iris



## link patient identifier ----


tmp <- read.table("data/GLASS_OD/Clinical_data/20230208_SentrixID_patient_link.txt", sep = "\t", header = T) |>
  dplyr::rename(GLASS_OD_patient = GLASS_OD_nr) |>
  dplyr::mutate(GLASS_OD_patient = paste0(GLASS_OD_patient, " [tmp-id]")) |>
  dplyr::mutate(Sentrix_ID = ifelse(Sentrix_ID == "", NA, Sentrix_ID)) |>
  dplyr::mutate(GLASS_OD_patient = ifelse(GLASS_OD_patient == "", NA, GLASS_OD_patient)) |> 
  dplyr::filter(grepl("NO_GLASS_OD", GLASS_OD_patient) == F) |> 
  dplyr::filter(!is.na(Sentrix_ID))


stopifnot(tmp$Sentrix_ID %in% glass_od.metadata.idat$Sentrix_ID)
stopifnot(glass_od.metadata.idat$Sentrix_ID %in% tmp$Sentrix_ID)
stopifnot(tmp$GLASS_OD_patient %in% glass_od.metadata.patients$GLASS_OD_patient)


glass_od.metadata.idat <- glass_od.metadata.idat |> 
  dplyr::left_join(tmp,
                   by=c('Sentrix_ID'='Sentrix_ID'), suffix=c('','')
  )



stopifnot(!is.na(glass_od.metadata.idat$GLASS_OD_patient))




## heidelberg reportBrain files ----


tmp <- c(
  list.files(path = "data/GLASS_OD/", pattern = "idat_reportBrain_*", recursive = TRUE),
  list.files(path = "data/GLASS_OD/", pattern = "report_website_mnp_brain_*", recursive = TRUE)
) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/GLASS_OD/", filename)) |>
  dplyr::filter(grepl("\\.pdf$", filename)) |>
  dplyr::mutate(basename = gsub("^.+/([^/]+)$", "\\1", filename)) |>
  dplyr::mutate(version = gsub("^.+_(v[0-9+]+[a-zA-Z0-9\\.]+)_.+$", "\\1", basename)) |>
  dplyr::filter(grepl("_sample_version_", basename) == F) |>
  dplyr::mutate(Sentrix_ID = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", filename)) |>
  dplyr::select(-basename) |>
  tidyr::pivot_wider(id_cols = Sentrix_ID, names_from = version, values_from = c(filename), names_prefix = "heidelberg_reportBrain_")


stopifnot(tmp$Sentrix_ID %in% glass_od.metadata.idat$Sentrix_ID)
stopifnot(duplicated(tmp$Sentrix_ID) == F)


stopifnot(nrow(glass_od.metadata.idat) == 202)


glass_od.metadata.idat <- glass_od.metadata.idat |> 
  dplyr::left_join(tmp, by=c('Sentrix_ID'='Sentrix_ID'), suffix=c('',''))
rm(tmp)


stopifnot(nrow(glass_od.metadata.idat) == 202)









## heidelberg CNV segment files ----


tmp <- list.files(path = "data/GLASS_OD/", pattern = "*.seg", recursive = TRUE) |> 
  data.frame(heidelberg_CNV_segments = _) |> 
  dplyr::mutate(heidelberg_CNV_segments = paste0("data/GLASS_OD/", heidelberg_CNV_segments)) |> 
  dplyr::mutate(heidelberg_cnvp_version = gsub("^.+/cnvp_([^/]+)/.+$","\\1", heidelberg_CNV_segments)) |> 
  dplyr::mutate(Sentrix_ID = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", heidelberg_CNV_segments)) |> 
  assertr::verify(!is.na(duplicated(Sentrix_ID)))


glass_od.metadata.idat <- glass_od.metadata.idat |> 
  dplyr::left_join(tmp, by=c('Sentrix_ID'='Sentrix_ID'), suffix=c('',''))
rm(tmp)


stopifnot(glass_od.metadata.idat |> 
  dplyr::filter(!is.na(heidelberg_reportBrain_v11b4) & is.na(heidelberg_CNV_segments)) |> 
  nrow() == 0)


## brain classifier 11 ----

x <- function(fn) {
  
  a = read.csv(fn) |> 
    tibble::column_to_rownames('X') |> 
    `colnames<-`('pval') |> 
    dplyr::arrange(-pval) |> 
    dplyr::mutate(pval = round(pval * 100,1)) 
  
  top <- a |> 
    tibble::rownames_to_column('class') |> 
    dplyr::slice(n=1) |> 
    dplyr::pull(class)
  
  a<- a |> 
    t() |> 
    as.data.frame() |> 
    dplyr::rename_with( ~ paste0("predictBrain_11_scores_cal_", .x)) |> 
    dplyr::mutate(predictBrain_11_scores_cal_class = top)
  
  return(a)
}


tmp <- list.files(path = "data/GLASS_OD/", pattern = "*_scores_cal.csv", recursive = TRUE) |> 
  data.frame(predictBrain_11 = _) |> 
  dplyr::mutate(predictBrain_11 = paste0("data/GLASS_OD/", predictBrain_11)) |> 
  dplyr::filter(grepl("predictBrain_v2\\.", predictBrain_11)) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = x(predictBrain_11)) |> 
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(Sentrix_ID = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", predictBrain_11)) 

stopifnot(duplicated(tmp$Sentrix_ID) == F)



stopifnot(nrow(glass_od.metadata.idat) == 202)


glass_od.metadata.idat <- glass_od.metadata.idat |> 
  dplyr::left_join(tmp, by=c('Sentrix_ID'='Sentrix_ID'), suffix=c('',''))
rm(tmp)


stopifnot(nrow(glass_od.metadata.idat) == 202)


rm(x)




## more metadata ----


md1 <- readxl::read_xlsx("data/GLASS_OD/Methylation_data/GLASS-NL_1p19qcodel/docs/Datasheet.xlsx")
#md1$Sample_ID

md2 <- readxl::read_xlsx("data/GLASS_OD/Methylation_data/GLASS-NL_1p19qcodel/docs/Datasheet EMC GLASS-OD+GLASS-NL.xlsx")
#md2$`Sentrix ID`

md3 <- readxl::read_xlsx("data/GLASS_OD/Methylation_data/MET2022-331-014/controls/MET2022-331-014_DATAControls.xlsx") |> 
  dplyr::mutate(Sentrix_ID = paste0(`Sentrix Barcode`, "_", `Sentrix Position`))
md3$Sentrix_ID

md4 <- read.csv("data/GLASS_OD/Methylation_data/MET2022-321-014/STS/MET2022-321-014.csv",skip=8) |> 
  dplyr::mutate(Sentrix_ID = paste0(Sentrix_ID, "_", Sentrix_Position))
md4$Sentrix_ID

md5 <- read.csv('data/GLASS_OD/Methylation_data/MET2020-240-014_plate2/STS/MET2020-240-014_plate2.csv',skip=7) |> 
  dplyr::mutate(Sentrix_ID = paste0(Sentrix_ID, "_", Sentrix_Position))


md6a <- read.csv('data/GLASS_OD/Methylation_data/MET2020-240-014/STS/MET2020-240-014.csv',skip=8) |> 
  dplyr::mutate(Sentrix_ID = paste0(Sentrix_ID, "_", Sentrix_Position))


md6d <- readxl::read_xlsx("data/GLASS_OD/Methylation_data/MET2020-240-014/docs/Datasheet MET2020-240-014.xlsx")
md6f <- readxl::read_xlsx("data/GLASS_OD/Methylation_data/MET2020-240-014/docs/MET2020-240-014.xlsx")






# 
# 
#   dplyr::mutate(in.md1 = Sentrix_ID %in% md1$Sample_ID) |> 
#   dplyr::mutate(in.md2 = Sentrix_ID %in% md2$`Sentrix ID`)|> 
#   dplyr::mutate(in.md3 = Sentrix_ID %in% md3$Sentrix_ID)|> 
#   dplyr::mutate(in.md4 = Sentrix_ID %in% md4$Sentrix_ID) |> 
#   dplyr::mutate(in.md5 = Sentrix_ID %in% md5$Sentrix_ID) |> 
#   dplyr::mutate(in.md6a = Sentrix_ID %in% md6a$Sentrix_ID)|> 
#   dplyr::mutate(in.md6d = Sentrix_ID %in% md6d$`Sentrix id`)|> 
#   dplyr::mutate(in.md6f = Sentrix_ID %in% md6f$`Sentrix id`)



##

# table(glass_od.metadata.idat$in.md1)
# table(glass_od.metadata.idat$in.md2)
# table(glass_od.metadata.idat$in.md3)
# table(glass_od.metadata.idat$in.md4)
# table(glass_od.metadata.idat$in.md5)
# table(glass_od.metadata.idat$in.md6a)
# table(glass_od.metadata.idat$in.md6d)
# table(glass_od.metadata.idat$in.md6f)


# 
# glass_od.metadata.idat |> 
#   dplyr::filter(!in.md1)  |> 
#   dplyr::filter(!in.md2)  |> 
#   dplyr::filter(!in.md3)  |> 
#   dplyr::filter(!in.md4)  |> 
#   #dplyr::filter(!in.md5)  |> 
#   dplyr::filter(!in.md6a)  |> 
#   dim()

stopifnot(nrow(glass_od.metadata.idat) == 202)

glass_od.metadata.idat <- glass_od.metadata.idat |> 
  dplyr::left_join(md1 ,  by=c('Sentrix_ID'='Sample_ID')) |> 
  dplyr::left_join(md2 , by=c('Sentrix_ID'='Sentrix ID'))|> 
  dplyr::left_join(md3 , by=c('Sentrix_ID'='Sentrix_ID'))|> 
  dplyr::left_join(md4 , by=c('Sentrix_ID'='Sentrix_ID')) |> 
  dplyr::left_join(md5 , by=c('Sentrix_ID'='Sentrix_ID')) |> 
  dplyr::left_join(md6a,  by=c('Sentrix_ID'='Sentrix_ID'))|> 
  dplyr::left_join(md6d, by=c('Sentrix_ID'='Sentrix id'))|> 
  dplyr::left_join(md6f, by=c('Sentrix_ID'='Sentrix id'))
rm(md1, md2, md3, md4, md5, md6a, md6d, md6f)

stopifnot(nrow(glass_od.metadata.idat) == 202)


# 
# glass_od.metadata.idat |> 
#   dplyr::filter(Sentrix_ID %in% c('204808700074_R07C01','204808700074_R06C01')) |> 
#   View()
# 
# 




