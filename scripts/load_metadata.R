#!/usr/bin/env R


# 1. idat level ----


glass_od.metadata.idat <- list.files(path = "data/GLASS_OD/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(filename = _) |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == 422)
    return(.)
  })() |> # equals in-pipe stopifnot(nrow(.) == 422)
  dplyr::mutate(idat_prefix = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", filename)) |>
  dplyr::mutate(channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", filename)) |>
  tidyr::pivot_wider(id_cols = idat_prefix, names_from = channel, values_from = c(filename)) |>
  dplyr::rename(channel_green = Grn) |>
  dplyr::rename(channel_red = Red) |>
  (function(.) {
    assertthat::assert_that(nrow(.) == (422 / 2))
    return(.)
  })() |>
  assertr::verify(!is.na(channel_green)) |>
  assertr::verify(!is.na(channel_red))




## heidelberg reportBrain files ----


tmp <- c(
  list.files(path = "data/GLASS_OD/", pattern = "idat_reportBrain_*", recursive = TRUE),
  list.files(path = "data/GLASS_OD/", pattern = "report_website_mnp_brain_*", recursive = TRUE)
) |> 
  data.frame(filename = _) |> 
  dplyr::mutate(filename = paste0("data/GLASS_OD/", filename)) |> 
  dplyr::filter(grepl("\\.pdf$", filename)) |> 
  dplyr::mutate(basename = gsub("^.+/([^/]+)$","\\1", filename)) |> 
  dplyr::mutate(version = gsub("^.+_(v[0-9+]+[a-zA-Z0-9\\.]+)_.+$","\\1", basename)) |> 
  dplyr::filter(grepl("_sample_version_", basename) == F) |> 
  dplyr::mutate(idat_prefix = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", filename)) |> 
  dplyr::select(-basename) |> 
  tidyr::pivot_wider(id_cols = idat_prefix, names_from = version, values_from = c(filename), names_prefix="heidelberg_reportBrain_")


stopifnot(tmp$idat_prefix %in% glass_od.metadata.idat$idat_prefix)
stopifnot(duplicated(tmp$idat_prefix)== F)


stopifnot(nrow(glass_od.metadata.idat) == 211)


glass_od.metadata.idat <- glass_od.metadata.idat |> 
  dplyr::left_join(tmp, by=c('idat_prefix'='idat_prefix'), suffix=c('',''))
rm(tmp)


stopifnot(nrow(glass_od.metadata.idat) == 211)


## heidelberg CNV segment files ----


tmp <- list.files(path = "data/GLASS_OD/", pattern = "*.seg", recursive = TRUE) |> 
  data.frame(heidelberg_CNV_segments = _) |> 
  dplyr::mutate(heidelberg_CNV_segments = paste0("data/GLASS_OD/", heidelberg_CNV_segments)) |> 
  dplyr::mutate(heidelberg_cnvp_version = gsub("^.+/cnvp_([^/]+)/.+$","\\1", heidelberg_CNV_segments)) |> 
  dplyr::mutate(idat_prefix = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", heidelberg_CNV_segments)) |> 
  assertr::verify(!is.na(duplicated(idat_prefix)))


glass_od.metadata.idat <- glass_od.metadata.idat |> 
  dplyr::left_join(tmp, by=c('idat_prefix'='idat_prefix'), suffix=c('',''))
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
  dplyr::mutate(idat_prefix = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", predictBrain_11)) 

stopifnot(duplicated(tmp$idat_prefix) == F)



stopifnot(nrow(glass_od.metadata.idat) == 211)


glass_od.metadata.idat <- glass_od.metadata.idat |> 
  dplyr::left_join(tmp, by=c('idat_prefix'='idat_prefix'), suffix=c('',''))

stopifnot(nrow(glass_od.metadata.idat) == 211)

#rm(tmp)





## more metadata ----


md1 <- readxl::read_xlsx("data/GLASS_OD/Methylation_data/Patients in GLASS-NL with 1p19q codel/Datasheet.xlsx")
md1$Sample_ID

md2 <- readxl::read_xlsx("data/GLASS_OD/Methylation_data/Patients in GLASS-NL with 1p19q codel/clinical data/Datasheet EMC GLASS-OD+GLASS-NL.xlsx")
md2$`Sentrix ID`

md3 <- readxl::read_xlsx("data/GLASS_OD/Methylation_data/MET2022-331-014_QPioNsjg8&cjud38d5>kj6&/controls/MET2022-331-014_DATAControls.xlsx") |> 
  dplyr::mutate(idat_prefix = paste0(`Sentrix Barcode`, "_", `Sentrix Position`))
md3$idat_prefix

md4 <- read.csv("data/GLASS_OD/Methylation_data/MET2022-321-014_78HBjiKHy>,86Cv#blPof1>*/STS/MET2022-321-014.csv",skip=8) |> 
  dplyr::mutate(idat_prefix = paste0(Sentrix_ID, "_", Sentrix_Position))
md4$idat_prefix

md5 <- read.csv('data/GLASS_OD/Methylation_data/MET2020-240-014_plate2_eSKLpmSd821/MET2020-240-014_plate2/STS/MET2020-240-014_plate2.csv',skip=7) |> 
  dplyr::mutate(idat_prefix = paste0(Sentrix_ID, "_", Sentrix_Position))


md6a <- read.csv('data/GLASS_OD/Methylation_data/MET2020-240-014_B@NMtyVaSQ843!@/MET2020-240-014/STS/MET2020-240-014.csv',skip=8) |> 
  dplyr::mutate(idat_prefix = paste0(Sentrix_ID, "_", Sentrix_Position))


#md6b <- readxl::read_xlsx("data/GLASS_OD/Methylation_data/MET2020-240-014_B@NMtyVaSQ843!@/20210120_GLASS Oligo_MartaPadovan.xlsx")
#md6c <- readxl::read_xlsx("data/GLASS_OD/Methylation_data/MET2020-240-014_B@NMtyVaSQ843!@/20210224_Methylation_GLASS-oligo_MartaPadovan_Redo.xlsx")
md6d <- readxl::read_xlsx("data/GLASS_OD/Methylation_data/MET2020-240-014_B@NMtyVaSQ843!@/Datasheet MET2020-240-014.xlsx")
#md6e <- readxl::read_xlsx("data/GLASS_OD/Methylation_data/MET2020-240-014_B@NMtyVaSQ843!@/GLASS Oligo blocks.xlsx")
md6f <- readxl::read_xlsx("data/GLASS_OD/Methylation_data/MET2020-240-014_B@NMtyVaSQ843!@/MET2020-240-014.xlsx")





  dplyr::mutate(in.md1 = idat_prefix %in% md1$Sample_ID) |> 
  dplyr::mutate(in.md2 = idat_prefix %in% md2$`Sentrix ID`)|> 
  dplyr::mutate(in.md3 = idat_prefix %in% md3$idat_prefix)|> 
  dplyr::mutate(in.md4 = idat_prefix %in% md4$idat_prefix) |> 
  dplyr::mutate(in.md5 = idat_prefix %in% md5$idat_prefix) |> 
  dplyr::mutate(in.md6a = idat_prefix %in% md6a$idat_prefix)|> 
  dplyr::mutate(in.md6d = idat_prefix %in% md6d$`Sentrix id`)|> 
  dplyr::mutate(in.md6f = idat_prefix %in% md6f$`Sentrix id`)



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

stopifnot(nrow(glass_od.metadata.idat) == 209)

glass_od.metadata.idat <- glass_od.metadata.idat |> 
  dplyr::left_join(md1 ,  by=c('idat_prefix'='Sample_ID')) |> 
  dplyr::left_join(md2 , by=c('idat_prefix'='Sentrix ID'))|> 
  dplyr::left_join(md3 , by=c('idat_prefix'='idat_prefix'))|> 
  dplyr::left_join(md4 , by=c('idat_prefix'='idat_prefix')) |> 
  dplyr::left_join(md5 , by=c('idat_prefix'='idat_prefix')) |> 
  dplyr::left_join(md6a,  by=c('idat_prefix'='idat_prefix'))|> 
  dplyr::left_join(md6d, by=c('idat_prefix'='Sentrix id'))|> 
  dplyr::left_join(md6f, by=c('idat_prefix'='Sentrix id'))
rm(md1, md2, md3, md4, md5, md6a, md6d, md6f)

stopifnot(nrow(glass_od.metadata.idat) == 209)


# 1. resection level ----

# no overall operations / resection table available, yet?
# some resections have 1 or more arrays



## per probe level? ----


# 3. patient level ----




# 
# readRDS('data/GLASS_OD/MET2022-321-014/meta_heidel_glassod.Rds')
#sample           batch chip.type     type                idat
#1  GLASS-oligo_86_R1 MET2022-321-014      EPIC DNA-FFPE 206467010069_R01C01
#2  GLASS-oligo_84_R2 MET2022-321-014      EPIC DNA-FFPE 206467010069_R02C01
#3  GLASS-oligo_84_R1 MET2022-321-014      EPIC DNA-FFPE 206467010069_R03C01
#4  GLASS-oligo_83_R2 MET2022-321-014      EPIC DNA-FFPE 206467010069_R04C01
#5  GLASS-oligo_83_R1 MET2022-321-014      EPIC DNA-FFPE 206467010069_R05C01
#6  GLASS-oligo_82_R2 MET2022-321-014      EPIC DNA-FFPE 206467010069_R06C01
#7  GLASS-oligo_82_R1 MET2022-321-014      EPIC DNA-FFPE 206467010069_R07C01
#8  GLASS-oligo_81_R2 MET2022-321-014      EPIC DNA-FFPE 206467010069_R08C01
#9  GLASS-oligo_77_R1 MET2022-321-014      EPIC DNA-FFPE 206467010089_R01C01

# glass_od.metadata.resections <- readxl::read_xlsx("data/GLASS_OD/Patients\ in\ GLASS-NL\ with\ 1p19q\ codel/Datasheet.xlsx")


