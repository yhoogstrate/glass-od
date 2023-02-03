#!/usr/bin/env R

# 1. idat level ----

metadata.glass_od.idat <- list.files(path="data/GLASS_OD/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |> 
  dplyr::tibble(.name_repair = ~c("filename")) |> 
  (function(.) {assertthat::assert_that(nrow(.) == 418); return(.)})() |>  # in-pipe stopifnot(nrow(metadata.glass_od.idat) == 418)
  dplyr::mutate(prefix = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$","\\1", filename)) |> 
  dplyr::mutate(channel = gsub("^.+_(Grn|Red)\\.idat$","\\1", filename)) |> 
  tidyr::pivot_wider(id_cols = prefix, names_from = channel, values_from = c(filename)) |> 
  dplyr::rename(channel_green = Grn) |> 
  dplyr::rename(channel_red = Red) |> 
  (function(.) {assertthat::assert_that(nrow(.) == (418/2)); return(.)})() |> 
  assertr::verify(!is.na(channel_green)) |> 
  assertr::verify(!is.na(channel_red))




## per probe level? ----

# 2. resection level ----


# 3. patient level ----


b = readRDS('data/GLASS_OD/MET2022-321-014/meta_heidel_glassod.Rds')

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

glass_od.metadata.resections <- readxl::read_xlsx("data/GLASS_OD/Patients\ in\ GLASS-NL\ with\ 1p19q\ codel/Datasheet.xlsx")


