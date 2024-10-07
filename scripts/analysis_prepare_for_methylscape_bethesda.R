#!/usr/bin/env R

# load ----


source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')


library(ggplot2)
library(patchwork)


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# GLASS-OD / OD ----


tab <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::select(resection_id, array_sentrix_id) |> 
  dplyr::rename(Sample_Name = resection_id) |> 
  dplyr::rename(Sample_IDAT = array_sentrix_id) |>
  dplyr::mutate(Sample_Plate = "") |> 
  dplyr::mutate(Sample_Group = "") |> 
  dplyr::mutate(Pool_ID = "") |> 
  dplyr::mutate(Sentrix_ID = gsub("^(.+)_(.+)$","\\1",Sample_IDAT)) |> 
  dplyr::mutate(Sentrix_Position = gsub("^(.+)_(.+)$","\\2",Sample_IDAT)) |> 
  dplyr::mutate(Material_Type = "") |> 
  dplyr::mutate(Gender = "") |> 
  dplyr::mutate(Surgical_Case = "") |> 
  dplyr::mutate(Diagnosis = "tumor") |> 
  dplyr::mutate(Age = "0") |> 
  dplyr::mutate(Notes = "") |> 
  dplyr::mutate(Tumor_site = "CNS") |> 
  dplyr::mutate(PI_Collaborator = "") |> 
  dplyr::mutate(Outside_ID = "") |> 
  dplyr::mutate(Surgery_date = "")

write.table(tab, file="/tmp/methylscape_data.txt", sep="\t", row.names=F, quote=F)



# clone files to "/tmp/idats/"
dir.create("/tmp/idats")


lst <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::select(array_channel_green, array_channel_red) |> 
  tidyr::pivot_longer(cols = c(array_channel_green, array_channel_red)) |> 
  dplyr::mutate(from = paste0(getwd(),"/",value)) |> 
  dplyr::select(from) |> 
  dplyr::mutate(to = paste0("/tmp/idats/",basename(from)))



file.symlink(lst$from, lst$to)

