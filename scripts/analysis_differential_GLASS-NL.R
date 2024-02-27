#!/usr/bin/env R


# load data ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')
#source('scripts/load_palette.R')
#source('scripts/load_themes.R')



if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}


if(!exists('data.intensities.combined.hq_samples')) {
  source('scripts/load_intensities_hq_samples.R')
}


if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}



# GLASS-NL ----
## prim rec & quality ----


metadata.glass_nl <- glass_nl.metadata.array_samples |> 
  dplyr::mutate(patiend_id = paste0("p", patient_id)) |>  # avoid numbers because regression models may go nuts
  filter_GLASS_NL_idats(218) |> 
  filter_primaries_and_last_recurrences(189)






## g2, g3/g4 & quality ----


