#!/usr/bin/env R


ad_bmc_clin_epi.metadata.array_samples <- readxl::read_xlsx('data/GLASS_OD/DNA Methylation - 450K arrays - Alzheimer/microarray_samples_AD.xlsx') |> 
  dplyr::mutate(DNAm_id = gsub("[ \\-]","_", `Sample ID`)) |> 
  dplyr::mutate(patient_id = as.numeric(gsub("^([0-9]+).+$","\\1", `Sample ID`))) |> 
  dplyr::mutate(reason_excluded = ifelse(grepl("internal control and replicate", Type),"technical control", as.character(NA)))


