#!/usr/bin/env R 


if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}



# keep original sids
glass_nl.proteomics <- read.csv('data/GLASS_NL/Proteomics/ProteinMatrix_30percentNA_cutoff_75percent_proteincutoff_MADnorm_MixedImputed_correct annotations_fixed-quotes_fixedspaces.csv',header=T) |> 
  dplyr::filter(X %in% c("HLA-B.1") == F) |>  # duplicated
  tibble::column_to_rownames('X')

  #dplyr::rename(!!! ( # rename to names as used in the raw data
  #  glass_nl.metadata.array_samples |> 
  #    dplyr::filter(!is.na(proteomics_sid)) |> 
  #    dplyr::pull(proteomics_sid, name=Sample_Name)
  #)) |> 
  #dplyr::rename(`153_R1` = GB_GIV_153_R1)


