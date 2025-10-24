#!/usr/bin/env R 


if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}


# data ----


# keep original sids
data.proteomics.glass_nl <- read.csv('data/GLASS_NL/Proteomics/ProteinMatrix_30percentNA_cutoff_75percent_proteincutoff_MADnorm_MixedImputed_correct annotations_fixed-quotes_fixedspaces.csv',header=T) |> 
  dplyr::filter(X %in% c("HLA-B.1") == F) |>  # duplicated
  tibble::column_to_rownames('X')

  #dplyr::rename(!!! ( # rename to names as used in the raw data
  #  glass_nl.metadata.array_samples |> 
  #    dplyr::filter(!is.na(proteomics_sid)) |> 
  #    dplyr::pull(proteomics_sid, name=Sample_Name)
  #)) |> 
  #dplyr::rename(`153_R1` = GB_GIV_153_R1)


# metadata / annotations ----

metadata.proteomics.glass_nl <- data.frame(gene_id = rownames(data.proteomics.glass_nl))


# load DMP x CGC ----


# DPA: CGC ----

fn <- "cache/analysis_differential_proteomics__GLASS-NL__stats.cgc.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(protein_id, logFC, t, P.Value, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DPA__GLASS_NL__CGC__", .x), .cols=!matches("^protein_id$", perl = T))
  
  
  metadata.proteomics.glass_nl <- metadata.proteomics.glass_nl |> 
    dplyr::left_join(tmp, by=c('gene_id'='protein_id'), suffix=c('','') )
  
  
  rm(tmp)
  
} else {
  warning("DPA GLASS-NL x CGC is missing")
}

rm(fn)




