#!/usr/bin/env R

#' https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-019-0672-7

ad_bmc_clin_epi.metadata.array_samples <- readxl::read_xlsx('data/GLASS_OD/DNA Methylation - 450K arrays - Alzheimer/microarray_samples_AD.xlsx') |> 
  dplyr::mutate(DNAm_id = gsub("[ \\-]","_", `Sample ID`)) |> 
  dplyr::mutate(patient_id = as.numeric(gsub("^([0-9]+).+$","\\1", `Sample ID`))) |> 
  dplyr::mutate(reason_excluded = ifelse(grepl("internal control", Type),"internal control and replicate", as.character(NA))) |> 
  (function(.) {
    
    print(dim(.))
    assertthat::assert_that(nrow(.) == (26 + 12))
    
    print(dim(.[is.na(.$reason_excluded),]))
    assertthat::assert_that(nrow(.[is.na(.$reason_excluded),]) == (26 + 12 - 1))
    
    return(.)
  })()


