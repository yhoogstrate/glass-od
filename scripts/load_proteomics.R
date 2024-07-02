#!/usr/bin/env R



tmp <- readxl::read_xlsx('data/GLASS_OD/Protein - Tobias Weiss/GLODprot_raw_protein_matrix.xlsx')


metadata.proteomics.proteins <- tmp |> 
  dplyr::select(Precursor.Id, Protein.Group, Protein.Ids,  Protein.Names, Genes,  First.Protein.Description, Proteotypic,  Stripped.Sequence, Modified.Sequence,  Precursor.Charge) |> 
  dplyr::mutate(peptide_id = Precursor.Id) |> 
  tibble::column_to_rownames("Precursor.Id")|> 
  tibble::tibble()


data.proteomics.glass_od <- tmp |> 
  dplyr::select(Precursor.Id, contains("/scratch/")) |> 
  tibble::column_to_rownames("Precursor.Id") |> 
  dplyr::rename_with(~ gsub("/scratch/DIANN_A314/WU300344/","",.x)) |> 
  dplyr::rename_with(~ gsub("\\.d$","",.x)) |> 
  dplyr::mutate(
    `20240215_003_S640489_GLODprot_pool` = NULL,
    `20240215_036_S640489_GLODprot_pool` = NULL,
    `20240215_069_S640489_GLODprot_pool` = NULL,
    `20240215_102__S640489_GLODprot_pool` = NULL,
    `20240215_135_S640489_GLODprot_pool` = NULL
  ) |> 
  dplyr::rename_with(~ gsub("^.+GLODprot_","GLODprot_",.x)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_GLASS_OD_ALL_PROTEOMICS)
    assertthat::assert_that(all(colnames(.) %in% glass_od.metadata.proteomics$proteomics_id))
    return(.)
  })() |> 
  dplyr::select(glass_od.metadata.proteomics$proteomics_id) # same order


rm(tmp)



