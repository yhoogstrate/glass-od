#!/usr/bin/env R

# load data ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')
#source('scripts/load_palette.R')


if(!exists('glass_od.metadata.proteomics')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# load proteomics raw per-peptide ----


#tmp <- readxl::read_xlsx('data/GLASS_OD/Protein - Tobias Weiss/GLODprot_raw_protein_matrix.xlsx')


# metadata.proteomics.proteins <- tmp |> 
#     dplyr::select(Precursor.Id, Protein.Group, Protein.Ids,  Protein.Names, Genes,  First.Protein.Description, Proteotypic,  Stripped.Sequence, Modified.Sequence,  Precursor.Charge) |> 
#     dplyr::mutate(peptide_id = Precursor.Id) |> 
#     tibble::column_to_rownames("Precursor.Id")|> 
#     tibble::tibble()


# data.proteomics.glass_od <- tmp |> 
#   dplyr::select(Precursor.Id, contains("/scratch/")) |> 
#   tibble::column_to_rownames("Precursor.Id") |> 
#   dplyr::rename_with(~ gsub("/scratch/DIANN_A314/WU300344/","",.x)) |> 
#   dplyr::rename_with(~ gsub("\\.d$","",.x)) |> 
#   dplyr::mutate(
#     `20240215_003_S640489_GLODprot_pool` = NULL,
#     `20240215_036_S640489_GLODprot_pool` = NULL,
#     `20240215_069_S640489_GLODprot_pool` = NULL,
#     `20240215_102__S640489_GLODprot_pool` = NULL,
#     `20240215_135_S640489_GLODprot_pool` = NULL
#   ) |> 
#   dplyr::rename_with(~ gsub("^.+GLODprot_","GLODprot_",.x)) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(ncol(.) == CONST_N_GLASS_OD_ALL_PROTEOMICS)
#     assertthat::assert_that(all(colnames(.) %in% glass_od.metadata.proteomics$proteomics_id))
#     return(.)
#   })() |> 
#   dplyr::select(glass_od.metadata.proteomics$proteomics_id) |># same order
#   log2()
# 
# 
# rm(tmp)


# load proteomics per-protein ----


tmp <- read.csv('data/GLASS_OD/Protein - Tobias Weiss/WU300344_report.pg_matrix.tsv',
                check.names=F,
                sep="\t") |> 
  dplyr::filter(!grepl("Y-FGCZCont", Protein.Ids)) |> 
  dplyr::filter(!is.na(Genes)) |> 
  dplyr::filter(Genes != "")


metadata.proteomics.proteins <- tmp |> 
  dplyr::select(Protein.Group,
                Protein.Ids,   
                Protein.Names,  
                Genes, 
                First.Protein.Description ) |> 
  dplyr::mutate(protein_id = Genes) |> 
  tibble::column_to_rownames("Genes")|> 
  tibble::tibble()



data.proteomics.glass_od <- tmp |> 
  dplyr::select(Genes, contains("scratch")) |> 
  tibble::column_to_rownames("Genes") |> 
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
  dplyr::select(glass_od.metadata.proteomics$proteomics_id) |>  # same order
  log2()


rm(tmp)


