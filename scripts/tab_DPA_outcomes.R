
# load ----



source('scripts/load_constants.R')
source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')
source('scripts/load_gene_annotations.R')



if(!exists('data.proteomics.glass_od')) {
  source('scripts/load_proteomics.R')
}


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}



# select and export ----

## glass-od outcomes ----


out <- metadata.proteomics.glass_od |> 
  dplyr::select(
    "Protein.Group",
    "Protein.Ids",
    "Genes",
    
    "DPA__GLASS_OD__CGC__logFC",
    "DPA__GLASS_OD__CGC__t",
    "DPA__GLASS_OD__CGC__P.Value",
    "DPA__GLASS_OD__CGC__adj.P.Val",
    
    "DPA__GLASS_OD__prim-rec__pat_corrected__logFC",
    "DPA__GLASS_OD__prim-rec__pat_corrected__t",
    "DPA__GLASS_OD__prim-rec__pat_corrected__P.Value",
    "DPA__GLASS_OD__prim-rec__pat_corrected__adj.P.Val",
    
    "DPA__GLASS_OD__grade__pat_corrected__logFC",
    "DPA__GLASS_OD__grade__pat_corrected__t",
    "DPA__GLASS_OD__grade__pat_corrected__P.Value",
    "DPA__GLASS_OD__grade__pat_corrected__adj.P.Val"
  )


write.table(out, "output/tables/tab_DPA_outcomes_GLASS-OD.txt", sep="\t", quote = F, row.names = F)



## glass-nl outcomes ----


out <- metadata.proteomics.glass_nl |> 
  dplyr::select(
    "gene_id",
    "DPA__GLASS_NL__CGC__logFC",
    "DPA__GLASS_NL__CGC__t",
    "DPA__GLASS_NL__CGC__P.Value",
    "DPA__GLASS_NL__CGC__adj.P.Val"
  )

write.table(out, "output/tables/tab_DPA_outcomes_GLASS-NL.txt", sep="\t", quote = F, row.names = F)




## glass-od normalised data ----




metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, starts_with("proteomics_")) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 137) 
    return(.)
  })()




data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id)
data$nas <- rowSums(is.na(data))
data <- data |> 
  dplyr::filter(nas < 60) |> 
  dplyr::mutate(nas=NULL)



stopifnot(metadata$proteomics_id == colnames(data))


write.table(
  rbind(
    c("#", metadata$resection_id),
    c("protein_id", metadata$proteomics_id)
    ),
   "output/tables/tab_protein_matrix_normalised.txt", sep="\t", quote = F, row.names = F, col.names = F)


write.table(
  data,
  "output/tables/tab_protein_matrix_normalised.txt", 
  append=T,
  sep="\t", quote = F, row.names = T, col.names = F)





