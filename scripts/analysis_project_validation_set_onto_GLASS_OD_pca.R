#!/usr/bin/env


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}

if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}



# select data ----


metadata <- glass_od.metadata.array_samples |> 
  filter_OD_validation_idats(CONST_N_OD_VALIDATION_INCLUDED_SAMPLES) 


data <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
    return(.)
  })()


# multiplication ---


data.pca.glass_od <- readRDS(file="cache/analysis_unsupervised_PCA_GLASS-OD_prcomp.Rds")


stats <- data |> 
  dplyr::select(metadata$array_sentrix_id)
  

pca.lifted_over <- predict(data.pca.glass_od, t(stats))|> 
  as.data.frame(stringsAsFactor=F) |>   # transform back from matrix to data.frame
  tibble::rownames_to_column('sentrix_id') |> 
  dplyr::rename_with( ~ paste0("array_", .x))



saveRDS(pca.lifted_over, file="cache/analysis_unsupervised_PCA_GLASS-OD__lifted_over_to_VALIDATION_x.Rds")


rm(data.pca.glass_od, stats, pca.lifted_over)
rm(data, metadata)



