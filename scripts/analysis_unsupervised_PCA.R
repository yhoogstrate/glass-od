#!/usr/bin/env R

# load data ----


if(!exists('glass_od.metadata.idats')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}



# GLASS-OD ----


metadata <- glass_od.metadata.idats |> 
  filter_GLASS_OD_idats(163)


data <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata$sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (694299))
    return(.)
  })()



# use all probes - no MAD filter, also no further detP filter
#(function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |> # this synthax, oh my
#  dplyr::arrange(mad) |> 
#  dplyr::mutate(mad = NULL)



data.pca.glass_od <- data |> 
  t() |> 
  prcomp()

data.pca.glass_od.x <- data.pca.glass_od |> 
  purrr::pluck('x') |> 
  as.data.frame(stringsAsFactors=F) |> 
  tibble::rownames_to_column('sentrix_id')


saveRDS(data.pca.glass_od, file="cache/analysis_unsupervised_PCA_GLASS-OD_prcomp.Rds")
saveRDS(data.pca.glass_od.x, file="cache/analysis_unsupervised_PCA_GLASS-OD_x.Rds")



rm(data, metadata, data.pca.glass_od)




# GLASS-NL ----





