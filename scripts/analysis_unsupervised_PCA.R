#!/usr/bin/env R

# load data ----


if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}



if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}

if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}

if(!exists('gsam.metadata.array_samples')) {
  source('scripts/load_GSAM_metadata.R')
}




ggplot(glass_nl.metadata.array_samples, aes(x=qc.pca.pc3purity.outlier, y=array_PC1)) +
  geom_boxplot()



plt <- glass_nl.metadata.array_samples |> 
  dplyr::select(array_percentage.detP.signi,
                array_median.overall.methylation,
                array_qc.pca.comp1,
                array_epiTOC2_tnsc,
                array_PC1,
                array_PC2,
                array_PC3,
                array_PC4,
                array_PC5, 
                array_A_IDH_HG__A_IDH_LG_lr_v12.8
                ) |> 
  dplyr::filter(!is.na(array_PC1))

corrplot::corrplot(abs(cor(plt, method="pearson")), order="hclust", tl.cex=0.75, tl.pos="l")



# GLASS-OD ----


metadata <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) # should be going toward 211


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
  dplyr::rename_with( ~ paste0("array_", .x)) |> 
  tibble::rownames_to_column('array_sentrix_id')


saveRDS(data.pca.glass_od, file="cache/analysis_unsupervised_PCA_GLASS-OD_prcomp.Rds")
saveRDS(data.pca.glass_od.x, file="cache/analysis_unsupervised_PCA_GLASS-OD_x.Rds")



# tmp check ~ R=0.98 as compared without new samples
# tmp <- data.pca.glass_od.x |> 
#   dplyr::rename_with( ~ paste0(.x, "_new"), .cols=!matches("^array_sentrix_id$",perl = T))
# 
# plt <- glass_od.metadata.array_samples |> 
#   filter_GLASS_OD_idats(215) |> 
#   dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))
# 
# ggplot(plt, aes(x=array_PC2, y=array_PC2_new)) +
#   geom_point()




rm(data, metadata, data.pca.glass_od)


# GLASS-OD + OD-validatie ----


metadata <- rbind(
  glass_od.metadata.array_samples |> 
    filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES),
  glass_od.metadata.array_samples |> 
    filter_OD_validation_idats(CONST_N_OD_VALIDATION_INCLUDED_SAMPLES)
)



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



data.pca.glass_od_and_validation <- data |> 
  t() |> 
  prcomp()

data.pca.glass_od_and_validation.x <- data.pca.glass_od_and_validation |> 
  purrr::pluck('x') |> 
  as.data.frame(stringsAsFactors=F) |> 
  dplyr::rename_with(~ paste0("array_GLASS_OD_VALIDATION_", .x)) |> 
  tibble::rownames_to_column('array_sentrix_id')


saveRDS(data.pca.glass_od_and_validation, file="cache/analysis_unsupervised_PCA_GLASS-OD_AND_VALIDATION_prcomp.Rds")
saveRDS(data.pca.glass_od_and_validation.x, file="cache/analysis_unsupervised_PCA_GLASS-OD_AND_VALIDATION_x.Rds")



# GLASS-NL [all] ----


# metadata <- glass_nl.metadata.array_samples |>
#   dplyr::select(array_sentrix_id) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == 235)
#     return(.)
#   })()
# 
# 
# 
# data <- data.mvalues.hq_samples |> 
#   tibble::rownames_to_column('probe_id') |> 
#   dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
#   tibble::column_to_rownames('probe_id') |> 
#   dplyr::select(metadata$array_sentrix_id) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
#     assertthat::assert_that(ncol(.) == 235)
#     return(.)
#   })()
# 
# 
# data.pca.glass_nl <- data |> 
#   t() |> 
#   prcomp()
# 
# 
# data.pca.glass_nl.x <- data.pca.glass_nl |> 
#   purrr::pluck('x') |> 
#   as.data.frame(stringsAsFactors=F) |> 
#   dplyr::rename_with( ~ paste0("array_", .x)) |> 
#   tibble::rownames_to_column('array_sentrix_id')
# 
# 
# saveRDS(data.pca.glass_nl, file="cache/analysis_unsupervised_PCA_GLASS-NL_all_prcomp.Rds")
# saveRDS(data.pca.glass_nl.x, file="cache/analysis_unsupervised_PCA_GLASS-NL_all_x.Rds")
# 
# 
# rm(data, metadata, data.pca.glass_nl)


# GLASS-NL [hq] ----


metadata <- glass_nl.metadata.array_samples |>
  filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES) |> 
  dplyr::select(array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_NL_INCLUDED_SAMPLES)
    return(.)
  })()


data <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
    assertthat::assert_that(ncol(.) == CONST_N_GLASS_NL_INCLUDED_SAMPLES)
    return(.)
  })()


data.pca.glass_nl <- data |> 
  t() |> 
  prcomp()


data.pca.glass_nl.x <- data.pca.glass_nl |> 
  purrr::pluck('x') |> 
  as.data.frame(stringsAsFactors=F) |> 
  dplyr::rename_with( ~ paste0("array_", .x)) |> 
  tibble::rownames_to_column('array_sentrix_id')


saveRDS(data.pca.glass_nl, file="cache/analysis_unsupervised_PCA_GLASS-NL_prcomp.Rds")
saveRDS(data.pca.glass_nl.x, file="cache/analysis_unsupervised_PCA_GLASS-NL_x.Rds")


rm(data, metadata, data.pca.glass_nl)





# G-SAM ----



metadata <- gsam.metadata.array_samples |>
  filter_GSAM_idats(CONST_N_GSAM_INCLUDED_SAMPLES) |> 
  dplyr::select(array_sentrix_id)


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


data.pca.gsam <- data |> 
  t() |> 
  prcomp()


data.pca.gsam.x <- data.pca.gsam |> 
  purrr::pluck('x') |> 
  as.data.frame(stringsAsFactors=F) |> 
  dplyr::rename_with( ~ paste0("array_", .x)) |> 
  tibble::rownames_to_column('array_sentrix_id')


saveRDS(data.pca.gsam, file="cache/analysis_unsupervised_PCA_GSAM_prcomp.Rds")
saveRDS(data.pca.gsam.x, file="cache/analysis_unsupervised_PCA_GSAM_x.Rds")


rm(data, metadata, data.pca.gsam)





# GLASS-[OD+NL] combi ----

# 
# 
# metadata <- rbind(
#   glass_od.metadata.array_samples |> filter_GLASS_OD_idats(163) |> dplyr::select(sentrix_id),
#   glass_nl.metadata.array_samples |> filter_GLASS_NL_idats(218) |> dplyr::select(sentrix_id)
# )
# 
# 
# data <- data.mvalues.hq_samples |> 
#   tibble::rownames_to_column('probe_id') |> 
#   dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
#   tibble::column_to_rownames('probe_id') |> 
#   dplyr::select(metadata$sentrix_id) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
#     return(.)
#   })()
# 
# 
# 
# data.pca.glass_od_nl <- data |> 
#   t() |> 
#   prcomp()
# 
# data.pca.glass_od_nl.x <- data.pca.glass_od_nl |> 
#   purrr::pluck('x') |> 
#   as.data.frame(stringsAsFactors=F) |> 
#   tibble::rownames_to_column('sentrix_id')
# 
# 
# saveRDS(data.pca.glass_od_nl.x, "cache/analysis_unsupervised_PCA_GLASS-OD_GLASS-NL_combined.Rds")
# 
# 
# 
# # GLASS-[OD+NL] combi excl 1P/19Q ----
# 
# 
# 
# 
# 
# metadata <- rbind(
#   glass_od.metadata.array_samples |> filter_GLASS_OD_idats(163) |> dplyr::select(sentrix_id),
#   glass_nl.metadata.array_samples |> filter_GLASS_NL_idats(218) |> dplyr::select(sentrix_id)
# )
# 
# 
# data <- data.mvalues.hq_samples |> 
#   tibble::rownames_to_column('probe_id') |> 
#   dplyr::filter(probe_id %in% (
#     data.mvalues.probes |> 
#       dplyr::filter(good_probe) |> 
#       dplyr::filter((is_1P | is_19Q) == F) |> 
#       dplyr::pull(probe_id)
#   )) |> 
#   tibble::column_to_rownames('probe_id') |> 
#   dplyr::select(metadata$sentrix_id) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == (639635))
#     return(.)
#   })()
# 
# 
# 
# data.pca.glass_od_nl_excl1P19Q <- data |> 
#   t() |> 
#   prcomp()
# 
# data.pca.glass_od_nl_excl1P19Q.x <- data.pca.glass_od_nl_excl1P19Q |> 
#   purrr::pluck('x') |> 
#   as.data.frame(stringsAsFactors=F) |> 
#   tibble::rownames_to_column('sentrix_id')
# 
# 
# saveRDS(data.pca.glass_od_nl_excl1P19Q.x, "cache/analysis_unsupervised_PCA_GLASS-OD_GLASS-NL_combined_no_1P19Q.Rds")
# 

