#!/usr/bin/env R


# load data ----


library(ggplot2)


if(!exists('youri_gg_theme')) {
  source('scripts/youri_gg_theme.R')
}


if(!exists('glass_od.metadata.idats')) {
  source('scripts/load_GLASS-OD_metadata.R')
}

if(!exists('glass_nl.metadata.idats')) {
  source('scripts/load_GLASS-NL_metadata.R')
}

if(!exists('gsam.metadata.idats')) {
  source('scripts/load_G-SAM_metadata.R')
}


source('scripts/load_functions.R')



# load all samples & data ----


metadata <- rbind(
  glass_od.metadata.idats |> 
    dplyr::select(sentrix_id, percentage.detP.signi, mnp_QC_predicted_sample_type, qc.pca.comp1, qc.pca.outlier) |> 
    dplyr::mutate(dataset = "GLASS-OD") # & CATNON sample -> codels & oligosarcoma's
  ,
  glass_nl.metadata.idats |> 
    dplyr::select(sentrix_id, percentage.detP.signi, mnp_QC_predicted_sample_type) |> 
    dplyr::mutate(qc.pca.comp1 = NA, qc.pca.outlier = NA) |> 
    dplyr::mutate(dataset = "GLASS-NL") # astro's
  ,
  gsam.metadata.idats |> 
    dplyr::select(sentrix_id, percentage.detP.signi, mnp_QC_predicted_sample_type) |> 
    dplyr::mutate(qc.pca.comp1 = NA, qc.pca.outlier = NA) |> 
    dplyr::mutate(dataset = "G-SAM")
  ) |> 
  assertr::verify(!is.na(sentrix_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (222 + 235 + 75))
    return(.)
  })()


data <- readRDS("cache/mvalues.all_samples.Rds") |> 
  dplyr::select(metadata$sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == (222 + 235 + 75))
    return(.)
  })()




# all samples ----

## pca ----


pca.raw <- data |> 
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |>  # this synthax, oh my
  dplyr::slice_head(n=200000) |> 
  dplyr::mutate(mad = NULL) |>  # remove the sd to obtain original vst matrix
  t() |> 
  prcomp()


tmp <- pca.raw |> 
  purrr::pluck('x')  |>   # take coordinates
  as.data.frame(stringsAsFactor=F) |>  # transform back from matrix to data.frame 
  dplyr::select(paste0('PC',1:20)) |> 
  tibble::rownames_to_column('sentrix_id')

metadata <- metadata |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'),suffix=c('',''))

rm(tmp)


## find the component w/ strongest cor w/ detP ----


metadata <- metadata |> 
  dplyr::mutate(qc.pca.detP.outlier = ifelse(PC1 >= 600 | percentage.detP.signi >= 2.5,T,F))


cor(metadata$percentage.detP.signi, metadata$PC1, method="spearman")
cor(metadata$percentage.detP.signi, metadata$PC2, method="spearman")
cor(metadata$percentage.detP.signi, metadata$PC3, method="spearman")
cor(metadata$percentage.detP.signi, metadata$PC4, method="spearman")


ggplot(metadata, aes(x=PC1, y=percentage.detP.signi, col=qc.pca.detP.outlier)) +
  geom_point() +
  coord_cartesian(ylim = c(0, 25)) +
  geom_hline(yintercept=2.5, col="red", lty=2, lwd=0.5) +
  geom_vline(xintercept=600, col="red", lty=2, lwd=0.5) +
  theme_bw()


ggplot(metadata, aes(x=PC3, y=percentage.detP.signi, col=qc.pca.detP.outlier)) +
  geom_point() +
  coord_cartesian(ylim = c(0, 25)) +
  geom_hline(yintercept=2.5, col="red", lty=2, lwd=0.5) +
  geom_vline(xintercept=600, col="red", lty=2, lwd=0.5) +
  theme_bw()


ggplot(metadata, aes(x=PC1, y=PC3, col=qc.pca.detP.outlier)) +
  geom_point() +
  theme_bw()



ggplot(metadata, aes(x=qc.pca.comp1, y=percentage.detP.signi, col=qc.pca.detP.outlier)) +
  geom_point() +
  coord_cartesian(ylim = c(0, 25))


ggplot(metadata, aes(x=qc.pca.comp1, y=PC1, col=qc.pca.detP.outlier)) +
  geom_point() 



ggplot(metadata, aes(x=PC3, y=percentage.detP.signi, col=qc.pca.outlier)) +
  geom_point() +
  coord_cartesian(ylim = c(0, 25))

ggplot(metadata, aes(x=PC1,y=PC3, col=qc.pca.outlier)) +
  geom_point()




out <- metadata |> 
  dplyr::select(sentrix_id, PC1, qc.pca.detP.outlier) |> 
  dplyr::rename(qc.pca.comp1 = PC1)


saveRDS(out, file="cache/unsupervised_qc_outliers_all.Rds")



# GLASS-NL low purity ----

plt <- glass_nl.metadata.idats |> 
  filter_GLASS_NL_idats(218) |> 
  dplyr::left_join(
    readRDS("cache/analysis_unsupervised_PCA_GLASS-NL_x.Rds"),
    by=c('sentrix_id'='sentrix_id'),
    suffix=c('','')
  )


