#!/usr/bin/env R


# load data ----


library(ggplot2)


if(!exists('youri_gg_theme')) {
  source('scripts/load_themes.R')
}


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}

if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}

if(!exists('gsam.metadata.array_samples')) {
  source('scripts/load_G-SAM_metadata.R')
}


source('scripts/load_functions.R')



# load all samples & data ----


metadata <- rbind(
  glass_od.metadata.array_samples |> 
    dplyr::filter(arraychip_version == 'EPICv1') |> 
    dplyr::select(array_sentrix_id, array_percentage.detP.signi, 
                  #array_mnp_QC_v12.8_predicted_sample_type,
                  array_qc.pca.comp1, array_qc.pca.detP.outlier) |> 
    dplyr::mutate(dataset = "GLASS-OD") # & CATNON sample -> codels & oligosarcoma's
  ,
  glass_nl.metadata.array_samples |> 
    dplyr::select(array_sentrix_id, array_percentage.detP.signi,
                  #array_mnp_QC_v12.8_predicted_sample_type,
                  array_qc.pca.comp1, array_qc.pca.detP.outlier) |> 
    dplyr::mutate(dataset = "GLASS-NL") # astro's
  ,
  gsam.metadata.array_samples |> 
    dplyr::select(array_sentrix_id, array_percentage.detP.signi, 
                  #array_mnp_QC_v12.8_predicted_sample_type,
                  array_qc.pca.comp1, array_qc.pca.detP.outlier) |> 
    dplyr::mutate(dataset = "G-SAM")
  ) |> 
  dplyr::rename(array_qc.pca.comp1.old = array_qc.pca.comp1) |> 
  dplyr::rename(array_qc.pca.detP.outlier.old = array_qc.pca.detP.outlier) |> 
  assertr::verify(!is.na(array_sentrix_id)) |> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_OD_VALIDATION_INCLUDED_SAMPLES + 235 + 79))
    return(.)
  })()


data <- readRDS("cache/mvalues.all_samples.Rds") |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == (CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_OD_VALIDATION_INCLUDED_SAMPLES + 235 + 79))
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
  tibble::rownames_to_column('array_sentrix_id') |> 
  assertr::verify(array_sentrix_id %in% metadata$array_sentrix_id) |> 
  assertr::verify(metadata$array_sentrix_id %in% array_sentrix_id)


metadata <- metadata |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')) |> 
  assertr::verify(!is.na(PC1))

rm(tmp)


## find the component w/ strongest cor w/ detP ----

#' is indeed former pc1
ggplot(metadata, aes(x=array_qc.pca.comp1.old, y=PC1)) + 
  geom_point()


metadata <- metadata |> 
  dplyr::mutate(array_qc.pca.detP.outlier.new = ifelse(PC1 >= 600 | array_percentage.detP.signi >= 2.5,T,F))


# those that changed status:
table(paste0(metadata$array_qc.pca.detP.outlier.old), paste0(metadata$array_qc.pca.detP.outlier.new))



cor(metadata$array_percentage.detP.signi, metadata$PC1, method="spearman")
cor(metadata$array_percentage.detP.signi, metadata$PC2, method="spearman")
cor(metadata$array_percentage.detP.signi, metadata$PC3, method="spearman")
cor(metadata$array_percentage.detP.signi, metadata$PC4, method="spearman")


ggplot(metadata, aes(x=PC1, y=array_percentage.detP.signi, col=array_qc.pca.detP.outlier.old)) +
  geom_point(data=subset(metadata, !is.na(array_qc.pca.detP.outlier.old))) +
  geom_point(data=subset(metadata, is.na(array_qc.pca.detP.outlier.old))) +
  coord_cartesian(ylim = c(0, 25)) +
  geom_hline(yintercept=2.5, col="red", lty=2, lwd=0.5) +
  geom_vline(xintercept=600, col="red", lty=2, lwd=0.5) +
  theme_bw()


ggplot(metadata, aes(x=PC3, y=array_percentage.detP.signi, col=qc.pca.detP.outlier)) +
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


# export ----

out <- metadata |> 
  dplyr::select(array_sentrix_id, PC1, array_qc.pca.detP.outlier.new) |> 
  dplyr::rename(array_qc.pca.detP.outlier = array_qc.pca.detP.outlier.new) |> 
  dplyr::rename(array_qc.pca.comp1 = PC1) |> 
  assertr::verify(!is.na(array_qc.pca.detP.outlier)) |> 
  assertr::verify(!is.na(array_qc.pca.comp1)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_GLASS_OD_ALL_SAMPLES + CONST_N_OD_VALIDATION_INCLUDED_SAMPLES + 235 + 79))
    return(.)
  })()



saveRDS(out, file="cache/unsupervised_qc_outliers_all.Rds")



# GLASS-NL low purity ----

plt <- glass_nl.metadata.array_samples |> 
  filter_GLASS_NL_idats(202) |> 
  dplyr::left_join(
    readRDS("cache/analysis_unsupervised_PCA_GLASS-NL_x.Rds"),
    by=c('sentrix_id'='sentrix_id'),
    suffix=c('','')
  ) |> 
  dplyr::mutate(label = dplyr::case_when(
    mnp_predictBrain_v12.8_cal_class == "CTRL_CBM" ~ "CTRL_...",
    mnp_predictBrain_v12.8_cal_class == "CTRL_CORPCAL" ~ "CTRL_...",
    mnp_predictBrain_v12.8_cal_class == "CTRL_HEMI" ~ "CTRL_...",
    mnp_predictBrain_v12.8_cal_class == "CTRL_HYPOTHAL" ~ "CTRL_...",
    mnp_predictBrain_v12.8_cal_class == "CTRL_REACTIVE" ~ "CTRL_...",
    T ~ mnp_predictBrain_v12.8_cal_class
  )) |> 
  dplyr::mutate(odd = ifelse(PC3 > 275 | RFpurity.absolute < 0.5,"yes","no")) |> 
  dplyr::mutate(cnv_profile = dplyr::case_when(
    Sample_Name %in% c("007_R3","007_R4","010_R2","018_R1","020_R3","029_R2","035_R2","039_R2","105_R2","128_R2","129_R2","135_R1","136_R1","141_R2","207_R1") ~ "poor",
    Sample_Name %in% c("103_R2","160_R1","178_R2") ~ "borderline",
    Sample_Name %in% c("025_R3","034_R2","126_R2","129_R3","206_R2","216_R3","122_R1") ~ "sufficient",
    T ~ "not checked"
  ))



ggplot(plt, aes(x=PC3, y=RFpurity.absolute, col=cnv_profile, label=Sample_Name)) +
  ggrepel::geom_text_repel(data = subset(plt, odd == "yes"),col="black",size=3,alpha=0.6,nudge_y = 0.005) +
  geom_point() +
  geom_vline(xintercept=300)



c = cor(plt |>  dplyr::select(RFpurity.estimate, RFpurity.absolute, mnp_predictBrain_v12.8_cal_CTRL_HEMI, paste0("PC",1:15)) |> as.matrix(), method="spearman")
corrplot::corrplot(c)

