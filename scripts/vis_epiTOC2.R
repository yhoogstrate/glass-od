#!/usr/bin/env R

# load ----


# deconvolution for cycling fraction?


source('scripts/load_functions.R')
source('scripts/load_themes.R')


library(ggplot2)


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# GLASS-OD / OD ----


metadata <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(163)


## cor ----


stat <- metadata |> 
  dplyr::select(array_A_IDH_HG__A_IDH_LG_lr_v12.8 , array_PC1, array_PC2, array_PC3, array_PC4, array_epiTOC2_tnsc, array_epiTOC2_hypoSC, array_percentage.detP.signi) |>
  dplyr::mutate(`-array_epiTOC2_hypoSC`= array_epiTOC2_hypoSC * -1, array_epiTOC2_hypoSC = NULL) |> 
  dplyr::mutate(`-array_PC2`= array_PC2 * -1,array_PC2 = NULL) |> 
  as.matrix()


corr <- cor(stat, method = "spearman")



corrplot::corrplot(corr, order="hclust", tl.cex=0.7) # , addgrid.col=NA




plt <- metadata |> 
  dplyr::select(array_A_IDH_HG__A_IDH_LG_lr_v12.8 , array_PC1, array_PC2, array_PC3, array_PC4, array_epiTOC2_tnsc, array_epiTOC2_hypoSC) 


ggplot(plt, aes(x=array_PC1, y=array_epiTOC2_hypoSC)) +
  geom_point() +
  theme_cellpress


ggplot(plt, aes(x=array_A_IDH_HG__A_IDH_LG_lr_v12.8, y=array_epiTOC2_tnsc)) +
  geom_point() +
  theme_cellpress


ggplot(plt, aes(x=array_A_IDH_HG__A_IDH_LG_lr_v12.8, y=array_epiTOC2_hypoSC)) +
  geom_point() +
  theme_cellpress



