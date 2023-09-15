#!/usr/bin/env R

source('scripts/load_functions.R')
source('scripts/load_themes.R')

library(ggplot2)

if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# check PC's ----

#
# plot(metadata$array_PC3 , metadata$array_percentage.detP.signi)
# load samples, x-check duplicates per patient

metadata <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(163) |> 
  dplyr::filter(array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH", "A_IDH_HG")) |> 
  assertr::verify(!duplicated(paste0(patient_id,".",array_mnp_predictBrain_v12.8_cal_class)))


for(pc in 1:55) {
  w = wilcox.test(
  metadata |> dplyr::filter(array_mnp_predictBrain_v12.8_cal_class == "OLIGOSARC_IDH") |> dplyr::pull(paste0("array_PC", pc)),
  metadata |> dplyr::filter(array_mnp_predictBrain_v12.8_cal_class == "A_IDH_HG") |> dplyr::pull(paste0("array_PC", pc))
  )
  
  if(w$p.value < 0.01) {
    print(paste0("PC",pc,": p=",w$p.value))  
  }
}



cor(plt$array_methylation_bins_1p19q_purity, plt$array_PC1, method="spearman")
cor(plt$array_methylation_bins_1p19q_purity, plt$array_PC2, method="spearman")
cor(plt$array_methylation_bins_1p19q_purity, plt$array_PC3, method="spearman")
cor(plt$array_methylation_bins_1p19q_purity, plt$array_PC4, method="spearman")
cor(plt$array_methylation_bins_1p19q_purity, plt$array_PC5, method="spearman")
cor(plt$array_methylation_bins_1p19q_purity, plt$array_PC6, method="spearman")




plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(163)

ggplot(plt, aes(x=-array_PC2, y=array_PC3, col = array_mnp_predictBrain_v12.8_cal_class)) +
  geom_point()


ggplot(plt, aes(x=-array_PC2, y=array_methylation_bins_1p19q_purity, col = array_mnp_predictBrain_v12.8_cal_class)) +
  geom_point()

ggplot(plt, aes(x=array_PC3, y=array_methylation_bins_1p19q_purity, col = array_mnp_predictBrain_v12.8_cal_class)) +
  geom_point()


ggplot(plt, aes(x=array_PC36, y=array_PC2, col = array_mnp_predictBrain_v12.8_cal_class)) +
  geom_point()



