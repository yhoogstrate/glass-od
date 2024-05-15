#!/usr/bin/env R

# load ----


source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')


library(ggplot2)
library(patchwork)


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# GLASS-OD / OD ----


metadata <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES)


# fig ----
# 


plt <- metadata |> 
  dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = ifelse(array_mnp_predictBrain_v12.8_cal_class %in% c(
    "O_IDH"    ,     "A_IDH_HG"   ,   "A_IDH_LG"  ,    "OLIGOSARC_IDH"
  ) == F, "Other", array_mnp_predictBrain_v12.8_cal_class))


ggplot(plt, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=array_methylation_bins_1p19q_purity, col=array_mnp_predictBrain_v12.8_cal_class, label=resection_id)) +
  geom_point(size=theme_nature_size/3) +

  labs(subtitle=format_subtitle("Cohort overview")) +
  scale_y_continuous(limits = c(0, 0.8)) +
  
  scale_color_manual(values = c(palette_mnp_12.8_6, `Other`="darkgray")) +
  labs(x="CGC[Ac] (Lasso)",y="Tumor cell fraction", col="MNP v12.8") +
  theme_nature


ggsave("output/figures/vis_oligosarcoma_and_purity__glass-od_scatter.pdf", width=8.5 * 0.975 / 3, height = 2.85)




# paths

exp <- metadata |> 
  dplyr::filter(resection_id %in% c("0017-R3","0008-R2","0121-R3","0054-R3","0047-R3","0069-R2","0024-R2"))


