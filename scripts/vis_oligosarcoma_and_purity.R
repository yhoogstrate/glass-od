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


# fig 1 ----



plt <- metadata |> 
  dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = ifelse(array_mnp_predictBrain_v12.8_cal_class %in% c("O_IDH", "A_IDH_HG", "A_IDH_LG","OLIGOSARC_IDH") == F, "Other", array_mnp_predictBrain_v12.8_cal_class)) |> 
  dplyr::mutate(selected_for_spike_in_with_non_tumor = resection_id %in% c("0017-R3","0008-R2","0121-R3","0054-R3"))


ggplot(plt, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=array_methylation_bins_1p19q_purity, col=array_mnp_predictBrain_v12.8_cal_class, label=resection_id)) +
  geom_point(size=theme_nature_size/3) +
  ggrepel::geom_text_repel(data=subset(plt, selected_for_spike_in_with_non_tumor), col="darkgray", nudge_x = 0.5, nudge_y = 0.03, size=theme_nature_size) +

  labs(subtitle=format_subtitle("Cohort overview")) +
  scale_y_continuous(limits = c(0, 0.8)) +
  
  scale_color_manual(values = c(palette_mnp_12.8_6, `Other`="darkgray")) +
  labs(x="CGC[Ac] (Lasso)",y="Tumor cell fraction", col="MNP v12.8") +
  theme_nature


ggsave("output/figures/vis_oligosarcoma_and_purity__glass-od_scatter.pdf", width=8.5 * 0.975 / 4, height = 2.4)



# fig 2 ----


plt <- metadata |> 
  dplyr::filter(array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG" , "OLIGOSARC_IDH")) |> 
  #dplyr::filter(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit > 1) |> 
  dplyr::mutate(lr_OLSC_A_HG = log(array_mnp_predictBrain_v12.8_cal_OLIGOSARC_IDH / array_mnp_predictBrain_v12.8_cal_A_IDH_HG))


ggplot(plt, aes(x= lr_OLSC_A_HG, y=array_methylation_bins_1p19q_purity, col=array_mnp_predictBrain_v12.8_cal_class)) +
  geom_point(size=theme_nature_size/3) +
  
  labs(subtitle=format_subtitle("Cohort overview")) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="black", cor.coef.name ="rho", size=theme_nature_size) +
  
  scale_color_manual(values = c(palette_mnp_12.8_6, `Other`="darkgray")) +
  #labs(x="CGC[Ac] (Lasso)",y="Tumor cell fraction", col="MNP v12.8") +
  theme_nature



# fig 3a ----


plt <- metadata |> 
  dplyr::filter(array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG" , "OLIGOSARC_IDH"))


ggplot(plt, aes(x=array_mnp_predictBrain_v12.8_cal_class, y=array_methylation_bins_1p19q_purity, col=array_mnp_predictBrain_v12.8_cal_class)) +
  geom_point(size=theme_nature_size/3) +
  
  labs(subtitle=format_subtitle("Purity oligosarc ~ A_IDH_HG")) +
  ggpubr::stat_compare_means( label.x.npc=0.25, method = "t.test", show_guide  = FALSE,  size=theme_nature_size) +
  
  scale_color_manual(values = c(palette_mnp_12.8_6, `Other`="darkgray")) +
  labs(x="",y="Tumor cell fraction", col="MNP v12.8") +
  theme_nature


ggsave("output/figures/vis_oligosarcoma_and_purity__glass-od_barplot.pdf", width=8.5 * 0.975 / 7, height = 2.4)


# fig 3b ----

plt <- metadata |> 
  dplyr::filter(array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG" , "OLIGOSARC_IDH")) |> 
  dplyr::mutate(resection_treatment_status_chemo = dplyr::case_when(
    resection_treatment_status_chemo == T ~ "Chemo +",
    resection_treatment_status_chemo == F ~ "Chemo -",
    T ~ as.character(NA)
  ))

ggplot(plt, aes(x=resection_treatment_status_chemo, y=array_methylation_bins_1p19q_purity, col=array_mnp_predictBrain_v12.8_cal_class)) +
  geom_point(size=theme_nature_size/3) +
  
  labs(subtitle=format_subtitle("Purity oligosarc ~ A_IDH_HG")) +
  ggpubr::stat_compare_means(comparisons=list(c('Chemo -','Chemo +')), label.x.npc=0.25, method = "t.test", show_guide  = FALSE,  size=theme_nature_size) +
  
  scale_color_manual(values = c(palette_mnp_12.8_6, `Other`="darkgray")) +
  labs(x="",y="Tumor cell fraction", col="MNP v12.8") +
  theme_nature

ggsave("output/figures/vis_oligosarcoma_and_purity__glass-od_barplot_chemo.pdf", width=8.5 * 0.975 / 7, height = 2.4)



# fig 3c ----


plt <- metadata |> 
  dplyr::filter(array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG" , "OLIGOSARC_IDH")) |> 
  dplyr::mutate(resection_treatment_status_radio = dplyr::case_when(
    resection_treatment_status_radio == T ~ "RT +",
    resection_treatment_status_radio == F ~ "RT -",
    T ~ as.character(NA)
  ))

ggplot(plt, aes(x=resection_treatment_status_radio, y=array_methylation_bins_1p19q_purity, col=array_mnp_predictBrain_v12.8_cal_class)) +
  geom_point(size=theme_nature_size/3) +
  
  labs(subtitle=format_subtitle("Purity oligosarc ~ A_IDH_HG")) +
  ggpubr::stat_compare_means(comparisons=list(c('RT -','RT +')), label.x.npc=0.25, method = "t.test", show_guide  = FALSE,  size=theme_nature_size) +
  
  scale_color_manual(values = c(palette_mnp_12.8_6, `Other`="darkgray")) +
  labs(x="",y="Tumor cell fraction", col="MNP v12.8") +
  theme_nature

ggsave("output/figures/vis_oligosarcoma_and_purity__glass-od_barplot_rt.pdf", width=8.5 * 0.975 / 7, height = 2.4)





# fig 4: spike-in paths ----


plt.tumor <- rbind(
  glass_od.metadata.array_samples |> 
    dplyr::filter(resection_id %in% c("0017-R3","0008-R2","0121-R3","0054-R3")) |> 
    dplyr::rename(array_tumor_id = resection_id) |> 
    dplyr::select(array_tumor_id, array_sentrix_id, array_mnp_predictBrain_v12.8_cal_class, array_methylation_bins_1p19q_purity) |> 
    dplyr::mutate(array_fraction_normal_brain = 0) |> 
    dplyr::mutate(array_control_id = "s107"),
  glass_od.metadata.array_samples |> 
    dplyr::filter(resection_id %in% c("0017-R3","0008-R2","0121-R3","0054-R3")) |> 
    dplyr::rename(array_tumor_id = resection_id) |> 
    dplyr::select(array_tumor_id, array_sentrix_id, array_mnp_predictBrain_v12.8_cal_class, array_methylation_bins_1p19q_purity) |> 
    dplyr::mutate(array_fraction_normal_brain = 0) |> 
    dplyr::mutate(array_control_id = "s108"),
  glass_od.metadata.array_samples |> 
    dplyr::filter(resection_id %in% c("0017-R3","0008-R2","0121-R3","0054-R3")) |> 
    dplyr::rename(array_tumor_id = resection_id) |> 
    dplyr::select(array_tumor_id, array_sentrix_id, array_mnp_predictBrain_v12.8_cal_class, array_methylation_bins_1p19q_purity) |> 
    dplyr::mutate(array_fraction_normal_brain = 0) |> 
    dplyr::mutate(array_control_id = "s110"),
  glass_od.metadata.array_samples |> 
    dplyr::filter(resection_id %in% c("0017-R3","0008-R2","0121-R3","0054-R3")) |> 
    dplyr::rename(array_tumor_id = resection_id) |> 
    dplyr::select(array_tumor_id, array_sentrix_id, array_mnp_predictBrain_v12.8_cal_class, array_methylation_bins_1p19q_purity) |> 
    dplyr::mutate(array_fraction_normal_brain = 0) |> 
    dplyr::mutate(array_control_id = "s112"))


plt.spiked <- cortex_dilution.metadata.array_samples |> 
    dplyr::select(array_tumor_id, array_sentrix_id, array_mnp_predictBrain_v12.8_cal_class, array_methylation_bins_1p19q_purity,
                  array_fraction_normal_brain, array_control_id)

plt.normal_cortex <- rbind(
  normal_cortex.metadata.array_samples |> 
    dplyr::filter(!is.na(array_control_id)) |> 
    dplyr::select(array_control_id, array_sentrix_id, array_mnp_predictBrain_v12.8_cal_class, array_methylation_bins_1p19q_purity) |> 
    dplyr::mutate(array_tumor_id = "0017-R3") |> 
    dplyr::mutate(array_fraction_normal_brain = 100),
  normal_cortex.metadata.array_samples |> 
    dplyr::filter(!is.na(array_control_id)) |> 
    dplyr::select(array_control_id, array_sentrix_id, array_mnp_predictBrain_v12.8_cal_class, array_methylation_bins_1p19q_purity) |> 
    dplyr::mutate(array_tumor_id = "0008-R2") |> 
    dplyr::mutate(array_fraction_normal_brain = 100),
  normal_cortex.metadata.array_samples |> 
    dplyr::filter(!is.na(array_control_id)) |> 
    dplyr::select(array_control_id, array_sentrix_id, array_mnp_predictBrain_v12.8_cal_class, array_methylation_bins_1p19q_purity) |> 
    dplyr::mutate(array_tumor_id = "0121-R3") |> 
    dplyr::mutate(array_fraction_normal_brain = 100),
  normal_cortex.metadata.array_samples |> 
    dplyr::filter(!is.na(array_control_id)) |> 
    dplyr::select(array_control_id, array_sentrix_id, array_mnp_predictBrain_v12.8_cal_class, array_methylation_bins_1p19q_purity) |> 
    dplyr::mutate(array_tumor_id = "0054-R3") |> 
    dplyr::mutate(array_fraction_normal_brain = 100)
)


plt <- rbind(plt.tumor, plt.spiked, plt.normal_cortex) |> 
  dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = ifelse(array_mnp_predictBrain_v12.8_cal_class %in% c("O_IDH", "A_IDH_HG", "A_IDH_LG","OLIGOSARC_IDH") == F, "Other", array_mnp_predictBrain_v12.8_cal_class))


ggplot(plt, aes(x=array_fraction_normal_brain, y=array_methylation_bins_1p19q_purity, col=array_mnp_predictBrain_v12.8_cal_class)) +
  facet_grid(cols = vars(array_tumor_id), rows = vars(array_control_id)) +
  
  geom_vline(xintercept=-5, lwd=theme_nature_lwd) +
  geom_hline(yintercept=0, lwd=theme_nature_lwd) +
  
  geom_point(size = theme_nature_size/3) +
  scale_color_manual(values = c(palette_mnp_12.8_6, `Other`="darkgray")) +
  scale_x_continuous(breaks=c(0, 20,  40, 60, 80, 100)) +
  coord_cartesian(xlim=c(0, 100), ylim=c(0, max(plt$array_methylation_bins_1p19q_purity))) +
  labs(x="Dilution fraction (normalised for array depth)",
       y = "Estimated purity (1P/19Q)", 
       col="",
       subtitle=format_subtitle("A_IDH_HG with spiked-in non malignant cortex"),
      ) +
  theme_nature


ggsave("output/figures/vis_oligosarcoma_and_purity__spike-in_nm_cortex.pdf", width=8.5 * 0.975 * 0.5 * 1.065, height = 3.5)




