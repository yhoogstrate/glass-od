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


if(!exists('cortex_dilution')) {
  source('scripts/load_nonmalignant_cortex_dilution_series.R')
}



# GLASS-OD / OD ----


metadata <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES)


# fig 1 ----
## our CNV bin purity fit ----


plt <- metadata |> 
  dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = ifelse(array_mnp_predictBrain_v12.8_cal_class %in% c("O_IDH", "A_IDH_HG", "A_IDH_LG","OLIGOSARC_IDH") == F, "Other", array_mnp_predictBrain_v12.8_cal_class)) |> 
  dplyr::mutate(selected_for_spike_in_with_non_tumor = resection_id %in% c("0017-R3","0008-R2","0121-R3","0054-R3"))


ggplot(plt, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=array_methylation_bins_1p19q_purity, col=array_mnp_predictBrain_v12.8_cal_class, label=resection_id)) +
  geom_point(size=theme_nature_size/3) +
  ggrepel::geom_text_repel(data=subset(plt, selected_for_spike_in_with_non_tumor), col="darkgray", nudge_x = 0.5, nudge_y = 0.03, size=theme_nature_size) +
  
  labs(subtitle=format_subtitle("Cohort overview")) +
  scale_y_continuous(limits = c(0, 0.8)) +
  
  scale_color_manual(values = c(palette_mnp_12.8_6, `Other`="darkgray")) +
  labs(x="CGC[Ac] (Lasso)",y="Tumor cell fraction (1p/19q logFC fit)", col="MNP v12.8") +
  theme_nature


ggsave("output/figures/vis_oligosarcoma_and_purity__glass-od_scatter.pdf", width=8.5 * 0.975 / 4, height = 2.4)


rm(plt)



## RFpurity validation ----


plt <- metadata |> 
  dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = ifelse(array_mnp_predictBrain_v12.8_cal_class %in% c("O_IDH", "A_IDH_HG", "A_IDH_LG","OLIGOSARC_IDH") == F, "Other", array_mnp_predictBrain_v12.8_cal_class)) |> 
  dplyr::mutate(selected_for_spike_in_with_non_tumor = resection_id %in% c("0017-R3","0008-R2","0121-R3","0054-R3"))


ggplot(plt, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=array_RFpurity_absolute, col=array_mnp_predictBrain_v12.8_cal_class, label=resection_id)) +
  geom_point(size=theme_nature_size/3) +
  ggrepel::geom_text_repel(data=subset(plt, selected_for_spike_in_with_non_tumor), col="darkgray", nudge_x = 0.5, nudge_y = 0.03, size=theme_nature_size) +
  
  labs(subtitle=format_subtitle("Cohort overview")) +
  #scale_y_continuous(limits = c(0, 0.8)) +
  
  scale_color_manual(values = c(palette_mnp_12.8_6, `Other`="darkgray")) +
  labs(x="CGC[Ac] (Lasso)",y="Tumor cell fraction (RFpurity: ABS)", col="MNP v12.8") +
  theme_nature


ggsave("output/figures/vis_oligosarcoma_and_purity__glass-od_scatter__RFpurity_absolute.pdf", width=8.5 * 0.975 / 4, height = 2.4)


rm(plt)



## RFpurity validation ----


plt <- metadata |> 
  dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = ifelse(array_mnp_predictBrain_v12.8_cal_class %in% c("O_IDH", "A_IDH_HG", "A_IDH_LG","OLIGOSARC_IDH") == F, "Other", array_mnp_predictBrain_v12.8_cal_class)) |> 
  dplyr::mutate(selected_for_spike_in_with_non_tumor = resection_id %in% c("0017-R3","0008-R2","0121-R3","0054-R3"))


ggplot(plt, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=array_RFpurity_estimate, col=array_mnp_predictBrain_v12.8_cal_class, label=resection_id)) +
  geom_point(size=theme_nature_size/3) +
  ggrepel::geom_text_repel(data=subset(plt, selected_for_spike_in_with_non_tumor), col="darkgray", nudge_x = 0.5, nudge_y = 0.03, size=theme_nature_size) +
  
  labs(subtitle=format_subtitle("Cohort overview")) +
  #scale_y_continuous(limits = c(0, 0.8)) +
  
  scale_color_manual(values = c(palette_mnp_12.8_6, `Other`="darkgray")) +
  labs(x="CGC[Ac] (Lasso)",y="Tumor cell fraction (RFpurity: EST)", col="MNP v12.8") +
  theme_nature


ggsave("output/figures/vis_oligosarcoma_and_purity__glass-od_scatter__RFpurity_estimate.pdf", width=8.5 * 0.975 / 4, height = 2.4)


rm(plt)




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





# fig 4g: spike-in paths ----


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
                  array_fraction_normal_brain_cli, array_control_id) |> 
  dplyr::rename(array_fraction_normal_brain = array_fraction_normal_brain_cli)

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




# Volcano Oligosarcoma - A_IDH_HG ----



plt.intervals <- 10^seq(log10(0.01), log10(1), by=2/6) # linear intervals in log scaled labels


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::mutate(significant = DMP__GLASS_OD__oligosarcoma__A_IDH_HG__PC1__adj.P.Val < 0.01)


n_oligosarcoma <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(array_mnp_predictBrain_v12.8_cal_class == "OLIGOSARC_IDH") |> 
  nrow()

n_a_idh_hg <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(array_mnp_predictBrain_v12.8_cal_class == "A_IDH_HG") |> 
  nrow()


ggplot(plt, aes(x=DMP__GLASS_OD__oligosarcoma__A_IDH_HG__PC1__logFC,
                y=DMP__GLASS_OD__oligosarcoma__A_IDH_HG__PC1__adj.P.Val)) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_hline(yintercept=plt.intervals[2:6], col="#DDDDDD", lwd=theme_nature_lwd) +
  geom_hline(yintercept=plt.intervals[1], col="red", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=plt.intervals[7], col="darkgray", lwd=theme_nature_lwd) +
  geom_point(pch=16, cex=0.01, alpha=0.05) +
  theme_nature +
  scale_y_continuous(trans=reverselog_trans(base=10),
                     breaks=plt.intervals,
                     labels=round(plt.intervals, 3)
                     #labels=trans_format("identity", function(x) -x)
  ) +
  labs(x = "log2FC",
       y = "-log10(adj. P)",
       subtitle=format_subtitle("VolcanoPlot DMP Oligsarcoma vs. A_IDH_HG"),
       caption=paste0("Oligodendrogliomas classified as: Oligosarcoma: n=",
                      n_oligosarcoma,
                      "  --  Astrocytoma: n=",
                      n_a_idh_hg)) +
  coord_cartesian(xlim=c(-2.5, 2.5)) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export


ggsave("output/figures/vis_oligosarcoma_and_purity__volcano_A_IDH_HG.png", width=(8.5 * 0.975)/2,height=2.375, dpi=300)


rm(plt, n_oligosarcoma, n_a_idh_hg, plt.intervals)



# sandbox ----


#!/usr/bin/env R

p1 = ggplot(metadata.hg_olsc, aes(x=array_mnp_predictBrain_v12.8_cal_class, y=array_methylation_bins_1p19q_purity, label=isolation_id)) +
  ggbeeswarm::geom_quasirandom(size=theme_cellpress_size/2) +
  theme_cellpress

p2 = ggplot(metadata.hg_olsc, aes(x=array_mnp_predictBrain_v12.8_cal_class, y=array_median.overall.methylation, label=isolation_id)) +
  ggbeeswarm::geom_quasirandom(size=theme_cellpress_size/2) +
  theme_cellpress


p3 = ggplot(metadata.hg_olsc, aes(x=array_mnp_predictBrain_v12.8_cal_class, y=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, label=isolation_id)) +
  ggbeeswarm::geom_quasirandom(size=theme_cellpress_size/2) +
  theme_cellpress


p4 = ggplot(metadata.hg_olsc, aes(x=array_mnp_predictBrain_v12.8_cal_class, y=array_qc.pca.comp1, label=isolation_id)) +
  ggbeeswarm::geom_quasirandom(size=theme_cellpress_size/2) +
  theme_cellpress


p1 + p2 + p3 + p4


