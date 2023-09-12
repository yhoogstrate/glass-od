#!/usr/bin/env R


# load ----


library(ggplot2)
library(patchwork)


source('scripts/load_functions.R')
source('scripts/load_themes.R')


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# plot A_IDH_rat ----



plt <-
  rbind(
    glass_od.metadata.array_samples |> 
      filter_GLASS_OD_idats(163) 
    ,
  glass_od.metadata.array_samples |> 
    dplyr::filter(study_name != "GLASS-OD")
  ) |> 
  dplyr::mutate(patient_id = as.factor(patient_id)) |> 
  dplyr::mutate(x = as.numeric(patient_id) + 1)|> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(last_recurrence = resection_number != 1 & resection_number == max(resection_number)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(panel = dplyr::recode(study_name, `CATNON` = 'non-canonical codels [CATNON, ...]', `GLASS-OD`='canonical codels [GLASS-OD]')) 



p1 <- ggplot(plt, aes(x = x, y=resection_number, fill=A_IDH_HG__A_IDH_LG_lr )) +
  facet_grid(cols = vars(panel)) +
  geom_rect(aes(xmin=x - 0.45,xmax=x + 0.45,
                ymin=resection_number + 0.5 - 0.45, ymax=resection_number + 0.5 + 0.45)) +
  theme_bw() +
  scale_fill_gradient2(low = "darkgreen", mid="gray90", high = "red", midpoint=0) 


p2 <- ggplot(plt, aes(x = x, y=A_IDH_HG__A_IDH_LG_lr, group=patient_id )) +
  facet_grid(cols = vars(panel)) +
  geom_line(alpha=0.3) +
  geom_point(aes(col=resection_number == 1, shape=resection_number == 1)) +
  theme_bw()


p1 / p2


# MNP brain classes ----



plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(180, exclude.suspected.noncodels = F) |> 
  dplyr::select(patient_id | patient_suspected_noncodel | resection_number | contains("_cal_class")) |> 
  tidyr::pivot_longer(cols = c(array_mnp_predictBrain_v2.0.1_cal_class,
                               array_mnp_predictBrain_v12.5_cal_class,
                               array_mnp_predictBrain_v12.8_cal_class), names_to = "classifier_version", values_to="class") |> 
  dplyr::mutate(col = as.factor(dplyr::case_when(
    #class %in% c("A_IDH_HG", "OLIGOSARC_IDH") ~ "HG (A_IDH_HG / Oligosarcoma)",
    class %in% c("A_IDH", "A_IDH_LG") ~ "A_IDH [_LG]",
    class %in% c("O_IDH") ~ "OLI",
    T ~ class
  ))) |> 
  dplyr::mutate(patient_suspected_noncodel = dplyr::case_when(
    patient_suspected_noncodel == T ~ "Non-codel",
    patient_suspected_noncodel == F ~ "Codel",
    T ~ as.character("-")
  )) |> 
  dplyr::mutate(classifier_version_txt = factor(gsub("^array_mnp_predictBrain_(.+)_cal_class$","\\1",classifier_version), levels=c("v2.0.1", "v12.5",  "v12.8" )))


plt.counts <- plt |>  dplyr::filter(classifier_version == "array_mnp_predictBrain_v12.8_cal_class") |> dplyr::pull("patient_suspected_noncodel") |>  table()
plt <- plt |> dplyr::mutate(patient_suspected_noncodel = paste0(patient_suspected_noncodel, " n=",plt.counts[plt$patient_suspected_noncodel]) )


cols = c('A_IDH [_LG]' = 'lightblue',
         'A_IDH_HG' = 'darkblue',
         
         'CTRL_CORPCAL' = 'orange',
         'GBM_RTK_I' = 'red',
         
         'OLI' = 'lightgreen',
         'OLIGOSARC_IDH' = 'darkgreen'
         )


ggplot(plt, aes(x = patient_id, y = resection_number, col=col)) +
  facet_grid(cols=vars(patient_suspected_noncodel),  rows=vars(classifier_version_txt), scales = "free", space="free") +
  geom_point(pch=15,size=1.75,alpha=0.65) +
  scale_color_manual(values=cols) +
  labs(x = "patient",y="resection #", col="",fill="") +
  ylim(0.5,5.5) +
  theme_cellpress_with_facet +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


ggsave("output/figures/vis_cohort_overview_MNP_classes.pdf", width=8.5 * 0.95, heigh = 2.5)



# misclassifications more common further away from Dx & R1 -
## no clear timing effect, but clear recurrence effect


## p ~ r ----


plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(163) |> 
  dplyr::filter(!is.na(resection_date) & !is.na(patient_diagnosis_date)) |> 
  #dplyr::select(resection_id, patient_diagnosis_date, resection_date) |> 
  dplyr::mutate(time_after_dx = resection_date - patient_diagnosis_date) |> 
  dplyr::filter(time_after_dx >= 0) |> 
  dplyr::mutate(misclass = dplyr::case_when(
    array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH","A_IDH_HG","A_IDH_LG") ~ "Yes: Astro",
    array_mnp_predictBrain_v12.8_cal_class %in% c("O_IDH","OLIGOSARC_IDH") ~ "No: Oligo",
    T ~ "Yes - other"
  )) |> 
  dplyr::mutate(stage_disease = ifelse(resection_number == 1, "primary tumours", "recurrent tumours")) |> 
  dplyr::select(misclass, stage_disease) |> 
  table() |> 
  as.data.frame()



ggplot(plt, aes(x=stage_disease, y = Freq, fill=misclass)) + 
  geom_bar(stat="identity") +
  labs(y= "frequency", x=NULL, col="Oligo misclassified") +
  theme_cellpress



