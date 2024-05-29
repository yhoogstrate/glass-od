#!/usr/bin/env R


# load ----


library(ggplot2)
library(patchwork)


source('scripts/load_functions.R')
source('scripts/load_themes.R')


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


tmp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(is.na(patient_last_follow_up_date)) |> 
  dplyr::select(patient_id, resection_id, patient_center_name, patient_last_follow_up_date, patient_last_follow_up_event) |> 
  View()


tmp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(is.na(resection_treatment_status_summary) | is.na(resection_treatment_status_radio) | is.na(resection_treatment_status_chemo)) |> 
  dplyr::select(patient_id, resection_id, patient_center_name, 
                resection_treatment_status_radio, 
                resection_treatment_status_chemo, 
                resection_treatment_status_summary) |> 
  View()




# Figure 1A: overview GLASS_OD ----


plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES + 18, exclude.suspected.noncodels = F) |> 
  dplyr::mutate(Source = dplyr::recode(isolation_material, `ffpe`="Source: FFPE", `tissue`="Source: fresh")) |> 
  
  dplyr::mutate(resection_treatment_status_chemo = dplyr::case_when(
    patient_suspected_noncodel ~ as.character(NA),
    as.character(resection_treatment_status_chemo) == "TRUE" ~ "Yes",
    as.character(resection_treatment_status_chemo) == "FALSE" ~ "No",
    T ~ as.character(NA))
  ) |> 
  dplyr::mutate(resection_treatment_status_radio = dplyr::case_when(
    patient_suspected_noncodel ~ as.character(NA),
    as.character(resection_treatment_status_radio) == "TRUE" ~ "Yes",
    as.character(resection_treatment_status_radio) == "FALSE" ~ "No",
    T ~ as.character(NA))
  ) |> 
  dplyr::mutate(too_low_purity = dplyr::case_when(
    patient_suspected_noncodel ~ as.character(NA),
    array_methylation_bins_1p19q_purity <= 0.1 ~ "Yes",
    array_methylation_bins_1p19q_purity > 0.1 ~ "No",
    T ~ as.character(NA)
  )) |> 
  #dplyr::mutate(too_low_purity = ifelse(array_methylation_bins_1p19q_purity <= 0.1, "Yes", "No")) |> 
  dplyr::select(patient_id, patient_suspected_noncodel, resection_number, resection_tumor_grade,
                resection_treatment_status_chemo, resection_treatment_status_radio,
                Source, contains("_cal_class"), too_low_purity) |> 
  dplyr::mutate(resection_tumor_grade = dplyr::case_when(
    is.na(resection_tumor_grade) | patient_suspected_noncodel ~ as.character(NA),
    !patient_suspected_noncodel ~ paste0("Grade ", resection_tumor_grade)
  )) |> 
  tidyr::pivot_longer(cols = c(resection_tumor_grade,
                               
                               array_mnp_predictBrain_v2.0.1_cal_class,
                               array_mnp_predictBrain_v12.5_cal_class,
                               array_mnp_predictBrain_v12.8_cal_class,
                               
                               resection_treatment_status_chemo,
                               resection_treatment_status_radio,
                               
                               #too_low_purity,
                               Source), names_to = "classifier_version", values_to="class") |> 
  dplyr::mutate(col = as.factor(dplyr::case_when(
    class %in% c("A_IDH", "A_IDH_LG") ~ "A_IDH [_LG]",
    class == "A_IDH_HG" ~ class,
    
    class %in% c("O_IDH") ~ "O_IDH",
    class %in% c("OLIGOSARC_IDH") ~ "OLIGOSARC_IDH",
    
    class %in% c("Grade 2", "Grade 3") ~ class,
    
    class %in% c("Source: FFPE", "Source: fresh") ~ class,
    
    class %in% c("Yes", "No") ~ class,
    
    is.na(class) ~ "N/A",
    
    T ~ "Other"
  ))) |> 
  dplyr::mutate(patient_suspected_noncodel = dplyr::case_when(
    patient_suspected_noncodel == T ~ "Non-codel",
    patient_suspected_noncodel == F ~ "Codel",
    T ~ as.character("-")
  )) |> 
  dplyr::mutate(classifier_version_txt = dplyr::case_when(
    classifier_version == "resection_tumor_grade" ~ "WHO grade",
    classifier_version == "resection_treatment_status_radio" ~ "Radio",
    classifier_version == "resection_treatment_status_chemo" ~ "Chemo",
    classifier_version == "too_low_purity" ~ "Insuff. purity",
    T ~ gsub("^array_mnp_predictBrain_(.+)_cal_class$","\\1",classifier_version)
  )) |> 
  dplyr::mutate(classifier_version_txt = factor(classifier_version_txt, levels=c("WHO grade","Radio", "Chemo", "v2.0.1", "v12.5", "v12.8", "Insuff. purity", "Source")))


plt.resection.counts <- plt |>  dplyr::filter(classifier_version == "array_mnp_predictBrain_v12.8_cal_class") |>
  dplyr::pull("patient_suspected_noncodel") |>
  table()
plt.patient.counts <- plt |>
  dplyr::filter(classifier_version == "array_mnp_predictBrain_v12.8_cal_class") |>
  dplyr::select(patient_id, patient_suspected_noncodel) |> 
  dplyr::distinct() |> 
  dplyr::pull("patient_suspected_noncodel") |>
  table()


plt <- plt |>
  dplyr::mutate(patient_suspected_noncodel = paste0(patient_suspected_noncodel, "\n",plt.resection.counts[plt$patient_suspected_noncodel], " resections\n",plt.patient.counts[plt$patient_suspected_noncodel], " patients" ) )


cols = c('A_IDH [_LG]' = 'lightblue',
         'A_IDH_HG' = 'darkblue',
         
         'O_IDH' = 'lightgreen',
         'OLIGOSARC_IDH' = 'darkgreen',
         
         'Grade 2' = mixcol( 'lightblue', 'lightgreen'),
         'Grade 3' = mixcol( 'darkblue', 'darkgreen'),
         
         'Yes' = c('#E49E27', '#59B2E6', '#009E74', '#CB75A4')[1],
         'No' = c('#E49E27', '#59B2E6', '#009E74', '#CB75A4')[2],
         
         'Source: FFPE' = mixcol('purple', "white", 0.2),
         'Source: fresh' = mixcol(mixcol('red','pink',0.3), "white", 0.2),
         
         'Other' = '#FFAC1C',
         'N/A' = 'darkgray'
)


ggplot(plt, aes(x = patient_id, y = resection_number, col=col)) +
  facet_grid(cols=vars(patient_suspected_noncodel),  rows=vars(classifier_version_txt), scales = "free", space="free") +
  geom_point(pch=15,size=1.4,alpha=0.65) +
  labs(x = "patient", y="resection #", col="",fill="") +
  labs(subtitle=format_subtitle("Cohort overview")) +
  ylim(0.5, 5.5) +
  scale_color_manual(values=cols) +
  theme_nature +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=theme_nature_lwd , linetype="solid"))


ggsave("output/figures/vis_cohort_overview_MNP_classes.pdf", width=8.5 * 0.975, height = 4.4)



# Figure S5A: GLASS_OD DMP selection ----

glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(patient_id %in% c("0018")) |>
  dplyr::select(patient_id, resection_id, resection_number, isolation_id, resection_tumor_grade) |> 
  
  dplyr::mutate(LG_HG_status = dplyr::case_when(
    resection_tumor_grade %in% c(2) ~ "LG",
    resection_tumor_grade %in% c(3) ~ "HG",
    T ~ as.character(NA),
  )) |> 
  
  dplyr::group_by(patient_id, LG_HG_status) |> 
  dplyr::mutate(first_LG = 
                  !is.na(LG_HG_status) & 
                  LG_HG_status == "LG" &
                  resection_number == min(resection_number)) |>
  dplyr::mutate(last_HG = 
                  !is.na(LG_HG_status) & 
                  LG_HG_status == "HG" &
                  resection_number == max(resection_number)) |>
  dplyr::ungroup() 



plt.selection.grade <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_first_G2_and_last_G3(156) |> 
  dplyr::mutate(DMP_grade_class = ifelse(resection_tumor_grade == 2, "first Grade 2", "last Grade 3")) |> 
  dplyr::select(array_sentrix_id, DMP_grade_class)


plt.selection.resection <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(180) |> 
  dplyr::mutate(DMP_resection_class = ifelse(resection_number == 1, "primary", "last recurrent")) |> 
  dplyr::select(array_sentrix_id, DMP_resection_class)


plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES, exclude.suspected.noncodels = T) |> 
  
  dplyr::left_join(plt.selection.grade, by=c('array_sentrix_id'='array_sentrix_id')) |> 
  dplyr::left_join(plt.selection.resection, by=c('array_sentrix_id'='array_sentrix_id')) |> 
  
  dplyr::select(patient_id, resection_number, resection_tumor_grade,
                DMP_grade_class, DMP_resection_class) |> 
  tidyr::pivot_longer(cols = c(DMP_grade_class,
                               DMP_resection_class
                               ), names_to = "classifier_version", values_to="class") |> 
  dplyr::mutate(col = as.factor(dplyr::case_when(
    class %in% c("first Grade 2", "last Grade 3") ~ class,
    
    class %in% c("primary", "last recurrent") ~ class,
    
    T ~ "Other"
  ))) |> 
  dplyr::mutate(classifier_version_txt = dplyr::case_when(
    classifier_version == "DMP_grade_class" ~ "WHO Grade\n(Pathology)",
    classifier_version == "DMP_resection_class" ~ "Resection\n(Time)",
    T ~ gsub("^array_mnp_predictBrain_(.+)_cal_class$","\\1",classifier_version)
  )) |> 
  dplyr::mutate(classifier_version_txt = factor(classifier_version_txt, levels=c("WHO Grade\n(Pathology)", "Resection\n(Time)")))




cols = c('primary' = 'lightblue',
         'last recurrent' = 'darkblue',
         
         'O_IDH' = 'lightgreen',
         'OLIGOSARC_IDH' = 'darkgreen',
         
         'first Grade 2' = mixcol( 'lightblue', 'lightgreen'),
         'last Grade 3' = mixcol( 'darkblue', 'darkgreen'),
         
         'Yes' = c('#E49E27', '#59B2E6', '#009E74', '#CB75A4')[1],
         'No' = c('#E49E27', '#59B2E6', '#009E74', '#CB75A4')[2],
         
         'Source: FFPE' = mixcol('purple', "white", 0.2),
         'Source: fresh' = mixcol(mixcol('red','pink',0.3), "white", 0.2),
         
         'Other' = 'darkgray',
         'N/A' = 'darkgray'
)


ggplot(plt, aes(x = patient_id, y = resection_number, col=col)) +
  facet_grid(rows=vars(classifier_version_txt), scales = "free", space="free") +
  geom_point(pch=15,size=1.4,alpha=0.65) +
  labs(x = "patient", y="resection #", col="",fill="") +
  labs(subtitle=format_subtitle("DMP sample selection")) +
  ylim(0.5, 5.5) +
  scale_color_manual(values=cols) +
  theme_nature +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=theme_nature_lwd , linetype="solid"))


ggsave("output/figures/vis_cohort_overview_MNP_classes__DMP_sample_selection.pdf", width=8.5 * 0.975 , height = 1.65)




# Figure 1B: overview OD-validation ----


plt <- glass_od.metadata.array_samples |> 
  dplyr::filter(patient_study_name == "OD-validation" & arraychip_version == "EPICv1") |>  # not only include those used, all available to track back to literature
  dplyr::mutate(resection_number = ifelse(resection_number == 0, 4, resection_number)) |> 
  dplyr::mutate(Source = dplyr::recode(isolation_material, `ffpe`="Source: FFPE", `tissue`="Source: fresh")) |> 
  dplyr::mutate(resection_treatment_status_chemo = dplyr::case_when(
    patient_suspected_noncodel ~ as.character(NA),
    as.character(resection_treatment_status_chemo) == "TRUE" ~ "Yes",
    as.character(resection_treatment_status_chemo) == "FALSE" ~ "No",
    T ~ as.character(NA))
    ) |> 
  dplyr::mutate(resection_treatment_status_radio = dplyr::case_when(
    patient_suspected_noncodel ~ as.character(NA),
    as.character(resection_treatment_status_radio) == "TRUE" ~ "Yes",
    as.character(resection_treatment_status_radio) == "FALSE" ~ "No",
    T ~ as.character(NA))
  ) |> 
  dplyr::mutate(insufficient_quality = dplyr::case_when(
    patient_suspected_noncodel ~ as.character(NA),
    array_qc.pca.detP.outlier == T ~ "Yes",
    array_qc.pca.detP.outlier == F ~ "No",
    T ~ as.character(NA)
  )) |> 
  dplyr::mutate(insufficient_purity = dplyr::case_when(
    patient_suspected_noncodel ~ as.character(NA),
    array_methylation_bins_1p19q_purity <= 0.1 ~ "Yes",
    array_methylation_bins_1p19q_purity > 0.1 ~ "No",
    T ~ as.character(NA)
  )) |> 
  dplyr::select(patient_id, patient_suspected_noncodel, resection_number, resection_tumor_grade,
                resection_treatment_status_chemo, resection_treatment_status_radio,
                
                array_notes,
                #array_qc.pca.detP.outlier,
                
                insufficient_quality,
                insufficient_purity,
                Source, contains("_cal_class")) |> 
  dplyr::mutate(resection_tumor_grade = dplyr::case_when(
    is.na(resection_tumor_grade) | patient_suspected_noncodel ~ as.character(NA),
    !patient_suspected_noncodel ~ paste0("Grade ", resection_tumor_grade)
  )) |> 
  dplyr::rename(dataset_source = array_notes) |> 
  dplyr::mutate(dataset_source = dplyr::case_when(
    patient_id %in% c("0106","0128","0129", "0130", "0132", "0133") ~ "internal",
    dataset_source == "from Oligosarcoma manuscript" ~ "Oligosarcoma study",
    dataset_source == "from GLASS-Methylome" ~ "GLASS-methylome",
    dataset_source == "PMID: 35998208 untreated data" ~ "PMID:35998208 naive",
    dataset_source == "PMID: 35998208 treated data" ~ "PMID:35998208 trt",
    dataset_source == "from Valor LGG x GBM" ~ "GSE147391 Valor",
    T ~ dataset_source
  )) |> 
  dplyr::mutate(dataset_source_pivot = dataset_source) |> 
  tidyr::pivot_longer(cols = c(resection_tumor_grade,
                               
                               array_mnp_predictBrain_v2.0.1_cal_class,
                               array_mnp_predictBrain_v12.5_cal_class,
                               array_mnp_predictBrain_v12.8_cal_class,
                               
                               resection_treatment_status_chemo,
                               resection_treatment_status_radio,
                               
                               dataset_source_pivot,
                               
                               insufficient_quality,
                               insufficient_purity,
                               
                               Source), names_to = "classifier_version", values_to="class") |> 
  dplyr::mutate(col = as.factor(dplyr::case_when(
    class %in% c("A_IDH", "A_IDH_LG") ~ "A_IDH [_LG]",
    class == "A_IDH_HG" ~ class,
    
    class %in% c("O_IDH") ~ "O_IDH",
    class %in% c("OLIGOSARC_IDH") ~ "OLIGOSARC_IDH",
    
    class %in% c("Grade 2", "Grade 3") ~ class,
    
    class %in% c("Source: FFPE", "Source: fresh") ~ class,
    
    class %in% c("Yes", "No") ~ class,
    
    class %in% c("Oligosarcoma study","internal", "GLASS-methylome", 'GSE147391 Valor',
                 "PMID:35998208 naive", "PMID:35998208 trt") ~ class,
    
    is.na(class) ~ "N/A",
    
    T ~ "Other"
  ))) |> 
  dplyr::mutate(patient_suspected_noncodel = dplyr::case_when(
    patient_suspected_noncodel == T ~ "Non-codel",
    patient_suspected_noncodel == F ~ "Codel",
    T ~ as.character("-")
  )) |> 
  dplyr::mutate(classifier_version_txt = dplyr::case_when(
    classifier_version == "resection_tumor_grade" ~ "WHO grade",
    classifier_version == "resection_treatment_status_radio" ~ "Radio",
    classifier_version == "resection_treatment_status_chemo" ~ "Chemo",
    classifier_version == "insufficient_quality" ~ "Insuff. quality",
    classifier_version == "insufficient_purity" ~ "Insuff. purity",
    classifier_version == "dataset_source_pivot" ~ "Dataset",
    T ~ gsub("^array_mnp_predictBrain_(.+)_cal_class$","\\1",classifier_version)
  )) |> 
  dplyr::mutate(classifier_version_txt = factor(classifier_version_txt, levels=c("WHO grade","Radio", "Chemo", "v2.0.1", "v12.5", "v12.8", "Dataset", "Insuff. quality", "Insuff. purity", "Source")))


plt.resection.counts <- plt |>  dplyr::filter(classifier_version == "array_mnp_predictBrain_v12.8_cal_class") |>
  dplyr::pull("patient_suspected_noncodel") |>
  table()
plt.patient.counts <- plt |>
  dplyr::filter(classifier_version == "array_mnp_predictBrain_v12.8_cal_class") |>
  dplyr::select(patient_id, patient_suspected_noncodel) |> 
  dplyr::distinct() |> 
  dplyr::pull("patient_suspected_noncodel") |>
  table()


plt <- plt |>
  dplyr::mutate(patient_suspected_noncodel = paste0(patient_suspected_noncodel, "\n", plt.resection.counts[plt$patient_suspected_noncodel], " resections\n",plt.patient.counts[plt$patient_suspected_noncodel], " patients" ) ) |>
  dplyr::arrange(dplyr::desc(dataset_source), patient_id) |> 
  dplyr::mutate(x = factor(patient_id, levels=unique(patient_id))) |> 
  dplyr::mutate()


cols <- c('A_IDH [_LG]' = 'lightblue',
         'A_IDH_HG' = 'darkblue',
         
         'O_IDH' = 'lightgreen',
         'OLIGOSARC_IDH' = 'darkgreen',
         
         'Grade 2' = mixcol( 'lightblue', 'lightgreen'),
         'Grade 3' = mixcol( 'darkblue', 'darkgreen'),
         
         'Yes' = c('#E49E27', '#59B2E6', '#009E74', '#CB75A4')[1],
         'No' = c('#E49E27', '#59B2E6', '#009E74', '#CB75A4')[2],
         
         'Source: FFPE' = mixcol('purple', "white", 0.2),
         'Source: fresh' = mixcol(mixcol('red','pink',0.3), "white", 0.2),
         
         'internal' = 'red',
         'GLASS-methylome' = '#686ae8', # #7a7be6
         'GSE147391 Valor' = '#58e0da', 
         'Oligosarcoma study' = '#5AE82C', #  #A9E82C
         'PMID:35998208 naive' = '#f5ce42', # #ffdb58
         'PMID:35998208 trt' = 'orange',
         
         'Other' = '#FFAC1C', # 
         'N/A' = 'darkgray'
)


ggplot(plt, aes(x = x, y = resection_number, col=col)) +
  facet_grid(cols=vars(patient_suspected_noncodel),  rows=vars(classifier_version_txt), scales = "free", space="free") +
  geom_point(pch=15,size=1.4,alpha=0.65) +
  geom_point(data=plt |> dplyr::filter(classifier_version_txt == "Insuff. purity" & col == "Yes"),pch=0,size=1.4 * 1.3, col="red") +
  labs(x = "patient", y="resection #", col="",fill="") +
  labs(subtitle=format_subtitle("Cohort overview")) +
  #ylim(0.5, 4.5) +
  scale_color_manual(values=cols) +
  scale_y_continuous(limits = c(0.5, 4.5),breaks = c(1,2,3,4), labels = c(1,2,3,"N/A")) +
  theme_nature +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=theme_nature_lwd , linetype="solid"))



ggsave("output/figures/vis_validation_cohort_overview_MNP_classes.pdf", width=8.5 * 0.975, height = 5.75)



# misclassifications more common further away from Dx & R1 -
## no clear timing effect, but clear recurrence effect


# Figure 5F [var1] ----

plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(!is.na(time_between_resection_and_array)) |> 
  dplyr::mutate(time_between_resection_and_array = as.numeric(time_between_resection_and_array))


plt.exp1 <- rbind(
  plt |> 
    dplyr::mutate(panel = "det-P") |> 
    dplyr::mutate(y = array_percentage.detP.signi)
  ,
  plt |>
    dplyr::mutate(panel = "det-P") |> 
    dplyr::mutate(y = 0),
  
  plt |> 
    dplyr::mutate(panel = "time") |> 
    dplyr::mutate(y = time_between_resection_and_array / 365.25)
  ,
  plt |>
    dplyr::mutate(panel = "time") |> 
    dplyr::mutate(y = 0),
  plt |> 
    dplyr::mutate(panel = "PC") |> 
    dplyr::mutate(y = array_PC1)
  ,
  plt |>
    dplyr::mutate(panel = "PC") |> 
    dplyr::mutate(y = 0)
) |> 
  dplyr::mutate(isolation_material = ifelse(is.na(isolation_material), "N/A", isolation_material))


p1 <- ggplot(plt.exp1, aes(x = reorder(array_sentrix_id, array_PC1),
                           y = y,
                           col = isolation_material)) +
  facet_grid(rows=vars(panel), cols=vars(isolation_material), scales = "free", space="free_x") +
  geom_line(lwd=1.05) +
  scale_color_manual(values=c('ffpe'=mixcol('purple', "white", 0.2),
                              'tissue' = mixcol(mixcol('red','pink',0.3), "white", 0.2),
                              'N/A'= 'darkgray')) +
  labs(x="Sample", y=NULL) +
  theme_nature +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())





cols <- c('Grade 2' = mixcol( 'lightblue', 'lightgreen'),
          'Grade 3' = mixcol( 'darkblue', 'darkgreen'),
          
          'Yes' = c('#E49E27', '#59B2E6', '#009E74', '#CB75A4')[1],
          'No' = c('#E49E27', '#59B2E6', '#009E74', '#CB75A4')[2],
          
          'Source: FFPE' = mixcol('purple', "white", 0.2),
          'Source: fresh' = mixcol(mixcol('red','pink',0.3), "white", 0.2),
          
          'internal' = 'red',
          'GLASS-methylome' = '#686ae8', # #7a7be6
          'GSE147391 Valor' = '#58e0da', 
          'Oligosarcoma study' = '#5AE82C', #  #A9E82C
          'PMID:35998208 naive' = '#f5ce42', # #ffdb58
          'PMID:35998208 trt' = 'orange',
          
          'Other' = '#FFAC1C', # 
          'N/A' = 'darkgray'
)




plt.exp2 <- rbind(
  plt |> 
    dplyr::mutate(value = ifelse(resection_treatment_status_chemo, "Yes", "No")) |> 
    dplyr::mutate(y = "chemo")
  ,
  plt |> 
    dplyr::mutate(value = ifelse(resection_treatment_status_radio, "Yes", "No")) |> 
    dplyr::mutate(y = "radio")
  ,
  plt |> 
    dplyr::mutate(value = paste0("Grade ",resection_tumor_grade)) |> 
    dplyr::mutate(y = "grade")
  ,
  plt |> 
    dplyr::mutate(value = NA) |> 
    dplyr::mutate(y = "HM")
) |> 
  dplyr::mutate(isolation_material = ifelse(is.na(isolation_material), "N/A", isolation_material)) |> 
  dplyr::mutate(value = ifelse(is.na(value), "N/A", value))


p2 <- ggplot(plt.exp2, aes(x = reorder(array_sentrix_id, array_PC1),
                           y = y,
                           col = value)) +
  facet_grid(cols=vars(isolation_material), scales = "free", space="free_x") +
  geom_point(pch=15,size=theme_nature_size / 3, alpha=0.65) +
  scale_color_manual(values=cols) +
  labs(x=NULL, y=NULL) +
  theme_nature +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


p2 / p1 + plot_layout(heights = c(1, 8))


ggsave("output/figures/vis_cohort_overview__overview_FFPE-time.pdf", width=8.5 * 0.975, height = 2.5)


rm(cols, p1, p2, plt, plt.exp1, plt.exp2)




# Figure 5F [var2] ----


plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(!is.na(time_between_resection_and_array)) |> 
  dplyr::mutate(time_between_resection_and_array = as.numeric(time_between_resection_and_array))


plt.exp1 <- rbind(
  plt |> 
    dplyr::mutate(panel = "det-P") |> 
    dplyr::mutate(y = array_percentage.detP.signi)
  ,
  plt |>
    dplyr::mutate(panel = "det-P") |> 
    dplyr::mutate(y = 0),
  
  plt |> 
    dplyr::mutate(panel = "time") |> 
    dplyr::mutate(y = time_between_resection_and_array / 365.25)
  ,
  plt |>
    dplyr::mutate(panel = "time") |> 
    dplyr::mutate(y = 0),
  plt |> 
    dplyr::mutate(panel = "PC") |> 
    dplyr::mutate(y = array_PC1)
  ,
  plt |>
    dplyr::mutate(panel = "PC") |> 
    dplyr::mutate(y = 0)
) |> 
  dplyr::mutate(isolation_material = ifelse(is.na(isolation_material), "N/A", isolation_material)) |> 
  dplyr::mutate(isolation_material = factor(isolation_material, levels=c("tissue", "ffpe", "N/A")))


ggplot(plt.exp1, aes(x = reorder(array_sentrix_id, array_PC1),
                           y = y,
                           col = isolation_material)) +
  facet_grid(rows=vars(panel), cols=vars(isolation_material), scales = "free", space="free_x") +
  geom_line(lwd=theme_nature_lwd * 2.4) +
  scale_color_manual(values=c('ffpe'=mixcol('purple', "white", 0.2),
                              'tissue' = mixcol(mixcol('red','pink',0.3), "white", 0.2),
                              'N/A'= 'darkgray')) +
  labs(x="Sample", y=NULL) +
  theme_nature +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



ggsave("output/figures/vis_cohort_overview__overview_FFPE-time.pdf", width=8.5 * 0.975 / 2 , height = 2.5)


rm(cols, p1, p2, plt, plt.exp1, plt.exp2)




# sandbox ----
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
  theme_nature


## A_IDH_ratio ----



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


