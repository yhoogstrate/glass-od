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


if(!exists('gsam.metadata.array_samples')) {
  source('scripts/load_G-SAM_metadata.R')
}


# GLASS-OD / OD ----


metadata <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::mutate(chemo = dplyr::case_when(
    is.na(resection_treatment_status_chemo) ~ as.character(NA),
    grepl("TMZ",resection_treatment_status_summary) & grepl("PCV", resection_treatment_status_summary) ~ "TMZ & PCV",
    grepl("TMZ",resection_treatment_status_summary) ~ "TMZ",
    grepl("PCV",resection_treatment_status_summary) ~ "PCV",
    resection_treatment_status_chemo ~ "other",
    T ~ "No"
  )) |> 
  dplyr::mutate(radio = dplyr::case_when(
    resection_treatment_status_radio == T ~ "Yes",
    resection_treatment_status_radio == F ~ "No",
    T ~ as.character(NA)
  )) |> 
  dplyr::mutate(CDKN2AB = `CDKN2A/B deletion status last recurrence`)




metadata.validation <- glass_od.metadata.array_samples |> 
  filter_OD_validation_idats(CONST_N_OD_VALIDATION_INCLUDED_SAMPLES) |> 
  dplyr::filter(array_methylation_bins_1p19q_purity > 0.1) |> # discard those with near flat CNV profile
  filter_OD_validation_idats(CONST_N_OD_VALIDATION_INCLUDED_SAMPLES) |> 
  dplyr::mutate(resection_number = ifelse(resection_number == 0, as.numeric(NA), resection_number))



## x-check batch effect isolation(s) ----


plt <- metadata |> 
  dplyr::mutate(col = ifelse(isolation_person_name == "USA / Duke", "Batch [US]", "Batch [EU]")) |> 
  dplyr::mutate(col2 = dplyr::case_when(
    isolation_id == "0100-R2-repA" ~ "US: DNA as is",
    isolation_id == "0100-R2-repB" ~ "US: DNA over column",
    
    isolation_id == "0101-R2-repA" ~ "US: DNA as is",
    isolation_id == "0101-R2-repB" ~ "US: DNA over column",
    
    isolation_id == "0102-R2-repA" ~ "US: DNA as is",
    isolation_id == "0102-R2-repB" ~ "US: DNA over column",
    
    isolation_id == "0104-R1-repA" ~ "US: DNA as is",
    isolation_id == "0104-R1-repB" ~ "US: DNA over column",
    
    resection_id == "101-R3" ~ "US: FFPE",
    resection_id == "102-R3" ~ "US: FFPE",
    
    resection_id == "0016-R2" ~ "EU: suspected frozen",
    
    T ~ col
  )) |> 
  dplyr::mutate(group = dplyr::case_when(
    isolation_id == "0100-R2-repA" ~ "0100-R2",
    isolation_id == "0100-R2-repB" ~ "0100-R2",
    
    isolation_id == "0101-R2-repA" ~ "0101-R2",
    isolation_id == "0101-R2-repB" ~ "0101-R2",
    
    isolation_id == "0102-R2-repA" ~ "0102-R2",
    isolation_id == "0102-R2-repB" ~ "0102-R2",
    
    isolation_id == "0104-R1-repA" ~ "0104-R1",
    isolation_id == "0104-R1-repB" ~ "0104-R1",
    
    T ~ resection_id
  ))




ggplot(plt, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=-array_PC2, col=col2, group=group)) +
  geom_line(lwd = 0.3, col="black") +
  geom_point(data = subset(plt, col2 %in% c("Batch [EU]", "Batch [US]")), alpha=0.26) +
  geom_point(data = subset(plt, col2 %in% c("Batch [EU]", "Batch [US]") == F)) +
  theme_nature





## Figure 2 ----

### new GLASS-NL ----


plt.split <- rbind(
  metadata |> 
    dplyr::mutate(col = dplyr::case_when(
      is.na(isolation_material) ~ "NA",
      isolation_material == "ffpe" ~ "FFPE",
      isolation_material == "tissue" ~ "Tissue"
    )) |> 
    dplyr::mutate(facet = "a. Material")  |> 
    dplyr::mutate(stat = col)
  ,
  metadata |> 
    dplyr::mutate(col = ifelse(resection_number == 1, "primary", "recurrent")) |> 
    dplyr::mutate(facet = "b. Resection")  |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata |>
    dplyr::mutate(col = as.factor(paste0("Grade ",resection_tumor_grade))) |> 
    dplyr::mutate(facet = "c. Histological grade") |> 
    dplyr::mutate(stat = 'NA')
  ,
  # array_mnp_CNVP_v12.8_v5.2_CNVP_logit_chr4hd
  #metadata |> 
  #  dplyr::mutate(col = ifelse(array_mnp_predictBrain_v2.0.1_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", array_mnp_predictBrain_v2.0.1_cal_class)) |> 
  #  dplyr::mutate(facet = "d. MNP CNS Classifier 114b/2.0.1")  |> 
  #  dplyr::mutate(stat = 'NA')
  metadata |> 
    dplyr::mutate(col = array_mnp_CNVP_v12.8_v5.2_CNVP_logit_chr4hd) |> 
    dplyr::mutate(facet = "d. loss of chr4")  |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata |> 
    dplyr::mutate(col = ifelse(array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", array_mnp_predictBrain_v12.8_cal_class)) |> 
    dplyr::mutate(facet = "e. MNP CNS Classifier 12.8")  |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata |> 
    dplyr::mutate(col = radio) |> 
    dplyr::mutate(facet = "f. Radio")  |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata |> 
    dplyr::mutate(col = chemo) |> 
    dplyr::mutate(facet = "g. Chemo")  |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata |> 
    dplyr::mutate(col = dplyr::case_when(
      CDKN2AB == "wt" ~ "CDKN2A/B wt",
      CDKN2AB == "SD, arm" ~ "CDKN2A/B SD arm/chr",
      CDKN2AB == "SD, chr" ~ "CDKN2A/B SD arm/chr",
      CDKN2AB == "SD, focal" ~ "CDKN2A/B SD focal",
      CDKN2AB == "HD" ~ "CDKN2A/B HD",
      T ~ "NA"
    )) |>
    dplyr::mutate(facet = "h. CDKN2A/B")  |> 
    dplyr::mutate(stat = 'NA')
)



# array_PC2 ~ 202 | array_GLASS_NL_g2_g3_sig ~ 218
ggplot(plt.split, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=array_GLASS_NL_g2_g3_sig, col=col)) + 
  facet_wrap(~facet, scales="free",ncol=5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label), col=stat), 
                   size=theme_nature_size,
                   family=theme_nature_font_family ,cor.coef.name ="rho", show.legend = FALSE, label.y.npc = "bottom") +
  geom_point(size=theme_nature_size/3) +
  scale_color_manual(values=c(
    "Grade 2"= col3(11)[10], "O_IDH"=col3(11)[10], "primary"=col3(11)[10], 
    "Grade 3"="red", "A_IDH_HG"="red", "recurrent"="red",
    
    "No"= col3(11)[10], "O_IDH"=col3(11)[10], "primary"=col3(11)[10], 
    "Yes"="red", "A_IDH_HG"="red", "recurrent"="red",
    
    "TMZ"="red",
    "PCV"="darkgreen",
    "TMZ & PCV" = "orange",
    
    "FFPE"="darkgreen",
    "OLIGOSARC_IDH" = "orange",
    "NA" = "gray80",
    "other" = "purple",
    "Tissue"="brown",
    
    "CDKN2A/B wt" = "darkgreen",
    "CDKN2A/B SD arm/chr" = "orange",
    "CDKN2A/B SD focal" = "red",
    "CDKN2A/B HD" = "brown"
    
  )) +
  #coord_equal() +
  labs(x="Astrocytoma CGC Lasso",y = "GLASS-NL signature (supervised primary - recurrent)", col="") +
  labs(subtitle=format_subtitle("Logit WHO grade")) +
  theme_nature


ggsave("output/figures/vis_LGC_x_PCA__scatter_A1.pdf",width=8.5 * 0.975, height = 4.42)




# array_PC2 ~ 202 | array_GLASS_NL_g2_g3_sig ~ 218
ggplot(plt.split, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=array_PC2, col=col)) + 
  facet_wrap(~facet, scales="free",ncol=5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label), col=stat),  size=theme_nature_size, cor.coef.name ="rho", show.legend = FALSE, label.y.npc = "bottom") +
  geom_point(size=theme_nature_size/3) +
  scale_color_manual(values=c(
    "Grade 2"= col3(11)[10], "O_IDH"=col3(11)[10], "primary"=col3(11)[10], 
    "Grade 3"="red", "A_IDH_HG"="red", "recurrent"="red",
    
    "FFPE"="darkgreen",
    "OLIGOSARC_IDH" = "orange",
    "NA" = "gray40",
    "other" = "purple",
    "Tissue"="brown"
  )) +
  labs(x="Astrocytoma CGC Lasso",y = "Unsupervised PC2", col="") +
  labs(subtitle=format_subtitle("Logit WHO grade")) +
  theme_nature


ggsave("output/figures/vis_LGC_x_PCA__scatter_A2.pdf",width=8.5 * 0.975, height = 2.72)





### new PC3 ----


#plot(metadata$array_PC3 , metadata$array_RFpurity_absolute, pch=19, cex=0.3)
#plot(metadata$array_PC3 , metadata$array_RFpurity_estimate, pch=19, cex=0.3)
#plot(-metadata$array_PC3 , metadata$array_dnaMethyAge__PCHorvathS2018, pch=19, cex=0.3)
#plot(metadata$array_RFpurity_estimate , metadata$array_dnaMethyAge__PCHorvathS2018, pch=19, cex=0.3)
#plot(metadata$array_RFpurity_absolute , metadata$array_dnaMethyAge__PCHorvathS2018, pch=19, cex=0.3)





plt.split <- rbind(
  metadata |> 
    dplyr::mutate(col = dplyr::case_when(
      is.na(isolation_material) ~ "NA",
      isolation_material == "ffpe" ~ "FFPE",
      isolation_material == "tissue" ~ "Tissue"
    )) |> 
    dplyr::mutate(facet = "a. Material")  |> 
    dplyr::mutate(stat = col)
  ,
  metadata |> 
    dplyr::mutate(col = ifelse(resection_number == 1, "primary", "recurrent")) |> 
    dplyr::mutate(facet = "b. Resection")  |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata |>
    dplyr::mutate(col = as.factor(paste0("Grade ",resection_tumor_grade))) |> 
    dplyr::mutate(facet = "c. Histological grade") |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata |> 
    dplyr::mutate(col = ifelse(array_mnp_predictBrain_v2.0.1_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", array_mnp_predictBrain_v2.0.1_cal_class)) |> 
    dplyr::mutate(facet = "d. MNP CNS Classifier 114b/2.0.1")  |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata |> 
    dplyr::mutate(col = ifelse(array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", array_mnp_predictBrain_v12.8_cal_class)) |> 
    dplyr::mutate(facet = "e. MNP CNS Classifier 12.8")  |> 
    dplyr::mutate(stat = 'NA')
  # ,
  # metadata |> 
  #   dplyr::mutate(col = radio) |> 
  #   dplyr::mutate(facet = "f. Radio")  |> 
  #   dplyr::mutate(stat = 'NA')
  # ,
  # metadata |> 
  #   dplyr::mutate(col = chemo) |> 
  #   dplyr::mutate(facet = "g. Chemo")  |> 
  #   dplyr::mutate(stat = 'NA')
)



#ggplot(plt.split, aes(x=array_dnaMethyAge__PCHorvathS2018, y=-array_PC3, col=col))
#ggplot(plt.split, aes(x=array_methylation_bins_1p19q_purity, y=-array_PC3, col=col)) + 
ggplot(plt.split, aes(x=array_dnaMethyAge__PCHorvathS2018, y=-array_PC3, col=col)) +
  facet_wrap(~facet, scales="free",ncol=5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label), col=stat),  size=theme_nature_size, cor.coef.name ="rho", show.legend = FALSE, label.y.npc = "top") +
  geom_point(size=theme_nature_size/3) +
  scale_color_manual(values=c(
    "Grade 2"= col3(11)[10], "O_IDH"=col3(11)[10], "primary"=col3(11)[10], 
    "Grade 3"="red", "A_IDH_HG"="red", "recurrent"="red",
    
    "No"= col3(11)[10], "O_IDH"=col3(11)[10], "primary"=col3(11)[10], 
    "Yes"="red", "A_IDH_HG"="red", "recurrent"="red",
    
    "TMZ"="red",
    "PCV"="darkgreen",
    "TMZ & PCV" = "orange",
    
    "FFPE"="darkgreen",
    "OLIGOSARC_IDH" = "orange",
    "NA" = "gray40",
    "other" = "purple",
    "Tissue"="brown"
  )) +
  labs(x="Horvath [S] 2018 (PC)",y = "(-) Unsupervised PC3", col="") +
  labs(subtitle=format_subtitle("Logit WHO grade")) +
  theme_nature


ggsave("output/figures/vis_LGC_x_PCA__scatter_A2_PC3.pdf",width=8.5 * 0.975, height = 2.72)





### new validationset ----


metadata.validation$patient_study_name
metadata.validation$array_notes


metadata.validation |>
  dplyr::filter(grepl("35998208 untrea", array_notes)) |> 
  dplyr::filter(resection_tumor_grade == 3) |> 
  dplyr::pull(resection_id)


# metadata.validation <- metadata.validation |> 
#   dplyr::filter(grepl("35998208", array_notes) )


plt.split.raw <- rbind(
  metadata.validation |> 
    dplyr::mutate(col = dplyr::case_when(
      patient_id %in% c("0106","0128","0129", "0130", "0132", "0133") ~ "internal",
      array_notes == "from Oligosarcoma manuscript" ~ "Oligosarcoma study",
      array_notes == "from GLASS-Methylome" ~ "GLASS-methylome",
      array_notes == "PMID: 35998208 untreated data" ~ "PMID:35998208 naive",
      array_notes == "PMID: 35998208 treated data" ~ "PMID:35998208 trt",
      array_notes == "from Valor LGG x GBM" ~ "GSE147391 Valor",
      T ~ array_notes
    )) |> 
    #dplyr::mutate(col = ifelse(col == "PMID:35998208 naive" & resection_tumor_grade == 3 & grepl("1378",resection_id), "PMID:35998208 naive" ,"other")) |> 
    #dplyr::mutate(col = ifelse(grepl("35998208", col), "PMID:35998208 naive" ,"other")) |> 
    dplyr::mutate(facet = "A) Dataset")  |> 
    dplyr::mutate(stat = col)
  ,
  metadata.validation |> 
    dplyr::mutate(col = ifelse(resection_number == 1, "primary", "recurrent")) |> 
    dplyr::mutate(facet = "B) Resection")  |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata.validation |>
    dplyr::mutate(col = as.factor(paste0("Grade ",resection_tumor_grade))) |> 
    dplyr::mutate(facet = "C) Histological grade") |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata.validation |> 
    dplyr::mutate(col = ifelse(array_mnp_predictBrain_v2.0.1_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", array_mnp_predictBrain_v2.0.1_cal_class)) |> 
    dplyr::mutate(facet = "D) MNP CNS Classifier 114b/2.0.1")  |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata.validation |> 
    dplyr::mutate(col = ifelse(array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", array_mnp_predictBrain_v12.8_cal_class)) |> 
    dplyr::mutate(facet = "E) MNP CNS Classifier 12.8")  |> 
    dplyr::mutate(stat = 'NA')
) 


#ids <- plt.split$array_sentrix_id

plt.split <- plt.split.raw |> 
  dplyr::filter(array_sentrix_id %in% ids == T)


# array_PC2
# array_PC2 ~ 202 | array_GLASS_NL_g2_g3_sig ~ 218
ggplot(plt.split, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=array_GLASS_NL_g2_g3_sig, col=col)) + 
  facet_wrap(~facet, scales="free",ncol=5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label), col=stat),  size=theme_nature_size, cor.coef.name ="rho", show.legend = FALSE, label.y.npc = "bottom") +
  geom_point(size=theme_nature_size/3) +
  scale_color_manual(values=c(
    "Grade 2"= col3(11)[10], "O_IDH"=col3(11)[10], "primary"=col3(11)[10], 
    "Grade 3"="red", "A_IDH_HG"="red", "recurrent"="red",
    
    "FFPE"="darkgreen",
    "OLIGOSARC_IDH" = "orange",
    "NA" = "gray40",
    "other" = "purple",
    "Tissue"="brown",
    
    'internal' = 'red',
    'GLASS-methylome' = '#686ae8', # #7a7be6
    'GSE147391 Valor' = '#58e0da', 
    'Oligosarcoma study' = '#5AE82C', #  #A9E82C
    'PMID:35998208 naive' = 'red', # 'f5ce42', # #ffdb58
    'PMID:35998208 trt' = 'darkgreen' #'orange'
    
  )) +
  labs(x="Astrocytoma CGC Lasso",y = "CGC[Ac]", col="") +
  labs(subtitle=format_subtitle("Logit WHO grade")) +
  theme_nature


ggsave("output/figures/vis_LGC_x_PCA__scatter_A1__validationset__mair.pdf",width=8.5 * 0.975, height = 2.72)
#ggsave("output/figures/vis_LGC_x_PCA__scatter_A1__validationset.pdf",width=8.5 * 0.975, height = 2.72)


### pick samples mair ----



plt.split.raw <- rbind(
  metadata.validation |> 
    dplyr::mutate(col = dplyr::case_when(
      patient_id %in% c("0106","0128","0129", "0130", "0132", "0133") ~ "internal",
      array_notes == "from Oligosarcoma manuscript" ~ "Oligosarcoma study",
      array_notes == "from GLASS-Methylome" ~ "GLASS-methylome",
      array_notes == "PMID: 35998208 untreated data" ~ "PMID:35998208 naive",
      array_notes == "PMID: 35998208 treated data" ~ "PMID:35998208 trt",
      array_notes == "from Valor LGG x GBM" ~ "GSE147391 Valor",
      T ~ array_notes
    )) |> 
    #dplyr::mutate(col = ifelse(col == "PMID:35998208 naive" & resection_tumor_grade == 3 & grepl("1378",resection_id), "PMID:35998208 naive" ,"other")) |> 
    #dplyr::mutate(col = ifelse(grepl("35998208", col), "PMID:35998208 naive" ,"other")) |> 
    dplyr::mutate(facet = "A) Dataset")  |> 
    dplyr::mutate(stat = col)
  ,
  metadata.validation |>
    dplyr::mutate(col = as.factor(paste0("Grade ",resection_tumor_grade))) |> 
    dplyr::mutate(facet = "C) Histological grade") |> 
    dplyr::mutate(stat = 'NA')
) 



plt.split <- plt.split.raw |> 
  dplyr::filter(array_sentrix_id %in% (
    metadata.validation |> dplyr::filter(grepl("35998208", array_notes)) |> dplyr::pull(array_sentrix_id)
  )) |> 
  dplyr::mutate(label = gsub("-R1","", resection_id)) |> 
  dplyr::mutate(label = gsub("PM35998208ut_","", label)) |> 
  dplyr::mutate(label = gsub("GSM6","G", label)) 


plt.split |> dplyr::filter(facet == "A) Dataset") |> dplyr::pull(array_notes) |> table()

# array_PC2
# array_PC2 ~ 202 | array_GLASS_NL_g2_g3_sig ~ 218
ggplot(plt.split, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=array_GLASS_NL_g2_g3_sig, col=col, label=label)) + 
  facet_wrap(~facet, scales="free",ncol=5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label), col=stat),  size=theme_nature_size, cor.coef.name ="rho", show.legend = FALSE, label.y.npc = "bottom") +
  geom_point(size=theme_nature_size/3) +
  scale_color_manual(values=c(
    "Grade 2"= col3(11)[10], "O_IDH"=col3(11)[10], "primary"=col3(11)[10], 
    "Grade 3"="red", "A_IDH_HG"="red", "recurrent"="red",
    
    "FFPE"="darkgreen",
    "OLIGOSARC_IDH" = "orange",
    "NA" = "gray40",
    "other" = "purple",
    "Tissue"="brown",
    
    'internal' = 'red',
    'GLASS-methylome' = '#686ae8', # #7a7be6
    'GSE147391 Valor' = '#58e0da', 
    'Oligosarcoma study' = '#5AE82C', #  #A9E82C
    'PMID:35998208 naive' = 'red', # 'f5ce42', # #ffdb58
    'PMID:35998208 trt' = 'darkgreen' #'orange'
    
  )) +
  labs(x="Astrocytoma CGC Lasso",y = "CGC[Ac]", col="") +
  labs(subtitle=format_subtitle("Logit WHO grade")) +
  theme_nature +
  ggrepel::geom_text_repel()





### old ----

plt.split <- rbind(
  metadata |>
    dplyr::mutate(col = as.factor(paste0("Grade ",resection_tumor_grade))) |> 
    dplyr::mutate(facet = "Histological grade") |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata |> 
    dplyr::mutate(col = ifelse(array_mnp_predictBrain_v2.0.1_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", array_mnp_predictBrain_v2.0.1_cal_class)) |> 
    dplyr::mutate(facet = "MNP CNS Classifier 114b/2.0.1")  |> 
    dplyr::mutate(stat = 'NA')
  # ,
  # metadata |> 
  #   dplyr::mutate(col = ifelse(array_mnp_predictBrain_v12.5_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", array_mnp_predictBrain_v12.8_cal_class)) |> 
  #   dplyr::mutate(facet = "MNP CNS Classifier 12.5")  |> 
  #   dplyr::mutate(stat = 'NA')
  ,
  metadata |> 
    dplyr::mutate(col = ifelse(array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", array_mnp_predictBrain_v12.8_cal_class)) |> 
    dplyr::mutate(facet = "MNP CNS Classifier 12.8")  |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata |> 
    dplyr::mutate(col = ifelse(isolation_person_name == "USA / Duke", "Batch [US]", "Batch [EU]")) |> 
    dplyr::mutate(facet = "Batch")  |> 
    dplyr::mutate(stat = col)
)




ggplot(plt.split, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=array_PC2, col=col)) + 
  #facet_grid(cols = vars(facet), scales = "free", space="free") +
  facet_wrap(~facet, scales="free",ncol=5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label), col=stat),  size=theme_nature_size, cor.coef.name ="rho", show_guide = FALSE) +
  geom_point(size=theme_nature_size/3) +
  scale_color_manual(values=c(
    "Grade 2"= col3(11)[10],"O_IDH"=col3(11)[10], 
    "Grade 3"="red","A_IDH_HG"="red",  
    
    "Batch [EU]"="darkgreen",
    "OLIGOSARC_IDH" = "orange",
    "NA" = "gray40",
    "other" = "purple",
    "Batch [US]"="brown"
  )) +
  labs(x="AcCGAP",y = "PC2", col="") +
  labs(subtitle=format_subtitle("Logit WHO grade")) +
  theme_cellpress





# pc loadings

pc <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_prcomp.Rds")


plot(sort(pc$rotation[,1]), type="l",col="gray")
lines(sort(pc$rotation[,2]), type="l", col="blue")
lines(sort(pc$rotation[,3]), type="l", col="gray")
lines(sort(pc$rotation[,4]), type="l", col="gray")
lines(sort(pc$rotation[,5]), type="l", col="gray")
lines(sort(pc$rotation[,6]), type="l", col="gray")
lines(sort(pc$rotation[,7]), type="l", col="gray")
lines(sort(pc$rotation[,8]), type="l", col="gray")
lines(sort(pc$rotation[,9]), type="l", col="gray")
lines(sort(pc$rotation[,10]), type="l", col="gray")
abline(h=0, col="red")


n = length(pc$rotation[,1])

s1 = sort(pc$rotation[,1])
fdiv1 = s1[2:n] - s1[1:(n-1)]

s2 = sort(pc$rotation[,2])
fdiv2 = s2[2:n] - s2[1:(n-1)]

s3 = sort(pc$rotation[,3])
fdiv3 = s3[2:n] - s3[1:(n-1)]

s4 = sort(pc$rotation[,4])
fdiv4 = s4[2:n] - s4[1:(n-1)]

s5 = sort(pc$rotation[,5])
fdiv5 = s5[2:n] - s5[1:(n-1)]


plot(sort(fdiv1, decreasing=T), type="l",col="gray",ylim=c(0,0.00001),xlim=c(0,2000))
lines(sort(fdiv2, decreasing=T), type="l",col="blue")
lines(sort(fdiv3, decreasing=T), type="l",col="gray")
lines(sort(fdiv4, decreasing=T), type="l",col="gray")
lines(sort(fdiv5, decreasing=T), type="l",col="gray")



#plot(sort(-pc$rotation[,2]), type="l")
#abline(h=0, col="red")

#plot(sort(pc$rotation[,3]), type="l")
#abline(h=0, col="red")

#plot(sort(pc$rotation[,4]), type="l")
#abline(h=0, col="red")

plot(sort(pc$rotation[,5]), type="l")
abline(h=0, col="red")





## Figure 2C: logit ----
### Grade ----
#### logistic AcCGAP x WHO grade unpaired ----



stats <- metadata |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(resection_tumor_grade__hg = ifelse(resection_tumor_grade == 3, 1 , 0)) |> 
  
  dplyr::mutate(resection = NULL) |> 
  dplyr::mutate(resection_recurrent = NULL) |> 
  
  dplyr::mutate(covar = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(covar_name = "Astrocytoma CGC Lasso")


model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_tumor_grade__hg - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_tumor_grade__hg + 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_tumor_grade__hg) |> 
    dplyr::mutate(col = resection_tumor_grade__hg)
) |> 
  dplyr::mutate(col = col  + 2)



p3 <- ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data"),
            col="white",
            lwd=theme_nature_lwd * 2, alpha=0.65
  ) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH","A_IDH_HG") == F),
            col="gray",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == T),
            col="darkblue",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH") == T),
            col="#3e9e8b",
            lwd=theme_nature_lwd) +  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 3), labels=c("Grade 2", "", "Grade 3"), oob = scales::squish) +
  theme_nature +
  annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[7], x = 0.3, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability", y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = paste0(unique(stats$covar_name)," x WHO Grade")) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("Grade 2", "Grade 3")) + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1))
p3




#### logistic AcCGAP x WHO grade paired ----


stats <- metadata |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(resection_tumor_grade__hg = ifelse(resection_tumor_grade == 3, 1 , 0)) |> 
  
  dplyr::mutate(resection = NULL) |> 
  dplyr::mutate(resection_recurrent = NULL) |> 
  
  dplyr::mutate(covar = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(covar_name = "Astrocytoma CGC Lasso") |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder"))))


model <- glm(resection_tumor_grade__hg ~ patient + covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


#expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(-15, 15), 500)) |> 
  dplyr::mutate(patient = "p_remainder")
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_tumor_grade__hg - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_tumor_grade__hg + 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(patient = NULL) |> 
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_tumor_grade__hg) |> 
    dplyr::mutate(col = resection_tumor_grade__hg)
) |> 
  dplyr::mutate(col = col  + 2)



p3b <- ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data"),
            col="white",
            lwd=theme_nature_lwd * 2, alpha=0.65
  ) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH","A_IDH_HG") == F),
            col="gray",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == T),
            col="darkblue",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH") == T),
            col="#3e9e8b",
            lwd=theme_nature_lwd) +  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 3), labels=c("Grade 2", "", "Grade 3"), oob = scales::squish) +
  theme_nature +
  annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[7], x = 0.3, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability", y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = paste0(unique(stats$covar_name)," x WHO Grade")) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("Grade 2", "Grade 3")) + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1)) +
  ylim(-15, 15)

p3b


### MNP CNS 12.8 ----
#### unpaired ----


stats <- metadata |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(high_grade_classes = ifelse(array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG", "OLIGOSARC_IDH"), 1 , 0)) |> 
  
  dplyr::mutate(resection = NULL) |> 
  dplyr::mutate(resection_recurrent = NULL) |> 
  
  dplyr::mutate(covar = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(covar_name = "predictBrain_v12.8")


model <- glm(high_grade_classes ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$high_grade_classes = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(high_grade_classes, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = high_grade_classes * 1.25 - 0.125) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = high_grade_classes ) |> 
    dplyr::mutate(col = high_grade_classes)
) |> 
  dplyr::mutate(col = col  + 2)



p <- ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  
  geom_vline(aes(xintercept=0), lty=3, lwd=theme_nature_lwd)  + 
  geom_vline(aes(xintercept=1), lty=3, lwd=theme_nature_lwd)  + 
  
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  
  #  A_IDH_HG      A_IDH_LG   GBM_MES_TYP      GBM_RTK1         O_IDH OLIGOSARC_IDH 
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == T),
                               col=palette_mnp_12.8_6['A_IDH_HG'],
                               size=theme_nature_size/3,
                               width=0.1) +
  
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_LG") == T),
                               col=palette_mnp_12.8_6['A_IDH_LG'],
                               size=theme_nature_size/3,
                               width=0.1) +
  
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("O_IDH") == T),
                               col=palette_mnp_12.8_6['O_IDH'],
                               size=theme_nature_size/3,
                               width=0.1) +
  
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH") == T),
                               col=palette_mnp_12.8_6['OLIGOSARC_IDH'],
                               size=theme_nature_size/3,
                               width=0.1) +
  
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("GBM_MES_TYP","GBM_RTK1","CTRL_CORPCAL") == T),
                               col='gray40',
                               size=theme_nature_size/3,
                               width=0.1) +
  
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 3),
                        labels=c("0", "", "1"), oob = scales::squish) +
  
  theme_nature +
  annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[7], x = 0.3, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability",
       y= "CGC''", fill=NULL, x=NULL, subtitle = paste0(unique(stats$covar_name)," x MNP")) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("p=0", "p=1")) + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1))
p


### NCI Methylscape ----
#### unpaired ----

# A_IDH      A_IDH_HG CONTR_CORPCAL   CONTR_REACT     GBM_RTK_I         MELAN         O_IDH    O_SARC_IDH 

stats <- metadata |> 
  dplyr::mutate(high_grade_classes = ifelse(array_methylscape_bethesda_class %in% c("A_IDH_HG", "O_SARC_IDH"), 1 , 0)) |> 
  
  dplyr::mutate(resection = NULL) |> 
  dplyr::mutate(resection_recurrent = NULL) |> 
  
  dplyr::mutate(covar = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(covar_name = "predictBrain_v12.8")


model <- glm(high_grade_classes ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$high_grade_classes = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(high_grade_classes, covar, array_methylscape_bethesda_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = high_grade_classes * 1.25 - 0.125) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_methylscape_bethesda_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = high_grade_classes ) |> 
    dplyr::mutate(col = high_grade_classes)
) |> 
  dplyr::mutate(col = col  + 2)



p <- ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  
  geom_vline(aes(xintercept=0), lty=3, lwd=theme_nature_lwd)  + 
  geom_vline(aes(xintercept=1), lty=3, lwd=theme_nature_lwd)  + 
  
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  
  #  metadata$array_methylscape_bethesda_class |> table()
  # A_IDH      A_IDH_HG CONTR_CORPCAL   CONTR_REACT     GBM_RTK_I         MELAN         O_IDH    O_SARC_IDH 
  #2             4             2             1             1             1           197             3 
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_methylscape_bethesda_class %in% c("A_IDH_HG") == T),
                               col=palette_mnp_12.8_6['A_IDH_HG'],
                               size=theme_nature_size/3,
                               width=0.1) +
  
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_methylscape_bethesda_class %in% c("A_IDH_LG", "A_IDH") == T),
                               col=palette_mnp_12.8_6['A_IDH_LG'],
                               size=theme_nature_size/3,
                               width=0.1) +
  
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_methylscape_bethesda_class %in% c("O_IDH") == T),
                               col=palette_mnp_12.8_6['O_IDH'],
                               size=theme_nature_size/3,
                               width=0.1) +
  
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_methylscape_bethesda_class %in% c("O_SARC_IDH") == T),
                               col=palette_mnp_12.8_6['OLIGOSARC_IDH'],
                               size=theme_nature_size/3,
                               width=0.1) +
  
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_methylscape_bethesda_class %in% c("CONTR_CORPCAL","CONTR_REACT","GBM_RTK_I","MELAN") == T),
                               col='gray40',
                               size=theme_nature_size/3,
                               width=0.1) +
  
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 3),
                        labels=c("0", "", "1"), oob = scales::squish) +
  
  theme_nature +
  annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[7], x = 0.3, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability",
       y= "CGC''", fill=NULL, x=NULL, subtitle = paste0(unique(stats$covar_name)," x MNP")) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("p=0", "p=1")) + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1))
p


### Radiotherapy ----
#### unpaired ----




stats <- metadata |> 
  dplyr::filter(!is.na(resection_treatment_status_radio)) |> 
  dplyr::mutate(trt_classes = resection_treatment_status_radio * 1) |> 
  dplyr::mutate(covar = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)

n_samples = nrow(stats)

model <- glm(trt_classes ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$trt_classes = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(trt_classes, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = trt_classes * 1.25 - 0.125) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = trt_classes ) |> 
    dplyr::mutate(col = trt_classes)
) |> 
  dplyr::mutate(col = col  + 2)



p <- ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  
  geom_vline(aes(xintercept=0), lty=3, lwd=theme_nature_lwd)  + 
  geom_vline(aes(xintercept=1), lty=3, lwd=theme_nature_lwd)  + 
  
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" ),
                               col='gray40',
                               size=theme_nature_size/3,
                               width=0.1) +
  
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 3),
                        labels=c("0", "", "1"), oob = scales::squish) +
  
  theme_nature +
  annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[7], x = 0.3, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability",
       y= "CGC''", fill=NULL, x=NULL, subtitle = paste0("MNP -- samples: ",n_samples)) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("p=0", "p=1")) + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1))
p


### Chemotherapy ----
#### unpaired ----




stats <- metadata |> 
  dplyr::filter(!is.na(resection_treatment_status_chemo)) |> 
  dplyr::mutate(trt_classes = resection_treatment_status_chemo * 1) |> 
  dplyr::mutate(covar = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(chemo = dplyr::case_when(
    is.na(resection_treatment_status_chemo) ~ as.character(NA),
    grepl("TMZ",resection_treatment_status_summary) & grepl("PCV", resection_treatment_status_summary) ~ "TMZ & PCV",
    grepl("TMZ",resection_treatment_status_summary) ~ "TMZ",
    grepl("PCV",resection_treatment_status_summary) ~ "PCV",
    resection_treatment_status_chemo ~ "other",
    T ~ "No"
  ))
  

n_samples = nrow(stats)

model <- glm(trt_classes ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$trt_classes = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(trt_classes, covar, chemo) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = trt_classes * 1.25 - 0.125) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(chemo = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = trt_classes ) |> 
    dplyr::mutate(col = trt_classes)
) |> 
  dplyr::mutate(col = col  + 2)



p <- ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  
  geom_vline(aes(xintercept=0), lty=3, lwd=theme_nature_lwd)  + 
  geom_vline(aes(xintercept=1), lty=3, lwd=theme_nature_lwd)  + 
  
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & chemo == "No"),
                               col='gray40',
                               size=theme_nature_size/3,
                               width=0.1) +
  
  
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & chemo == "TMZ"),
                               col='red',
                               size=theme_nature_size/3,
                               width=0.1) +
  
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & chemo == "PCV"),
                               col=col3(11)[10],
                               size=theme_nature_size/3,
                               width=0.1) +
  
  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & chemo == "TMZ & PCV"),
                               col="purple",
                               size=theme_nature_size/3,
                               width=0.1) +

  ggbeeswarm::geom_quasirandom(data = plt.logit.simplistic |> dplyr::filter(type == "data" & chemo == "other"),
                               col="orange",
                               size=theme_nature_size/3,
                               width=0.1) +
  
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 3),
                        labels=c("0", "", "1"), oob = scales::squish) +
  
  theme_nature +
  annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[7], x = 0.3, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability",
       y= "CGC''", fill=NULL, x=NULL, subtitle = paste0("MNP -- samples: ",n_samples)) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("p=0", "p=1")) + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1))
p




### logistic GLASS-NL MM sig x WHO grade validation ----
 


stats <- metadata.validation |> 
#  dplyr::filter(!grepl("PMID: 35998208 treated", metadata.validation$array_notes)) |> # s.t.h. is odd about this data
  
  dplyr::filter(!is.na(resection_tumor_grade)) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(resection_tumor_grade__hg = ifelse(resection_tumor_grade == 3, 1, 0)) |> 
  
  dplyr::mutate(resection = NULL) |> 
  dplyr::mutate(resection_recurrent = NULL) |> 
  
  dplyr::mutate(covar = array_GLASS_NL_g2_g3_sig) |> 
  dplyr::mutate(covar_name = "GLASS-NL signature")


model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.175
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_tumor_grade__hg - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_tumor_grade__hg + 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_tumor_grade__hg) |> 
    dplyr::mutate(col = resection_tumor_grade__hg)
) |> 
  dplyr::mutate(col = col  + 2)



ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data"),
            col="white",
            lwd=theme_nature_lwd * 2, alpha=0.65
  ) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH","A_IDH_HG") == F),
            col="gray",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == T),
            col="darkblue",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH") == T),
            col="#3e9e8b",
            lwd=theme_nature_lwd) +
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 2.50, 3), labels=c("Grade 2","","", "Grade 3"), oob = scales::squish) +
  theme_nature +
annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[2], x = 0.3, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability", y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = paste0(unique(stats$covar_name)," x WHO Grade")) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("Grade 2", "Grade 3")) + 
  scale_y_reverse() + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1))




### logistic GLASS-NL MM sig x resection ----


stats <- metadata |> 
  dplyr::mutate(resection_tumor_grade = NULL) |> 
  dplyr::mutate(resection_tumor_grade__hg = NULL) |> 
  
  dplyr::mutate(resection = ifelse(resection_number == 1, "primary", "recurrent")) |> 
  dplyr::mutate(resection_recurrent = ifelse(resection == "primary", 0 , 1)) |> 

  dplyr::mutate(covar = array_GLASS_NL_g2_g3_sig) |> 
  dplyr::mutate(covar_name = "GLASS-NL signature")


model <- glm(resection_recurrent ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$resection_recurrent = predict(model, Predicted_data, type="response")


plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_recurrent, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1 ) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_recurrent - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_recurrent, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1 ) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_recurrent + 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_recurrent) |> 
    dplyr::mutate(col = resection_recurrent)
)




p2 <- ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd= theme_nature_lwd * 5) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data"),
            col="white",
            lwd=theme_nature_lwd * 2, alpha=0.65
  ) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH","A_IDH_HG") == F),
            col="gray",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == T),
            col="darkblue",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH") == T),
            col="#3e9e8b",
            lwd=theme_nature_lwd) +  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(0, 1), breaks=c(0, 0.50, 1), labels=c("Primary","", "Recurrent"), oob = scales::squish) +
  theme_nature +
annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[2], x = 0.3, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability", y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = paste0(unique(stats$covar_name)," x Resection")) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("Primary", "Recurrent")) + 
  scale_y_reverse() + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1))
p2



# complex color variant:
# plt.logit.restyled <- rbind(
#   stats |>  # left point line:
#     dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class, resection_recurrent) |> 
#     dplyr::mutate(col = resection_tumor_grade__hg + 2) |> 
#     dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
#     dplyr::mutate(type = "data") |> 
#     dplyr::mutate(x = resection_recurrent - 0.06 ) |> 
#     dplyr::rename(y = covar)
#   ,
#   stats |>  # right point line:
#     dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class, resection_recurrent) |> 
#     dplyr::mutate(col = resection_tumor_grade__hg + 2) |> 
#     dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
#     dplyr::mutate(type = "data") |> 
#     dplyr::mutate(x = resection_recurrent + 0.06 ) |> 
#     dplyr::rename(y = covar)
#   ,
#   Predicted_data |> # logit fit
#     dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
#     dplyr::mutate(resection_tumor_grade__hg = 0) |> 
#     dplyr::mutate(group = "logit fit") |> 
#     dplyr::mutate(type = "fit") |> 
#     dplyr::rename(y = covar) |> 
#     dplyr::mutate(x = resection_recurrent) |> 
#     dplyr::mutate(col = resection_recurrent)
# )
# 
# 
# p2 <- ggplot(plt.logit.restyled, aes(x=x, y=y, group=group, col=col)) +
#   geom_line(data = plt.logit.restyled |> dplyr::filter(type == "fit") ,
#             aes(col=col),
#             lwd=2) +
#   scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(0, 1), breaks=c(0, 1), labels=c("Primary", "Recurrent"), oob = scales::squish) +
#   labs(col=NULL) +
#   ggnewscale::new_scale_colour() +
#   geom_line(data = plt.logit.restyled |> dplyr::filter(type == "data"),
#             col="white",
#             lwd=theme_nature_lwd * 2, alpha=0.65
#   ) +
#   geom_line(data = plt.logit.restyled |> dplyr::filter(type == "data"),
#             aes(col=col),
#             lwd=theme_nature_lwd) +
#   scale_color_gradientn(colours = c("aquamarine3","aquamarine3","#d34394ff","#d34394ff"), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 2.50, 3), labels=c("Grade 2","","", "Grade 3"), oob = scales::squish) +
#   theme_nature +
# annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[2], x = 0.2, label = paste0("p = ",format.pval(pval)), size=theme_nature_size) +
#   labs(col=NULL, y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = "GLASS-NL MedMeth x Resection") +
#   scale_x_continuous(breaks = c(0,1),
#                      labels=c("Primary", "Recurrent")) + 
#   theme(legend.box = "vertical") + # space dependent
#   theme(legend.key.size = unit(0.6, 'lines'))
# p2




### logistic AcCGAP x WHO grade [validation] ----


stats <- metadata.validation |> 
  dplyr::filter(!is.na(resection_tumor_grade)) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(resection_tumor_grade__hg = ifelse(resection_tumor_grade == 3, 1 , 0)) |> 
  
  dplyr::mutate(resection = NULL) |> 
  dplyr::mutate(resection_recurrent = NULL) |> 
  
  dplyr::mutate(covar = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(covar_name = "Astrocytoma CGC Lasso")


model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_tumor_grade__hg - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_tumor_grade__hg + 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_tumor_grade__hg) |> 
    dplyr::mutate(col = resection_tumor_grade__hg)
) |> 
  dplyr::mutate(col = col  + 2)



ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data"),
            col="white",
            lwd=theme_nature_lwd * 2, alpha=0.65
  ) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH","A_IDH_HG") == F),
            col="gray",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == T),
            col="darkblue",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH") == T),
            col="#3e9e8b",
            lwd=theme_nature_lwd) +  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 3), labels=c("Grade 2", "", "Grade 3"), oob = scales::squish) +
  theme_nature +
  annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[7], x = 0.3, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability", y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = paste0(unique(stats$covar_name)," x WHO Grade")) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("Grade 2", "Grade 3")) + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1))




# complex color variant:
# plt.logit.restyled <- rbind(
#   stats |>  # left point line:
#     dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class, resection_recurrent) |> 
#     dplyr::mutate(col = resection_recurrent + 1) |> 
#     dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
#     dplyr::mutate(type = "data") |> 
#     dplyr::mutate(x = ifelse(resection_tumor_grade__hg == 0, resection_tumor_grade__hg, resection_tumor_grade__hg - 0.175)) |> 
#     dplyr::rename(y = covar),
#   
#   stats |>  # right point line:
#     dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class, resection_recurrent) |> 
#     dplyr::mutate(col = resection_recurrent + 1) |> 
#     dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
#     dplyr::mutate(type = "data") |> 
#     dplyr::mutate(x = ifelse(resection_tumor_grade__hg == 1, resection_tumor_grade__hg, resection_tumor_grade__hg + 0.175)) |> 
#     dplyr::rename(y = covar),
#   
#   Predicted_data |> # logit fit
#     dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
#     dplyr::mutate(resection_recurrent = 0) |> 
#     dplyr::mutate(group = "logit fit") |> 
#     dplyr::mutate(type = "fit") |> 
#     dplyr::rename(y = covar) |> 
#     dplyr::mutate(x = resection_tumor_grade__hg) |> 
#     dplyr::mutate(col = resection_tumor_grade__hg)
# )
# 
# 
# plt.logit.restyled.grade <- plt.logit.restyled
# 
# 
# p3 <- ggplot(plt.logit.restyled, aes(x=x, y=y, group=group, col=col)) +
#   geom_line(data = plt.logit.restyled |> dplyr::filter(type == "data"),
#             aes(col=col),
#             lwd=theme_nature_lwd) +
#   scale_color_gradientn(colours = c("aquamarine3","aquamarine3","#d34394ff","#d34394ff"), na.value = "grey50", limits = c(1, 2), breaks=c(1, 1.50, 1.50, 2), labels=c("primary","","", "recurrent"), oob = scales::squish) +
#   labs(col=NULL) +
#   ggnewscale::new_scale_colour() +
#   geom_line(data = plt.logit.restyled |> dplyr::filter(type == "fit") ,
#             aes(col=col),
#             lwd=1.0) +
#   theme_nature +
#   scale_fill_manual(values = palette_mnp_12.8_6) +
# annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 16)[15], x = 0.375, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
#   scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(0, 1), breaks=c(0, 1), labels=c("Grade 2", "Grade 3"), oob = scales::squish) +
#   labs(col=NULL, y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = "CGC Lasso x WHO Grade") +
#   scale_x_continuous(breaks = c(0,1), labels=c("Grade 2", "Grade 3")) + 
#   theme(legend.box = "vertical") + # space dependent
#   theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
# p3



### paired ----



#### logistic AcCGAP x resection paired ----


stats <- metadata |> 
  dplyr::mutate(resection_tumor_grade = NULL) |> 
  dplyr::mutate(resection_tumor_grade__hg = NULL) |> 
  
  dplyr::mutate(resection = ifelse(resection_number == 1, "primary", "recurrent")) |> 
  dplyr::mutate(resection_recurrent = ifelse(resection == "primary", 0, 1)) |> 
  
  dplyr::mutate(covar = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(covar_name = "Astrocytoma CGC Lasso") |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p", ifelse(is.paired,patient_id, "_remainder"))))



model <- glm(resection_recurrent ~ patient + covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


#expnd <- (max(stats$covar)-min(stats$covar)) * 0.085
Predicted_data <- data.frame(covar=modelr::seq_range(c(-15, 15), 500)) |> 
  dplyr::mutate(patient = "p_remainder")
Predicted_data$resection_recurrent = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_recurrent, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1 ) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_recurrent - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_recurrent, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1 ) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_recurrent + 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(patient = NULL) |> 
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_recurrent) |> 
    dplyr::mutate(col = resection_recurrent)
)


p4b <- ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data"),
            col="white",
            lwd=theme_nature_lwd * 2, alpha=0.65
  ) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH","A_IDH_HG") == F),
            col="gray",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == T),
            col="darkblue",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH") == T),
            col="#3e9e8b",
            lwd=theme_nature_lwd) +  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(0, 1),
                                                           breaks=c(0, 0.50, 1), labels=c("Primary","","Recurrent"),
                                                           oob = scales::squish) +
  theme_nature +
  annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[7], x = 0.3, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability", y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = paste0(unique(stats$covar_name)," x Resection")) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("Primary", "Recurrent")) + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1)) +
  ylim(-15, 15)
p4b




#### logistic AcCGAP x WHO grade paired [validation] ----


stats <- metadata.validation |> 
  dplyr::filter(!is.na(resection_tumor_grade)) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(resection_tumor_grade__hg = ifelse(resection_tumor_grade == 3, 1 , 0)) |> 
  
  dplyr::mutate(resection = NULL) |> 
  dplyr::mutate(resection_recurrent = NULL) |> 
  
  dplyr::mutate(covar = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(covar_name = "Astrocytoma CGC Lasso") |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder"))))


model <- glm(resection_tumor_grade__hg ~ patient + covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


#expnd <- (max(stats$covar)-min(stats$covar)) * 0.085
Predicted_data <- data.frame(covar=modelr::seq_range(c(-15, 15), 500)) |> 
  dplyr::mutate(patient = "p_remainder")
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_tumor_grade__hg - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_tumor_grade__hg + 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(patient = NULL) |> 
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_tumor_grade__hg) |> 
    dplyr::mutate(col = resection_tumor_grade__hg)
) |> 
  dplyr::mutate(col = col  + 2)



p5b <- ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data"),
            col="white",
            lwd=theme_nature_lwd * 2, alpha=0.65
  ) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH","A_IDH_HG") == F),
            col="gray",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == T),
            col="darkblue",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH") == T),
            col="#3e9e8b",
            lwd=theme_nature_lwd) +  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 3), labels=c("Grade 2", "", "Grade 3"), oob = scales::squish) +
  theme_nature +
annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[7], x = 0.3, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability", y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = paste0(unique(stats$covar_name)," x WHO Grade")) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("Grade 2", "Grade 3")) + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1)) +
  ylim(-15, 15)
p5b



#### export ----



library(patchwork)
p3b + p4b + p5b + plot_layout(ncol=3) 


ggsave("output/figures/vis_LGC_PCA_logistic__paired.pdf", width=(8.5*0.975)*(0.5), height=3.0)



### stats ----


stats <- metadata |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(resection_tumor_grade__hg = ifelse(resection_tumor_grade == 3, 1 , 0)) |> 
  
  dplyr::mutate(resection = NULL) |> 
  dplyr::mutate(resection_recurrent = NULL) |> 
  
  dplyr::mutate(covar = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(covar_name = "Astrocytoma CGC Lasso") |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder")))) |> 
  
  dplyr::mutate(heidelberg_A_IDH_HG_OLIGOSARC = array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG", "OLIGOSARC_IDH"))

dim(stats)

summary(lm(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit ~ patient + heidelberg_A_IDH_HG_OLIGOSARC, data=stats))

glm(as.factor(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) ~ patient + heidelberg_A_IDH_HG_OLIGOSARC, data=stats, family = binomial)



plt <- stats |> 
  dplyr::mutate(label = ifelse(heidelberg_A_IDH_HG_OLIGOSARC, "A_IDH_HG /\nOLIGOSARC_IDH", "other")) |> 
  dplyr::mutate(label = factor(label, levels=c("other", "A_IDH_HG /\nOLIGOSARC_IDH")))

ggplot(plt, aes(x = label, y=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, col=array_mnp_predictBrain_v12.8_cal_class, group=label)) +
  ggbeeswarm::geom_quasirandom(size = theme_nature_size/3, width=0.3) +
  scale_color_manual(values = c(palette_mnp_12.8_6, `GBM_RTK1`="gray40", `GBM_MES_TYP` = "gray40")) +
  labs(x = NULL, y="CGC") +
  #annotate("text", family = theme_nature_font_family, 
  #         size=theme_nature_size,
  #         x=1.9, y=-10, 
  #         label=paste0("p < 2e-16 ***\nlm; pat. corrected")) +
  ggpubr::stat_compare_means(method = "t.test", show_guide  = FALSE, size=theme_nature_size, family=theme_nature_font_family) +
  theme_nature

ggsave("output/figures/vis_LGC_x_PCA__lm_OLIGOSARC_A_IDH_HG__Other.pdf", width=1.5, height=3 * 0.9305)



### logistic AcCGAP x resection ----


stats <- metadata |> 
  dplyr::mutate(resection_tumor_grade = NULL) |> 
  dplyr::mutate(resection_tumor_grade__hg = NULL) |> 
  
  dplyr::mutate(resection = ifelse(resection_number == 1, "primary", "recurrent")) |> 
  dplyr::mutate(resection_recurrent = ifelse(resection == "primary", 0 , 1)) |> 
  
  dplyr::mutate(covar = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(covar_name = "Astrocytoma CGC Lasso")


model <- glm(resection_recurrent ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$resection_recurrent = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_recurrent, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1 ) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_recurrent - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_recurrent, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1 ) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_recurrent + 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_recurrent) |> 
    dplyr::mutate(col = resection_recurrent)
)


p4 <- ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data"),
            col="white",
            lwd=theme_nature_lwd * 2, alpha=0.65
  ) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH","A_IDH_HG") == F),
            col="gray",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == T),
            col="darkblue",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH") == T),
            col="#3e9e8b",
            lwd=theme_nature_lwd) +  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(0, 1),
                        breaks=c(0, 0.50, 1), labels=c("Primary","","Recurrent"),
                        oob = scales::squish) +
  theme_nature +
annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[7], x = 0.3, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability", y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = paste0(unique(stats$covar_name)," x Resection")) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("Primary", "Recurrent")) + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1))
p4



# complex color variant:
# plt.logit.restyled <- rbind(
#   stats |>  # left point line:
#     dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class, resection_recurrent) |> 
#     dplyr::mutate(col = resection_tumor_grade__hg + 2) |> 
#     dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
#     dplyr::mutate(type = "data") |> 
#     dplyr::mutate(x = ifelse(resection_recurrent == 0, resection_recurrent, resection_recurrent - 0.175)) |> 
#     #dplyr::mutate(x = resection_recurrent - 0.06 ) |> 
#     dplyr::rename(y = covar)
#   ,
#   stats |>  # right point line:
#     dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class, resection_recurrent) |> 
#     dplyr::mutate(col = resection_tumor_grade__hg + 2) |> 
#     dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
#     dplyr::mutate(type = "data") |> 
#     #dplyr::mutate(x = resection_recurrent + 0.06 ) |> 
#     dplyr::mutate(x = ifelse(resection_recurrent == 1, resection_recurrent, resection_recurrent + 0.175)) |>
#     dplyr::rename(y = covar)
#   ,
#   Predicted_data |> # logit fit
#     dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
#     dplyr::mutate(resection_tumor_grade__hg = 0) |> 
#     dplyr::mutate(group = "logit fit") |> 
#     dplyr::mutate(type = "fit") |> 
#     dplyr::rename(y = covar) |> 
#     dplyr::mutate(x = resection_recurrent) |> 
#     dplyr::mutate(col = resection_recurrent)
# )

# plt.logit.restyled.resection <- plt.logit.restyled

# p4 <- ggplot(plt.logit.restyled, aes(x=x, y=y, group=group, col=col)) +
#   # geom_line(data = plt.logit.restyled |> dplyr::filter(type == "data"),
#   #           col="white",
#   #           lwd=theme_nature_lwd * 2, alpha=0.65
#   # ) +
#   geom_line(data = plt.logit.restyled |> dplyr::filter(type == "data"),
#             aes(col=col),
#             lwd=theme_nature_lwd) +
#   scale_color_gradientn(colours = c("aquamarine3","aquamarine3","#d34394ff","#d34394ff"), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 2.50, 3), labels=c("Grade 2","","", "Grade 3"), oob = scales::squish) +
#   labs(col=NULL) +
#   ggnewscale::new_scale_colour() +
#   geom_line(data = plt.logit.restyled |> dplyr::filter(type == "fit") ,
#             aes(col=col),
#             lwd=1.0) +
#   scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(0, 1), breaks=c(0, 1), labels=c("Primary", "Recurrent"), oob = scales::squish) +
#   theme_nature +
# annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 16)[15], x = 0.375, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
#   labs(col=NULL, y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = "CGC Lasso x Resection") +
#   scale_x_continuous(breaks = c(0,1), labels=c("Primary", "Recurrent")) + 
#   theme(legend.box = "vertical") + # space dependent
#   theme(legend.key.size = unit(0.6, 'lines'))
# p4




### logistic PC2 x WHO grade ----


stats <- metadata |> 
  assertr::verify(resection_tumor_grade %in% c(2, 3)) |> 
  dplyr::mutate(resection_tumor_grade__hg = ifelse(resection_tumor_grade == 3, 1 , 0)) |> 
  
  dplyr::mutate(resection = NULL) |> 
  dplyr::mutate(resection_recurrent = NULL) |> 
  
  dplyr::mutate(covar = array_PC2) |> 
  dplyr::mutate(covar_name = "PC2 (unsupervised)")



model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_tumor_grade__hg - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_tumor_grade__hg + 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_tumor_grade__hg) |> 
    dplyr::mutate(col = resection_tumor_grade__hg)
) |> 
  dplyr::mutate(col = col  + 2)



p5 <- ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data"),
            col="white",
            lwd=theme_nature_lwd * 2, alpha=0.65
  ) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH","A_IDH_HG") == F),
            col="gray",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == T),
            col="darkblue",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH") == T),
            col="#3e9e8b",
            lwd=theme_nature_lwd) +  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 2.50, 3), labels=c("Grade 2","","", "Grade 3"), oob = scales::squish) +
  theme_nature +
annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[2], x = 0.3, label = paste0("p = ", format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability", y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = paste0(unique(stats$covar_name)," x WHO Grade")) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("Grade 2", "Grade 3")) + 
  scale_y_reverse() + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1))
p5



t.test(
  stats |>
    dplyr::filter(resection_tumor_grade == 2) |> 
    dplyr::pull(covar),
  stats |>
    dplyr::filter(resection_tumor_grade == 3) |> 
    dplyr::pull(covar)
)


### logistic PC2 x resection ----


stats <- metadata |> 
  dplyr::mutate(resection_tumor_grade = NULL) |> 
  dplyr::mutate(resection_tumor_grade__hg = NULL) |> 
  
  dplyr::mutate(resection = ifelse(resection_number == 1, "primary", "recurrent")) |> 
  dplyr::mutate(resection_recurrent = ifelse(resection == "primary", 0 , 1)) |> 
  
  dplyr::mutate(covar = array_PC2) |> 
  dplyr::mutate(covar_name = "PC2 (unsupervised)")



model <- glm(resection_recurrent ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$resection_recurrent = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_recurrent, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1 ) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_recurrent - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_recurrent, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1 ) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_recurrent + 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_recurrent) |> 
    dplyr::mutate(col = resection_recurrent)
)


p6 <- ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data"),
            col="white",
            lwd=theme_nature_lwd * 2, alpha=0.65
  ) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH","A_IDH_HG") == F),
            col="gray",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == T),
            col="darkblue",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH") == T),
            col="#3e9e8b",
            lwd=theme_nature_lwd) +  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(0, 1),
                                                           breaks=c(0, 0.50, 1), labels=c("Primary","","Recurrent"),
                                                           oob = scales::squish) +
  theme_nature +
annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[2], x = 0.3, label = paste0("p = ", format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability", y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = paste0(unique(stats$covar_name)," x Resection")) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("Primary", "Recurrent")) + 
  scale_y_reverse() + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1))
p6




### logistic PC3 x WHO grade ----


stats <- metadata |> 
  assertr::verify(resection_tumor_grade %in% c(2, 3)) |> 
  dplyr::mutate(resection_tumor_grade__hg = ifelse(resection_tumor_grade == 3, 1 , 0)) |> 
  
  dplyr::mutate(resection = NULL) |> 
  dplyr::mutate(resection_recurrent = NULL) |> 
  
  dplyr::mutate(covar = array_PC3) |> 
  dplyr::mutate(covar_name = "PC3 (unsupervised)")



model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_tumor_grade__hg - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_tumor_grade__hg + 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_tumor_grade__hg) |> 
    dplyr::mutate(col = resection_tumor_grade__hg)
) |> 
  dplyr::mutate(col = col  + 2)



p7 <- ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data"),
            col="white",
            lwd=theme_nature_lwd * 2, alpha=0.65
  ) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH","A_IDH_HG") == F),
            col="gray",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == T),
            col="darkblue",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH") == T),
            col="#3e9e8b",
            lwd=theme_nature_lwd) +  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 2.50, 3), labels=c("Grade 2","","", "Grade 3"), oob = scales::squish) +
  theme_nature +
annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[2], x = 0.3, label = paste0("p = ", format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability", y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = paste0(unique(stats$covar_name)," x WHO Grade")) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("Grade 2", "Grade 3")) + 
  scale_y_reverse() + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1))
p7



t.test(
  stats |>
    dplyr::filter(resection_tumor_grade == 2) |> 
    dplyr::pull(covar),
  stats |>
    dplyr::filter(resection_tumor_grade == 3) |> 
    dplyr::pull(covar)
)


### logistic PC3 x resection ----


stats <- metadata |> 
  dplyr::mutate(resection_tumor_grade = NULL) |> 
  dplyr::mutate(resection_tumor_grade__hg = NULL) |> 
  
  dplyr::mutate(resection = ifelse(resection_number == 1, "primary", "recurrent")) |> 
  dplyr::mutate(resection_recurrent = ifelse(resection == "primary", 0 , 1)) |> 
  
  dplyr::mutate(covar = array_PC3) |> 
  dplyr::mutate(covar_name = "PC3 (unsupervised)")



model <- glm(resection_recurrent ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$resection_recurrent = predict(model, Predicted_data, type="response")



plt.logit.simplistic <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_recurrent, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1 ) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_recurrent - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_recurrent, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(col = -1 ) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_recurrent + 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_recurrent) |> 
    dplyr::mutate(col = resection_recurrent)
)


p8 <- ggplot(plt.logit.simplistic, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=theme_nature_lwd * 5) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data"),
            col="white",
            lwd=theme_nature_lwd * 2, alpha=0.65
  ) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH","A_IDH_HG") == F),
            col="gray",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == T),
            col="darkblue",
            lwd=theme_nature_lwd) +
  geom_line(data = plt.logit.simplistic |> dplyr::filter(type == "data" & array_mnp_predictBrain_v12.8_cal_class %in% c("OLIGOSARC_IDH") == T),
            col="#3e9e8b",
            lwd=theme_nature_lwd) +  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(0, 1),
                                                           breaks=c(0, 0.50, 1), labels=c("Primary","","Recurrent"),
                                                           oob = scales::squish) +
  theme_nature +
annotate("text", family = theme_nature_font_family,  y = modelr::seq_range(stats$covar, 8)[2], x = 0.3, label = paste0("p = ", format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col="Probability", y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = paste0(unique(stats$covar_name)," x Resection")) +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("Primary", "Recurrent")) + 
  scale_y_reverse() + 
  theme(legend.box = "vertical",
        legend.key.size = unit(theme_nature_lwd * 1.5, 'lines'),
        legend.title =  element_text(vjust=1))
p8


t.test(
  stats |>
    dplyr::filter(resection_number == 1) |> 
    dplyr::pull(covar),
  stats |>
    dplyr::filter(resection_number > 1) |> 
    dplyr::pull(covar)
)






### export ----

p1 + p3 + p5 + p2 + p4 + p6 + plot_layout(ncol=3) + plot_annotation( subtitle = format_subtitle("PC2 & CGC x Grade & Resection timing") )


#ggsave("output/figures/vis_LGC_PCA_logistic.pdf", width=(8.5*0.975)*(2/3), height=4.5)



library(patchwork)
p3 + p5 + p4 + p6 + plot_layout(ncol=4) + plot_annotation(subtitle = format_subtitle("PC2 & CGC x Grade & Resection timing"), theme=theme_nature )

ggsave("output/figures/vis_LGC_PCA_logistic.pdf", width=(8.5*0.975)*(4/5), height=2.5)





### logistic Zhang2017 Mortality risk clock ----




scaling_factor_x <- 0.002*5*(0.5/0.83)
stats <- stats |> 
  dplyr::mutate(covar = array_dnaMethyAge__ZhangY2017) |> 
  dplyr::mutate(covar_name = "Zhang 2017 Mortality risk clock")



model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


Predicted_data <- data.frame(covar=modelr::seq_range(stats$covar, 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  stats |> 
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data"),
  Predicted_data |> 
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2) |> 
  dplyr::mutate(xmin=covar - scaling_factor_x,
                xmax=covar + scaling_factor_x,
                xmin_wash = covar - (2.25 * scaling_factor_x), 
                xmax_wash = covar + (2.25 * scaling_factor_x), 
                ymin=`resection tumor grade` - 0.06,
                ymax=`resection tumor grade` + 0.06)



p3 <- ggplot(plt.logit, aes(x = covar, y=`resection tumor grade`)) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), aes(col=`resection tumor grade`),lwd=2) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin_wash,xmax=xmax_wash,ymin=ymin,ymax=ymax), fill=alpha("white",0.65)) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=`array_mnp_predictBrain_v12.8_cal_class`)) +
annotate("text", family = theme_nature_font_family,  x = modelr::seq_range(stats$covar, 8)[2], y = 2.5, label = paste0("p = ",format.pval(pval)), size=theme_nature_size) + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col=NULL, x = stats |> dplyr::pull(covar_name) |> unique(), fill=NULL) +
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 3), oob = scales::squish) +
  scale_fill_manual(values = palette_mnp_12.8_6) +
  theme_cellpress +
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
p3



#### t-test ----


t.test(
  stats |>
    dplyr::filter(resection_tumor_grade == 2) |> 
    dplyr::pull(covar),
  stats |>
    dplyr::filter(resection_tumor_grade == 3) |> 
    dplyr::pull(covar)
)




### logistic epiTOC2 ----



scaling_factor_x <- 16 * (0.5 / 0.58)
stats <- stats |> 
  dplyr::mutate(covar = array_epiTOC2_tnsc) |> 
  dplyr::mutate(covar_name = "epiTOC2")



model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


Predicted_data <- data.frame(covar=modelr::seq_range(stats$covar, 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  stats |> 
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data"),
  Predicted_data |> 
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2) |> 
  dplyr::mutate(xmin=covar - scaling_factor_x,
                xmax=covar + scaling_factor_x,
                xmin_wash = covar - (2.25 * scaling_factor_x), 
                xmax_wash = covar + (2.25 * scaling_factor_x), 
                ymin=`resection tumor grade` - 0.06,
                ymax=`resection tumor grade` + 0.06)



p4 <- ggplot(plt.logit, aes(x = covar, y=`resection tumor grade`)) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), aes(col=`resection tumor grade`),lwd=2) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin_wash,xmax=xmax_wash,ymin=ymin,ymax=ymax), fill=alpha("white",0.65)) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=`array_mnp_predictBrain_v12.8_cal_class`)) +
annotate("text", family = theme_nature_font_family,  x = modelr::seq_range(stats$covar, 8)[2], y = 2.5, label = paste0("p = ",format.pval(pval)), size=theme_nature_size) + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col=NULL, x = stats |> dplyr::pull(covar_name) |> unique(), fill=NULL) +
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 3), oob = scales::squish) +
  scale_fill_manual(values = palette_mnp_12.8_6) +
  theme_cellpress +
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
p4



#### t-test ----


t.test(
  stats |>
    dplyr::filter(resection_tumor_grade == 2) |> 
    dplyr::pull(covar),
  stats |>
    dplyr::filter(resection_tumor_grade == 3) |> 
    dplyr::pull(covar)
)




### logistic median.overall.methylation ----



scaling_factor_x <- 0.0185 * 0.16 * (0.5/0.48)
stats <- stats |> 
  dplyr::mutate(covar = array_median.overall.methylation) |> 
  dplyr::mutate(covar_name = "Median methylation")



model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


Predicted_data <- data.frame(covar=modelr::seq_range(stats$covar, 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  stats |> 
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data"),
  Predicted_data |> 
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2) |> 
  dplyr::mutate(xmin=covar - scaling_factor_x,
                xmax=covar + scaling_factor_x,
                xmin_wash = covar - (2.25 * scaling_factor_x), 
                xmax_wash = covar + (2.25 * scaling_factor_x), 
                ymin=`resection tumor grade` - 0.06,
                ymax=`resection tumor grade` + 0.06)



p3 <- ggplot(plt.logit, aes(x = covar, y=`resection tumor grade`)) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), aes(col=`resection tumor grade`),lwd=2) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin-(5 * scaling_factor_x),xmax=xmax+(5 * scaling_factor_x),ymin=ymin,ymax=ymax), fill=alpha("white",0.65)) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=`array_mnp_predictBrain_v12.8_cal_class`)) +
annotate("text", family = theme_nature_font_family,  x = modelr::seq_range(stats$covar, 8)[2], y = 2.5, label = paste0("p = ",format.pval(pval)), size=theme_nature_size) + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col=NULL, x = stats |> dplyr::pull(covar_name) |> unique(), fill=NULL) +
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 3), oob = scales::squish) +
  scale_fill_manual(values = palette_mnp_12.8_6) +
  theme_cellpress +
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
p3



### logistic median.glass_nl_supervised.methylation ----



scaling_factor_x <- 0.0185 * 0.22 * (0.5/0.39)
stats <- stats |> 
  dplyr::mutate(covar = array_median.glass_nl_supervised.methylation) |> 
  dplyr::mutate(covar_name = "Median methylation (AC primary - recurrence DMP's)")



model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


Predicted_data <- data.frame(covar=modelr::seq_range(stats$covar, 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  stats |> 
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data"),
  Predicted_data |> 
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2) |> 
  dplyr::mutate(xmin=covar - scaling_factor_x,
                xmax=covar + scaling_factor_x,
                xmin_wash = covar - (2.25 * scaling_factor_x), 
                xmax_wash = covar + (2.25 * scaling_factor_x), 
                ymin=`resection tumor grade` - 0.06,
                ymax=`resection tumor grade` + 0.06)



p4 <- ggplot(plt.logit, aes(x = covar, y=`resection tumor grade`)) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), aes(col=`resection tumor grade`),lwd=2) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin_wash,xmax=xmax_wash,ymin=ymin,ymax=ymax), fill=alpha("white",0.65)) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=`array_mnp_predictBrain_v12.8_cal_class`)) +
annotate("text", family = theme_nature_font_family,  x = modelr::seq_range(stats$covar, 8)[2], y = 2.5, label = paste0("p = ",format.pval(pval)), size=theme_nature_size) + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col=NULL, x = stats |> dplyr::pull(covar_name) |> unique(), fill=NULL) +
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 3), oob = scales::squish) +
  scale_fill_manual(values = palette_mnp_12.8_6) +
  theme_cellpress +
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
p4





### logistic hypoSC solo-WCGW ----



scaling_factor_x <- 0.00035 * (0.5 / 0.38)
stats <- stats |> 
  dplyr::mutate(covar = array_epiTOC2_hypoSC) |> 
  dplyr::mutate(covar_name = "solo-WCGW HypoClock score")



model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


Predicted_data <- data.frame(covar=modelr::seq_range(stats$covar, 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  stats |> 
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data"),
  Predicted_data |> 
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2) |> 
  dplyr::mutate(xmin=covar - scaling_factor_x,
                xmax=covar + scaling_factor_x,
                xmin_wash = covar - (2.25 * scaling_factor_x), 
                xmax_wash = covar + (2.25 * scaling_factor_x),
                ymin=`resection tumor grade` - 0.06,
                ymax=`resection tumor grade` + 0.06)




p6 <- ggplot(plt.logit, aes(x = covar, y=`resection tumor grade`)) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), aes(col=`resection tumor grade`),lwd=2) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin_wash,xmax=xmax_wash,ymin=ymin,ymax=ymax), fill=alpha("white",0.65)) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=`array_mnp_predictBrain_v12.8_cal_class`)) +
annotate("text", family = theme_nature_font_family,  x = modelr::seq_range(stats$covar, 8)[2], y = 2.5, label = paste0("p = ",format.pval(pval)), size=theme_nature_size) + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col=NULL, x = stats |> dplyr::pull(covar_name) |> unique(), fill=NULL) +
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 3), oob = scales::squish) +
  scale_fill_manual(values = palette_mnp_12.8_6) +
  theme_cellpress +
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
p6



### combined export ----

library(patchwork)

(p1 + p2) /
(p3 + p4) +
  plot_annotation(theme = theme_cellpress,
                  subtitle=format_subtitle("Logit WHO grade"))


ggsave("output/figures/vis_LGC_x_PCA__logits.pdf", width = (8.5 * 0.975), height = 4)



### multi LM ----




model <- lm(resection_tumor_grade__hg ~ PC2 + A_IDH_HG__A_IDH_LG_lr__lasso_fit, data = stats)
model.g <- glm(resection_tumor_grade__hg ~ PC2 + A_IDH_HG__A_IDH_LG_lr__lasso_fit, data = stats,family = binomial)
summary(model)
summary(model.g)
sjPlot::plot_model(model)


model2 <- lm(A_IDH_HG__A_IDH_LG_lr__lasso_fit ~ resection_tumor_grade__hg + batch_us, data = stats)
summary(model2)
sjPlot::plot_model(model2)

model3 <- lm(-PC2 ~ resection_tumor_grade__hg + batch_us, data = stats)
summary(model3)
sjPlot::plot_model(model3)


model3 <- lm(PC2 ~ resection_tumor_grade__hg, data = stats)
summary(model3)





# G-SAM ----


metadata <- gsam.metadata.array_samples |> 
  filter_GSAM_idats(CONST_N_GSAM_INCLUDED_SAMPLES) |> 
  dplyr::mutate(resection__rec = ifelse(resection == "R1",0,1))



## logistic LR & recurrence ----


data <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
    return(.)
  })() |> 
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |> # this syntax, oh my
  dplyr::arrange(mad) |> 
  dplyr::mutate(mad = NULL)


ggplot(metadata, aes(x=resection, y=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) +
  ggbeeswarm::geom_quasirandom() +
  ggpubr::stat_compare_means(method = "t.test") +
  theme_nature


model <- glm(resection__rec ~  A_IDH_HG__A_IDH_LG_lr, 
             data = metadata,
             family = binomial)
summary(model)


pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "A_IDH_HG__A_IDH_LG_lr") |> 
  dplyr::pull(`Pr(>|z|)`)


Predicted_data <- data.frame(A_IDH_HG__A_IDH_LG_lr=modelr::seq_range(metadata$A_IDH_HG__A_IDH_LG_lr, 500))
Predicted_data$resection__rec = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  metadata |> 
    dplyr::select(resection__rec, A_IDH_HG__A_IDH_LG_lr, mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data")
  ,
  Predicted_data |> 
    dplyr::mutate(mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
)


ggplot(plt.logit, aes(x = A_IDH_HG__A_IDH_LG_lr, y=`resection__rec`, col= mnp_predictBrain_v12.8_cal_class)) +
  geom_point(data = plt.logit |> dplyr::filter(type == "data"),pch="|",size=5) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), col="darkgreen") +
  theme_bw() + 
annotate("text", family = theme_nature_font_family,  x = 2.5, y = 0.5, label = paste0("p = ",format.pval(pval))) + 
  theme(legend.position = 'bottom') + 
  scale_y_continuous(breaks = c(0,1)) + 
  labs(col="mnp v12.8", y = "Resection (0 primary, 1 recurrent)")


## lm/glm/limma ----


tmp <- metadata |>
  dplyr::select(resection_id, patient_id, resection, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"remainder")))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection == "R1","primary","recurrence"),levels=c("primary","recurrence")))

fit = lm(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit ~ patient + resection, data=tmp)
summary(fit)





coef_summary <- summary(fit)$coefficients
conf_int <- confint(fit)


p_x <- coef_summary |>
  as.data.frame() |> 
  tibble::rownames_to_column("Term") |>
  dplyr::filter(Term == "resectionR2") |>
  dplyr::select(Estimate) |>
  as.numeric()

p <- coef_summary |>
  as.data.frame() |> 
  tibble::rownames_to_column("Term") |>
  dplyr::filter(Term == "resectionR2") |>
  dplyr::select(`Pr(>|t|)`) |>
  as.numeric()

# Combine coefficients and confidence intervals into one data frame
forest_data <- data.frame(
  Term = rownames(coef_summary),
  Estimate = coef_summary[, 1],
  Lower = conf_int[, 1],
  Upper = conf_int[, 2]
) |> 
  dplyr::filter(Term == "resectionR2") |> 
  dplyr::mutate(Term = "recurrent resection")

ggplot(forest_data, aes(x = Estimate, y = Term)) +
  geom_point(size=theme_nature_size/3) + 
  geom_vline(xintercept = 0, linetype="dotted",  color = "gray", lwd=theme_nature_lwd) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.52, lwd=theme_nature_lwd) + 
  labs(x="CGC: LM Estimate (glioblastoma)", y=NULL) +
  annotate("text", x = p_x, y = 1.4, label = paste0("p = ",round(p, 3)), size=theme_nature_size, family  = theme_nature_font_family) +
  theme_nature

ggsave("output/figures/vis_LGC_x_PCA__GSAM__lm.pdf", width=8.5*0.975*1/3, height=0.6)



## PCA loadings comparison ----


metadata <- gsam.metadata.array_samples |> 
  filter_GSAM_idats(CONST_N_GSAM_INCLUDED_SAMPLES - 2 ) |> 
  dplyr::select(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, 
                array_percentage.detP.signi,
                #array_dnaMethyAge__PCHorvathS2013, 
                array_dnaMethyAge__PCHorvathS2018,
                
                paste0("array_PC",1:5))

ggcorrplot(metadata)


ggplot(metadata, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=array_PC1)) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)),  size=theme_nature_size, cor.coef.name ="rho", show.legend = FALSE, label.y.npc = "top", label.x.npc = "right", hjust=1) + # spearman bc of  outlier at bottom
  geom_point(size=theme_nature_size/3) +
  theme_nature



# get loadings of original

loadings_glod <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_prcomp.Rds")
loadings_gsam <- readRDS("cache/analysis_unsupervised_PCA_GSAM_prcomp.Rds")

a = data.frame(PC2_glod = loadings_glod$rotation[,"PC2"]) |> 
  tibble::rownames_to_column("probe_id") 
b = as.data.frame(loadings_gsam$rotation[,paste0("PC",1:20)]) |> 
  tibble::rownames_to_column("probe_id") 

m <- a |>
  dplyr::left_join(b, by=c("probe_id"="probe_id")) |> 
  tibble::column_to_rownames("probe_id")


ggcorrplot(m)




# GLASS-NL / AC ----


metadata <- glass_nl.metadata.array_samples |> 
  dplyr::filter(!qc.pca.detP.outlier) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (218))
    return(.)
  })()


## PCA ----


masked <- readRDS('cache/mvalues.HQ_samples.detP_mask.Rds') |> 
  is.na() |> 
  rowSums() |> 
  as.data.frame() |> 
  dplyr::rename(count = 1) |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(count > 0) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (66116))
    return(.)
  })() |> 
  dplyr::pull('probe_id')




data <- readRDS('cache/mvalues.HQ_samples.Rds') 

data <- data |> 
  dplyr::select(metadata$sentrix_id)

data <- data |> 
  tibble::rownames_to_column('probe_id')

data <- data |> 
  dplyr::filter(probe_id %in% (masked) == F) |> 
  tibble::column_to_rownames('probe_id')

dim(data)


data <- data |> 
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |> # this synthax, oh my
  dplyr::arrange(mad) |> 
  dplyr::mutate(mad = NULL)


data.all.pca.obj <- data |> 
  #dplyr::slice_head(n=100000) |> 
  t() |> 
  prcomp()

data.all.pca <- data.all.pca.obj |> 
  purrr::pluck('x') |> 
  as.data.frame(stringsAsFactors=F) |> 
  dplyr::select(paste0("PC",1:40)) |> 
  tibble::rownames_to_column('sentrix_id')


plt <- metadata |> 
  dplyr::left_join(data.all.pca, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))


plt <- plt |> 
  dplyr::mutate(mnp_predictBrain_v12.8_cal_A_IDH_LG_odds = (mnp_predictBrain_v12.8_cal_A_IDH_LG/100) / (1 - (mnp_predictBrain_v12.8_cal_A_IDH_LG/100))) |> 
  dplyr::mutate(mnp_predictBrain_v12.8_cal_A_IDH_HG_odds = (mnp_predictBrain_v12.8_cal_A_IDH_HG/100) / (1 - (mnp_predictBrain_v12.8_cal_A_IDH_HG/100))) |> 
  dplyr::mutate(A_IDH_HG__A_IDH_LG_lr = log(mnp_predictBrain_v12.8_cal_A_IDH_HG_odds / mnp_predictBrain_v12.8_cal_A_IDH_LG_odds))


ggplot(plt, aes(x=A_IDH_HG__A_IDH_LG_lr, y=PC1)) +
  geom_point()




# detP mask filter
# order on MAD
# pick top 200.000
# prcomp


# plot(plt$PC2, plt$median.overall.methylation)
# plot(plt$PC2, plt$A_IDH_HG__A_IDH_LG_lr)
# plot(plt$PC2, plt$A_IDH_HG__A_IDH_LG_lr__lasso_fit)
# plot(plt$median.overall.methylation, plt$A_IDH_HG__A_IDH_LG_lr__lasso_fit)
# plt$A_IDH_HG__O_IDH_lr, plt$A_IDH_HG__A_IDH_LG_lr)


# c <- plt |>
#   dplyr::mutate(mnp_predictBrain_v12.8_cal_A_IDH_HG = -log(mnp_predictBrain_v12.8_cal_A_IDH_HG)) |> 
#   dplyr::select(
#   PC2,
#   A_IDH_HG__A_IDH_LG_lr,
#   A_IDH_HG__A_IDH_LG_lr__lasso_fit,
#   A_IDH_HGoligsarc__O_IDH_lr,
#   median.overall.methylation,
#   median.glass_nl_supervised.methylation,
#   A_IDH_HG__O_IDH_lr
#   #mnp_predictBrain_v12.8_cal_A_IDH_HG
# ) |> 
#   as.matrix()
# 
# 
corrplot::corrplot(abs(cor(c, method="spearman")), order="hclust")
corrplot::corrplot(abs(cor(c, method="kendall")), order="hclust")
# 
# 
# plot(plt$A_IDH_HGoligsarc__O_IDH_lr, plt$PC2)
# plot(plt$A_IDH_HG__O_IDH_lr, plt$PC2)
# 
# plot(plt$A_IDH_HGoligsarc__O_IDH_lr, plt$PC2)
# plot(plt$PC2)
# 
