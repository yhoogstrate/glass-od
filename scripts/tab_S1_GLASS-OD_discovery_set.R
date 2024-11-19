#!/usr/bin/env R

# load data ----


library(ggplot2)

source('scripts/load_constants.R')
source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}



# export ----


out_meth <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::select(array_sentrix_id,
                resection_id,
                patient_id,
                resection_number,
                resection_tumor_grade,
                
                array_A_IDH_HG__A_IDH_LG_lr__lasso_fit,
                
                resection_treatment_status_radio,
                resection_treatment_status_chemo,
                resection_treatment_status_summary,
                
                time_between_resection_and_last_event,
                patient_last_follow_up_event,
                
                array_mnp_predictBrain_v2.0.1_cal_class,
                array_mnp_predictBrain_v12.5_cal_class,
                array_mnp_predictBrain_v12.8_cal_class
                
                )


write.table(out_meth, file="/tmp/exp.txt", sep="\t", quote = F, row.names = F)


out_ki67 <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::select(sentrix_id)


out_proteomics <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::select(sentrix_id)



