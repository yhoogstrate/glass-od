#!/usr/bin/env R


# load data ----


source('scripts/load_functions.R')


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


if(!exists('cortex_dilution')) {
  source('scripts/load_nonmalignant_cortex_dilution_series.R')
}



# export to table ----


cortex_dilution.metadata.array_samples |> 
  dplyr::mutate(array_tumor_id = factor(array_tumor_id, levels=c("0017-R3","0008-R2","0121-R3","0054-R3"))) |> 
  dplyr::arrange(array_fraction_normal_brain, array_tumor_id, array_control_id) |> 
  dplyr::select(array_sentrix_id, array_tumor_id, array_control_id ,array_fraction_normal_brain,
                
                array_mnp_predictBrain_v12.8_cal_A_IDH_HG,
                array_mnp_predictBrain_v12.8_cal_OLIGOSARC_IDH,
                array_mnp_predictBrain_v12.8_cal_CTRL_HEMI,
  ) 

cortex_dilution.metadata.array_samples |> 
  dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = NULL) |> 
  dplyr::select(contains("edictBrain_v12.8_cal_")) |> 
  tidyr::pivot_longer(cols = dplyr::everything()) |> 
  dplyr::group_by(name) |> 
  dplyr::summarise(maxp = max(value)) |> 
  dplyr::ungroup() |> 
  dplyr::arrange(-maxp) |> 
  View()




## func ----

#cortex_dilution.metadata.array_samples$array_tumor_id |> unique()



exporttable <- function(tumor, normal) {
  
  # diluted
  exp.diluted <- cortex_dilution.metadata.array_samples |> 
    dplyr::filter(array_tumor_id == tumor & array_control_id == normal) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) %in% c(10))
      return(.)
    })() |> 
    dplyr::arrange(array_fraction_normal_brain) |> 
    dplyr::select(array_fraction_normal_brain,

                  
                  array_mnp_predictBrain_v12.8_cal_A_IDH_HG,
                  array_mnp_predictBrain_v12.8_cal_OLIGOSARC_IDH,
                  array_mnp_predictBrain_v12.8_cal_O_IDH,
                  array_mnp_predictBrain_v12.8_cal_A_IDH_LG,
                  
                  array_mnp_predictBrain_v12.8_cal_CTRL_HEMI,
                  array_mnp_predictBrain_v12.8_cal_CTRL_HYPOTHAL,
                  
                  array_mnp_predictBrain_v12.8_cal_PA_CORT,
                  array_mnp_predictBrain_v12.8_cal_CTRL_CORPCAL,
                  array_mnp_predictBrain_v12.8_cal_DNET,
                  #array_mnp_predictBrain_v12.8_cal_GBM_MES_TYP,
                  array_mnp_predictBrain_v12.8_cal_GG,
                  array_mnp_predictBrain_v12.8_cal_PGNT,
                  array_mnp_predictBrain_v12.8_cal_HGG_F
    )
  
  
  # control
  exp.ctrl <- normal_cortex.metadata.array_samples |> 
    dplyr::filter(!is.na(array_control_id) & array_control_id == normal) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == 1)
      return(.)
    })() |> 
    dplyr::select(
      array_mnp_predictBrain_v12.8_cal_A_IDH_HG,
      array_mnp_predictBrain_v12.8_cal_OLIGOSARC_IDH,
      array_mnp_predictBrain_v12.8_cal_O_IDH,
      array_mnp_predictBrain_v12.8_cal_A_IDH_LG,
      
      array_mnp_predictBrain_v12.8_cal_CTRL_HEMI,
      array_mnp_predictBrain_v12.8_cal_CTRL_HYPOTHAL,
      
      array_mnp_predictBrain_v12.8_cal_PA_CORT,
      array_mnp_predictBrain_v12.8_cal_CTRL_CORPCAL,
      array_mnp_predictBrain_v12.8_cal_DNET,
      #array_mnp_predictBrain_v12.8_cal_GBM_MES_TYP,
      array_mnp_predictBrain_v12.8_cal_GG,
      array_mnp_predictBrain_v12.8_cal_PGNT,
      array_mnp_predictBrain_v12.8_cal_HGG_F
    ) |> 
    round(3) |> 
    dplyr::mutate(array_fraction_normal_brain = "100") |> 
    dplyr::select(array_fraction_normal_brain,
                  
                  
                  array_mnp_predictBrain_v12.8_cal_A_IDH_HG,
                  array_mnp_predictBrain_v12.8_cal_OLIGOSARC_IDH,
                  array_mnp_predictBrain_v12.8_cal_O_IDH,
                  array_mnp_predictBrain_v12.8_cal_A_IDH_LG,
                  
                  array_mnp_predictBrain_v12.8_cal_CTRL_HEMI,
                  array_mnp_predictBrain_v12.8_cal_CTRL_HYPOTHAL,
                  
                  array_mnp_predictBrain_v12.8_cal_PA_CORT,
                  array_mnp_predictBrain_v12.8_cal_CTRL_CORPCAL,
                  array_mnp_predictBrain_v12.8_cal_DNET,
                  #array_mnp_predictBrain_v12.8_cal_GBM_MES_TYP,
                  array_mnp_predictBrain_v12.8_cal_GG,
                  array_mnp_predictBrain_v12.8_cal_PGNT,
                  array_mnp_predictBrain_v12.8_cal_HGG_F
    )
  
  
  
  # tumor
  exp.tumor <- glass_od.metadata.array_samples |> 
    filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
    dplyr::filter(resection_id == tumor) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == 1)
      return(.)
    })() |> 
    dplyr::select(
      array_mnp_predictBrain_v12.8_cal_A_IDH_HG,
      array_mnp_predictBrain_v12.8_cal_OLIGOSARC_IDH,
      array_mnp_predictBrain_v12.8_cal_O_IDH,
      array_mnp_predictBrain_v12.8_cal_A_IDH_LG,
      
      array_mnp_predictBrain_v12.8_cal_CTRL_HEMI,
      array_mnp_predictBrain_v12.8_cal_CTRL_HYPOTHAL,
      
      array_mnp_predictBrain_v12.8_cal_PA_CORT,
      array_mnp_predictBrain_v12.8_cal_CTRL_CORPCAL,
      array_mnp_predictBrain_v12.8_cal_DNET,
      #array_mnp_predictBrain_v12.8_cal_GBM_MES_TYP,
      array_mnp_predictBrain_v12.8_cal_GG,
      array_mnp_predictBrain_v12.8_cal_PGNT,
      array_mnp_predictBrain_v12.8_cal_HGG_F
    ) |> 
    round(3) |> 
    dplyr::mutate(array_fraction_normal_brain = "0") |> 
    dplyr::select(array_fraction_normal_brain,
                  
                  
                  array_mnp_predictBrain_v12.8_cal_A_IDH_HG,
                  array_mnp_predictBrain_v12.8_cal_OLIGOSARC_IDH,
                  array_mnp_predictBrain_v12.8_cal_O_IDH,
                  array_mnp_predictBrain_v12.8_cal_A_IDH_LG,
                  
                  array_mnp_predictBrain_v12.8_cal_CTRL_HEMI,
                  array_mnp_predictBrain_v12.8_cal_CTRL_HYPOTHAL,
                  
                  array_mnp_predictBrain_v12.8_cal_PA_CORT,
                  array_mnp_predictBrain_v12.8_cal_CTRL_CORPCAL,
                  array_mnp_predictBrain_v12.8_cal_DNET,
                  #array_mnp_predictBrain_v12.8_cal_GBM_MES_TYP,
                  array_mnp_predictBrain_v12.8_cal_GG,
                  array_mnp_predictBrain_v12.8_cal_PGNT,
                  array_mnp_predictBrain_v12.8_cal_HGG_F
    )
  
  
  stopifnot(colnames(exp.ctrl) == colnames(exp.diluted))
  stopifnot(colnames(exp.ctrl) == colnames(exp.tumor))
  
  exp <- rbind(exp.ctrl, exp.diluted, exp.tumor)
  rm(exp.ctrl, exp.diluted, exp.tumor)
  
  
  exp <- exp |> 
    dplyr::rename_with(~ gsub("array_mnp_predictBrain_v12.8_cal_","cal. P ",.x)) |> 
    tibble::column_to_rownames('array_fraction_normal_brain') |> 
    round(3) |> 
    tibble::rownames_to_column('fraction_normal') |> 
    dplyr::mutate(fraction_normal = as.numeric(fraction_normal)) |> 
    dplyr::arrange(fraction_normal)

  write.table(exp, file=paste0("/tmp/dilution_series_",tumor,"_",normal,".txt"), sep="\t", quote = F, row.names = F)
  
}



## export them all ----



exporttable('0017-R3', 's107')
exporttable('0017-R3', 's108')
exporttable('0017-R3', 's110')
exporttable('0017-R3', 's112')

exporttable('0008-R2', 's107')
exporttable('0008-R2', 's108')
exporttable('0008-R2', 's110')
exporttable('0008-R2', 's112')

exporttable('0121-R3', 's107')
exporttable('0121-R3', 's108')
exporttable('0121-R3', 's110')
exporttable('0121-R3', 's112')

exporttable('0054-R3', 's107')
exporttable('0054-R3', 's108')
exporttable('0054-R3', 's110')
exporttable('0054-R3', 's112')


