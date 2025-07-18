#!/usr/bin/env R 

# load ----



source('scripts/load_constants.R')
source('scripts/load_functions.R')
#source('scripts/load_palette.R')
#source('scripts/load_themes.R')



if(!exists('metadata.cg_probes.epic')) {
  print("First loading probe annotations")
  source('scripts/load_probe_annotations.R')
  print("Done")
}



# EPIC: all hq ----


data.mvalues.hq_samples <- readRDS("cache/mvalues/mvalues_hq.Rds") |>  # former "cache/mvalues.HQ_samples.Rds"
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_SAMPLES)
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
    return(.)
  })()

#data.mvalues.mask.hq_samples <- readRDS("cache/masks_hq/masks_hq.Rds") |> 
#  (function(.) {
#    print(dim(.))
#    assertthat::assert_that(ncol(.) == CONST_N_SAMPLES)
#    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
#    return(.)
#  })()

#stopifnot(colnames(data.mvalues.hq_samples) == colnames(data.mvalues.mask.hq_samples))
#stopifnot(rownames(data.mvalues.hq_samples) == rownames(data.mvalues.mask.hq_samples))




## probe table - mask / detP counts ----


fn <- "cache/data.mvalues.probes.init.Rds"
#file.remove(fn)

if(file.exists(fn)) {
  
  message(paste0("loading '",fn,"' from cache"))
  tmp <- readRDS(fn)
  
} else {
  
  tmp <- load_masks_hq() |> 
    is.na() |> 
    rowSums() |> 
    as.data.frame() |> 
    dplyr::rename(n_na = 1) |> 
    dplyr::mutate(detP_good_probe = n_na == 0) |> 
    tibble::rownames_to_column('probe_id') |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
      return(.)
    })()
  
  
  saveRDS(tmp, file=fn)
  
}


data.mvalues.probes <- tmp |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
    return(.)
  })()

rm(tmp, fn)








## probe table - probe annotation ----


data.mvalues.probes <- data.mvalues.probes |> 
  dplyr::left_join(metadata.cg_probes.epic, by=c('probe_id'='probe_id'), suffix=c('','')) |>  # + annotation of ALL probes
  assertr::verify(MASK_general == F) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
    return(.)
  })()



## filter for good probes ----


data.mvalues.good_probes <- data.mvalues.probes |> 
  dplyr::filter(detP_good_probe) |> 
  dplyr::filter(grepl("^cg", probe_id)) |> # non CG probes https://knowledge.illumina.com/microarray/general/microarray-general-troubleshooting-list/000005501
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
    return(.)
  })() |> 
  assertr::verify(MASK_general == F) |> 
  dplyr::pull(probe_id)



## medians across hq samples ----


tmp <- readRDS(file="cache/analysis_probe_median_statistics.Rds") |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
    return(.)
  })()


data.mvalues.probes <- data.mvalues.probes |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(median.beta.primary.glass_od))


rm(tmp)



## DMP outcomes ----


### prim rec ----


fn <- "cache/analysis_differential__primary_recurrence__partial_paired_nc__stats.Rds"

if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, t, adj.P.Val, logFC) |> 
    dplyr::rename_with(~paste0("DMP__primary_recurrence__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DMP result primary - recurrence is missing")
}

rm(fn)



### prim rec + quality (PC1) ----


fn <- "cache/analysis_differential__primary_recurrence__PC1__partial_paired_nc__stats.Rds"

if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, t, logFC, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__primary_recurrence__pp_nc_PC1__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DMP result primary - recurrence is missing")
}

rm(fn)




### prim rec INTENSITIES ----

# 
# fn <- "cache/analysis_differential_intensities__primary_recurrence__partial_paired_nc__stats.Rds"
# if(file.exists(fn)) {
#   
#   tmp <- readRDS(fn) |> 
#     dplyr::select(probe_id, t, adj.P.Val) |> 
#     dplyr::rename_with(~paste0("DPI__primary_recurrence__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
#     (function(.) {
#       print(dim(.))
#       assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
#       return(.)
#     })()
#   
#   
#   data.mvalues.probes <- data.mvalues.probes |> 
#     dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
#   
#   rm(tmp)
#   
# } else {
#   warning("DMP result primary - recurrence is missing")
# }
# 
# rm(fn)
# 
# 
# ### prim rec INTENSITIES non-log ----
# 
# 
# fn <- "cache/analysis_differential_intensities__primary_recurrence__partial_paired_nc__non_log__stats.Rds"
# if(file.exists(fn)) {
#   
#   tmp <- readRDS(fn) |> 
#     dplyr::select(probe_id, t, adj.P.Val) |> 
#     dplyr::rename_with(~paste0("DPI__primary_recurrence__pp_nc__non_log__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
#     (function(.) {
#       print(dim(.))
#       assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
#       return(.)
#     })()
#   
#   
#   data.mvalues.probes <- data.mvalues.probes |> 
#     dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
#   
#   rm(tmp)
#   
# } else {
#   warning("DMP result primary - recurrence is missing")
# }
# 
# rm(fn)
# 



### g3 g2 ----


fn <- "cache/analysis_differential__g2_g3__partial_paired_nc__stats.Rds"
if(file.exists(fn)) {
  
  tmp <-  readRDS(fn) |> 
    dplyr::select(probe_id, t, adj.P.Val, logFC) |> 
    dplyr::rename_with(~paste0("DMP__g2_g3__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))
  
  rm(tmp)
  
} else {
  warning("DMP result g2 - g3 is missing")
}

rm(fn)



### g3 g2 + quality (PC1) ----
#### pct detP

fn <- "cache/analysis_differential__g2_g3__PC1__partial_paired_nc__stats.Rds"

if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, t, adj.P.Val, logFC) |> 
    dplyr::rename_with(~paste0("DMP__g2_g3__pp_nc_PC1__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DMP result primary - recurrence is missing")
}

rm(fn)



### g3 g2 INTENSITIES ----

# 
# fn <- "cache/analysis_differential_intensities__g2_g3__partial_paired_nc__stats.Rds"
# if(file.exists(fn)) {
#   
#   tmp <-  readRDS(fn) |> 
#     dplyr::select(probe_id, t, adj.P.Val) |> 
#     dplyr::rename_with(~paste0("DPI__g2_g3__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
#     (function(.) {
#       print(dim(.))
#       assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
#       return(.)
#     })()
#   
#   data.mvalues.probes <- data.mvalues.probes |> 
#     dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))
#   
#   rm(tmp)
#   
# } else {
#   warning("DMP result g2 - g3 is missing")
# }
# 
# rm(fn)
# 


### g3 g2 INTENSITIES non log ----

# 
# fn <- "cache/analysis_differential_intensities__g2_g3__partial_paired_nc__non_log__stats.Rds"
# if(file.exists(fn)) {
#   
#   tmp <-  readRDS(fn) |> 
#     dplyr::select(probe_id, t, adj.P.Val) |> 
#     dplyr::rename_with(~paste0("DPI__g2_g3__pp_nc__non_log__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
#     (function(.) {
#       print(dim(.))
#       assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
#       return(.)
#     })()
#   
#   data.mvalues.probes <- data.mvalues.probes |> 
#     dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))
#   
#   rm(tmp)
#   
# } else {
#   warning("DMP result g2 - g3 is missing")
# }
# 
# rm(fn)



### PC1 + PC2 + PC3 ----


fn <- "cache/analysis_differential__PC1_PC2_PC3.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, starts_with("t"), starts_with("adj.P.Val")) |> 
    dplyr::rename_with(~paste0("DMP__PC1_PC2_PC3_multivariate__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DMP result primary - recurrence is missing")
}

rm(fn)



### Oligosarcoma A_IDH_HG ----



fn <- "cache/analysis_differential__Oligosarcoma_A_IDH_HG__PC1__stats.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    #dplyr::select(probe_id, t, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__GLASS_OD__oligosarcoma__A_IDH_HG__PC1__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DMP result primary - recurrence is missing")
}

rm(fn)



### GLASS-NL ----
#### prim rec naive GLASS-NL ----


fn <- "cache/analysis_differential__GLASS-NL__primary_recurrence__partial_paired_nc__stats.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, t, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__GLASS_NL__prim_rec__pp_nc_naive__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DMP result primary - recurrence is missing")
}

rm(fn)



#### prim rec + quality (intrinsic. PC1) GLASS-NL ----


fn <- "cache/analysis_differential__GLASS-NL__primary_recurrence__PC1__partial_paired_nc__stats.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, t, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__GLASS_NL__prim_rec__pp_nc_PC1__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DMP result primary - recurrence is missing")
}

rm(fn)






#### g2 vs. g3,g4 + naive GLASS-NL ----


fn <- "cache/analysis_differential__GLASS-NL__g2_g3_4__partial_paired_nc__stats.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, t, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__GLASS_NL__g2_g3.4__pp_nc_naive__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DMP result primary - recurrence is missing")
}

rm(fn)



#### g2 vs. g3,g4 + quality (PC1) GLASS-NL ----


fn <- "cache/analysis_differential__GLASS-NL__g2_g3_4__PC1__partial_paired_nc__stats.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, t, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__GLASS_NL__g2_g3.4__pp_nc_PC1__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DMP result primary - recurrence is missing")
}

rm(fn)






#### g2,g3 vs. g4 + naive GLASS-NL ----


fn <- "cache/analysis_differential__GLASS-NL__g2_3_g4__partial_paired_nc__stats.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, t, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__GLASS_NL__g2.3_g4__pp_nc_naive__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DMP result primary - recurrence is missing")
}

rm(fn)



#### g2,g3 vs. g4 + quality (PC1) GLASS-NL ----


fn <- "cache/analysis_differential__GLASS-NL__g2_3_g4__PC1__partial_paired_nc__stats.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, t, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__GLASS_NL__g2.3_g4__pp_nc_PC1__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DMP result primary - recurrence is missing")
}

rm(fn)




## exp1 - exp8 ----



for(i in 1:9) {
  fn <- paste0("cache/analysis_differential__GLASS-NL__exp",i,"__stats.Rds")
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, t, adj.P.Val, logFC) |> 
    dplyr::rename_with(~paste0("DMP__GLASS_NL__exp",i,"__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
}









### PCs multivariate 1st 8 components ----


fn <- "cache/analysis_differential__PCs__partial_paired_nc__stats.Rds"

if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, ends_with("_t") | ends_with("_adj.P.Val")) |> 
    dplyr::rename_with(~paste0("DMP__PCs__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DMP result PCs is missing")
}

rm(fn)




### QC ----


fns <- c("analysis_differential__pct_detP_signi__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_STAINING_DNP_High_Red_smaller_9000_3000__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_STAINING_DNP_Bkg_Red_larger_750_1500__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5__partial_paired_nc__stats.Rds",
         "analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5__partial_paired_nc__stats.Rds"
         # "analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5__partial_paired_nc__stats.Rds", 
         # "analysis_differential__array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4__partial_paired_nc__stats.Rds",
         # "analysis_differential__array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6__partial_paired_nc__stats.Rds"
         # 
         )


for (fn in fns) {
  txt <- gsub("array_qc_","mnp_qc_",gsub("^.+_differential__(.+)__partial_pai.+","\\1",fn))
  print(txt)
  
  fn <- paste0("cache/", fn)
  
  if(file.exists(fn)) {
     
     tmp <- readRDS(fn) |> 
       dplyr::select(probe_id, t, adj.P.Val) |> 
       dplyr::rename_with(~paste0("DMP__",txt,"__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
       (function(.) {
         print(dim(.))
         assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
         return(.)
       })()
     
     data.mvalues.probes <- data.mvalues.probes |> 
       dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
     
    rm(tmp)
   
  } else {
   warning(paste0("DMP result qc: ",txt," is missing"))
  }
  
  rm(fn)
}





### ewastools ----


fns <- c(
         #"analysis_differential__ewastools_Restoration__partial_paired_nc__stats.Rds",
         #"analysis_differential__ewastools_Staining.Green__partial_paired_nc__stats.Rds",
         #"analysis_differential__ewastools_Staining.Red__partial_paired_nc__stats.Rds",
         #"analysis_differential__ewastools_Extension.Green__partial_paired_nc__stats.Rds",
         #"analysis_differential__ewastools_Extension.Red__partial_paired_nc__stats.Rds",
         #"analysis_differential__ewastools_Hybridization.High.Medium__partial_paired_nc__stats.Rds",
         #"analysis_differential__ewastools_Hybridization.Medium.Low__partial_paired_nc__stats.Rds",
         #"analysis_differential__ewastools_Bisulfite.Conversion.I.Green__partial_paired_nc__stats.Rds",
         #"analysis_differential__ewastools_Bisulfite.Conversion.I.Red__partial_paired_nc__stats.Rds",
         #"analysis_differential__ewastools_Bisulfite.Conversion.II__partial_paired_nc__stats.Rds",
         "analysis_differential__ewastools_Specificity.I.Green__partial_paired_nc__stats.Rds",
         "analysis_differential__ewastools_Specificity.I.Red__partial_paired_nc__stats.Rds",
         "analysis_differential__ewastools_Non.polymorphic.Green__partial_paired_nc__stats.Rds",
         #"analysis_differential__ewastools_Specificity.II__partial_paired_nc__stats.Rds",
         "analysis_differential__ewastools_Non.polymorphic.Red__partial_paired_nc__stats.Rds"
         #"analysis_differential__ewastools_Target.Removal.1__partial_paired_nc__stats.Rds",
         #"analysis_differential__ewastools_Target.Removal.2__partial_paired_nc__stats.Rds"
         )

for (fn in fns) {
  txt <- gsub("^.+l__ewastools_(.+)__partial_pai.+","\\1",fn)
  print(txt)
  
  fn <- paste0("cache/", fn)
  
  if(file.exists(fn)) {
    
    tmp <- readRDS(fn) |> 
      dplyr::select(probe_id, t) |> 
      dplyr::rename_with(~paste0("DMP__ewastools_",txt,"__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
      (function(.) {
        print(dim(.))
        assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
        return(.)
      })()
    
    data.mvalues.probes <- data.mvalues.probes |> 
      dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
    
    rm(tmp)
    
  } else {
    warning(paste0("DMP result ewastools:",txt," is missing"))
  }
  
  rm(fn)
}

rm(txt, fns)



### CGC aka AcCGAP ----


fn <- "cache/analysis_differential__AcCGAP__partial_paired_nc__stats.Rds"

if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, t, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__AcCGAP__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
    
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DMP result AcCGAP is missing")
}

rm(fn)




### GLASS-NL signature ----


fn <- "cache/analysis_differential__GLASS_NL_sig__partial_paired_nc__stats.Rds"

if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, t, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__GLASS_NL_signature__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DMP result GLASS-NL signature is missing")
}

rm(fn)




### Tissue or FFPE ----


fn <- "cache/analysis_differential__ffpe_or_ff__partial_paired_nc__stats.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, t, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__FFPE_or_FF__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    dplyr::select(probe_id, ends_with("__logFC"), ends_with("__t") ,  ends_with("__adj.P.Val"), ends_with("__AveExpr")) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
    
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))
  
  rm(tmp)
  
} else {
  warning("DMP result FFPE | FF PP is missing")
}

rm(fn)





### FFPE & Freezer Decay time multivar PP ----



fn <- "cache/analysis_differential__ffpe_freezer_decay-time_multivar__partial_paired_nc__stats.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(probe_id, starts_with("t_")) |> 
    dplyr::rename_with(~paste0("DMP__FFPE_and_freezer_decay_time__multivar__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))
  
  rm(tmp)
  
} else {
  warning("DMP result FFP & Freezer decay multivariate PP time is missing")
}

rm(fn)





### Freezer Decay time PP ----


fn <- "cache/analysis_differential__freezer-decay-time__partial_paired_nc__stats.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |>
    #dplyr::select(probe_id, t, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__freezer_decay_time__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |>
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  data.mvalues.probes <- data.mvalues.probes |>
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))
  
  rm(tmp)
  
} else {
  warning("DMP result Freezer decay PP time is missing")
}

rm(fn)






### FFPE Decay time PP ----


fn <- "cache/analysis_differential__ffpe-decay-time__partial_paired_nc__stats.Rds"

if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |>
    dplyr::select(probe_id, t, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__FFPE_decay_time__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |>
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  data.mvalues.probes <- data.mvalues.probes |>
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))
  
  rm(tmp)
  
} else {
  warning("DMP result FFPE decay PP time is missing")
}

rm(fn)



### log1p FFPE Decay time PP ----


fn <- "cache/analysis_differential__log_ffpe_freezer_decay-time__partial_paired_nc__stats.Rds"

if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |>
    dplyr::select(probe_id, t, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__FFPE_log1p__decay_time__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |>
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  data.mvalues.probes <- data.mvalues.probes |>
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))
  
  rm(tmp)
  
} else {
  warning("DMP result log1p FFPE decay PP time is missing")
}

rm(fn)







### FFPE Decay time PP intensities combi ----

# 
# fn <- "cache/analysis_differential_intensities__ffpe-decay-time__partial_paired_nc__stats.Rds"
# if(file.exists(fn)) {
#   
#   tmp <- readRDS(fn) |> 
#     dplyr::select(probe_id, t, adj.P.Val) |> 
#     dplyr::rename_with(~paste0("DPI__FFPE_decay_time__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
#     (function(.) {
#       print(dim(.))
#       assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
#       return(.)
#     })()
#   
#   data.mvalues.probes <- data.mvalues.probes |> 
#     dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))
#   
#   rm(tmp)
#   
# } else {
#   warning("DMP result FFPE decay PP time is missing")
# }
# 
# rm(fn)




### FFPE Decay time PP intensities methylated ----

# 
# fn <- "cache/analysis_differential_methylated_intensities__ffpe-decay-time__partial_paired_nc__stats.Rds"
# if(file.exists(fn)) {
#   
#   tmp <- readRDS(fn) |> 
#     dplyr::select(probe_id, t, adj.P.Val) |> 
#     dplyr::rename_with(~paste0("DMPI__FFPE_decay_time__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
#     (function(.) {
#       print(dim(.))
#       assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
#       return(.)
#     })()
#   
#   data.mvalues.probes <- data.mvalues.probes |> 
#     dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))
#   
#   rm(tmp)
#   
# } else {
#   warning("DMP result FFPE decay PP time is missing")
# }
# 
# rm(fn)




### FFPE Decay time PP intensities unmethylated ----

# 
# fn <- "cache/analysis_differential_unmethylated_intensities__ffpe-decay-time__partial_paired_nc__stats.Rds"
# if(file.exists(fn)) {
#   
#   tmp <- readRDS(fn) |> 
#     dplyr::select(probe_id, t, adj.P.Val) |> 
#     dplyr::rename_with(~paste0("DUPI__FFPE_decay_time__pp_nc__", .x), .cols=!matches("^probe_id$",perl = T)) |> 
#     (function(.) {
#       print(dim(.))
#       assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
#       return(.)
#     })()
#   
#   data.mvalues.probes <- data.mvalues.probes |> 
#     dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))
#   
#   rm(tmp)
#   
# } else {
#   warning("DMP result FFPE decay PP time is missing")
# }
# 
# rm(fn)
# 

# plot(c(0, nrow(data.mvalues.probes)), y=c(-10,5), type="n")
# lines(sort(data.mvalues.probes$DMPI__FFPE_decay_time__pp_nc__t))
# lines(sort(data.mvalues.probes$DUPI__FFPE_decay_time__pp_nc__t))
# lines(sort(data.mvalues.probes$DPI__FFPE_decay_time__pp_nc__t))
# lines(sort(data.mvalues.probes$DMP__FFPE_decay_time__pp_nc__t))
# 
# 
# sum(is.na(data.mvalues.probes$DMPI__FFPE_decay_time__pp_nc__t))
# sum(is.na(data.mvalues.probes$DUPI__FFPE_decay_time__pp_nc__t))
# sum(is.na(data.mvalues.probes$DPI__FFPE_decay_time__pp_nc__t))
# sum(is.na(data.mvalues.probes$DMP__FFPE_decay_time__pp_nc__t))

# 
# plt <- data.mvalues.probes |> 
#   dplyr::filter(!is.na(DMPI__FFPE_decay_time__pp_nc__t)) |> 
#   dplyr::select(DMPI__FFPE_decay_time__pp_nc__t, DUPI__FFPE_decay_time__pp_nc__t, DPI__FFPE_decay_time__pp_nc__t, DMP__FFPE_decay_time__pp_nc__t) |> 
#   dplyr::rename(`Intensity (M)` = DMPI__FFPE_decay_time__pp_nc__t) |> 
#   dplyr::rename(`Intensity (UM)` = DUPI__FFPE_decay_time__pp_nc__t) |> 
#   dplyr::rename(`Intensity (M + UM)` = DPI__FFPE_decay_time__pp_nc__t) |> 
#   dplyr::rename(`M UM ratio (M-value)` = DMP__FFPE_decay_time__pp_nc__t)
# 
# c <- cor(plt)
# corrplot::corrplot(c, order="hclust")



### epiTOC2 & dnaMethyAge epiGenetic clocks ----


clocks <- list.files(path = "cache/", pattern = "*.Rds$", recursive = TRUE) |>
  data.frame(filename = _) |> 
  dplyr::mutate(filename = paste0("cache/", filename)) |> 
  dplyr::filter(grepl("stats.Rds$", filename)) |> 
  dplyr::filter(grepl("(epiTOC2|dnaMethyAge|RepliTali)", filename)) |> 
  dplyr::pull(filename)


for(clock in clocks) {
  clock_name <- gsub("__stats.Rds","",gsub("unpaired","up",gsub("cache/analysis_differential__","", clock)))
  print(clock_name)
  
  tmp <- readRDS(clock) |> 
    dplyr::select(probe_id, t, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DMP__",clock_name,"__", .x), .cols=!matches("^probe_id$",perl = T))  |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
      return(.)
    })()
  
  data.mvalues.probes <- data.mvalues.probes |> 
    dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))
  
  rm(tmp, clock_name)
  
}

rm(clocks)




## s[c/n]ATAC-seq bins ----


tmp <- readRDS(file="cache/analysis_ATAC_oligo.Rds")
data.mvalues.probes <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  
  dplyr::mutate(center = (Start_hg38 + End_hg38) / 2) |> 
  dplyr::left_join(tmp,
                   by=dplyr::join_by(CHR_hg38 == ATAC_oligo_chr,
                                     dplyr::closest(center >= ATAC_oligo_start),
                                     dplyr::closest(center <= ATAC_oligo_end)
                   ))  |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::mutate(center = NULL)

rm(tmp)





tmp <- readRDS(file="cache/analysis_ATAC_astro.Rds")
data.mvalues.probes <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  
  dplyr::mutate(center = (Start_hg38 + End_hg38) / 2) |> 
  dplyr::left_join(tmp,
                   by=dplyr::join_by(CHR_hg38 == ATAC_astro_chr,
                                     dplyr::closest(center >= ATAC_astro_start),
                                     dplyr::closest(center <= ATAC_astro_end)
                   ))  |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::mutate(center = NULL)

rm(tmp)





# data.mvalues.probes |> 
#   dplyr::mutate(exp = !is.na(ATAC_oligo_data_per_bin_per_base)) |> 
#   dplyr::pull(exp) |> 
#   table()
# 
# data.mvalues.probes |> 
#   dplyr::mutate(exp = !is.na(ATAC_astro_data_per_bin_per_base)) |> 
#   dplyr::pull(exp) |> 
#   table()
# 
# data.mvalues.probes |> 
#   dplyr::mutate(exp = !is.na(ATAC_astro_data_per_bin_per_base) &  !is.na(ATAC_oligo_data_per_bin_per_base)) |> 
#   dplyr::pull(exp) |> 
#   table()






# 450K ----
## AD: from Beta-value exported files ----
#' https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-019-0672-7


if(!exists('ad_bmc_clin_epi.metadata.array_samples')) {
  source('scripts/load_AD_BMC_Clin_Epi.R')
}



tmp.1 <- read.delim("data/GLASS_OD/DNA Methylation - 450K arrays - Alzheimer/SampleMethFinalReport_ctrl-bkg.txt", skip=8, header=T) |> 
  tibble::column_to_rownames("TargetID") |> 
  dplyr::select(contains("_Beta")) |> 
  dplyr::mutate_all(function(x){return (  log2( `x` / (1 - `x`))  )}) |>  # Beta to M-values
  tibble::rownames_to_column("probe_id")



tmp.2 <- read.delim("data/GLASS_OD/DNA Methylation - 450K arrays - Alzheimer/SampleMethFinalReport_nonorm.txt", skip=8, header=T) |> 
  tibble::column_to_rownames("TargetID") |> 
  dplyr::select(contains("_Beta")) |> 
  dplyr::mutate_all(function(x){return (  log2( `x` / (1 - `x`))  )}) |>  # Beta to M-values
  tibble::rownames_to_column("probe_id")


stopifnot(tmp.1$probe_id == tmp.2$probe_id)
stopifnot(nrow(tmp.1) == nrow(tmp.2))


data.mvalues.alzheimer.dirty <- tmp.1 |> 
  dplyr::left_join(tmp.2, by=c('probe_id'='probe_id'), suffix=c('','')) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::rename_with(~ gsub("\\.AVG_Beta$","",.x)) |> 
  dplyr::rename_with(~ gsub("^X","",.x)) |> 
  dplyr::rename_with(~ gsub("\\.","_",.x)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == (26 + 12))
    assertthat::assert_that(all(colnames(.) %in% ad_bmc_clin_epi.metadata.array_samples$DNAm_id))

    return(.)
  })() |> 
  # dplyr::select(ad_bmc_clin_epi.metadata.array_samples |> dplyr::filter(is.na(reason_excluded)) |>  dplyr::pull(DNAm_id)) |> 
  # (function(.) {
  #   print(dim(.))
  #   assertthat::assert_that(ncol(.) == (26 + 12 - 1))
  #   
  #   return(.)
  # })() |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (485577))

    .$NAs <- rowSums(is.na(.))
    . <- .[.$NAs == 0,]
    .$NAs <- NULL

    print(dim(.))
    assertthat::assert_that(sum(is.na(.)) == 0)
    assertthat::assert_that(nrow(.) == (470392))  #470691
    return(.)
  })()



rm(tmp.1, tmp.2)




