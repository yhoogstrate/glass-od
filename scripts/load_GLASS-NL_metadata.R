#!/usr/bin/env R

# GLASS-NL ----


source('scripts/load_functions.R')
source('scripts/load_constants.R')


## idats ----


glass_nl.metadata.array_samples <-  list.files(path = "data/GLASS_NL/Methylation/Methylation Array Data/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("data/GLASS_NL/Methylation/Methylation Array Data/", filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (658)) # no idea how there can be 658 idats in here?
    return(.)
  })() |>
  assertr::verify(file.exists(filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", filename)) |>
  dplyr::mutate(channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", filename)) |>
  tidyr::pivot_wider(id_cols = sentrix_id, names_from = channel, values_from = c(filename)) |>
  dplyr::rename(channel_green = Grn) |>
  dplyr::rename(channel_red = Red) |> 
  assertr::verify(!is.na(sentrix_id)) |> 
  assertr::verify(!is.na(channel_green)) |> 
  assertr::verify(!is.na(channel_red)) |> 
  dplyr::rename_with( ~ paste0("array_", .x)) 



## study identifier: Sample_Name & resection metadata ----


tmp <- read.csv("data/GLASS_NL/Metadata/Samples/Master Datasheet_ALL METHODS_27012023.csv") |>
  dplyr::select( Surgery_ID, GLASS_ID,  Sample_Name, Sample_Sex, Sample_Type, Resectie, Sample_ID, Recurrent_Select.Meth, Matched_Pair.Meth) |> 
  dplyr::rename(array_sentrix_id = Sample_ID) |> 
  dplyr::rename(resection_number = Resectie) |> 
  dplyr::filter(!is.na(array_sentrix_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (235)) 
    return(.)
  })() |> 
  dplyr::mutate(patient_id = as.factor(gsub("^.+_[0]*([0-9]+)$", "\\1", GLASS_ID)))


glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |>
  dplyr::filter(!is.na(Sample_Name)) |>  # exclude hundreds present in directory but absent in metadata
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (235)) 
    return(.)
  })()


rm(tmp)




## WHO classes ----


tmp <- read.csv("data/GLASS_NL/Metadata/Samples/WHOclassification_03052022.csv")  |> 
  dplyr::select(Surgery_ID, WHO_Classification2021) |> # join on Surgery_ID
  dplyr::mutate(WHO_Classification2021 = ifelse(WHO_Classification2021 == "Therapy Effects", as.character(NA), WHO_Classification2021)) |> 
  assertr::verify(!is.na(Surgery_ID)) |> 
  dplyr::filter(Surgery_ID %in% glass_nl.metadata.array_samples$Surgery_ID) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (235)) 
    return(.)
  })()



glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('Surgery_ID'='Surgery_ID'), suffix=c('','')) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (235)) 
    return(.)
  })()


rm(tmp)


## Percentage detP probes ----
#' from: scripts/analysis_percentage_detP_probes.R


tmp <- read.table("output/tables/percentage_detP_probes.txt") |> 
  assertr::verify(!is.na(array_percentage.detP.signi) & is.numeric(array_percentage.detP.signi)) |> 
  assertr::verify(glass_nl.metadata.array_samples$array_sentrix_id %in% array_sentrix_id)


glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_percentage.detP.signi)) 

rm(tmp)


## QC PCA outlier ----

tmp <- readRDS('cache/unsupervised_qc_outliers_all.Rds') |> 
  assertr::verify(!is.na(array_qc.pca.comp1)) |> 
  assertr::verify(is.numeric(array_qc.pca.comp1)) |> 
  assertr::verify(!is.na(array_qc.pca.detP.outlier)) |> 
  assertr::verify(is.logical(array_qc.pca.detP.outlier)) |> 
  assertr::verify(glass_nl.metadata.array_samples$array_sentrix_id %in% array_sentrix_id) |> 
  dplyr::filter(array_sentrix_id %in% glass_nl.metadata.array_samples$array_sentrix_id)


glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_qc.pca.comp1)) |> 
  assertr::verify(!is.na(array_qc.pca.detP.outlier))

rm(tmp)



## QC Purity PCA outlier 2 ----
#' PCA was performed on only 202 HQ arrays
#' these outliers in PC3 (val >300) were CONTR_ and sometimes difficultly classifiable

PC3_low_purity_samples <- c("203175700013_R08C01", "203189480016_R03C01", "203189480016_R05C01", "203986510092_R04C01", "203986510125_R04C01", "203989100024_R08C01",
                            "203989100035_R02C01", "203989100096_R03C01", "203989100096_R05C01", "203989100142_R02C01", "203991400003_R01C01", "203991400003_R06C01",
                            "204073520032_R07C01", "204073570005_R02C01", "203519500055_R03C01", "203519500055_R05C01")

glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |>
  dplyr::mutate(qc.pca.pc3purity.outlier = array_sentrix_id %in% PC3_low_purity_samples) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 235)
    return(.)
  })()

rm(PC3_low_purity_samples)





## ewastools ----


tmp <- read.table("output/tables/ewastools.txt") |> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(glass_nl.metadata.array_samples$array_sentrix_id %in% array_sentrix_id)


glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))


rm(tmp)


## Heidelberg 12.8 reportBrain files ----


tmp <- c(
  list.files(path = "data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/",
             pattern = "*_scores_cal.csv", recursive = TRUE)) |>
  data.frame(array_mnp_predictBrain_v12.8_filename = _) |>
  dplyr::mutate(array_mnp_predictBrain_v12.8_filename = paste0("data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/", array_mnp_predictBrain_v12.8_filename)) |>
  dplyr::mutate(array_mnp_predictBrain_v12.8_version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", array_mnp_predictBrain_v12.8_filename)) |> 
  assertr::verify(array_mnp_predictBrain_v12.8_version == "v12.8") |> 
  dplyr::mutate(array_mnp_predictBrain_v12.8_version = NULL) |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+/([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", array_mnp_predictBrain_v12.8_filename)) |>
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |>  # only one version per sentrix_id desired
  assertr::verify(array_sentrix_id %in% glass_nl.metadata.array_samples$array_sentrix_id) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_reportBrain_csv(array_mnp_predictBrain_v12.8_filename, paste0("array_mnp_predictBrain_v12.8_"))) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(array_mnp_predictBrain_v12.8_filename = NULL) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 235)
    return(.)
  })() |> 
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_class)) |>  # version is hardcoded here
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_A_IDH_LG)) |>
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_A_IDH_HG)) |>
  dplyr::mutate(array_A_IDH_HG__A_IDH_LG_lr_v12.8 = log(array_mnp_predictBrain_v12.8_cal_A_IDH_HG / array_mnp_predictBrain_v12.8_cal_A_IDH_LG))
  #dplyr::mutate(A_IDH_HG__A_IDH_LG_lr_neat = log(    (mnp_predictBrain_v12.8_cal_A_IDH_HG / (1-mnp_predictBrain_v12.8_cal_A_IDH_HG))    /     (mnp_predictBrain_v12.8_cal_A_IDH_LG / (1-mnp_predictBrain_v12.8_cal_A_IDH_LG))   ))



glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_predictBrain_v12.8_cal_class)) |> 
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_LG_lr_v12.8)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 235)
    return(.)
  })()

rm(tmp)



## Heidelberg 12.8 Frozen ~ FFPE status ----


tmp <- list.files(path = "data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/",  pattern = "_ffpe_frozen.txt", recursive = TRUE) |> 
  data.frame(array_mnp_QC_v12.8_FrozenFFPEstatus_table = _) |> 
  dplyr::mutate(array_mnp_QC_v12.8_FrozenFFPEstatus_table = paste0("data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/", array_mnp_QC_v12.8_FrozenFFPEstatus_table)) |>
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_FrozenFFPEstatus_table(array_mnp_QC_v12.8_FrozenFFPEstatus_table, "array_mnp_QC_v12.8_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(array_mnp_QC_v12.8_FrozenFFPEstatus_table = NULL) |> 
  assertr::verify(!is.na(array_sentrix_id)) |> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(array_sentrix_id %in% glass_nl.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 235)
    return(.)
  })()



glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_QC_v12.8_predicted_array_type)) |> 
  assertr::verify(!is.na(array_mnp_QC_v12.8_predicted_sample_type)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 235)
    return(.)
  })() 

rm(tmp)



## Heidelberg 12.8 CNVP segment files ----


tmp <- query_mnp_12.8_CNVP_segment_csv(
  "data/GLASS_NL/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/",
  235,
  glass_nl.metadata.array_samples$array_sentrix_id,
  "array_mnp_CNVP_v12.8_v5.2_")


glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_CNVP_v12.8_v5.2_CNVP_segments))  |> 
  assertr::verify(!is.na(array_mnp_CNVP_v12.8_v5.2_CNVP_version)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 235)
    return(.)
  })() 

rm(tmp)






## RF purity calls ----


tmp <- read.table("data/GLASS_NL/Methylation/Analysis/RFpurity/purities_RFpurity.txt") |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+/([0-9]{10,}_R[0-9]+C[0-9]+).+\\.idat$","\\1",fn)) |> 
  assertr::verify(!is.na(array_sentrix_id)) |> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  dplyr::select(array_sentrix_id, absolute,  estimate) |> 
  assertr::verify(glass_nl.metadata.array_samples$array_sentrix_id %in% array_sentrix_id) |> 
  dplyr::rename(array_RFpurity.absolute = absolute) |> 
  dplyr::rename(array_RFpurity.estimate = estimate)


glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_RFpurity.absolute)) |> 
  assertr::verify(!is.na(array_RFpurity.estimate))

rm(tmp)


## RNA signature ----
### cell-cycling ----
#' let wel - deze is supervised na primary - recurrent diff


tmp <- dplyr::inner_join(
  read.csv("data/GLASS_NL/Metadata/Cleaned_clinical/metadata_2022/Samplesheet_GLASS_RNAseq__ALL.csv") |> 
    dplyr::rename(rnaseq_sample_id = GS_ID) |> 
    dplyr::select(rnaseq_sample_id, Sample_Name),
  read.csv('data/GLASS_NL/Methylation/Metadata/Datasheet4.csv') |> 
    dplyr::rename(array_sentrix_id = Sample_ID) |> 
    dplyr::select(array_sentrix_id, Sample_Name),
  by=c('Sample_Name'='Sample_Name')
) |> 
  dplyr::mutate(Sample_Name = NULL) |> 
  dplyr::inner_join(
    readRDS("cache/glass-nl_transcriptional.signatures.Rds") |> 
      dplyr::rename(rnaseq_cell.cycling.signature = lts.up1) |>  # x-checked, up1
      dplyr::rename(rnaseq_ECM.signature = lts.up2) |>  # x-checked, up2
      dplyr::rename(rnaseq_ACdown.signature = lts.down) |>  # x-checked, up2
      dplyr::rename(rnaseq_sample_id = genomescan.sid) |> 
      dplyr::select(rnaseq_sample_id, rnaseq_cell.cycling.signature, rnaseq_ECM.signature, lts.up3, rnaseq_ACdown.signature),
    by=c('rnaseq_sample_id'='rnaseq_sample_id')
  ) |> 
  dplyr::mutate(rnaseq_sample_id = NULL)



glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(sum(!is.na(.$rnaseq_cell.cycling.signature)) == 167)
    return(.)
  })()


rm(tmp)



## ++ below: re-build because mvalue normalisation ++ ----

## Median methylation levels ----


tmp <- readRDS("cache/analysis_median_methylation.Rds") |> 
  dplyr::rename_with( ~ paste0("array_", .x)) |> 
  assertr::verify(!is.na(array_median.overall.methylation)) |> 
  assertr::verify(!is.na(array_median.glass_nl_supervised.methylation)) |> 
  dplyr::filter(array_sentrix_id %in% glass_nl.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_NL_INCLUDED_SAMPLES)
    return(.)
  })()


glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(qc.pca.pc3purity.outlier == T | array_qc.pca.detP.outlier == T | !is.na(array_median.overall.methylation)) 



rm(tmp)


## A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV ----


tmp <- readRDS(file="cache/analysis_A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV.Rds") |> 
  dplyr::rename(`array_A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV_v12.8` = `array_A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV`) |> 
  assertr::verify(!is.na(`array_A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV_v12.8`)) |> 
  dplyr::filter(array_sentrix_id %in% glass_nl.metadata.array_samples$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_NL_INCLUDED_SAMPLES)
    return(.)
  })()


glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))



## unsupervised PCA ----


tmp <- readRDS(file="cache/analysis_unsupervised_PCA_GLASS-NL_x.Rds") |>
  assertr::verify(!is.na(array_PC1)) |>
  assertr::verify(!is.na(array_PC2)) |>
  assertr::verify(!is.na(array_PC3)) |>
  assertr::verify(!is.na(array_PC202)) |>
  dplyr::filter(array_sentrix_id %in% glass_nl.metadata.array_samples$array_sentrix_id) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_NL_INCLUDED_SAMPLES)
    return(.)
  })()


glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))



# tmp <- glass_nl.metadata.array_samples |>
#   dplyr::filter(!is.na(array_PC1))

# cor(tmp$array_percentage.detP.signi , tmp$array_PC1, method="spearman")
# plot(log1p(tmp$array_percentage.detP.signi) , tmp$array_PC1, cex=0.1)
# plot(log1p(tmp$array_percentage.detP.signi) , tmp$array_PC2, cex=0.1)
# plot(log1p(tmp$array_percentage.detP.signi) , tmp$array_PC3, cex=0.1)
# 
# 
# tmp <- glass_nl.metadata.array_samples |> 
#    dplyr::filter(!is.na(array_PC1)) |> 
#    dplyr::filter(!is.na(WHO_Classification2021))
# 
# plot(as.factor(tmp$WHO_Classification2021) , tmp$array_PC1, cex=0.1)
# plot(as.factor(tmp$WHO_Classification2021) , tmp$array_PC2, cex=0.1)
# #plot(as.factor(tmp$WHO_Classification2021) , tmp$array_PC3, cex=0.1)
# #plot(as.factor(tmp$WHO_Classification2021) , tmp$array_PC4, cex=0.1)
# 
# cor(log1p(tmp$array_percentage.detP.signi) , tmp$array_PC1, cex=0.1)
# cor(log1p(tmp$array_percentage.detP.signi) , tmp$array_PC2, cex=0.1)
# cor(log1p(tmp$array_percentage.detP.signi) , tmp$array_PC3, cex=0.1)


## unsupervised PCA [GLASS-OD + GLASS-NL combi] ----


# tmp <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_GLASS-NL_combined.Rds") |> 
#   dplyr::rename_with(~ gsub("^PC","PC.GLASS_OD_NL_combined.",.x), .cols = matches("^PC[0-9]", perl = T)) |> 
#   dplyr::filter(sentrix_id %in% glass_nl.metadata.array_samples$array_sentrix_id) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == CONST_N_GLASS_NL_INCLUDED_SAMPLES)
#     return(.)
#   })()
# 
# 
# glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |> 
#   dplyr::left_join(tmp, by=c('sentrix_id'='array_sentrix_id'), suffix=c('',''))
# rm(tmp)



tmp <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_GLASS-NL_combined_no_1P19Q.Rds") |>
  dplyr::rename_with( ~ paste0("array_", .x)) |>
  dplyr::rename_with(~ gsub("^array_PC","array_PC.GLASS_OD_NL_combined_excl_1P19Q.",.x), .cols = matches("^array_PC[0-9]", perl = T)) |>
  dplyr::filter(array_sentrix_id %in% glass_nl.metadata.array_samples$array_sentrix_id) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == (CONST_N_GLASS_OD_INCLUDED_SAMPLES + CONST_N_GLASS_NL_INCLUDED_SAMPLES + 1))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_NL_INCLUDED_SAMPLES)
    return(.)
  })()


glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)




## lift-over PCA from GLASS-OD ----

fn <- "cache/load_GLASS-NL_PCA_liftover_from_GLASS-OD.Rds"
#file.remove(fn)
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    
    (function(.) {
      print(dim(.))
      assertthat::assert_that(ncol(.) == 1 + CONST_N_GLASS_OD_INCLUDED_SAMPLES) 
      return(.)
    })()

} else {
  
  data <- data.mvalues.hq_samples |> 
    tibble::rownames_to_column('probe_id') |> 
    dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
    tibble::column_to_rownames('probe_id') |> 
    dplyr::select(all_of(intersect(glass_nl.metadata.array_samples$array_sentrix_id, colnames(data.mvalues.hq_samples))))
  
  pc <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_prcomp.Rds")
  tmp <- predict(pc, t(data)) |> 
    as.data.frame() |> 
    
    (function(.) {
      print(dim(.))
      assertthat::assert_that(ncol(.) == CONST_N_GLASS_OD_INCLUDED_SAMPLES) 
      return(.)
    })() |> 
    dplyr::rename_with(~ gsub("^PC","PC__from_OD__",.x)) |> 
    tibble::rownames_to_column('array_sentrix_id')

  rm(data, pc)
  
  saveRDS(tmp, file="cache/load_GLASS-NL_PCA_liftover_from_GLASS-OD.Rds")
  
  #rm(tmp)
}


glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))

rm(tmp)



## epiTOC2 ----


tmp <- readRDS("cache/analysis_EPITOC2.Rds") |> 
  dplyr::filter(array_sentrix_id %in% glass_nl.metadata.array_samples$array_sentrix_id) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_GLASS_NL_INCLUDED_SAMPLES)
    return(.)
  })()


glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))
rm(tmp)


## RNA + proteomics identifiers ----

tmp <- read.csv("data/GLASS-NL_identifiers.csv") |> 
  dplyr::mutate(X = NULL) |> 
  dplyr::filter(!is.na(methylation.sid)) |> 
  dplyr::rename(array_sentrix_id = methylation.sid) |> 
  dplyr::rename(RNA_seq_sid = genomescan.sid) |> 
  dplyr::rename(proteomics_sid = ProtID)


glass_nl.metadata.array_samples <- glass_nl.metadata.array_samples |>
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))
rm(tmp)





