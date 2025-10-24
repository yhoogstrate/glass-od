#!/usr/bin/env

# libs ----


library(limma)
library(ggplot2)
library(patchwork)
#library(EnhancedVolcano)
#library(recursiveCorPlot)


# load data ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')



if(!exists('metadata.proteomics.glass_od')) {
  source('scripts/load_GLASS-OD_proteomics.R')
}


if(!exists('metadata.proteomics.glass_nl')) {
  source('scripts/load_GLASS-NL_proteomics.R')
}



# limma N/A ----



pp <- function(measurements, patient_ids, conditions, coeff) {
  stopifnot(length(patient_ids) == length(conditions))
  #stopifnot(is.factor(conditions))
  stopifnot(length(patient_ids) == length(measurements))
  
  metadata <- data.frame(
    patient_id = patient_ids[!is.na(measurements)],
    condition  = conditions[  !is.na(measurements)]
  ) |> 
    dplyr::mutate(patient_id = as.character(patient_id)) |> 
    dplyr::group_by(patient_id) |> 
    dplyr::mutate(is.paired = dplyr::n() == 2) |> 
    dplyr::ungroup() |> 
    
    dplyr::mutate(patient = as.factor(ifelse(is.paired, paste0("p",patient_id), "a_remainder"))) |> 
    dplyr::mutate(patient_id = NULL)
  
  measurements <- tibble::tibble(t(tibble::tibble(measurements[!is.na(measurements)])))

  tryCatch({
    design <- model.matrix(~factor(patient) + condition, data=metadata)
    fit <- limma::lmFit(measurements, design)
    fit <- limma::eBayes(fit, trend=T)
    
    #print(colnames(fit$coefficients))
    
    stats <- limma::topTable(fit,
                             n=nrow(measurements),
                             coef=coeff,
                             sort.by = "none",
                             adjust.method="fdr") |> 
      tibble::rownames_to_column('probe_id') |> 
      dplyr::mutate(probe_id = NULL) |> 
      dplyr::mutate(n = nrow(metadata))
    
    return(stats)
  },
  error = function(w) { return(NA)}, # in case there are not sufficient samples per condition (too many NA's)
  warning = function(w) { return(NA)}
  )
  
}





## correlation w/ KI67 ----





ki67.proteomics <-  data.proteomics.glass_od  |>
  as.data.frame() |> 
  tibble::rownames_to_column('gene') |> 
  dplyr::filter(grepl("KI67", gene)) |> 
  tibble::column_to_rownames('gene') |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('proteomics_id') |> 
  dplyr::filter(proteomics_id != "GLODprot_03_1521") |> # doesn't matter which one
  dplyr::filter(proteomics_id != "GLODprot_09_1522") |> # no value
  dplyr::filter(proteomics_id != "GLODprot_11_1524")
dim(ki67.proteomics)


ki67.proteomics <-  ki67.proteomics |> 
  dplyr::left_join(
    glass_od.metadata.proteomics |>  dplyr::select(proteomics_id, resection_id, proteomics_notes),
    by=c('proteomics_id'='proteomics_id')
  )
dim(ki67.proteomics)
stopifnot(duplicated(ki67.proteomics$resection_id) == 0)








ki67.if <- glass_od.metadata.array_samples |>
  dplyr::select(resection_id,isolation_id, contains("67")) |>
  dplyr::filter(!is.na(staining_KI67_md5sum) & isolation_id %in% c("0047-RX", "0047-RY") == F)

sum(duplicated(ki67.if$resection_id))

dups <- ki67.if |>
  dplyr::filter(duplicated(resection_id)) |> 
  dplyr::pull(resection_id)

ki67.if <- ki67.if |> dplyr::filter((resection_id %in% dups == F) | (resection_id %in% dups == T & !is.na(array_PC67)))

sum(duplicated(ki67.if$resection_id))




ki67.joined = ki67.proteomics |> dplyr::inner_join(ki67.if, by=c('resection_id'='resection_id'))
sum(duplicated(ki67.joined$resection_id))


ggplot(ki67.joined, aes(x=MKI67, y=staining_KI67_pos_per_detected_cells)) +
  ggpubr::stat_cor(method = "spearman") +
  geom_point() +
  theme_nature



ggplot(ki67.joined, aes(x=MKI67, y=staining_KI67_lr_pos_neg_cells)) +
  ggpubr::stat_cor(method = "spearman") +
  geom_point() +
  theme_nature




ggplot(ki67.joined, aes(x=MKI67, y=staining_KI67_pos_neg_cell_density)) +
  ggpubr::stat_cor(method = "spearman") +
  geom_point() +
  theme_nature




ggplot(ki67.joined, aes(x=MKI67, y=staining_KI67_pos_per_area_um2)) +
  ggpubr::stat_cor(method = "spearman") +
  geom_point() +
  theme_nature




cor(exp(ki67.joined$MKI67),
    ki67.joined$staining_KI67_pos_per_area_um2,
    use="complete.obs", method="spearman")



# 0.4589544



## WHO Grade ----
### naive ----


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::filter(!is.na(resection_tumor_grade)) |>  # 0053-R2 pathology was not decisve
  dplyr::select(resection_id, patient_id, starts_with("proteomics_"), starts_with("resection_"), contains("_excluded")) |> 
  dplyr::mutate(condition = factor(paste0("g", resection_tumor_grade), levels=c("g2", "g3"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 139) 
    return(.)
  })()




data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id)
data$nas <- rowSums(is.na(data))
data <- data |> 
  dplyr::filter(nas < 60) |> 
  dplyr::mutate(nas=NULL)


stopifnot(metadata$proteomics_id == colnames(data))



design <- model.matrix(~factor(condition), data=metadata)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit, trend=T)
stats.grade <- limma::topTable(fit,
                               n=nrow(data),
                               coef="factor(condition)g3",
                               sort.by = "none",
                               adjust.method="fdr") |> 
  as.data.frame() |> 
  tibble::rownames_to_column('protein_id')

sum(stats.grade$adj.P.Val < 0.01)


rm(design, fit, data, metadata)



saveRDS(stats.grade, file="cache/analysis_differential_proteomics__GLASS-OD__stats.grade.naive.Rds")


rm(stats.grade)




### pat corrected ----



metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::filter(!is.na(resection_tumor_grade)) |>  # 0053-R2 pathology was not decisve
  
  dplyr::select(proteomics_id, patient_id, resection_tumor_grade) |> 
  dplyr::filter(!is.na(resection_tumor_grade)) |> 
  dplyr::mutate(resection_tumor_grade = factor(ifelse(resection_tumor_grade == 2, "Grade2", "Grade3"), levels=c("Grade2", "Grade3"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 139) 
    return(.)
  })()


data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id) 



stats.grade.pat.corrected <- data |> 
  tibble::rownames_to_column("protein_id") |> 
  dplyr::select("protein_id") |> 
  dplyr::mutate(tmp = pbapply::pbapply(data, 1, pp, metadata$patient_id, metadata$resection_tumor_grade, "conditionGrade3")) |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(adj.P.Val = NULL) |> 
  dplyr::mutate(adj.P.Val = p.adjust(P.Value, "fdr"))
  #tibble::column_to_rownames("protein_id")


sum((stats.grade.pat.corrected |> dplyr::filter(!is.na(P.Value)) |> dplyr::pull(P.Value) ) < 0.01)
sum((stats.grade.pat.corrected |> dplyr::filter(!is.na(adj.P.Val)) |> dplyr::pull(adj.P.Val) ) < 0.01)


rm(metadata, data)




saveRDS(stats.grade.pat.corrected, file="cache/analysis_differential_proteomics__GLASS-OD__stats.grade.pat-corrected.Rds")


rm(stats.grade.pat.corrected)




## prim-rec ----
### naive ----


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, patient_id, starts_with("proteomics_"), starts_with("resection_"), contains("_excluded")) |> 
  dplyr::mutate(condition = factor(ifelse(resection_number == 1,"primary","recurrence"), levels=c("primary","recurrence"))) |> 

  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 140) 
    return(.)
  })()




data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id)
data$nas <- rowSums(is.na(data))
data <- data |> 
  dplyr::filter(nas < 60) |> 
  dplyr::mutate(nas=NULL)

stopifnot(metadata$proteomics_id == colnames(data))




design <- model.matrix(~ factor(condition), data=metadata)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit, trend=T)
stats.time <- limma::topTable(fit,
                         n=nrow(data),
                         coef="factor(condition)recurrence",
                         sort.by = "none",
                         adjust.method="fdr") |> 
  as.data.frame() |> 
  tibble::rownames_to_column('protein_id')


sum(stats.time$adj.P.Val < 0.01)


rm(design, fit, data, metadata)


saveRDS(stats.time, file="cache/analysis_differential_proteomics__GLASS-OD__stats.prim-rec.naive.Rds")


rm(stats.time)



### pat corrected ----


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 

  dplyr::select(proteomics_id, patient_id, resection_number) |> 
  dplyr::mutate(primrec = factor(ifelse(resection_number == 1,"primary","recurrence"), levels=c("primary","recurrence"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 140) 
    return(.)
  })()



data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id) 



stats.time.pat.corrected <- data |> 
  tibble::rownames_to_column("protein_id") |> 
  dplyr::select("protein_id") |> 
  dplyr::mutate(tmp = pbapply::pbapply(data, 1, pp, metadata$patient_id, metadata$primrec, "conditionrecurrence")) |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(adj.P.Val = NULL) |> 
  dplyr::mutate(adj.P.Val = p.adjust(P.Value, "fdr"))
  #tibble::column_to_rownames("protein_id")


sum((stats.time.pat.corrected |> dplyr::filter(!is.na(P.Value)) |> dplyr::pull(P.Value) ) < 0.01)
sum((stats.time.pat.corrected |> dplyr::filter(!is.na(adj.P.Val)) |> dplyr::pull(adj.P.Val) ) < 0.01)



rm(metadata, data)



saveRDS(stats.time.pat.corrected, file="cache/analysis_differential_proteomics__GLASS-OD__stats.prim-rec.pat-corrected.Rds")


rm(stats.time.pat.corrected)




## CGC - regular ----


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, starts_with("proteomics_")) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 140) 
    return(.)
  })()

metadata <- metadata |> 
  dplyr::inner_join(
    glass_od.metadata.array_samples |> 
      filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
      dplyr::filter(resection_id %in% metadata$resection_id) |> 
      assertr::verify(!duplicated(resection_id)),
    by=c('resection_id'='resection_id')
  ) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 118) 
    return(.)
  })() |> 

  # dplyr::group_by(patient_id) |>
  # dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  # dplyr::ungroup() |> 
  # dplyr::mutate(patient = as.factor(ifelse(is.paired,paste0("p",patient_id),"a_remainder"))) |> 
   
  dplyr::mutate(CGC = scale(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit))




data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id)
data$nas <- rowSums(is.na(data))
data <- data |> 
  dplyr::filter(nas < 60) |> 
  dplyr::mutate(nas=NULL)

stopifnot(metadata$proteomics_id == colnames(data))





design <- model.matrix(~CGC, data=metadata)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit, trend=T)
stats.cgc <- limma::topTable(fit,
                              n=nrow(data),
                              coef="CGC",
                              sort.by = "none",
                              adjust.method="fdr") |> 
  tibble::rownames_to_column('protein_id')


rm(design, fit, data, metadata)


saveRDS(stats.cgc, file="cache/analysis_differential_proteomics__GLASS-OD__stats.cgc.Rds")

rm(stats.cgc)




# GLASS-NL ----
## CGC - regular ----


metadata <- glass_nl.metadata.array_samples |> 
  filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES) |> 
  dplyr::filter(!is.na(array_A_IDH_HG__A_IDH_LG_lr_v12.8 & !is.na(proteomics_sid))) |> 
  dplyr::select(Sample_Name, array_A_IDH_HG__A_IDH_LG_lr_v12.8, proteomics_sid) |> 
  dplyr::filter(!is.na(proteomics_sid)) |> 
  dplyr::rename(CGC = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
  dplyr::mutate(CGC = scale(CGC, center=T))



data <- data.proteomics.glass_nl |>
  dplyr::select(metadata$proteomics_sid)



design <- model.matrix(~CGC, data=metadata)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit, trend=T)
stats.cgc.glass_nl <- limma::topTable(fit,
                                      n=nrow(data),
                                      coef="CGC",
                                      sort.by = "none",
                                      adjust.method="fdr") |> 
  as.data.frame() |> 
  tibble::rownames_to_column('protein_id') 



saveRDS(stats.cgc.glass_nl, file="cache/analysis_differential_proteomics__GLASS-NL__stats.cgc.Rds")

rm(stats.cgc.glass_nl)





