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
source('scripts/load_themes.R')
source('scripts/load_palette.R')



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


pp2 <- function(measurements, conditions, coeff) {
  stopifnot(length(conditions) == length(measurements))
  
  metadata <- data.frame(
#    patient_id = patient_ids[!is.na(measurements)],
    condition  = conditions[  !is.na(measurements)]
  )
  
  measurements <- tibble::tibble(t(tibble::tibble(measurements[!is.na(measurements)])))
  
  #print(dim(metadata))
  #print(dim(measurements))
  
  tryCatch({
    design <- model.matrix(~condition, data=metadata)
    fit <- limma::lmFit(measurements, design)
    fit <- limma::eBayes(fit, trend=T)
    
    #print(fit)
    
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
  dplyr::filter(proteomics_id != "GLODprot_11_1524") |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F)  # no pass qc
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
  theme_cellpress



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



dim(ki67.joined |> dplyr::filter(!is.na(MKI67)))
cor(exp(ki67.joined$MKI67),
    ki67.joined$staining_KI67_pos_per_area_um2,
    use="complete.obs", method="spearman")



# 0.4589544



## WHO Grade ----
### naive ----


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc

  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::filter(!is.na(resection_tumor_grade)) |>  # 0053-R2 pathology was not decisve
  dplyr::select(resection_id, patient_id, starts_with("proteomics_"), starts_with("resection_"), contains("_excluded")) |> 
  dplyr::mutate(condition = factor(paste0("g", resection_tumor_grade), levels=c("g2", "g3"))) |> 
  assertr::verify(!duplicated(resection_id)) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 136) 
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
# 28 zonder slechte samples erin opgenomen
# 43 met slechte erin opgenomen


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::filter(!is.na(resection_tumor_grade)) |>  # 0053-R2 pathology was not decisve
  
  assertr::verify(!duplicated(resection_id)) |> 
  
  
  dplyr::select(proteomics_id, patient_id, resection_tumor_grade) |> 
  dplyr::filter(!is.na(resection_tumor_grade)) |> 
  dplyr::mutate(resection_tumor_grade = factor(ifelse(resection_tumor_grade == 2, "Grade2", "Grade3"), levels=c("Grade2", "Grade3"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 136) 
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


readRDS(file="cache/analysis_differential_proteomics__GLASS-OD__stats.grade.pat-corrected.Rds") |> 
  dplyr::filter(!is.na(adj.P.Val) & adj.P.Val < 0.01) 




## prim-rec ----
### naive ----


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, patient_id, starts_with("proteomics_"), starts_with("resection_"), contains("_excluded")) |> 
  dplyr::mutate(condition = factor(ifelse(resection_number == 1,"primary","recurrence"), levels=c("primary","recurrence"))) |> 
  assertr::verify(!duplicated(resection_id)) |> 

  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 137) 
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
stats.time |> dplyr::filter(adj.P.Val < 0.01)

rm(design, fit, data, metadata)


saveRDS(stats.time, file="cache/analysis_differential_proteomics__GLASS-OD__stats.prim-rec.naive.Rds")


rm(stats.time)



### pat corrected ----


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  assertr::verify(!duplicated(resection_id)) |> 
  
  dplyr::select(proteomics_id, patient_id, resection_number) |> 
  dplyr::mutate(primrec = factor(ifelse(resection_number == 1,"primary","recurrence"), levels=c("primary","recurrence"))) |> 
  

  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 137) 
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


stats.time.pat.corrected |> 
  dplyr::filter(!is.na(adj.P.Val) & adj.P.Val < 0.01) 


rm(metadata, data)



saveRDS(stats.time.pat.corrected, file="cache/analysis_differential_proteomics__GLASS-OD__stats.prim-rec.pat-corrected.Rds")


rm(stats.time.pat.corrected)


# 1. FTL
readRDS(file="cache/analysis_differential_proteomics__GLASS-OD__stats.prim-rec.pat-corrected.Rds") |> 
  dplyr::filter(!is.na(adj.P.Val) & adj.P.Val < 0.01) 



## CGC - regular ----


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, starts_with("proteomics_")) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 137) 
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


stats.cgc |> 
  dplyr::filter(abs(logFC)>0.5 & stats.cgc$adj.P.Val < 0.01) |> 
  dim()


rm(design, fit, data, metadata)


saveRDS(stats.cgc, file="cache/analysis_differential_proteomics__GLASS-OD__stats.Rds")

rm(stats.cgc)



## CGC - pp style ----


metadata <- glass_od.metadata.proteomics |> 
  #dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
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
  dplyr::mutate(CGC = scale(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit))




data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id)
data$nas <- rowSums(is.na(data))
data <- data |> 
  dplyr::filter(nas < 60) |> 
  dplyr::mutate(nas=NULL)

stopifnot(metadata$proteomics_id == colnames(data))



stats.cgc.pp2 <- data |> 
  tibble::rownames_to_column("protein_id") |> 
  dplyr::select("protein_id") |> 
  dplyr::mutate(tmp = pbapply::pbapply(data, 1, pp2, metadata$CGC, "condition")) |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(adj.P.Val = NULL) |> 
  dplyr::mutate(adj.P.Val = p.adjust(P.Value, "fdr"))
#tibble::column_to_rownames("protein_id")



rm(design, fit, data, metadata)


saveRDS(stats.cgc.pp2, file="cache/analysis_differential_proteomics__GLASS-OD__stats.cgc-pp2-style.Rds")

rm(stats.cgc.pp2)

# using pp2 differs,but slighly, cor ~0.98, an virtually no corr between n samples and t stat

plt <- dplyr::left_join(
  stats.cgc,
  stats.cgc.pp2,
  by=c('protein_id'='protein_id'),
  suffix = c('.naive','.pp2')
)

ggplot(plt, aes(x=t.pp2, y=t.naive)) +
  geom_point(pch=19, cex=0.5)


ggplot(plt, aes(x=t.pp2, y=n)) +
  geom_point(pch=19, cex=0.5)



## M - F ----



metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::filter(!is.na(resection_tumor_grade)) |>  # 0053-R2 pathology was not decisve
  dplyr::select(resection_id, patient_sex, patient_id, starts_with("proteomics_"), starts_with("resection_"), contains("_excluded")) |> 
  dplyr::mutate(condition = factor(patient_sex, levels=c("female", "male"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 139) 
    return(.)
  })()



plt <- data.proteomics.glass_od |>
  tibble::rownames_to_column('protein') |>
  dplyr::filter(protein %in% c("XIST", "FOXL2", "RSPO1", "CYP19A1", "ZP3", "FIGLA", "SRY", "SOX9", "DAZ", "PRM1", "PRM2", "TSPY", "AMH")) |> 
  tibble::column_to_rownames('protein') |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('proteomics_id') |> 
  dplyr::left_join(
    metadata, by=c('proteomics_id'='proteomics_id')
  )

ggplot(plt, aes(x=condition,y=SOX9)) +
  geom_point()





data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id)
data$nas <- rowSums(is.na(data))
data <- data |> 
  dplyr::filter(nas < 60) |> 
  dplyr::mutate(nas=NULL)




data.proteomics.glass_od |>
  tibble::rownames_to_column('protein') |> 
  dplyr::filter(grepl("^X", protein)) |>
  dplyr::select(protein)



stopifnot(metadata$proteomics_id == colnames(data))



design <- model.matrix(~factor(condition), data=metadata)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit, trend=T)
stats.grade <- limma::topTable(fit,
                               n=nrow(data),
                               coef="factor(condition)male",
                               sort.by = "none",
                               adjust.method="fdr") |> 
  as.data.frame() |> 
  tibble::rownames_to_column('protein_id')

sum(stats.grade$adj.P.Val < 0.01)


rm(design, fit, data, metadata)







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




# per gene intersection ----

## TMPO ----

metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, starts_with("proteomics_")) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 137) 
    return(.)
  })() |> 
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
  })()



data.p <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id) |> 
  tibble::rownames_to_column('protein_id') |> 
  dplyr::filter(grepl("^TMPO$", protein_id)) |> 
  tibble::column_to_rownames('protein_id') |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('proteomics_id')


data.m.features <- data.mvalues.probes |> dplyr::filter(grepl("TMPO", UCSC_RefGene_Name) |
                                                          grepl("TMPO", GencodeCompV12_NAME)) |> 
  dplyr::filter(grepl("^cg", probe_id)) |> 
  dplyr::select(probe_id, UCSC_RefGene_Name, GencodeCompV12_NAME, CpG_chrm,   CpG_beg) |> 
  dplyr::arrange( CpG_chrm, CpG_beg) |> 
  dplyr::pull(probe_id)


data.m <- data.mvalues.hq_samples |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(all_of(data.m.features)) |> 
  tibble::rownames_to_column('array_sentrix_id')


data <- metadata |> 
  dplyr::select(resection_id, proteomics_id, array_sentrix_id) |> 
  dplyr::left_join(data.p, by=c('proteomics_id'='proteomics_id')) |> 
  dplyr::left_join(data.m, by=c('array_sentrix_id'='array_sentrix_id')) |> 
  dplyr::mutate(proteomics_id =  NULL, array_sentrix_id = NULL) |> 
  tibble::column_to_rownames('resection_id')


ggcorrplot(data, reorder=F)


## COL1A1 ----


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, starts_with("proteomics_")) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 137) 
    return(.)
  })() |> 
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
  })()



data.p <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id) |> 
  tibble::rownames_to_column('protein_id') |> 
  dplyr::filter(grepl("^COL1A1$", protein_id)) |> 
  tibble::column_to_rownames('protein_id') |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('proteomics_id')


data.m.features <- data.mvalues.probes |> dplyr::filter(grepl("COL1A1", UCSC_RefGene_Name) |
                                                          grepl("COL1A1", GencodeCompV12_NAME)) |> 
  dplyr::filter(grepl("^cg", probe_id)) |> 
  dplyr::select(probe_id, UCSC_RefGene_Name, GencodeCompV12_NAME, CpG_chrm,   CpG_beg) |> 
  dplyr::arrange( CpG_chrm,   CpG_beg) |> 
  dplyr::pull(probe_id)


data.m <- data.mvalues.hq_samples |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(all_of(data.m.features)) |> 
  tibble::rownames_to_column('array_sentrix_id')


data <- metadata |> 
  dplyr::select(resection_id, proteomics_id, array_sentrix_id) |> 
  dplyr::left_join(data.p, by=c('proteomics_id'='proteomics_id')) |> 
  dplyr::left_join(data.m, by=c('array_sentrix_id'='array_sentrix_id')) |> 
  dplyr::mutate(proteomics_id =  NULL, array_sentrix_id = NULL) |> 
  tibble::column_to_rownames('resection_id')


ggcorrplot(data, reorder=F)




## COL1A2 ----

metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, starts_with("proteomics_")) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 137) 
    return(.)
  })() |> 
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
  })()



data.p <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id) |> 
  tibble::rownames_to_column('protein_id') |> 
  dplyr::filter(grepl("^COL1A2$", protein_id)) |> 
  tibble::column_to_rownames('protein_id') |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('proteomics_id')


data.m.features <- data.mvalues.probes |> dplyr::filter(grepl("COL1A2", UCSC_RefGene_Name) |
                                                          grepl("COL1A2", GencodeCompV12_NAME)) |> 
  dplyr::select(probe_id, UCSC_RefGene_Name, GencodeCompV12_NAME, CpG_chrm,   CpG_beg) |> 
  dplyr::filter(grepl("^cg", probe_id)) |> 
  dplyr::arrange( CpG_chrm,   CpG_beg) |> 
  dplyr::pull(probe_id)


data.m <- data.mvalues.hq_samples |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(all_of(data.m.features)) |> 
  tibble::rownames_to_column('array_sentrix_id')


data <- metadata |> 
  dplyr::select(resection_id, proteomics_id, array_sentrix_id, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::rename(CGC = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::left_join(data.p, by=c('proteomics_id'='proteomics_id')) |> 
  dplyr::left_join(data.m, by=c('array_sentrix_id'='array_sentrix_id')) |> 
  dplyr::mutate(proteomics_id =  NULL, array_sentrix_id = NULL) |> 
  tibble::column_to_rownames('resection_id')


ggcorrplot(data, reorder=F)




## PLP1 ----

metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, starts_with("proteomics_")) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 137) 
    return(.)
  })() |> 
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
  })()



data.p <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id) |> 
  tibble::rownames_to_column('protein_id') |> 
  dplyr::filter(grepl("^PLP1$", protein_id)) |> 
  tibble::column_to_rownames('protein_id') |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('proteomics_id')


data.m.features <- data.mvalues.probes |> 
  dplyr::filter(grepl("PLP1", UCSC_RefGene_Name) | grepl("PLP1", GencodeCompV12_NAME)) |> 
  dplyr::filter(!grepl("RPLP1", UCSC_RefGene_Name)) |> 
  dplyr::filter(!grepl("RPLP1", GencodeCompV12_NAME)) |> 
  dplyr::filter(!grepl("APLP1", UCSC_RefGene_Name)) |> 
  dplyr::filter(!grepl("APLP1", GencodeCompV12_NAME)) |> 
  dplyr::select(probe_id, UCSC_RefGene_Name, GencodeCompV12_NAME, CpG_chrm,   CpG_beg) |> 
  dplyr::filter(grepl("^cg", probe_id)) |> 
  dplyr::arrange( CpG_chrm,   CpG_beg) |> 
  dplyr::pull(probe_id)


data.m <- data.mvalues.hq_samples |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(all_of(data.m.features)) |> 
  tibble::rownames_to_column('array_sentrix_id')


data <- metadata |> 
  dplyr::select(resection_id, proteomics_id, array_sentrix_id, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::rename(CGC = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(CGC = NULL) |> 
  dplyr::left_join(data.p, by=c('proteomics_id'='proteomics_id')) |> 
  dplyr::left_join(data.m, by=c('array_sentrix_id'='array_sentrix_id')) |> 
  dplyr::mutate(proteomics_id =  NULL, array_sentrix_id = NULL) |> 
  tibble::column_to_rownames('resection_id')


ggcorrplot(data, reorder=F)



## FN1 ----

metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, starts_with("proteomics_")) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 137) 
    return(.)
  })() |> 
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
  })()



data.p <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id) |> 
  tibble::rownames_to_column('protein_id') |> 
  dplyr::filter(grepl("^FN1|FINC$", protein_id)) |> 
  tibble::column_to_rownames('protein_id') |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('proteomics_id')


data.m.features <- data.mvalues.probes |> 
  dplyr::filter(grepl("FN1|FINC", UCSC_RefGene_Name) | grepl("FN1|FINC", GencodeCompV12_NAME)) |> 
  dplyr::filter(!grepl("[A-Z]FN1", UCSC_RefGene_Name)) |> 
  dplyr::filter(!grepl("[A-Z]FN1", GencodeCompV12_NAME)) |> 
  dplyr::select(probe_id, UCSC_RefGene_Name, GencodeCompV12_NAME, CpG_chrm,   CpG_beg) |> 
  dplyr::filter(grepl("^cg", probe_id)) |> 
  dplyr::arrange( CpG_chrm,   CpG_beg) |> 
  dplyr::pull(probe_id)


data.m <- data.mvalues.hq_samples |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(all_of(data.m.features)) |> 
  tibble::rownames_to_column('array_sentrix_id')


data <- metadata |> 
  dplyr::select(resection_id, proteomics_id, array_sentrix_id, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::rename(CGC = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(CGC = NULL) |> 
  dplyr::left_join(data.p, by=c('proteomics_id'='proteomics_id')) |> 
  dplyr::left_join(data.m, by=c('array_sentrix_id'='array_sentrix_id')) |> 
  dplyr::mutate(proteomics_id =  NULL, array_sentrix_id = NULL) |> 
  tibble::column_to_rownames('resection_id')


ggcorrplot(data, reorder=F)


## ANXA1 ----


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, starts_with("proteomics_")) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 137) 
    return(.)
  })() |> 
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
  })()



data.p <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id) |> 
  tibble::rownames_to_column('protein_id') |> 
  dplyr::filter(grepl("^ANXA1$", protein_id)) |> 
  tibble::column_to_rownames('protein_id') |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('proteomics_id')


data.m.features <- data.mvalues.probes |> 
  dplyr::filter(grepl("ANXA1", UCSC_RefGene_Name) | grepl("ANXA1", GencodeCompV12_NAME)) |> 
  dplyr::filter(!grepl("[A-Z]ANXA1", UCSC_RefGene_Name)) |> 
  dplyr::filter(!grepl("[A-Z]ANXA1", GencodeCompV12_NAME)) |> 
  dplyr::select(probe_id, UCSC_RefGene_Name, GencodeCompV12_NAME, CpG_chrm,   CpG_beg) |> 
  dplyr::filter(grepl("^cg", probe_id)) |> 
  dplyr::arrange( CpG_chrm,   CpG_beg) |> 
  dplyr::pull(probe_id)


data.m <- data.mvalues.hq_samples |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(all_of(data.m.features)) |> 
  tibble::rownames_to_column('array_sentrix_id')


data <- metadata |> 
  dplyr::select(resection_id, proteomics_id, array_sentrix_id, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::rename(CGC = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(CGC = NULL) |> 
  dplyr::left_join(data.p, by=c('proteomics_id'='proteomics_id')) |> 
  dplyr::left_join(data.m, by=c('array_sentrix_id'='array_sentrix_id')) |> 
  dplyr::mutate(proteomics_id =  NULL, array_sentrix_id = NULL) |> 
  tibble::column_to_rownames('resection_id')


ggcorrplot(data, reorder=F)


## PCNA ----


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, starts_with("proteomics_")) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 137) 
    return(.)
  })() |> 
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
  })()



data.p <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id) |> 
  tibble::rownames_to_column('protein_id') |> 
  dplyr::filter(grepl("^PCNA$", protein_id)) |> 
  tibble::column_to_rownames('protein_id') |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('proteomics_id')


data.m.features <- data.mvalues.probes |> 
  dplyr::filter(grepl("PCNA", UCSC_RefGene_Name) | grepl("PCNA", GencodeCompV12_NAME)) |> 
  dplyr::filter(!grepl("[A-Z]PCNA", UCSC_RefGene_Name)) |> 
  dplyr::filter(!grepl("[A-Z]PCNA", GencodeCompV12_NAME)) |> 
  dplyr::select(probe_id, UCSC_RefGene_Name, GencodeCompV12_NAME, CpG_chrm, CpG_beg) |> 
  dplyr::filter(grepl("^cg", probe_id)) |> 
  dplyr::arrange( CpG_chrm,   CpG_beg) |> 
  dplyr::pull(probe_id)


data.m <- data.mvalues.hq_samples |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(all_of(data.m.features)) |> 
  tibble::rownames_to_column('array_sentrix_id')


data <- metadata |> 
  dplyr::select(resection_id, proteomics_id, array_sentrix_id, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::rename(CGC = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(CGC = NULL) |> 
  dplyr::left_join(data.p, by=c('proteomics_id'='proteomics_id')) |> 
  dplyr::left_join(data.m, by=c('array_sentrix_id'='array_sentrix_id')) |> 
  dplyr::mutate(proteomics_id =  NULL, array_sentrix_id = NULL) |> 
  tibble::column_to_rownames('resection_id') |> 
  (\(x) dplyr::filter(x, complete.cases(x)))()


ggcorrplot(data, reorder=F)




## SLC1A3 ----

metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, starts_with("proteomics_")) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 137) 
    return(.)
  })() |> 
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
  })()



data.p <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id) |> 
  tibble::rownames_to_column('protein_id') |> 
  dplyr::filter(grepl("^SLC1A3$", protein_id)) |> 
  tibble::column_to_rownames('protein_id') |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('proteomics_id')


data.m.features <- data.mvalues.probes |> 
  dplyr::filter(grepl("SLC1A3", UCSC_RefGene_Name) | grepl("SLC1A3", GencodeCompV12_NAME)) |> 
  dplyr::filter(!grepl("[A-Z]SLC1A3", UCSC_RefGene_Name)) |> 
  dplyr::filter(!grepl("[A-Z]SLC1A3", GencodeCompV12_NAME)) |> 
  dplyr::select(probe_id, UCSC_RefGene_Name, GencodeCompV12_NAME, CpG_chrm, CpG_beg) |> 
  dplyr::filter(grepl("^cg", probe_id)) |> 
  dplyr::arrange( CpG_chrm,   CpG_beg) |> 
  dplyr::pull(probe_id)


data.m <- data.mvalues.hq_samples |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(all_of(data.m.features)) |> 
  tibble::rownames_to_column('array_sentrix_id')


data <- metadata |> 
  dplyr::select(resection_id, proteomics_id, array_sentrix_id, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::rename(CGC = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(CGC = NULL) |> 
  dplyr::left_join(data.p, by=c('proteomics_id'='proteomics_id')) |> 
  dplyr::left_join(data.m, by=c('array_sentrix_id'='array_sentrix_id')) |> 
  dplyr::mutate(proteomics_id =  NULL, array_sentrix_id = NULL) |> 
  tibble::column_to_rownames('resection_id') |> 
  (\(x) dplyr::filter(x, complete.cases(x)))()


ggcorrplot(data, reorder=F)






# checks ----




data.proteomics.glass_od |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column("proteomics_id") |> 
  dplyr::left_join(glass_od.metadata.proteomics |> dplyr::select(proteomics_id, resection_id), by=c('proteomics_id'='proteomics_id')) |> 
  dplyr::filter(resection_id %in% c("0056-R1", "0056-R3", "0042-R1", "0012-R1") )  |> 
  dplyr::pull(proteomics_id)



expr = data.proteomics.glass_od |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column("proteomics_id") |> 
  dplyr::left_join(glass_od.metadata.proteomics |> dplyr::select(proteomics_id, resection_id), by=c('proteomics_id'='proteomics_id')) |> 
  dplyr::group_by(resection_id) |> 
  dplyr::mutate(resection_id = if(dplyr::n() > 1) paste0(resection_id, "-rep", dplyr::row_number()) else resection_id) |> 
  dplyr::ungroup() |> 
  #dplyr::mutate(proteomics_id = NULL) |> 
  dplyr::filter(resection_id %in% c("0056-R1", "0056-R3", "0042-R1", "0012-R1") == F) |> # odd outliers
  tibble::column_to_rownames("resection_id") |> 
  t() |> 
  as.matrix()



ggcorrplot(na.omit(expr), abs=T)



library(matrixStats)
sample_totals <- colSums(2^expr, na.rm=TRUE)   # or colSums(expr, na.rm=TRUE) if already linear
missing_by_sample <- colMeans(is.na(expr))
missing_by_protein <- rowMeans(is.na(expr))
# replicate correlation (pearson)
cor_mat <- cor(expr, use="pairwise.complete.obs", method="pearson")
median_rep_corr <- median(cor_mat[lower.tri(cor_mat)])


corrplot::corrplot(cor_mat,  order="hclust", tl.cex=0.35)



pca <- prcomp(t(na.omit(expr)), scale.=F)  # or use imputed data
plt <- pca$x |> 
  as.data.frame() |> 
  tibble::rownames_to_column('resection_id') |> 
  dplyr::mutate(rep = grepl("rep",resection_id)) |> 
  dplyr::mutate(resection_id = gsub("-rep.$","", resection_id)) |> 
  dplyr::left_join(glass_od.metadata.resections, by=c('resection_id'='resection_id')) |> 
  dplyr::left_join(glass_od.metadata.array_samples |>
                     dplyr::filter(!is.na(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |>
                     dplyr::select(resection_id, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit,
                                   array_PC1, array_PC2, array_PC3), by=c('resection_id'='resection_id')) |> 
  dplyr::rename(CGC = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)



corrplot::corrplot(
  plt |> 
    dplyr::filter(!is.na(CGC)) |> 
    dplyr::select(CGC,paste0("PC",1:6), paste0("array_PC",2:3)) |> 
    cor(method="spearman")
)


ggplot(plt, aes(x=PC1, y=PC4, label=resection_id, col=as.factor(resection_tumor_grade))) +
  geom_point() +
  #ggrepel::geom_text_repel(size = 3.6) +
  theme_cellpress



ggplot(plt, aes(x=PC1, y=PC2, label=resection_id, col=resection_tumor_grade)) +
  geom_point() +
  #ggrepel::geom_text_repel(size = 3.6) +
  theme_cellpress


ggplot(plt, aes(x=PC1, y=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, label=resection_id)) +
  geom_point() +
  #ggrepel::geom_text_repel(size = 3.6) +
  theme_cellpress

ggplot(plt, aes(x=PC10, y=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, label=resection_id)) +
  geom_point() +
  #ggrepel::geom_text_repel(size = 3.6) +
  theme_cellpress





ggplot(plt, aes(x=PC2, y=PC5, label=resection_id, col=resection_tumor_grade)) +
  geom_point() +
  #ggrepel::geom_text_repel(size = 3.6) +
  theme_cellpress

ggplot(plt, aes(x=PC1, y=PC4, label=resection_id, col=resection_tumor_grade)) +
  geom_point() +
  #ggrepel::geom_text_repel(size = 3.6) +
  theme_cellpress

ggplot(plt, aes(x=CGC, y=PC5, label=resection_id, col=resection_tumor_grade)) +
  geom_point() +
  #ggrepel::geom_text_repel(size = 3.6) +
  theme_cellpress

ggplot(plt, aes(x=PC1, y=PC6, label=resection_id, col=resection_tumor_grade)) +
  geom_point() +
  #ggrepel::geom_text_repel(size = 3.6) +
  theme_cellpress





# 1. Get loadings (variable contributions)
loadings <- pca$rotation  # matrix: variables (proteins) × PCs

# 2. Extract absolute contributions to PC1
pc1_load <- abs(loadings[, 1])

# 3. Sort and take top N (e.g. top 20)
top_pc1 <- sort(pc1_load, decreasing = TRUE)[1:40]

# 4. View or get names
top_pc1_genes <- names(top_pc1)
top_pc1_genes



