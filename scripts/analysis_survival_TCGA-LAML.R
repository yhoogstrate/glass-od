#!/usr/bin/env R
#' combined reading and normalization of GLASS-NL and GLASS-OD data

# load data ----


source('scripts/load_functions.R')


if(!exists('tcga_laml.metadata.array_samples')) {
  source('scripts/load_TCGA-LAML_metadata.R')
}



# svvl analysis ----


library(survival, survminer)


## all samples ----


stats.meta.all <- tcga_laml.metadata.array_samples |> 
  dplyr::filter(array_type == "450k") |>
  dplyr::filter(!is.na(last_event_status)) |>  # must have survival info
  dplyr::select(day_to_last_event, last_event_status, array_sentrix_id)


stats.mat.all <- readRDS(file="cache/mvalues.TCGA-LAML.Rds") |> 
  dplyr::select(all_of(stats.meta.all |> dplyr::pull(array_sentrix_id))) |> 
  #head(n=5) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::left_join(stats.meta.all, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  tibble::column_to_rownames('array_sentrix_id')



stats.cox.all <- RegParallel::RegParallel(
  data = stats.mat.all,
  formula = 'Surv(day_to_last_event, last_event_status) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(stats.mat.all) |>  purrr::keep(function(x) grepl("^cg",x)),
  blocksize = 5,
  cores = 16,
  nestedParallel = FALSE,
  conflevel = 95)


saveRDS(stats.cox.all, file="cache/analysis_suvival_TCGA-LAML.all.R")




## IDH-pos samples ----

stats.meta.IDHpositive <- tcga_laml.metadata.array_samples |> 
  dplyr::filter(array_type == "450k" & IDH) |>
  dplyr::filter(!is.na(last_event_status)) |>  # must have survival info
  dplyr::select(day_to_last_event, last_event_status, array_sentrix_id)


stats.mat.IDHpositive <- readRDS(file="cache/mvalues.TCGA-LAML.Rds") |> 
  dplyr::select(all_of(stats.meta.IDHpositive |> dplyr::pull(array_sentrix_id))) |> 
  #head(n=5) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::left_join(stats.meta.IDHpositive, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  tibble::column_to_rownames('array_sentrix_id')


stats.cox.IDHpositive <- RegParallel::RegParallel(
  data = stats.mat.IDHpositive,
  formula = 'Surv(day_to_last_event, last_event_status) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(stats.mat.IDHpositive) |>  purrr::keep(function(x) grepl("^cg",x)),
  blocksize = 5,
  cores = 16,
  nestedParallel = FALSE,
  conflevel = 95)


saveRDS(stats.cox.all, file="cache/analysis_suvival_TCGA-LAML.IDHpositive.R")





## IDH-neg samples ----


stats.meta.IDHnegative <- tcga_laml.metadata.array_samples |> 
  dplyr::filter(array_type == "450k" & !IDH) |>
  dplyr::filter(!is.na(last_event_status)) |>  # must have survival info
  dplyr::select(day_to_last_event, last_event_status, array_sentrix_id)


stats.mat.IDHnegative <- readRDS(file="cache/mvalues.TCGA-LAML.Rds") |> 
  dplyr::select(all_of(stats.meta.IDHnegative |> dplyr::pull(array_sentrix_id))) |> 
  #head(n=5) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::left_join(stats.meta.IDHnegative, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  tibble::column_to_rownames('array_sentrix_id')


stats.cox.IDHnegative <- RegParallel::RegParallel(
  data = stats.mat.IDHnegative,
  formula = 'Surv(day_to_last_event, last_event_status) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(stats.mat.IDHnegative) |>  purrr::keep(function(x) grepl("^cg",x)),
  blocksize = 5,
  cores = 16,
  nestedParallel = FALSE,
  conflevel = 95)


saveRDS(stats.cox.all, file="cache/analysis_suvival_TCGA-LAML.IDHnegative.R")






