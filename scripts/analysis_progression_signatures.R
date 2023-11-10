#!/usr/bin/env R

# load ----


source('scripts/load_mvalues_hq_samples.R')


# progression signature according to GLASS-NL paper ----

#' 1. load all mvalues
#' 2. load n=~1000 probes
#' 3. take per sampel median of ~1000 probes


stats <- read.csv("data/GLASS_NL/Metadata/(Epi)genetic_data/ProbeSelection_IvR_FDR 1e-9 Delta 1_06072022.csv") |> 
  dplyr::rename(probe_id = Probe_ID) |> 
  dplyr::mutate(deep_significant = T) |>
  dplyr::select(probe_id) |> 
  dplyr::left_join(
    data.mvalues.hq_samples |> tibble::rownames_to_column('probe_id'),
    by=c('probe_id'='probe_id'), suffix=c('','')
  ) |> 
  tibble::column_to_rownames('probe_id') |> 
  as.matrix() |> 
  matrixStats::colMedians() |> 
  as.data.frame() |> 
  dplyr::rename(array_GLASS_NL_g2_g3_sig = 1) |> 
  tibble::rownames_to_column('array_sentrix_id')



saveRDS(stats, file="cache/analysis_progression_signatures__GLASS-NL_prim_rec_signature.Rds")

rm(stats)


# glass-od g2/g3 ----
## median type 1: all ----


stats1 <- data.mvalues.probes |>
  dplyr::filter(DMP__g2_g3__pp_nc__adj.P.Val < 0.01 & abs(DMP__g2_g3__pp_nc__logFC) >= 0.5) |> 
  dplyr::select(probe_id) |> 
  dplyr::left_join(
    data.mvalues.hq_samples |> tibble::rownames_to_column('probe_id'),
    by=c('probe_id'='probe_id'), suffix=c('','')
  ) |> 
  tibble::column_to_rownames('probe_id') |> 
  as.matrix() |> 
  matrixStats::colMedians() |> 
  as.data.frame() |> 
  dplyr::rename(array_GLASS_OD_g2_g3_sig1 = 1) |> 
  tibble::rownames_to_column('array_sentrix_id')




## median type 2 ----


stats2 <- data.mvalues.probes |>
  dplyr::filter(DMP__g2_g3__pp_nc__adj.P.Val < 0.01 & abs(DMP__g2_g3__pp_nc__logFC) >= 0.5) |> 
  dplyr::arrange(desc(abs(DMP__g2_g3__pp_nc__t))) |> 
  head(n=10000) |> 
  dplyr::select(probe_id) |> 
  dplyr::left_join(
    data.mvalues.hq_samples |> tibble::rownames_to_column('probe_id'),
    by=c('probe_id'='probe_id'), suffix=c('','')
  ) |> 
  tibble::column_to_rownames('probe_id') |> 
  as.matrix() |> 
  matrixStats::colMedians() |> 
  as.data.frame() |> 
  dplyr::rename(array_GLASS_OD_g2_g3_sig2 = 1) |> 
  tibble::rownames_to_column('array_sentrix_id')



## median type 3 ----


stats3 <- data.mvalues.probes |>
  dplyr::filter(DMP__g2_g3__pp_nc__adj.P.Val < 0.01 & DMP__g2_g3__pp_nc__logFC <= -0.5) |> 
  dplyr::arrange(desc(abs(DMP__g2_g3__pp_nc__t))) |> 
  head(n=10000) |> 
  dplyr::select(probe_id) |> 
  dplyr::left_join(
    data.mvalues.hq_samples |> tibble::rownames_to_column('probe_id'),
    by=c('probe_id'='probe_id'), suffix=c('','')
  ) |> 
  tibble::column_to_rownames('probe_id') |> 
  as.matrix() |> 
  matrixStats::colMedians() |> 
  as.data.frame() |> 
  dplyr::rename(array_GLASS_OD_g2_g3_sig3 = 1) |> 
  tibble::rownames_to_column('array_sentrix_id')


## median type 4 ----


stats4 <- data.mvalues.probes |>
  dplyr::filter(DMP__g2_g3__pp_nc__adj.P.Val < 0.01 & DMP__g2_g3__pp_nc__logFC >= 0.5) |> 
  dplyr::arrange(desc(abs(DMP__g2_g3__pp_nc__t))) |> 
  head(n=10000) |> 
  dplyr::select(probe_id) |> 
  dplyr::left_join(
    data.mvalues.hq_samples |> tibble::rownames_to_column('probe_id'),
    by=c('probe_id'='probe_id'), suffix=c('','')
  ) |> 
  tibble::column_to_rownames('probe_id') |> 
  as.matrix() |> 
  matrixStats::colMedians() |> 
  as.data.frame() |> 
  dplyr::rename(array_GLASS_OD_g2_g3_sig4 = 1) |> 
  tibble::rownames_to_column('array_sentrix_id')



## supervised PCA ----


tmp <- data.mvalues.probes |>
  dplyr::filter(DMP__g2_g3__pp_nc__adj.P.Val < 0.01 & abs(DMP__g2_g3__pp_nc__logFC) >= 0.5) |> 
  dplyr::select(probe_id) |> 
  dplyr::left_join(
    data.mvalues.hq_samples |> tibble::rownames_to_column('probe_id'),
    by=c('probe_id'='probe_id'), suffix=c('','')
  ) |> 
  tibble::column_to_rownames('probe_id') 


pc <- tmp |> 
  t() |> 
  prcomp() |> 
  purrr::pluck('x') |> 
  as.data.frame() |> 
  dplyr::select(PC1, PC2, PC3, PC4) |>
  dplyr::rename(array_GLASS_OD_g2_g3_sig5 = PC1) |> 
  dplyr::rename(array_GLASS_OD_g2_g3_sig6 = PC2) |> 
  dplyr::rename(array_GLASS_OD_g2_g3_sig7 = PC3) |> 
  dplyr::rename(array_GLASS_OD_g2_g3_sig8 = PC4) |> 
  tibble::rownames_to_column('array_sentrix_id')



## exp ----


exp <- stats1 |> 
  dplyr::left_join(stats2, by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')) |> 
  dplyr::left_join(stats3, by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')) |> 
  dplyr::left_join(stats4, by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')) |> 
  dplyr::left_join(pc, by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('',''))


saveRDS(exp, file="cache/analysis_progression_signatures__GLASS-OD__g2_g3.Rds")


rm(stats1, stats2, stats3, stats4, tmp, pc)


# glass-od primary/recurrence ----
## median type 1: all ----


stats1 <- data.mvalues.probes |>
  dplyr::filter(DMP__primary_recurrence__pp_nc__adj.P.Val < 0.01 & abs(DMP__primary_recurrence__pp_nc__logFC) >= 0.5) |> 
  dplyr::select(probe_id) |> 
  dplyr::left_join(
    data.mvalues.hq_samples |> tibble::rownames_to_column('probe_id'),
    by=c('probe_id'='probe_id'), suffix=c('','')
  ) |> 
  tibble::column_to_rownames('probe_id') |> 
  as.matrix() |> 
  matrixStats::colMedians() |> 
  as.data.frame() |> 
  dplyr::rename(array_GLASS_OD_prim_rec_sig1 = 1) |> 
  tibble::rownames_to_column('array_sentrix_id')




## median type 2 ----


stats2 <- data.mvalues.probes |>
  dplyr::filter(DMP__primary_recurrence__pp_nc__adj.P.Val < 0.01 & abs(DMP__primary_recurrence__pp_nc__logFC) >= 0.5) |> 
  dplyr::arrange(desc(abs(DMP__primary_recurrence__pp_nc__t))) |> 
  head(n=10000) |> 
  dplyr::select(probe_id) |> 
  dplyr::left_join(
    data.mvalues.hq_samples |> tibble::rownames_to_column('probe_id'),
    by=c('probe_id'='probe_id'), suffix=c('','')
  ) |> 
  tibble::column_to_rownames('probe_id') |> 
  as.matrix() |> 
  matrixStats::colMedians() |> 
  as.data.frame() |> 
  dplyr::rename(array_GLASS_OD_prim_rec_sig2 = 1) |> 
  tibble::rownames_to_column('array_sentrix_id')



## median type 3 ----


stats3 <- data.mvalues.probes |>
  dplyr::filter(DMP__primary_recurrence__pp_nc__adj.P.Val < 0.01 & DMP__primary_recurrence__pp_nc__logFC <= -0.5) |> 
  dplyr::arrange(desc(abs(DMP__primary_recurrence__pp_nc__t))) |> 
  head(n=10000) |> 
  dplyr::select(probe_id) |> 
  dplyr::left_join(
    data.mvalues.hq_samples |> tibble::rownames_to_column('probe_id'),
    by=c('probe_id'='probe_id'), suffix=c('','')
  ) |> 
  tibble::column_to_rownames('probe_id') |> 
  as.matrix() |> 
  matrixStats::colMedians() |> 
  as.data.frame() |> 
  dplyr::rename(array_GLASS_OD_prim_rec_sig3 = 1) |> 
  tibble::rownames_to_column('array_sentrix_id')


## median type 4 ----


stats4 <- data.mvalues.probes |>
  dplyr::filter(DMP__primary_recurrence__pp_nc__adj.P.Val < 0.01 & DMP__primary_recurrence__pp_nc__logFC >= 0.5) |> 
  dplyr::arrange(desc(abs(DMP__primary_recurrence__pp_nc__t))) |> 
  head(n=10000) |> 
  dplyr::select(probe_id) |> 
  dplyr::left_join(
    data.mvalues.hq_samples |> tibble::rownames_to_column('probe_id'),
    by=c('probe_id'='probe_id'), suffix=c('','')
  ) |> 
  tibble::column_to_rownames('probe_id') |> 
  as.matrix() |> 
  matrixStats::colMedians() |> 
  as.data.frame() |> 
  dplyr::rename(array_GLASS_OD_prim_rec_sig4 = 1) |> 
  tibble::rownames_to_column('array_sentrix_id')



## supervised PCA ----


tmp <- data.mvalues.probes |>
  dplyr::filter(DMP__primary_recurrence__pp_nc__adj.P.Val < 0.01 & abs(DMP__primary_recurrence__pp_nc__logFC) >= 0.5) |> 
  dplyr::select(probe_id) |> 
  dplyr::left_join(
    data.mvalues.hq_samples |> tibble::rownames_to_column('probe_id'),
    by=c('probe_id'='probe_id'), suffix=c('','')
  ) |> 
  tibble::column_to_rownames('probe_id') 


pc <- tmp |> 
  t() |> 
  prcomp() |> 
  purrr::pluck('x') |> 
  as.data.frame() |> 
  dplyr::select(PC1, PC2, PC3, PC4) |>
  dplyr::rename(array_GLASS_OD_prim_rec_sig5 = PC1) |> 
  dplyr::rename(array_GLASS_OD_prim_rec_sig6 = PC2) |> 
  dplyr::rename(array_GLASS_OD_prim_rec_sig7 = PC3) |> 
  dplyr::rename(array_GLASS_OD_prim_rec_sig8 = PC4) |> 
  tibble::rownames_to_column('array_sentrix_id')



## exp ----


exp <- stats1 |> 
  dplyr::left_join(stats2, by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')) |> 
  dplyr::left_join(stats3, by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')) |> 
  dplyr::left_join(stats4, by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')) |> 
  dplyr::left_join(pc, by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('',''))


saveRDS(exp, file="cache/analysis_progression_signatures__GLASS-OD__primary_recurrence.Rds")



