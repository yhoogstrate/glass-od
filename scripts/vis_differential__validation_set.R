#!/usr/bin/env R



# match both back to grade in main dataset ----


readRDS("cache/analysis_differential__g2_g3__partial_paired_nc__validationset__stats.Rds") |>
  dplyr::mutate(signi = adj.P.Val < 0.01) |>
  dplyr::pull(signi) |>
  table()


readRDS("cache/analysis_differential__p_r__PC1__partial_paired_nc__validationset__stats.Rds") |>
  dplyr::mutate(signi = adj.P.Val < 0.01) |>
  dplyr::pull(signi) |>
  table()


readRDS("cache/analysis_differential__g2_g3__partial_paired_nc__validationset__stats.Rds") |>
  dplyr::mutate(signi = adj.P.Val < 0.01) |>
  dplyr::pull(signi) |>
  table()


readRDS("cache/analysis_differential__p_r__PC1__partial_paired_nc__validationset__stats.Rds") |>
  dplyr::mutate(signi = adj.P.Val < 0.01) |>
  dplyr::pull(signi) |>
  table()





plt <- data.mvalues.probes |> 
  dplyr::select(probe_id,
                DMP__primary_recurrence__pp_nc__t,
                DMP__primary_recurrence__pp_nc_PC1__t,
                DMP__g2_g3__pp_nc__t,
                DMP__g2_g3__pp_nc_PC1__t) |> 
  dplyr::left_join(
    readRDS("cache/analysis_differential__g2_g3__partial_paired_nc__validationset__stats.Rds") |> 
      dplyr::select(probe_id, t) |> 
      dplyr::rename( t_validation_g2_g3 = t),
    by=c('probe_id'='probe_id')
  )|> 
  dplyr::left_join(
    readRDS("cache/analysis_differential__p_r__partial_paired_nc__validationset__stats.Rds") |> 
      dplyr::select(probe_id, t) |> 
      dplyr::rename(t_validation_p_r = t),
    by=c('probe_id'='probe_id')
  ) |> 
  dplyr::left_join(
    readRDS("cache/analysis_differential__g2_g3__PC1__partial_paired_nc__validationset__stats.Rds") |> 
      dplyr::select(probe_id, t) |> 
      dplyr::rename( t_validation_g2_g3__PC1 = t),
    by=c('probe_id'='probe_id')
  )|> 
  dplyr::left_join(
    readRDS("cache/analysis_differential__p_r__PC1__partial_paired_nc__validationset__stats.Rds") |> 
      dplyr::select(probe_id, t) |> 
      dplyr::rename(t_validation_p_r__PC1 = t),
    by=c('probe_id'='probe_id')
  ) |> 
    dplyr::filter(!is.na(DMP__primary_recurrence__pp_nc__t))


c = cor(plt |> tibble::column_to_rownames('probe_id'))

c = cor(plt |> dplyr::select(probe_id,
                             DMP__g2_g3__pp_nc__t, DMP__g2_g3__pp_nc_PC1__t,
                             t_validation_p_r,t_validation_p_r__PC1) |> tibble::column_to_rownames('probe_id'))


corrplot::corrplot(c, order="hclust", shade.lwd=0.5, tl.cex=0.4, cl.cex=0.4, tl.pos="l" )



# scatter GLASS-OD x validatie PC1 ----



