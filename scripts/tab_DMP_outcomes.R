
# load ----



source('scripts/load_constants.R')
source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')
source('scripts/load_gene_annotations.R')



if(!exists('data.mvalues.probes')) {
  source('scripts/load_mvalues_hq_samples.R')
}


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}



# select and export ----


out <- data.mvalues.probes |> 
  dplyr::mutate(PCA_strongest_fit = dplyr::case_when(
    (abs(DMP__PC1_PC2_PC3_multivariate__t_PC1) > abs(DMP__PC1_PC2_PC3_multivariate__t_PC2)) & (abs(DMP__PC1_PC2_PC3_multivariate__t_PC1) > abs(DMP__PC1_PC2_PC3_multivariate__t_PC3))  ~ "PC1",
    (abs(DMP__PC1_PC2_PC3_multivariate__t_PC2) > abs(DMP__PC1_PC2_PC3_multivariate__t_PC1)) & (abs(DMP__PC1_PC2_PC3_multivariate__t_PC2) > abs(DMP__PC1_PC2_PC3_multivariate__t_PC3))  ~ "PC2",
    (abs(DMP__PC1_PC2_PC3_multivariate__t_PC3) > abs(DMP__PC1_PC2_PC3_multivariate__t_PC1)) & (abs(DMP__PC1_PC2_PC3_multivariate__t_PC3) > abs(DMP__PC1_PC2_PC3_multivariate__t_PC2))  ~ "PC3",
    T ~ "err"
  )) |> 
  dplyr::select(probe_id, gc_sequence_context_2_new,
                DMP__g2_g3__pp_nc__t,
                DMP__g2_g3__pp_nc_PC1__t,
                DMP__primary_recurrence__pp_nc__t,
                DMP__primary_recurrence__pp_nc_PC1__t,
                DMP__FFPE_log1p__decay_time__pp_nc__t,
                
                DMP__g2_g3__pp_nc__adj.P.Val,
                DMP__g2_g3__pp_nc_PC1__adj.P.Val,
                DMP__primary_recurrence__pp_nc__adj.P.Val,
                DMP__primary_recurrence__pp_nc_PC1__adj.P.Val,
                DMP__FFPE_log1p__decay_time__pp_nc__adj.P.Val,
                
                DMP__g2_g3__pp_nc__P.Value,
                DMP__g2_g3__pp_nc_PC1__P.Value,
                DMP__primary_recurrence__pp_nc__P.Value,
                DMP__primary_recurrence__pp_nc_PC1__P.Value,
                DMP__FFPE_log1p__decay_time__pp_nc__P.Value,
                
                DMP__g2_g3__pp_nc__logFC,
                DMP__g2_g3__pp_nc_PC1__logFC,
                DMP__primary_recurrence__pp_nc__logFC,
                DMP__primary_recurrence__pp_nc_PC1__logFC,
                DMP__FFPE_log1p__decay_time__pp_nc__logFC,
                
                PCA_strongest_fit
                )



write.table(out, "output/tables/tab_DMP_outcomes.txt", sep="\t", quote = F, row.names = F)



