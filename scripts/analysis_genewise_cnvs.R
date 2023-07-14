#!/usr/bin/env R


library(ggplot2)


# load data ----


if(!exists('glass_od.metadata.idats')) {
  source('scripts/load_metadata.R')
}


# CDKN2A/B ----


plt <- glass_od.metadata.idats |> 
  dplyr::filter(is.na(reason_excluded_patient)) |> 
  dplyr::filter(is.na(reason_excluded_resection)) |> 
  dplyr::filter(is.na(reason_excluded_resection_isolation)) |> 
  dplyr::filter(is.na(reason_excluded_array_sample)) |>  # those are typically the replicates
  assertr::verify(!duplicated(resection_id)) |> 
  assertr::verify(sentrix_id != "204808700074_R04C01") |>  # 0003-R3-repA
  dplyr::filter(study_name == "GLASS-NL")


ggplot(plt, aes(x=`heidelberg_cnvp_C19MC`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_CCND1`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_CCND2`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_CDK4`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_CDK6`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_CDKN2A/B`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_EGFR`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_FGFR1/TACC1`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_FGFR3/TACC3`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_GLI2`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_KIAA1549/BRAF`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_MDM2`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_MDM4`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_MET`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_MGMT`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_MYB`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_MYBL1`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_MYC`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_MYCN`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_NF1`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_NF2 `, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_PDGFRA`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_PPM1D`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_PTCH1`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_PTEN`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_RB1`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_SMARCB1`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_TERT`, y=methylation_bins_1p19q_purity)) + geom_point()
ggplot(plt, aes(x=`heidelberg_cnvp_TP53`, y=methylation_bins_1p19q_purity)) + geom_point()





