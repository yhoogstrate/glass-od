#!/usr/bin/env R 


# load data ----


if(!exists('glass_od.metadata.idats')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


if(!exists('glass_nl.metadata.idats')) {
  source('scripts/load_GLASS-NL_metadata.R')
}



# including this may have some odd biological implications
# if(!exists('gsam.metadata.idats')) {
#   source('scripts/load_G-SAM_metadata.R')
# }

# plot ----

plt <- rbind(
  glass_od.metadata.idats |>
    filter_GLASS_OD_idats(163) |>
    dplyr::select(
      sentrix_id,
      A_IDH_HG__A_IDH_LG_lr__lasso_fit,
      mnp_predictBrain_v12.8_cal_class
    ) |>
    dplyr::rename(AcCGAP_score = A_IDH_HG__A_IDH_LG_lr__lasso_fit) |>
    dplyr::mutate(
      type = "predictions by full AcCGAP model",
      dataset = "GLASS-OD"
    ),
  glass_nl.metadata.idats |>
    filter_GLASS_NL_idats(218) |>
    dplyr::select(
      sentrix_id,
      A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV,
      mnp_predictBrain_v12.8_cal_class
    ) |>
    dplyr::rename(AcCGAP_score = A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV) |>
    dplyr::mutate(
      type = "predictions by 10xCV AcCGAP model",
      dataset = "GLASS-NL"
    )
)


ggplot(plt, aes(x=dataset, y=AcCGAP_score, col=type)) +
  ggbeeswarm::geom_quasirandom() +
  ggpubr::stat_compare_means(label.x.npc=0.4) +
  theme_bw()


