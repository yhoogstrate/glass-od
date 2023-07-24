#!/usr/bin/env R


# load ----

library(ggplot2)
library(patchwork)


if(!exists('youri_gg_theme')) {
  source('scripts/youri_gg_theme.R')
}


if(!exists('glass_od.metadata.idats')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


source('scripts/load_functions.R')


# plot ----



plt <-
  rbind(
    glass_od.metadata.idats |> 
      filter_GLASS_OD_idats(163) 
    ,
  glass_od.metadata.idats |> 
    dplyr::filter(study_name != "GLASS-OD")
  ) |> 
  dplyr::mutate(patient_id = as.factor(patient_id)) |> 
  dplyr::mutate(x = as.numeric(patient_id) + 1)|> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(last_recurrence = resection_number != 1 & resection_number == max(resection_number)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(panel = dplyr::recode(study_name, `CATNON` = 'non-canonical codels [CATNON, ...]', `GLASS-OD`='canonical codels [GLASS-OD]')) 



p1 <- ggplot(plt, aes(x = x, y=resection_number, fill=A_IDH_HG__A_IDH_LG_lr )) +
  facet_grid(cols = vars(panel)) +
  geom_rect(aes(xmin=x - 0.45,xmax=x + 0.45,
                ymin=resection_number + 0.5 - 0.45, ymax=resection_number + 0.5 + 0.45)) +
  theme_bw() +
  scale_fill_gradient2(low = "darkgreen", mid="gray90", high = "red", midpoint=0) 


p2 <- ggplot(plt, aes(x = x, y=A_IDH_HG__A_IDH_LG_lr, group=patient_id )) +
  facet_grid(cols = vars(panel)) +
  geom_line(alpha=0.3) +
  geom_point(aes(col=resection_number == 1, shape=resection_number == 1)) +
  theme_bw()


p1 / p2


