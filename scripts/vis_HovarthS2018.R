#!/usr/bin/env R


# load ----


source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')


library(ggplot2)
library(patchwork)


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# GLASS-OD / OD ----


metadata <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES)



# tab ----


# metadata |>
#   dplyr::select(patient_id, time_between_birth_and_resection) |> 
#   dplyr::mutate(as_year = as.difftime(time_between_birth_and_resection, "years") )


plt <- metadata |>
  dplyr::select(patient_id, resection_id, resection_number, contains("HorvathS2018"), time_between_birth_and_resection, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, resection_tumor_grade) |> 
  dplyr::mutate(time_between_birth_and_resection = as.numeric(time_between_birth_and_resection / 365.25) ) |> 
  dplyr::mutate(col = cut_number(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, 4)) |> 
  dplyr::mutate(col2 = cut(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, 5)) |> 
  dplyr::mutate(prim_rec = ifelse(resection_number == 1, "primary", "recurrent")) |> 
  dplyr::mutate(resection_tumor_grade = as.factor(paste0("Grade ",resection_tumor_grade))) |> 
  dplyr::mutate(delta_t = array_dnaMethyAge__PCHorvathS2018 -  time_between_birth_and_resection)


p1 = ggplot(plt, aes(x = time_between_birth_and_resection, y=array_dnaMethyAge__PCHorvathS2018, col=col)) +
  geom_abline(intercept = 1, lwd=theme_nature_lwd, lty=2) +
  geom_smooth(method="lm", se=F,  lwd=theme_nature_lwd,  formula= y~x) +
  geom_point(size=theme_nature_size/3) +
  
  labs(x = "Age at resection", y = "PCHorvathS2018", col="CGC") +
  scale_color_manual(values = rev(col4(4))) +
  theme_nature



p2 = ggplot(plt, aes(x = time_between_birth_and_resection, y=array_dnaMethyAge__PCHorvathS2018, col=resection_tumor_grade)) +
  geom_abline(intercept = 1, lwd=theme_nature_lwd, lty=2) +
  geom_smooth(method="lm", se=F,  lwd=theme_nature_lwd,  formula= y~x) +
  geom_point(size=theme_nature_size/3) +
  
  labs(x = "Age at resection", y = "PCHorvathS2018 (Epi-genetic age)", col=NULL) +
  scale_color_manual(values=c('#009E74', '#CB75A4')) +
  theme_nature


p3 = ggplot(plt, aes(x = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=(array_dnaMethyAge__PCHorvathS2018 -  time_between_birth_and_resection), col=resection_tumor_grade)) +
  geom_hline(yintercept=0, lwd=theme_nature_lwd, lty=2) +
  geom_point(size=theme_nature_size/3) +
  scale_color_manual(values=c('#009E74', '#CB75A4')) +
  labs(x = "CGC", y = "Epi-genetic age  beyond  actual age", col=NULL) +
  theme_nature


library(patchwork)
p1 + p2 + p3


ggsave("output/figures/vis_HovarthS2018__x_tissue_age.pdf", width=8.5*0.975*2/3, height=2.2)



