#!/usr/bin/env R 


# load data ----


source('scripts/load_functions.R')
source('scripts/load_themes.R')
source('scripts/load_palette.R')


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}



# including this may have some odd biological implications
# if(!exists('gsam.metadata.array_samples')) {
#   source('scripts/load_G-SAM_metadata.R')
# }


# plot LASSO 10xFC training ----


cv_model_probe_based <- readRDS("cache/analysis_A_IDH_HG__A_IDH_LG_lr__lasso_fit__probe_based__train_paramters.Rds")
#cv_model_probe_based  <- readRDS("cache/LGC_predictor_probe_based_lm.Rds")

plot(cv_model_probe_based)
# @todo export & ggplot



plt <- data.frame(
  lambda = cv_model_probe_based$lambda,
  cvm = cv_model_probe_based$cvm,
  cvlo = cv_model_probe_based$cvlo,
  cvup = cv_model_probe_based$cvup,
  nzero = cv_model_probe_based$nzero
) |> 
  dplyr::mutate(i = 1:dplyr::n()) |> 
  dplyr::mutate(show_label = i %in% c(1,10,20,30,40,50,60,70,80,90,100))


ggplot(plt, aes(x=log(lambda), y=cvm, label=nzero)) +
  geom_ribbon(aes(ymin=cvlo, ymax=cvup),alpha=0.1) +
  geom_vline(xintercept=log(cv_model_probe_based$lambda.min), col="gray50", lty=2, lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=log(cv_model_probe_based$lambda.1se), col="gray50", lty=2, lwd=theme_cellpress_lwd) +
  geom_point(cex=0.4, col="red") +
  geom_text(data=subset(plt, show_label == T),y=50, size=theme_cellpress_size) +
  scale_y_continuous(limits = c(0, 50)) +
  labs(x = expression(Log(lambda)), y="Mean-Squared Error",
       subtitle=format_subtitle("LASSO parameter fit")
       ) +
  theme_nature


ggsave("output/figures/vis_AcCGAP-OD_x_AcCGAP-10xCV-AC_LASSO_10xFC_training.pdf", width=(8.5*0.95)/3, height=2.5)



# plot training vs predicted scatter ----


#"array_A_IDH_HG__A_IDH_lr"               
#"array_A_IDH_HG__A_IDH_LG_lr_v12.8"   
#"array_A_IDH_HG__A_IDH_LG_lr__lasso_fit"


plt <- glass_nl.metadata.array_samples |>
  filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES) |> 
  dplyr::mutate(col = dplyr::case_when(
    array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_LG", "A_IDH_HG") ~ array_mnp_predictBrain_v12.8_cal_class,
    T ~ "other"
  ))


ggplot(plt, aes(
  y=`array_A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV_v12.8`,
  x=array_A_IDH_HG__A_IDH_LG_lr_v12.8,
  col=col,
  #label=sentrix_id
)) +
  geom_point(size=theme_nature_size/3) +
  #ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", fill="1", size=theme_cellpress_size) +
  labs(subtitle = "CGC[phi] in GLASS-NL",
       x= "Actual CGC",
       y="CGC[Phi] (predicted)", col="", fill="") +
  scale_color_manual(values = c(`A_IDH_LG`='lightblue',`A_IDH_HG`='darkblue', `other`=mixcol("#444444","red",0.75))) +
  theme_nature



ggsave("output/figures/vis_AcCGAP-OD_x_AcCGAP-10xCV-AC_LASSO_10xFC_scatter_observed_predicted.pdf",
       width=2.3725, height=2.65)




# plot training (AC) vs predicted (OLI) ----


plt <- rbind(
  glass_od.metadata.array_samples |>
    filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |>
    dplyr::select(
      array_sentrix_id,
      array_A_IDH_HG__A_IDH_LG_lr__lasso_fit,
      array_mnp_predictBrain_v12.8_cal_class
    ) |>
    dplyr::rename(AcCGAP_score = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |>
    dplyr::mutate(
      type = "LASSO full",
      dataset = "GLASS-OD"
    )
  ,
  glass_nl.metadata.array_samples |>
    filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES) |>
    dplyr::select(
      array_sentrix_id,
      `array_A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV_v12.8`,
      array_mnp_predictBrain_v12.8_cal_class
    ) |>
    dplyr::rename(AcCGAP_score = `array_A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV_v12.8`) |>
    dplyr::mutate(
      type = "LASSO 10xCV",
      dataset = "GLASS-NL"
    ) ,
  glass_nl.metadata.array_samples |>
    filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES) |>
    dplyr::select(
      array_sentrix_id,
      `array_A_IDH_HG__A_IDH_LG_lr_v12.8`,
      array_mnp_predictBrain_v12.8_cal_class
    ) |>
    dplyr::rename(AcCGAP_score = `array_A_IDH_HG__A_IDH_LG_lr_v12.8`) |>
    dplyr::mutate(
      type = "CGC",
      dataset = "GLASS-NL"
    ) 
  )|> 
  dplyr::mutate(class = paste0(dataset, ":\n", type))


ggplot(plt, aes(x=class, y=AcCGAP_score, fill=class, col=class)) +
  ggbeeswarm::geom_quasirandom(size=theme_nature_size/3, show.legend = FALSE) +
  ggpubr::stat_compare_means(comparisons=list(c('GLASS-NL:\nLASSO 10xCV',
                                                'GLASS-OD:\nLASSO full')), method="wilcox.test",
    label.x.npc=theme_nature_lwd, size=theme_nature_size,
    family=theme_nature_font_family, 
    show_guide = FALSE) +
  labs(col="", fill="",x = "", y="AcCGAP score", subtitle=format_subtitle("CGC LASSO"),caption=paste0(
    "GLASS-NL: n=",CONST_N_GLASS_NL_INCLUDED_SAMPLES, "  --  GLASS-OD: n=",CONST_N_GLASS_OD_INCLUDED_SAMPLES, "  --  p: wilcox test"
  )) +
  scale_color_manual(values = c('GLASS-NL:\nLASSO 10xCV'='darkblue',
                                'GLASS-OD:\nLASSO full'='darkgreen',
                                'GLASS-NL:\nCGC' = '#888888')) +
  theme_nature


ggsave("output/figures/vis_AcCGAP-OD_x_AcCGAP-10xCV-AC_LASSO_10xFC_AC_OLI_wilcox.pdf",
       width=(8.5*0.95)/3, 
       height=2.69)


