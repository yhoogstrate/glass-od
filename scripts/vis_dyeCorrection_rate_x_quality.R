#!/usr/bin/env R


# load ----


library(ggplot2)
library(patchwork)



source('scripts/load_constants.R')
source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')



if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}



# plt ----



plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |>
  dplyr::select(contains("detP") | contains("minfi") | contains("sentrix")) |> 
  dplyr::mutate(lr_red_green = log(array_minfi_dyeCorrection_Red.avg / array_minfi_dyeCorrection_Green.avg)) |> 
  dplyr::mutate(A_ACTG = log(array_minfi_dyeCorrection_A.avg / ((array_minfi_dyeCorrection_A.avg + array_minfi_dyeCorrection_C.avg + array_minfi_dyeCorrection_T.avg + array_minfi_dyeCorrection_G.avg)/4))) |> 
  dplyr::mutate(C_ACTG = log(array_minfi_dyeCorrection_C.avg / ((array_minfi_dyeCorrection_A.avg + array_minfi_dyeCorrection_C.avg + array_minfi_dyeCorrection_T.avg + array_minfi_dyeCorrection_G.avg)/4))) |> 
  dplyr::mutate(T_ACTG = log(array_minfi_dyeCorrection_T.avg / ((array_minfi_dyeCorrection_A.avg + array_minfi_dyeCorrection_C.avg + array_minfi_dyeCorrection_T.avg + array_minfi_dyeCorrection_G.avg)/4))) |> 
  dplyr::mutate(G_ACTG = log(array_minfi_dyeCorrection_G.avg / ((array_minfi_dyeCorrection_A.avg + array_minfi_dyeCorrection_C.avg + array_minfi_dyeCorrection_T.avg + array_minfi_dyeCorrection_G.avg)/4)))



ggplot(plt, aes(
  x = log(array_percentage.detP.signi/100)/(1-(array_percentage.detP.signi/100)),
  #x=log1p(array_percentage.detP.signi),
                y=lr_red_green
)) +
  geom_hline(yintercept=0, col="red", lty=2, lwd=theme_nature_lwd) +
  geom_point(size=theme_nature_size/3) +
  
  #geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.225), lwd=theme_nature_lwd * 7) +
  #geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.225), lwd=theme_nature_lwd * 6) +
  #geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.225), lwd=theme_nature_lwd * 5) +
  geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.3), lwd=theme_nature_lwd * 4) +
  geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.3), lwd=theme_nature_lwd * 3) +
  geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.3), lwd=theme_nature_lwd * 2) +
  geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.3), lwd=theme_nature_lwd * 1) +
  geom_smooth(method='lm', formula= y~x, se = F, lty=1, col="#6ba6e5", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)),
                   #col="#6ba6e5",
                   family=theme_nature_font_family,
                   size=theme_nature_size) +
  
  scale_y_continuous(limits=c(-2,2)) +
  labs(x = "sample quality: log(detP failed / detP passed)",
       y = "log(Red intensity / Green intensity) ctrl probes\n<- [more green]                         [more red] ->",
       subtitle=format_subtitle("det-P x dyeCorrection"),
       caption = paste0("n=",nrow(plt), " GLASS-OD samples")
       ) +
  theme_nature


ggsave("output/figures/vis_dyeCorrection_rate_x_quality.pdf", width = 8.5 * 0.975 * 1/4, height = 2.2)






plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES)  |> 
  dplyr::mutate(LR_detP = log(array_percentage.detP.signi/100)/(1-(array_percentage.detP.signi/100))) |> 
  dplyr::mutate(resection_type = ifelse(resection_number == 1, "primary", "recurrent"))


ggplot(plt, aes(x=resection_type, y=LR_detP)) +
  ggbeeswarm::geom_quasirandom() +
  ggpubr::stat_compare_means(label.x.npc=0.25, method = "wilcox.test", show_guide  = FALSE,  size=theme_nature_size) + 
  theme_nature



# ggplot(plt, aes(x=array_percentage.detP.signi,
#                 y=A_ACTG)) +
#   geom_point() +
#   coord_trans(x = "log1p") +
#   theme_nature
# 
# ggplot(plt, aes(x=array_percentage.detP.signi,
#                 y=T_ACTG)) +
#   geom_point() +
#   coord_trans(x = "log1p") +
#   theme_nature
# 
# 
# ggplot(plt, aes(x=array_percentage.detP.signi,
#                 y=G_ACTG)) +
#   geom_point() +
#   coord_trans(x = "log1p") +
#   theme_nature
# 
# 
# ggplot(plt, aes(x=array_percentage.detP.signi,
#                 y=C_ACTG)) +
#   geom_point() +
#   coord_trans(x = "log1p") +
#   theme_nature
# 


