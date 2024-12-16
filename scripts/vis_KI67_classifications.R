
# load ----


library(ggplot2)
library(patchwork)


source('scripts/load_palette.R')
source('scripts/load_functions.R')
source('scripts/load_themes.R')


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


## figs ----


plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |>
  dplyr::filter(!is.na(staining_KI67_pos_per_area_um2)) |> 
  dplyr::mutate(col = ifelse(
    array_mnp_predictBrain_v12.8_cal_class %in% c("GBM_MES_TYP", "GBM_RTK1"), "Other", array_mnp_predictBrain_v12.8_cal_class
  )) |> 
  dplyr::mutate(highlighted = staining_KI67_filename %in% c("KI67_09-5718 I.ndpi","KI67_20-9179 B.ndpi","KI67_22-3789 I-C.ndpi","KI67_H04-10641 I.ndpi"))



plt.facet <- rbind(
  plt |>
    dplyr::mutate(facet = "PC2") |> 
    dplyr::mutate(y = array_PC2),
  plt |>
    dplyr::mutate(facet = "PC3") |> 
    dplyr::mutate(y = array_PC3),
  plt |>
    dplyr::mutate(facet = "PC HorvathS2018 * -1") |> 
    dplyr::mutate(y = -1 * array_dnaMethyAge__PCHorvathS2018)
) |> 
  dplyr::mutate(facet = factor(facet, levels=c("PC2", "PC3", "PC HorvathS2018 * -1")))


plt.intervals <- 10^seq(log10(max(plt$staining_KI67_pos_per_area_um2)), log10(min(plt$staining_KI67_pos_per_area_um2)), length.out = 5) # linear intervals in log scaled labels


ggplot(plt.facet, aes(x=staining_KI67_pos_per_area_um2,
                y=y,
                col=col)) +
  facet_wrap(~facet, scales="free") +
  geom_point(size=theme_nature_size/3) + 
  ggpubr::stat_cor(method = "pearson",
                   aes(label = after_stat(r.label)),
                   label.x.npc = "right", hjust=1,
                   col="black",
                   size=theme_nature_size,
                   family = theme_nature_font_family) +
  scale_color_manual(values = c(palette_mnp_12.8_6, `Other`="darkgray")) +
  scale_x_log10() +
  labs(x = expression("KI67+ cells per um^2"),
       #y="PC HorvathS2018 * -1",
       col=NULL,
       subtitle=format_subtitle("KI67 correlations"),
       caption=paste0("n=",nrow(plt)," samples with KI67 staining and EPIC array data")) + 
  theme_nature




ggsave("output/figures/vis_KI67_classifications__scatterplots.pdf", width = 8.5*0.975 * (3/5), height=2.34)







# ggplot(glass_od.metadata.array_samples, aes(x=staining_KI67_lr_pos_neg_cells, 
#                                             y=array_dnaMethyAge__PCHorvathS2018)) +
#   geom_point()



# ggplot(glass_od.metadata.array_samples, aes(x=staining_KI67_lr_pos_neg_cells, y=array_PC2)) +
#   geom_point()


# ggplot(glass_od.metadata.array_samples, aes(x=staining_KI67_lr_pos_neg_cells, y=array_PC3)) +
#   geom_point()




ggplot(glass_od.metadata.array_samples, aes(x=log(staining_KI67_pos_per_detected_cells), y=array_PC2)) +
  geom_point()


ggplot(glass_od.metadata.array_samples, aes(x=log(staining_KI67_pos_per_detected_cells), y=array_PC3)) +
  geom_point()



ggplot(glass_od.metadata.array_samples, aes(x=staining_KI67_lr_pos_neg_cells, y=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) +
  geom_point()


ggplot(glass_od.metadata.array_samples, aes(x=log(staining_KI67_pos_per_detected_cells), y=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) +
  geom_point()


ggplot(glass_od.metadata.array_samples, aes(x=log(staining_KI67_pos_per_area_um2), y=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) +
  geom_point()


ggplot(glass_od.metadata.array_samples, aes(
  x=log(staining_KI67_Num.KI67pos / (staining_KI67_Area_um_2)),
  y=array_PC3)) +
  geom_point()



## corr ----



dat <- glass_od.metadata.array_samples |> 
  dplyr::mutate(log_staining_KI67_pos_per_area_um2 = log(staining_KI67_pos_per_area_um2)) |> 
  #dplyr::mutate(log_staining_KI67_pos_neg_cell_density = log(staining_KI67_pos_neg_cell_density)) |> 
  dplyr::select(
                staining_KI67_lr_pos_neg_cells,
                log_staining_KI67_pos_per_area_um2,
                #staining_KI67_pos_neg_cell_density,
                #log_staining_KI67_pos_neg_cell_density,
                
                array_PC2, array_PC3,
                array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, array_GLASS_NL_g2_g3_sig,
                array_dnaMethyAge__PCHorvathS2018) |> 
  dplyr::filter(!is.na(staining_KI67_lr_pos_neg_cells) & !is.na(array_PC3)) |> 
  dplyr::mutate(`-1 * array_PC2` = -1 * array_PC2, array_PC2=NULL) |> 
  dplyr::mutate(`-1 * array_PC3` = -1 * array_PC3, array_PC3=NULL) |> 
  dplyr::mutate(`-1 * array_GLASS_NL_g2_g3_sig` = -1 * array_GLASS_NL_g2_g3_sig, array_GLASS_NL_g2_g3_sig=NULL)




pdf("output/figures/vis_KI67_classifications__corrplot.pdf",width=2.5, height=2.15)



corrplot::corrplot(cor(dat), order="hclust", tl.cex=0.52, tl.pos="l", cl.cex=0.52,
                   

                   
                   title = paste0("\nsample: n=", nrow(dat))
                   )



dev.off()





## cohort overview ----


plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(!is.na(staining_KI67_lr_pos_neg_cells)) |> 
  dplyr::mutate(resection_grade = paste0("WHO Grade ", resection_tumor_grade)) |> 
  dplyr::mutate(resection_type = ifelse(resection_number == 1, "primary", "recurrent")) |> 
  dplyr::mutate(rank = order(order(staining_KI67_pos_per_detected_cells))) |> 
  dplyr::mutate(col = ifelse(
    array_mnp_predictBrain_v12.8_cal_class %in% c("GBM_MES_TYP", "GBM_RTK1"), "Other", array_mnp_predictBrain_v12.8_cal_class
  )) |> 
  dplyr::filter(!is.na(staining_KI67_pos_per_detected_cells)) |> 
  dplyr::mutate(highlighted = staining_KI67_filename %in% c("KI67_09-5718 I.ndpi","KI67_20-9179 B.ndpi","KI67_22-3789 I-C.ndpi","KI67_H04-10641 I.ndpi"))


# plt |> dplyr::select(resection_id, resection_number, staining_KI67_filename, staining_KI67_pos_per_detected_cells) |> 
#   dplyr::mutate(pct = staining_KI67_pos_per_detected_cells * 100) |> 
#   View()



ggplot(plt, aes(x=reorder(resection_id, rank),
                y=staining_KI67_pos_per_detected_cells,
                col=col,
                label=paste0(round(staining_KI67_pos_per_detected_cells * 100,1),"%"))) +
  geom_hline(yintercept=0.05, col="red", lty=2,lwd=theme_nature_lwd) +
  geom_hline(yintercept=0.1, col="red", lty=2,lwd=theme_nature_lwd) +
  scale_y_continuous(labels = scales::percent) +
  geom_point(size=theme_nature_size/3) +
  ggrepel::geom_text_repel(data=subset(plt, highlighted),
                           col="black",
                           size=theme_nature_size,
                           family=theme_nature_font_family,
                           segment.size=theme_nature_lwd,
                           nudge_y = 0.075) +
  scale_color_manual(values = c(palette_mnp_12.8_6, `Other`="darkgray")) +
  theme_nature +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x='resection', y="Percentage cells classified Ki67+", col=NULL)

ggsave("output/figures/vis_KI67_classifications.pdf", width=0.975 * 8.5, height=2.25)





## who grade ----




plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::mutate(resection_grade = paste0("WHO Grade ", resection_tumor_grade)) |> 
  dplyr::filter(!is.na(staining_KI67_lr_pos_neg_cells))


ggplot(plt, aes(x=resection_grade, y=staining_KI67_lr_pos_neg_cells)) +
  ggbeeswarm::geom_quasirandom(size=theme_nature_size/3,width=0.175) + 
  labs(subtitle=format_subtitle("KI67 x WHO Grade"),
       caption=paste0("samples  --  ",
                      "Grade 2: n=", (plt |> dplyr::filter(resection_grade == "WHO Grade 2") |>  nrow()), "  --  ",
                      "Grade 3: n=", (plt |> dplyr::filter(resection_grade == "WHO Grade 3") |>  nrow()
       ))) +
  ggpubr::stat_compare_means(label.x.npc=0.5, method = "t.test", show_guide  = FALSE,  size=theme_nature_size, col="darkgray") +
  labs(y="log(Ki-67 positive cells / Ki-67 negative cells)", x=NULL) +
  theme_nature



ggsave("output/figures/vis_KI67_classifications__who_grade.pdf",width=1.5, height=2.25)





## prim rec ----



plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::mutate(resection_type = ifelse(resection_number == 1, "primary", "recurrent")) |> 
  dplyr::filter(!is.na(staining_KI67_lr_pos_neg_cells))


ggplot(plt, aes(x=resection_type, y=staining_KI67_lr_pos_neg_cells)) +
  ggbeeswarm::geom_quasirandom(size=theme_nature_size/3,width=0.175) + 
  labs(subtitle=format_subtitle("KI67 x resection type"),
       caption=paste0("samples  --  ",
                      "primary: n=",   (plt |> dplyr::filter(resection_type == "primary") |>    nrow()),"  --  ",
                      "recurrent: n=", (plt |> dplyr::filter(resection_type == "recurrent") |>  nrow()
       ))) +
  ggpubr::stat_compare_means( label.x.npc=0.5, method = "t.test", show_guide  = FALSE,  size=theme_nature_size, col="darkgray") +
  labs(y="log(Ki-67 positive cells / Ki-67 negative cells)", x=NULL) +
  theme_nature


ggsave("output/figures/vis_KI67_classifications__prim_rec.pdf",width=1.5, height=2.25)





