#!/usr/bin/env R


# load ----


library(ggplot2)
library(patchwork)



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





if(!exists('data.intensities.probes')) {
  source('scripts/load_intensities_hq_samples.R')
}




# plots ----

## 00 statistical power plot ----


data.mvalues.probes |> 
  dplyr::filter(!is.na(DMP__primary_recurrence__pp_nc__adj.P.Val)) |> 
  dplyr::mutate(significant = DMP__primary_recurrence__pp_nc__adj.P.Val < 0.01) |> 
  dplyr::pull(significant) |> 
  table()

data.mvalues.probes |> 
  dplyr::filter(!is.na(DMP__g2_g3__pp_nc__P.Value)) |> 
  dplyr::mutate(significant = DMP__g2_g3__pp_nc__P.Value < 0.01) |> 
  dplyr::pull(significant) |> 
  table()

(193478 / (193478  + 491687))   /    (44375 / ( 640790 + 44375  ))




plt <- rbind(
  data.mvalues.probes |> 
    dplyr::filter(!is.na(DMP__g2_g3__pp_nc__P.Value)) |> 
    dplyr::mutate(significant = DMP__g2_g3__pp_nc__P.Value < 0.01) |> 
    dplyr::select(significant) |> 
    dplyr::mutate(type = "WHO Grade")
  ,
  data.mvalues.probes |> 
    dplyr::filter(!is.na(DMP__primary_recurrence__pp_nc__adj.P.Val)) |> 
    dplyr::mutate(significant = DMP__primary_recurrence__pp_nc__adj.P.Val < 0.01) |> 
    dplyr::select(significant) |> 
    dplyr::mutate(type = "Resection")
) |> 
  dplyr::mutate(significant = ifelse(significant, "< 0.01" , ">= 0.01"))



ggplot(plt, aes(fill=significant, x=type)) +
  geom_bar() + 
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) +
  scale_fill_manual(values=c("< 0.01"= mixcol('darkblue',  'darkgreen'),
                             ">= 0.01"=mixcol('lightblue', 'lightgreen'))) +
  theme_nature +
  labs(fill="adj. P-value", y = "tested probes", x="", format_subtitle("DMP power"))

ggsave("output/figures/vis_differential__DMP_power.pdf", width=8.5 * 0.975 / 6, height = 1.4)


## A: overall density ----

### 1. gray scale ----

n_samples_grade <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_first_G2_and_last_G3(156) |> 
  nrow()

n_samples_prim_rec <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(180) |> 
  nrow()


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t)) +
  #geom_vline(xintercept=0, col="red") +
  #geom_hline(yintercept=0, col="red") +
  #geom_vline(xintercept=0, col="red", alpha=0.1) + # do in illustrator
  #geom_hline(yintercept=0, col="red", alpha=0.1) + # do in illustrator
  
  geom_point(pch=16, cex=0.001, alpha=0.0085, col="black") +
  
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)), col="black", size=theme_nature_size, label.x=5) +
  
  
  labs(x="Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent",
                caption=paste0("Included samples per test: n=",n_samples_grade, " (grade), n=",n_samples_prim_rec," (resection type)")) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))))+
  coord_cartesian(ylim=c(-max(abs(plt$DMP__primary_recurrence__pp_nc__t)), max(abs(plt$DMP__primary_recurrence__pp_nc__t))))+
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


ggsave(paste0("output/figures/vis_differential__overall_density.png"), width=(8.5 * 0.975 * (2/5)), height=3.48, dpi=300)


rm(n_samples_grade, n_samples_prim_rec)




### 2. EFF1 & EFF2 ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t,
                y=DMP__primary_recurrence__pp_nc__t,
                col = abs(DMP__PCs__pp_nc__PC1_t) - abs(DMP__PCs__pp_nc__PC2_t)
                )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3", y="Per probe t-score Primary ~ Recurrent",col = "Fit to PC1 vs. PC2") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col5(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


ggsave(paste0("output/figures/vis_differential__EFF1-PC1_EFF2-PC2.png"), width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))



rm(n_samples_grade, n_samples_prim_rec)



### 3. PC1 / qual corrected ----


n_samples_grade <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_first_G2_and_last_G3(156) |> 
  nrow()

n_samples_prim_rec <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(180) |> 
  nrow()



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                y=DMP__primary_recurrence__pp_nc_PC1__t)) +
  
  #geom_vline(xintercept=0, col="red") +
  #geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.0085, col="black") +
  
  #geom_vline(xintercept=0, col="red", alpha=0.1) +
  #geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  #geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.5), lwd=theme_nature_lwd) +
  #geom_smooth(method='lm', formula= y~x, se = F, lty=2, col="#6ba6e5", lwd=theme_nature_lwd) +
  
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)), col="black", size=theme_nature_size, label.x=5) +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent",
       caption=paste0("Included samples per test: n=",n_samples_grade, " (grade), n=",n_samples_prim_rec," (resection type)")
  ) +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))))+
  coord_cartesian(ylim=c(-max(abs(plt$DMP__primary_recurrence__pp_nc__t)), max(abs(plt$DMP__primary_recurrence__pp_nc__t)))) +
  
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(0, 1), oob = scales::squish)
  #theme(plot.background = element_rect(fill="white", colour=NA))  # png export




ggsave(paste0("output/figures/vis_differential__overall_density__quality_corrected.png"), width=(8.5 * 0.975 * (2/5)), height=3.48, dpi=300)



rm(n_samples_grade, n_samples_prim_rec)



### 4. PC1 / qual corrected HOX CC etc. ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()




## B: detP fraction ----



n_samples_grade <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_first_G2_and_last_G3(156) |> 
  nrow()

n_samples_prim_rec <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(180) |> 
  nrow()


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::mutate(col = DMP__pct_detP_signi__pp_nc__t)



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t,
                color=col)) +
  #geom_vline(xintercept=0, col="red", alpha=0.1) + # do in illustrator
  #geom_hline(yintercept=0, col="red", alpha=0.1) + # do in illustrator
  
  geom_point(pch=16, cex=0.001 , alpha=0.15) + 
  
  labs(x="Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent",
       col="",
       caption=paste0("Included samples per test: n=",n_samples_grade, " (grade), n=",n_samples_prim_rec," (resection type)")) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))))+
  coord_cartesian(ylim=c(-max(abs(plt$DMP__primary_recurrence__pp_nc__t)), max(abs(plt$DMP__primary_recurrence__pp_nc__t))))+
  
  #theme(plot.background = element_rect(fill="white", colour = NA)) + # png export
  
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
  
  theme(axis.line =     element_line(linewidth = theme_nature_lwd / 2)) # somehow w/ png export




# 871 x 872
# 436 x 436
ggsave(paste0("output/figures/vis_differential__detP.png"), width=(8.5 * 0.975 * (1/5) * 1.12425), height=2.428, dpi=600)





## B: PC1 - PC8 multivariate ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


plt.expanded <- rbind(
  plt |> 
    dplyr::mutate(col = -DMP__PCs__pp_nc__PC1_t) |> 
    dplyr::mutate(grid_x = 1) |>
    dplyr::mutate(grid_y = 1)
  ,
  plt |> 
    dplyr::mutate(col = DMP__PCs__pp_nc__PC2_t) |> 
    dplyr::mutate(grid_x = 2) |>
    dplyr::mutate(grid_y = 1)
  ,
  plt |> 
    dplyr::mutate(col = DMP__PCs__pp_nc__PC3_t) |> 
    dplyr::mutate(grid_x = 3) |>
    dplyr::mutate(grid_y = 1)
  ,
  plt |> 
    dplyr::mutate(col = -DMP__PCs__pp_nc__PC4_t) |> 
    dplyr::mutate(grid_x = 4) |>
    dplyr::mutate(grid_y = 1)
  ,
  plt |> 
    dplyr::mutate(col = -DMP__PCs__pp_nc__PC5_t) |> 
    dplyr::mutate(grid_x = 1) |>
    dplyr::mutate(grid_y = 2)
  ,
  plt |> 
    dplyr::mutate(col = -DMP__PCs__pp_nc__PC6_t) |> 
    dplyr::mutate(grid_x = 2) |>
    dplyr::mutate(grid_y = 2)
  ,
  plt |> 
    dplyr::mutate(col = -DMP__PCs__pp_nc__PC7_t) |> 
    dplyr::mutate(grid_x = 3) |>
    dplyr::mutate(grid_y = 2)
  ,
  plt |> 
    dplyr::mutate(col = -DMP__PCs__pp_nc__PC8_t) |> 
    dplyr::mutate(grid_x = 4) |>
    dplyr::mutate(grid_y = 2)
)


p_PCs <- ggplot(plt.expanded, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t,
                         color=col)) +
  facet_grid(cols = vars(grid_x), rows = vars(grid_y)) +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent",
       col="Association PC1") +
  
  #theme_nature +
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


p_PCs

corrplot::corrplot(abs(cor(
  glass_od.metadata.array_samples |> 
    dplyr::filter(!is.na(array_PC1)) |> 
    dplyr::filter(!is.na(array_methylation_bins_1p19q_purity)) |> 
    dplyr::select(array_PC1, array_PC2, array_PC3, array_PC4,
                  array_A_IDH_HG__A_IDH_LG_lr__lasso_fit,
                  array_percentage.detP.signi,
                  array_methylation_bins_1p19q_purity,
                  array_epiTOC2_tnsc
                  )
  
,method="spearman")))


plot(
  glass_od.metadata.array_samples$array_methylation_bins_1p19q_purity,
  glass_od.metadata.array_samples$array_PC1
)

plot(
  glass_od.metadata.array_samples$array_methylation_bins_1p19q_purity,
  glass_od.metadata.array_samples$array_PC2 # purity?
)

plot(
  glass_od.metadata.array_samples$array_methylation_bins_1p19q_purity,
  glass_od.metadata.array_samples$array_PC3
)





plot(
  glass_od.metadata.array_samples$array_A_IDH_HG__A_IDH_LG_lr__lasso_fit,
  glass_od.metadata.array_samples$array_PC1
)

plot(
  glass_od.metadata.array_samples$array_A_IDH_HG__A_IDH_LG_lr__lasso_fit,
  glass_od.metadata.array_samples$array_PC2 # CGC
)

plot(
  glass_od.metadata.array_samples$array_A_IDH_HG__A_IDH_LG_lr__lasso_fit,
  glass_od.metadata.array_samples$array_PC3
)




## C: AcCGAP in OD ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


p_B <- ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t,
                       col=DMP__AcCGAP__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3", y="Per probe t-score Primary ~ Recurrent", col="Association A_IDH_LG ~ A_IDH_HG") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)




## D: GLASS-NL probes ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


p_C <- ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=glass_nl_prim_rec__deep_significant)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == F),pch=16, cex=0.001, alpha=0.15, show.legend=F) + 
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == T),pch=16, cex=0.1, alpha=0.35, show.legend=F) + 
  geom_point(data = head(plt, n=0),pch=16, show.legend=T) + 

  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3", y="Per probe t-score Primary ~ Recurrent", col="Significant GLASS-NL\nPrimary ~ Recurrent") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_manual(values=c('TRUE'= col3(3)[1], 'FALSE'='gray80'))


p_A + p_B + p_C


#ggsave("/tmp/papbpc.png",width = 8.5*0.975, height=3.2, dpi=600)


## E: GLASS-OD x GLASS-NL + PC1 ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                y=DMP__GLASS_NL__g2.3_g4__pp_nc_PC1__t)) +
  geom_vline(xintercept=0, col="red", lwd=theme_nature_lwd) +
  geom_hline(yintercept=0, col="red", lwd=theme_nature_lwd) +
  
  geom_point(pch=16, cex=0.001, alpha=0.01, col="black") +  
  
  geom_vline(xintercept=0, col="red", alpha=0.1, lwd=theme_nature_lwd) +
  geom_hline(yintercept=0, col="red", alpha=0.1, lwd=theme_nature_lwd) +
  
  #geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.1), lwd=theme_nature_lwd * 5) +
  geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.1), lwd=theme_nature_lwd * 4) +
  #geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.1), lwd=theme_nature_lwd * 3) +
  geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.1), lwd=theme_nature_lwd * 2) +
  geom_smooth(method='lm', formula= y~x, se = F, lty=2, col="#6ba6e5", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)), col="#6ba6e5", size=theme_nature_size, family = theme_nature_font_family) +
  
  labs(x = "Per probe t-score GLASS-OD Grade 2 ~ Grade 3",
       y = "Per probe t-score GLASS-NL Grade 2 & 3 ~ Grade 4",
       subtitle=format_subtitle("Overlap outcome oligodendroglioma vs. astrocytoma"),
       caption="Overlap Grade associated differences oligodendroglioma & astrocytoma (patient and PC1 corrected)"
  ) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  #scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(0, 1), oob = scales::squish) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export

ggsave("output/figures/vis_differential__GLASS_NL__x__GLASS_OD__grade__PC1.png", width=(8.5 * 0.97 / 4), height=(8.5 * 0.97 / 4))





## G: chromosomal ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |>
  (function(.) {
   print(dim(.))
   assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
   return(.)
  })() |>
  
  dplyr::filter(CHR_hg38 %in% c("chrX", "chrY") == F ) |> 
  dplyr::mutate(chr_arm = dplyr::case_when(
    CHR_hg38 == "chr1" & Start_hg38 <= 123400000 ~ "p",
    CHR_hg38 == "chr1" & Start_hg38 >  123400000 ~ "q",
    
    CHR_hg38 == "chr10" & Start_hg38 <= 39800000 ~ "p", 
    CHR_hg38 == "chr10" & Start_hg38 > 39800000 ~ "q", 
    
    CHR_hg38 == "chr11" & Start_hg38 <= 53400000 ~ "p", 
    CHR_hg38 == "chr11" & Start_hg38 > 53400000 ~ "q", 
    
    CHR_hg38 == "chr12" & Start_hg38 <= 35500000 ~ "p", 
    CHR_hg38 == "chr12" & Start_hg38 > 35500000 ~ "q", 
    
    CHR_hg38 == "chr13" & Start_hg38 <= 17700000 ~ "p", 
    CHR_hg38 == "chr13" & Start_hg38 > 17700000 ~ "q", 
    
    CHR_hg38 == "chr14" & Start_hg38 <= 17200000 ~ "p", 
    CHR_hg38 == "chr14" & Start_hg38 > 17200000 ~ "q", 
    
    CHR_hg38 == "chr15" & Start_hg38 <= 19000000 ~ "p", 
    CHR_hg38 == "chr15" & Start_hg38 > 19000000 ~ "q", 
    
    CHR_hg38 == "chr16" & Start_hg38 <= 36800000 ~ "p", 
    CHR_hg38 == "chr16" & Start_hg38 > 36800000 ~ "q", 
    
    CHR_hg38 == "chr17" & Start_hg38 <= 25100000 ~ "p", 
    CHR_hg38 == "chr17" & Start_hg38 > 25100000 ~ "q", 
    
    CHR_hg38 == "chr18" & Start_hg38 <= 18500000 ~ "p", 
    CHR_hg38 == "chr18" & Start_hg38 > 18500000 ~ "q", 
    
    CHR_hg38 == "chr19" & Start_hg38 <= 26200000 ~ "p",
    CHR_hg38 == "chr19" & Start_hg38 > 26200000 ~ "q",
    
    CHR_hg38 == "chr2" & Start_hg38 <= 93900000 ~ "p", 
    CHR_hg38 == "chr2" & Start_hg38 > 93900000 ~ "q", 
    
    CHR_hg38 == "chr20" & Start_hg38 <= 28100000 ~ "p", 
    CHR_hg38 == "chr20" & Start_hg38 > 28100000 ~ "q", 
    
    CHR_hg38 == "chr21" & Start_hg38 <= 12000000 ~ "p", 
    CHR_hg38 == "chr21" & Start_hg38 > 12000000 ~ "q", 
    
    CHR_hg38 == "chr22" & Start_hg38 <= 15000000 ~ "p", 
    CHR_hg38 == "chr22" & Start_hg38 > 15000000 ~ "q", 
    
    CHR_hg38 == "chr3" & Start_hg38 <= 90900000 ~ "p", 
    CHR_hg38 == "chr3" & Start_hg38 > 90900000 ~ "q", 
    
    CHR_hg38 == "chr4" & Start_hg38 <= 50000000 ~ "p", 
    CHR_hg38 == "chr4" & Start_hg38 > 50000000 ~ "q", 
    
    CHR_hg38 == "chr5" & Start_hg38 <= 48800000 ~ "p", 
    CHR_hg38 == "chr5" & Start_hg38 > 48800000 ~ "q", 
    
    CHR_hg38 == "chr6" & Start_hg38 <= 59800000 ~ "p", 
    CHR_hg38 == "chr6" & Start_hg38 > 59800000 ~ "q", 
    
    CHR_hg38 == "chr7" & Start_hg38 <= 60100000 ~ "p", 
    CHR_hg38 == "chr7" & Start_hg38 > 60100000 ~ "q", 
    
    CHR_hg38 == "chr8" & Start_hg38 <= 45200000 ~ "p", 
    CHR_hg38 == "chr8" & Start_hg38 > 45200000 ~ "q", 
    
    CHR_hg38 == "chr9" & Start_hg38 <= 43000000 ~ "p", 
    CHR_hg38 == "chr9" & Start_hg38 > 43000000 ~ "q", 
    
    T ~ "?"
  )) |> 
  dplyr::mutate(chr_arm_full = paste0(CHR_hg38, chr_arm))


# only take the arms with sufficient probes
plt <- plt |> 
  dplyr::filter(
    chr_arm_full %in% c(
      
      plt.oligo |> 
        dplyr::pull(chr_arm_full) |> 
        table() |> 
        as.data.frame() |> 
        dplyr::rename(arm = Var1, probes = Freq) |> 
        dplyr::filter(probes > 10) |> 
        dplyr::pull(arm) |> 
        as.character()
      
    ))




plt.oligo <- plt  |> 
  dplyr::select(CHR_hg38, chr_arm, chr_arm_full, DMP__g2_g3__pp_nc_PC1__t)


  
# add NA to those missing
plt.oligo <- rbind(plt.oligo, 
                     data.frame(
                       CHR_hg38 = c("chr13", "chr14", "chr15", "chr21", "chr22"),
                       chr_arm = rep("p", 5)
                     ) |> 
                       dplyr::mutate(chr_arm_full = paste0(CHR_hg38, chr_arm)) |> 
                       dplyr::mutate(DMP__g2_g3__pp_nc_PC1__t = NA)
                   ) |> 
  dplyr::mutate(tumortype = "Oligodendroglioma")



# reorder lexicographically
plt.oligo <- plt.oligo |>
  dplyr::mutate(col = ifelse(chr_arm_full %in% c("chr1p", "chr19q"), as.character(chr_arm_full),"other")) |> 
  dplyr::mutate(chr_arm_full = factor(chr_arm_full, levels=gtools::mixedsort(unique(as.character(chr_arm_full))))) |> 
  dplyr::mutate(CHR_hg38 = factor(CHR_hg38, levels=gtools::mixedsort(unique(as.character(CHR_hg38)))))




ggplot(plt.oligo) +
  geom_hline(yintercept=0, col="red", lwd=theme_nature_lwd) +
  geom_half_violin( # buggy library, needs AES within function itself otherwise it places the violins very odd
    aes(x = CHR_hg38,
        y = DMP__g2_g3__pp_nc_PC1__t,
        split = as.factor(chr_arm),
        fill = col),
    position = "identity",
    draw_quantiles = c(0.5),
    linewidth=theme_nature_lwd,
    width = 1.05
    #adjust = 2.95
  ) + 
  scale_fill_manual(values=c('chr1p'='lightblue',
                             #'chr20p'='lightpink',
                             #'chr20q'='lightpink',
                             'chr19q'='lightgreen',
                             'other'='lightgray')) + 
  theme_nature +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  labs(x=NULL, fill=NULL , y="per-probe t-score GLASS-OD") +
  coord_cartesian(ylim=c(-6.5, 4.75))









plt.astro <- plt  |> 
  dplyr::select(CHR_hg38, chr_arm, chr_arm_full, DMP__GLASS_NL__g2.3_g4__pp_nc_PC1__t)


# add NA to those missing
plt.astro <- rbind(plt.astro, 
                   data.frame(
                     CHR_hg38 = c("chr13", "chr14", "chr15", "chr21", "chr22"),
                     chr_arm = rep("p", 5)
                   ) |> 
                     dplyr::mutate(chr_arm_full = paste0(CHR_hg38, chr_arm)) |> 
                     dplyr::mutate(DMP__GLASS_NL__g2.3_g4__pp_nc_PC1__t = NA)
) |> 
  dplyr::mutate(tumortype = "Astrocytoma")


# reorder lexicographically
plt.astro <- plt.astro |>
  dplyr::mutate(col = ifelse(chr_arm_full %in% c(), as.character(chr_arm_full),"other")) |> 
  dplyr::mutate(chr_arm_full = factor(chr_arm_full, levels=gtools::mixedsort(unique(as.character(chr_arm_full))))) |> 
  dplyr::mutate(CHR_hg38 = factor(CHR_hg38, levels=gtools::mixedsort(unique(as.character(CHR_hg38)))))


ggplot(plt.astro) +
  geom_hline(yintercept=0, col="red", lwd=theme_nature_lwd) +
  geom_half_violin( # buggy library, needs AES within function itself otherwise it places the violins very odd
    aes(x = CHR_hg38,
        y = DMP__GLASS_NL__g2.3_g4__pp_nc_PC1__t,
        split = as.factor(chr_arm),
        fill = col),
    position = "identity",
    draw_quantiles = c( 0.5),
    linewidth=theme_nature_lwd,
    trim=T,
    nudge = 0,
    width = 1.05
  ) + 
  scale_fill_manual(values=c(#'chr20p'='lightpink',
                             #'chr20q'='lightpink',
                             'other'='lightgray')) + 
  theme_nature +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x=NULL, fill=NULL , y="per-probe t-score GLASS-NL") +
  coord_cartesian(ylim=c(-9, 6))




library(patchwork)
p1 / p2

ggsave("output/figures/vis_differential__GLASS_OD__GLASS_NL__grade__PC1__chromosomal_arm.pdf", width=8.5*0.975 * (6/8), height = 6)






## E: FFPE decay time ----
### 1. FFPE decay time ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t,
                col=DMP__FFPE_decay_time__pp_nc__t)) +
  #geom_vline(xintercept=0, col="red") +
  #geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001 , alpha=0.15) + 
  
  labs(x="Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent",
       col="",
       caption=paste0("todo")) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))))+
  coord_cartesian(ylim=c(-max(abs(plt$DMP__primary_recurrence__pp_nc__t)), max(abs(plt$DMP__primary_recurrence__pp_nc__t))))+
  
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
  
  theme(axis.line = element_line(linewidth = theme_nature_lwd))



ggsave(paste0("output/figures/vis_differential__FFPE_decay_time.png"), width=(8.5 * 0.975 * (1/5) * 1.12425), height=2.428, dpi=600)




ggplot(plt |> dplyr::filter(!grepl(" \\? ", probe_type_orientation)),
       aes(x=DMP__g2_g3__pp_nc__t,
           y=DMP__primary_recurrence__pp_nc__t,
           col=DMP__FFPE_decay_time__pp_nc__t)) +
  facet_wrap(~probe_type_orientation, ncol=2) + 
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


rm(plt)






### 2. Freezer decay time ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t,
                y=DMP__primary_recurrence__pp_nc__t,
              col=DMP__freezer_decay_time__pp_nc__t)) +
  
  #geom_vline(xintercept=0, col="red") +
  #geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001 , alpha=0.15) + 
  
  labs(x="Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent",
       col="",
       caption=paste0("todo")) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))))+
  coord_cartesian(ylim=c(-max(abs(plt$DMP__primary_recurrence__pp_nc__t)), max(abs(plt$DMP__primary_recurrence__pp_nc__t))))+
  
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
  
  theme(axis.line = element_line(linewidth = theme_nature_lwd )) # somehow w/ high-res png export



ggsave(paste0("output/figures/vis_differential__freezer_decay_time.png"), width=(8.5 * 0.975 * (1/5) * 1.12425), height=2.428, dpi=600)



rm(plt)





### 3. FFPE & freezer time multivar ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() 



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t,
                y=DMP__primary_recurrence__pp_nc__t,
              col=DMP__FFPE_and_freezer_decay_time__multivar__pp_nc__t__ffpe)) +
  
  #geom_vline(xintercept=0, col="red") +
  #geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001 , alpha=0.15) + 
  
  labs(x="Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent",
       col="",
       caption=paste0("todo")) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))))+
  coord_cartesian(ylim=c(-max(abs(plt$DMP__primary_recurrence__pp_nc__t)), max(abs(plt$DMP__primary_recurrence__pp_nc__t))))+
  
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
  
  theme(axis.line = element_line(linewidth = theme_nature_lwd)) # somehow w/ high-res png export


ggsave(paste0("output/figures/vis_differential__FFPE_and_freezer_decay_time__multivar__ffpe.png"), width=(8.5 * 0.975 * (1/5) * 1.12425), height=2.428, dpi=600)




ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t,
                y=DMP__primary_recurrence__pp_nc__t,
                col=DMP__FFPE_and_freezer_decay_time__multivar__pp_nc__t__freezer)) +
  
  #geom_vline(xintercept=0, col="red") +
  #geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001 , alpha=0.15) + 
  
  labs(x="Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent",
       col="",
       caption=paste0("todo")) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))))+
  coord_cartesian(ylim=c(-max(abs(plt$DMP__primary_recurrence__pp_nc__t)), max(abs(plt$DMP__primary_recurrence__pp_nc__t))))+
  
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
  
  theme(axis.line = element_line(linewidth = theme_nature_lwd)) # somehow w/ high-res png export


ggsave(paste0("output/figures/vis_differential__FFPE_and_freezer_decay_time__multivar__freezer.png"), width=(8.5 * 0.975 * (1/5) * 1.12425), height=2.428, dpi=600)







## F: Treatment effect(s) ----





## G: MNP-QC & detP ----


varss <- data.mvalues.probes |> 
  dplyr::select(contains("DMP__pct_detP_signi__pp_nc__t") | contains("_mnp_qc_")) |> 
  colnames()


for(var in varss) {
  txt <- gsub("^DMP__(.*?)__.+$","\\1", var)
  print(txt)
  
  plt <- data.mvalues.probes |> 
    data.table::copy() |> # odd hack needed, because setnames also affects the former "glass_od.metadata.array_samples" object...
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
      return(.)
    })() |> 
    dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
      return(.)
    })() |> 
    dplyr::filter(!grepl(" \\? ", probe_type_orientation)) |> 
    data.table::setnames(old = c(var), new = c("ewas_covar"))
  
  
  ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t,
                  y=DMP__primary_recurrence__pp_nc__t,
                  color=ewas_covar)
  ) +
    
    geom_vline(xintercept=0, col="red") +
    geom_hline(yintercept=0, col="red") +
    
    geom_point(pch=16, cex=0.001, alpha=0.15) + 
    
    geom_vline(xintercept=0, col="red", alpha=0.1) +
    geom_hline(yintercept=0, col="red", alpha=0.1) +
    
    labs(x = "Per probe t-score Grade 2 ~ Grade 3",
         y="Per probe t-score Primary ~ Recurrent",
         col="Association PC1",
         subtitle = paste0("qc: ", txt)) +
    
    #theme_nature +
    theme_nature +
    theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
    scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
    theme(plot.background = element_rect(fill="white", colour=NA)) # png export
    
    
    ggsave(paste0("output/figures/vis_differential__",txt,".png"), width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))
  
}


## H: ewastools ----


varss <- data.mvalues.probes |> 
  dplyr::select(contains("DMP__ewastools_")) |> 
  colnames()


for(var in varss) {
  txt <- gsub("^DMP__ewastools_(.*?)__.+$","\\1", var)
  print(txt)
  
  plt <- data.mvalues.probes |> 
    data.table::copy() |> # odd hack needed, because setnames also affects the former "glass_od.metadata.array_samples" object...
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
      return(.)
    })() |> 
    dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
      return(.)
    })() |> 
    dplyr::filter(!grepl(" \\? ", probe_type_orientation)) |> 
    data.table::setnames(old = c(var), new = c("ewas_covar"))
  
  
  ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t,
                  y=DMP__primary_recurrence__pp_nc__t,
                  color=ewas_covar)
         ) +

    geom_vline(xintercept=0, col="red") +
    geom_hline(yintercept=0, col="red") +
    
    geom_point(pch=16, cex=0.001, alpha=0.15) + 
    
    geom_vline(xintercept=0, col="red", alpha=0.1) +
    geom_hline(yintercept=0, col="red", alpha=0.1) +
    
    labs(x = "Per probe t-score Grade 2 ~ Grade 3",
         y="Per probe t-score Primary ~ Recurrent",
         col="Association PC1",
         subtitle = paste0("ewastools: ", txt)) +
    
    #theme_nature +
    theme_nature +
    theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
    scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
    theme(plot.background = element_rect(fill="white", colour=NA))  # png export
  
  
  ggsave(paste0("output/figures/vis_differential__ewastools_",txt,".png"), width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))
  
}

## I: HOX and PAX ----


### EMB-UP ----

emb_up <- c("cg12634591", "cg01128482", "cg24704177", "cg14666564", "cg25938806", "cg13316854", "cg01946191", "cg10159630", "cg23697546", "cg23697546",
            "cg00893471", "cg04316624", "cg24080247", "cg26551975", "cg02000318", "cg25287257", "cg01581084", "cg18702197", "cg25942990", "cg06145336", "cg24040595",
            "cg20540209", "cg05021643", "cg21778348", "cg00014998", "cg13019491", "cg08938793", "cg24308654", "cg09173768", "cg24541426", "cg01152019", "cg26521404",
            "cg27160395", "cg03255182", "cg03255182", "cg10753836", "cg16642791", "cg20501518", "cg25829490", "cg16304215", "cg12127282", "cg10682155", "cg22151644",
            "cg17104824", "cg01014615", "cg23981871", "cg04771946", "cg04822748", "cg14391419", "cg11836236", "cg26019295", "cg26019295", "cg10624122", "cg17466857",
            "cg06911354", "cg10126205", "cg16038003", "cg05065989", "cg04904385", "cg01901262", "cg05864326", "cg16168668", "cg15726154", "cg09495769", "cg18978493",
            "cg11357746", "cg15244786", "cg13352750", "cg16104915", "cg08431536", "cg17097782", "cg03048654", "cg16001495", "cg25062797", "cg27117509", "cg25953239",
            "cg04415798", "cg05060949")


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


plt <- plt |> 
  dplyr::mutate(label = hg38.manifest.gencode.v36.genesUniq) |> 
  dplyr::mutate(label = dplyr::case_when(
    (probe_id %in% emb_up) & grepl("(^|;)HOX", hg38.manifest.gencode.v36.genesUniq) ~ "EMB-UP (HOX)",
    (probe_id %in% emb_up) & (grepl("(^|;)HOX", hg38.manifest.gencode.v36.genesUniq) == F) ~ "EMB-UP (not HOX)",
    T ~ "-"
  ))


plt |> 
  dplyr::pull(label) |>
  table()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                y=DMP__primary_recurrence__pp_nc_PC1__t,
                color=DMP__g2_g3__pp_nc_PC1__t)) +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, label == "-"), pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data=subset(plt, label == "EMB-UP (HOX)"), pch=16, cex=0.001, col="black") + 
  geom_point(data=subset(plt, label == "EMB-UP (not HOX)"), pch=16, cex=0.001, col="purple") + 
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent",
       col="Association PC1",
       subtitle = paste0("EMP-UP: ", txt)) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export




### HOX ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



plt <- plt |> 
  dplyr::mutate(label = gsub("PHOX","", hg38.manifest.gencode.v36.genesUniq)) |> 
  dplyr::mutate(label = gsub("PHOX","", label)) |> 
  dplyr::mutate(label = gsub("SHOX","", label)) |>
  dplyr::mutate(label = gsub("RHOXF","", label)) |>
  dplyr::mutate(label = ifelse(grepl("HOX",label), "HOX+", "-"))



ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                  y=DMP__primary_recurrence__pp_nc_PC1__t,
                  color=DMP__g2_g3__pp_nc_PC1__t)) +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, label == "-"), pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data=subset(plt, label == "HOX+"), pch=16, cex=0.001, col="black") + 
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent",
       col="Association PC1",
       subtitle = paste0("HOX associated probes: ", txt)) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export






### PAX ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



plt <- plt |> 
  dplyr::mutate(label = hg38.manifest.gencode.v36.genesUniq) |> 
  
  dplyr::mutate(label = gsub("PAXIP1","", label)) |> 
  dplyr::mutate(label = gsub("PAXBP1","", label)) |>
  dplyr::mutate(label = gsub("PAXX","", label)) |>
  
  dplyr::mutate(label = ifelse(grepl("PAX9",label), "PAX+", "-"))




# plt |> 
#   dplyr::mutate(label = hg38.manifest.gencode.v36.genesUniq) |> 
#   
#   dplyr::mutate(label = gsub("PAXIP1","", label)) |> 
#   dplyr::mutate(label = gsub("PAXBP1","", label)) |>
#   dplyr::mutate(label = gsub("PAXX","", label)) |>
#   
#   dplyr::filter(grepl("PAX", label)) |> 
#   dplyr::pull(label) |> 
#   unique()



ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                y=DMP__primary_recurrence__pp_nc_PC1__t,
                color=DMP__g2_g3__pp_nc_PC1__t)) +
  
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, label == "-"), pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data=subset(plt, label == "PAX+"), pch=16, cex=0.001, col="black") + 
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent",
       col="Association PC1",
       subtitle = paste0("HOX associated probes: ", txt)) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export






### TET ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


plt |>
  dplyr::mutate(label = hg38.manifest.gencode.v36.genesUniq) |>

  dplyr::mutate(label = gsub("PAXIP1","", label)) |>

  dplyr::filter(grepl("IDH", label)) |>
  dplyr::pull(label) |>
  unique()



plt <- plt |> 
  dplyr::mutate(label = hg38.manifest.gencode.v36.genesUniq) |> 
  
  dplyr::mutate(label = dplyr::case_when(
    grepl("IDH1", label) ~ "IDH1",
    grepl("IDH2", label) ~ "IDH2",
    grepl("IDH3", label) ~ "IDH3",
    T ~ "-"
  ))



ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                y=DMP__primary_recurrence__pp_nc_PC1__t,
                color=DMP__g2_g3__pp_nc_PC1__t)) +
  
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, label == "-"), pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data=subset(plt, label == "IDH2"), pch=16, cex=0.001, col="black") + 
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent",
       col="Association PC1",
       subtitle = paste0("HOX associated probes: ", txt)) +
  
  #theme_nature +
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export




### FN1 ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



plt <- plt |> 
  dplyr::mutate(label = hg38.manifest.gencode.v36.genesUniq) |> 
  
  dplyr::mutate(label = gsub("[A-Z]FN1","", label)) |>
  
  dplyr::mutate(label = ifelse(grepl("NODAL",label), "FN1", "-"))



plt |>
  dplyr::mutate(label = hg38.manifest.gencode.v36.genesUniq) |>
  dplyr::mutate(label = gsub("[A-Z]FN1","", label)) |>
  dplyr::filter(grepl("FN1", label)) |>
  dplyr::pull(label) |>
  unique()



ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                y=DMP__primary_recurrence__pp_nc_PC1__t,
                color=DMP__g2_g3__pp_nc_PC1__t)) +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, label == "-"), pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data=subset(plt, label == "FN1"), pch=16, cex=0.001, col="black") + 
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent",
       col="Association PC1",
       subtitle = paste0("HOX associated probes: ", txt)) +
  
  #theme_nature +
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export




## J: median at grade 2 ----
### m-val ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t,
                y=DMP__primary_recurrence__pp_nc__t,
                col = median.mvalue.primary
)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3", y="Per probe t-score Primary ~ Recurrent",col = "Fit to PC1 vs. PC2") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)




### beta ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t,
                y=DMP__primary_recurrence__pp_nc__t,
                col = median.beta.recurrent
)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3", y="Per probe t-score Primary ~ Recurrent",col = "Median beta-value in primary tumours") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(0, 1), oob = scales::squish) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export


ggsave(paste0("output/figures/vis_differential__median_beta-value_primary.png"), width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))


## K: quality corrected ----
### dupli check ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::select(DMP__g2_g3__pp_nc_PC1__t,
                DMP__primary_recurrence__pp_nc_PC1__t,
                DMP__GLASS_NL__g2_g3.4__pp_nc_PC1__t,
                DMP__GLASS_NL__prim_rec__pp_nc_PC1__t
                )

corrplot::corrplot(cor(plt), order="hclust", tl.cex=0.75, tl.pos="l")




### GLASS-OD ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                y=DMP__primary_recurrence__pp_nc_PC1__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.0085, col="black") +  
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.5), lwd=theme_nature_lwd) +
  geom_smooth(method='lm', formula= y~x, se = F, lty=2, col="#6ba6e5", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)), col="#6ba6e5") +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
         y="Per probe t-score Primary ~ Recurrent"
  ) +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(0, 1), oob = scales::squish) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export



ggsave(paste0("output/figures/vis_differential__overall_density__quality_corrected.png"), width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))






### GLASS-OD x GLASS-NL x ATAC-seq ? ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(ATAC_astro_counts_per_bin_per_base) & !is.na(ATAC_oligo_counts_per_bin_per_base)) |> 
  dplyr::mutate(col = ATAC_astro_counts_per_bin_per_base - ATAC_oligo_counts_per_bin_per_base) |> 
  dplyr::mutate(rank = order(order(col))) |> 
  dplyr::mutate(rank_fraction = (rank - 1) / max(rank - 1)) 
#   dplyr::mutate(col = dplyr::case_when(
#     ATAC_astro_counts_per_bin_per_base > 2 ~ "strong",
#     T ~ "weak"
# )) 
  




ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                y=DMP__GLASS_NL__g2_g3.4__pp_nc_PC1__t,
                col=col)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  #geom_point(pch=16, cex=0.35) +   # , col="black"
  geom_point(pch=16, cex=0.3, alpha=0.1125 * 4) + # , col="black"
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.5), lwd=theme_nature_lwd) +
  geom_smooth(method='lm', formula= y~x, se = F, lty=2, col="#6ba6e5", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)), col="#6ba6e5") +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3 OD",
       y = "Per probe t-score Grade 2 ~ Grade 3 & 4 AC"
  ) +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-1.5, 1.5), oob = scales::squish) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export


ggsave(paste0("output/figures/vis_differential__GLASS_OD__x__GLASS_NL__overall_density_grade__quality_corrected.png"), width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))





ggplot(plt, aes(x=DMP__primary_recurrence__pp_nc_PC1__t,
                y=DMP__GLASS_NL__prim_rec__pp_nc_naive__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.0125*2, col="black") +  
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.5), lwd=theme_nature_lwd) +
  geom_smooth(method='lm', formula= y~x, se = F, lty=2, col="#6ba6e5", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)), col="#6ba6e5") +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3 OD",
       y = "Per probe t-score Grade 2 ~ Grade 3 & 4 AC"
  ) +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(0, 1), oob = scales::squish) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export



ggsave(paste0("output/figures/vis_differential__GLASS_OD__x__GLASS_NL__overall_density_prim_rec_corrected.png"), width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))


#### col chrs ----


plt <- plt |> 
  dplyr::mutate(col = CpG_chrm == "chr19" & is_19Q) |>  # chr2? chr13? chr18? [chr1 vs chr19]
  dplyr::mutate(col = is_1P | is_19Q) |> 
  dplyr::mutate(col = gc_sequence_context_2_new == "AA[CG]TC")


ggplot(plt, aes(x=DMP__primary_recurrence__pp_nc_PC1__t,
                y=DMP__GLASS_NL__prim_rec__pp_nc_naive__t,
                col=col)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, col==F), pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data=subset(plt, col==T), pch=16, cex=0.015, alpha=0.25) + 
  #geom_point(pch=16, cex=0.001, alpha=0.0125, col="black") +  
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.5), lwd=theme_nature_lwd) +
  geom_smooth(method='lm', formula= y~x, se = F, lty=2, col="#6ba6e5", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)), col="#6ba6e5") +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3 OD",
       y = "Per probe t-score Grade 2 ~ Grade 3 & 4 AC"
  ) +
  
  scale_color_manual(values=c('TRUE'='darkgreen','FALSE'='orange')) +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export


ggsave(paste0("output/figures/vis_differential__epiGenetic_clock__PC_corr_axes__PCHorvathS2018.png"), width = (8.5*0.975/2), height = 4.3)




### CGC ac / AcCGAP ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                y=DMP__primary_recurrence__pp_nc_PC1__t,
                col=DMP__AcCGAP__pp_nc__t
)) +
  
  geom_point(pch=16, cex=0.001 , alpha=0.15) + 
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent"
       #caption=paste0("Included samples per test: n=",n_samples_grade, " (grade), n=",n_samples_prim_rec," (resection type)")
  ) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))))+
  coord_cartesian(ylim=c(-max(abs(plt$DMP__primary_recurrence__pp_nc__t)), max(abs(plt$DMP__primary_recurrence__pp_nc__t)))) +
  
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


ggsave(paste0("output/figures/vis_differential__epiGenetic_clock__PC_corr_axes__AcCGAP.png"), width = (8.5*0.975/2), height = 4.3)



### PC2 ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                y=DMP__primary_recurrence__pp_nc_PC1__t,
                col=DMP__PCs__pp_nc__PC2_t
)) +
  
  geom_point(pch=16, cex=0.001 , alpha=0.15) + 
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent"
       #caption=paste0("Included samples per test: n=",n_samples_grade, " (grade), n=",n_samples_prim_rec," (resection type)")
  ) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))))+
  coord_cartesian(ylim=c(-max(abs(plt$DMP__primary_recurrence__pp_nc__t)), max(abs(plt$DMP__primary_recurrence__pp_nc__t)))) +
  
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


ggsave(paste0("output/figures/vis_differential__epiGenetic_clock__PC_corr_axes__PC2.png"), width = (8.5*0.975/2), height = 4.3)


### PC3 ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                y=DMP__primary_recurrence__pp_nc_PC1__t,
                col=DMP__PCs__pp_nc__PC3_t
)) +
  
  geom_point(pch=16, cex=0.001 , alpha=0.15) + 
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent"
       #caption=paste0("Included samples per test: n=",n_samples_grade, " (grade), n=",n_samples_prim_rec," (resection type)")
  ) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))))+
  coord_cartesian(ylim=c(-max(abs(plt$DMP__primary_recurrence__pp_nc__t)), max(abs(plt$DMP__primary_recurrence__pp_nc__t)))) +
  
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


ggsave(paste0("output/figures/vis_differential__epiGenetic_clock__PC_corr_axes__PC3.png"), width = (8.5*0.975/2), height = 4.3)




### PCHorvathS2018 qual cor ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()




ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                y=DMP__primary_recurrence__pp_nc_PC1__t,
                col=DMP__dnaMethyAge__PCHorvathS2018__up_nc__t
                )) +

  geom_point(pch=16, cex=0.001 , alpha=0.15) + 
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent"
       #caption=paste0("Included samples per test: n=",n_samples_grade, " (grade), n=",n_samples_prim_rec," (resection type)")
  ) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))))+
  coord_cartesian(ylim=c(-max(abs(plt$DMP__primary_recurrence__pp_nc__t)), max(abs(plt$DMP__primary_recurrence__pp_nc__t)))) +
  
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


ggsave(paste0("output/figures/vis_differential__epiGenetic_clock__PC_corr_axes__PCHorvathS2018.png"), width = (8.5*0.975/2), height = 4.3)





## C: G-CIMP ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


### a: IDHmut ----


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=-10*GCIMP_IDHmut_probe)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, GCIMP_IDHmut_probe == F),pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, GCIMP_IDHmut_probe == T),pch=16, cex=0.1, alpha=0.35) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


### b: IDHwt ----



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=-10*GCIMP_IDHwt_probe)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, GCIMP_IDHwt_probe == F),pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, GCIMP_IDHwt_probe == T),pch=16, cex=0.1, alpha=0.35) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)



### c: PanGlioma ----



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=-10*GCIMP_PanGlioma_probe)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, GCIMP_PanGlioma_probe == F),pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, GCIMP_PanGlioma_probe == T),pch=16, cex=0.1, alpha=0.35) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


### d: G-CIMP low 90 ----



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=-10*GCIMP_low_90)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, GCIMP_low_90 == F),pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, GCIMP_low_90 == T),pch=16, cex=0.1, alpha=0.35) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


## X: RepliTali predictors ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=-10*RepliTali_coef)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, RepliTali_coef == F | is.na(RepliTali_coef)),pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, RepliTali_coef == T),pch=16, cex=1, alpha=0.35) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(RepliTali_coef) & RepliTali_coef != F)



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=RepliTali_coef)) +
  #geom_vline(xintercept=0, col="red") +
  #geom_hline(yintercept=0, col="red") +
  
  #geom_point(data = subset(plt, RepliTali_coef == F | is.na(RepliTali_coef)),pch=16, cex=0.001, alpha=0.15) + 
  #geom_point(data = subset(plt, RepliTali_coef == T),pch=16, cex=1, alpha=0.35) + 
  geom_point(pch=16, cex=1, alpha=0.35) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox


ggplot(plt, aes(x=DMP__primary_recurrence__pp_nc__t, y=RepliTali_coef)) +
  #geom_vline(xintercept=0, col="red") +
  #geom_hline(yintercept=0, col="red") +
  
  #geom_point(data = subset(plt, RepliTali_coef == F | is.na(RepliTali_coef)),pch=16, cex=0.001, alpha=0.15) + 
  #geom_point(data = subset(plt, RepliTali_coef == T),pch=16, cex=1, alpha=0.35) + 
  geom_point(pch=16, cex=1, alpha=0.35) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox




### death clock ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::mutate(death_clock_probe = probe_id %in% c("cg03725309","cg25763716","cg13854219","cg25189904","cg15459165","cg19266329","cg24397007","cg23079012","cg27241845","cg06905155","cg16503724","cg19859270","cg02657160","cg14975410","cg05575921","cg14817490","cg21161138","cg12513616","cg20732076","cg06126421","cg15342087","cg01612140","cg25983901","cg12510708","cg00285394","cg01140244","cg23190089","cg07123182","cg26963277","cg18550212","cg10321156","cg25193885","cg07986378","cg23665802","cg04987734","cg19459791","cg00310412","cg23842572","cg19572487","cg01572694","cg08546016","cg18181703","cg03636183","cg24704287","cg11341610","cg14085840","cg26470501","cg05492306","cg25607249","cg01406381","cg07626482","cg03707168","cg08362785"))



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=-10*death_clock_probe, label=probe_id)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, death_clock_probe == F | is.na(death_clock_probe)),pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, death_clock_probe == T),pch=16, cex=1, alpha=0.35) + 
  
  #ggrepel::geom_text_repel(data = subset(plt, probe_id %in% 
  #                                         c("cg06905155","cg18181703","cg03636183","cg24704287","cg26470501")),nudge_y = 1) +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)






## D: FFPE | Tissue ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__FFPE_or_FF__pp_nc__t)) + # 
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)






### sandbox ----
#' uitzoeken wat 


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(abs(DMP__FFPE_decay_time__pp_nc__t) > 4.5) |> 
  dplyr::filter(grepl("?", probe_type_orientation))


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=c_content)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-8, 8), oob = scales::squish)



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=g_content)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(0, 1), oob = scales::squish)



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_content)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(0, 1), oob = scales::squish)



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=probeCpGcnt)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", trans = "log1p",
                                 breaks = 0:11,
                                 labels = c(0,"",2,"","","",6,"","","","",11)
  )  # , oob = scales::squish





ggplot(plt, aes(x=DMP__FFPE_decay_time__pp_nc__t, y=gc_content, col=probeCpGcnt)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", trans = "log1p",
                                 breaks = 0:11,
                                 labels = c(0,"",2,"","","",6,"","","","",11)
  )  # , oob = scales::squish





ggplot(plt, aes(x=DMP__FFPE_decay_time__pp_nc__t, y=gc_content, col=probeCpGcnt)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", trans = "log1p",
                                 breaks = 0:11,
                                 labels = c(0,"",2,"","","",6,"","","","",11)
  )  # , oob = scales::squish



ggplot(plt, aes(x=DMP__FFPE_decay_time__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=probeCpGcnt)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", trans = "log1p",
                                 breaks = 0:11,
                                 labels = c(0,"",2,"","","",6,"","","","",11)
  )  # , oob = scales::squish





## probe type




ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__FFPE_decay_time__pp_nc__t)) +
  facet_wrap(~type, scales="free",ncol=2) +

  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.5, alpha=0.75) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-7, 7), oob = scales::squish)




ggplot(plt, aes(x=type, y=probeCpGcnt)) +
  ggbeeswarm::geom_quasirandom(size=0.01)


ggplot(plt, aes(x=as.factor(probeCpGcnt), y=DMP__FFPE_decay_time__pp_nc__t)) +
  facet_wrap(~type, scales="free",ncol=2) +
  geom_violin(draw_quantiles = c(0.5), col="red")

#' wanneer een andere CpG wordt aangepast, wordt deze de-methylated?


## uberhaupt aantal C's



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(abs(DMP__epiTOC2_tnsc__up_nc__t) > 4.5 & !is.na(gc_content))


ggplot(plt, aes(x=probeCpGcnt, y=DMP__primary_recurrence__pp_nc__t, col=gc_content)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(0, 1), oob = scales::squish)




## uberhaupt aantal G's
## overall GC content
## overall intensity
## median b-value





## F: epiTOC2 / tsnc ----
#' normal epiTOC


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__epiTOC2_tnsc__up_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 

  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)




## G: epiTOC2 / [hypoSC!!] ----

plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()

ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__epiTOC2_hypoSC__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == F),pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == T),pch=16, cex=0.1, alpha=0.35) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)



## C content ----
### on probe ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(sequence_pre_48_c_content)) |> 
  dplyr::mutate(sequence_pre_48_c_content = scale(sequence_pre_48_c_content))


plt <- plt |> 
  dplyr::mutate(
    sequence_pre_48_c_content_f = cut(sequence_pre_48_c_content,
                      breaks=quantile(plt$sequence_pre_48_c_content,
                      probs = seq(0,1, by=0.2),
                      na.rm=TRUE, 
                      names=TRUE,
                      include.lowest=TRUE, 
                      right = TRUE)))



p1 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=sequence_pre_48_c_content)) +
  #facet_grid(cols = vars(Infinium_Design_Type))  +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.1) + 
  
  geom_point(data=subset(plt, sequence_pre_48_c_content_f == levels(plt$sequence_pre_48_c_content_f)[5]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_pre_48_c_content_f == levels(plt$sequence_pre_48_c_content_f)[4]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_pre_48_c_content_f == levels(plt$sequence_pre_48_c_content_f)[3]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_pre_48_c_content_f == levels(plt$sequence_pre_48_c_content_f)[2]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_pre_48_c_content_f == levels(plt$sequence_pre_48_c_content_f)[1]), pch=16, cex=0.001, alpha=0.1) + 

  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", 
                                 limits = c(-1.85, 1.85), oob = scales::squish ) + 
  
  labs(col = "C content on probe (scaled)") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))




### opposite to probe ----



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(sequence_post_48_c_content)) |> 
  dplyr::mutate(sequence_post_48_c_content = scale(sequence_post_48_c_content))



plt <- plt |> 
  dplyr::mutate(
    sequence_post_48_c_content_f = cut(sequence_post_48_c_content,
                      breaks=quantile(plt$sequence_post_48_c_content,
                                      probs = seq(0,1, by=0.2),
                                      na.rm=TRUE, 
                                      names=TRUE,
                                      include.lowest=TRUE, 
                                      right = TRUE)))



p2 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=sequence_post_48_c_content)) +
  #facet_grid(cols = vars(Infinium_Design_Type))  +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.1) + 
  
  geom_point(data=subset(plt, sequence_post_48_c_content_f == levels(plt$sequence_post_48_c_content_f)[5]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_post_48_c_content_f == levels(plt$sequence_post_48_c_content_f)[4]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_post_48_c_content_f == levels(plt$sequence_post_48_c_content_f)[3]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_post_48_c_content_f == levels(plt$sequence_post_48_c_content_f)[2]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_post_48_c_content_f == levels(plt$sequence_post_48_c_content_f)[1]), pch=16, cex=0.001, alpha=0.1) + 

  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", 
                                 limits = c(-1.85, 1.85), oob = scales::squish ) + 
  
  labs(col = "C content after probe (scaled)") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))




## G content ----


### on probe ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(sequence_pre_48_g_content)) |> 
  dplyr::mutate(sequence_pre_48_g_content = scale(sequence_pre_48_g_content))


plt <- plt |> 
  dplyr::mutate(
    sequence_pre_48_g_content_f = cut(sequence_pre_48_g_content,
                                      breaks=quantile(plt$sequence_pre_48_g_content,
                                                      probs = seq(0,1, by=0.2),
                                                      na.rm=TRUE, 
                                                      names=TRUE,
                                                      include.lowest=TRUE, 
                                                      right = TRUE)))



p3 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=sequence_pre_48_g_content)) +
  #facet_grid(cols = vars(Infinium_Design_Type))  +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.1) + 
  
  geom_point(data=subset(plt, sequence_pre_48_g_content_f == levels(plt$sequence_pre_48_g_content_f)[5]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_pre_48_g_content_f == levels(plt$sequence_pre_48_g_content_f)[4]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_pre_48_g_content_f == levels(plt$sequence_pre_48_g_content_f)[3]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_pre_48_g_content_f == levels(plt$sequence_pre_48_g_content_f)[2]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_pre_48_g_content_f == levels(plt$sequence_pre_48_g_content_f)[1]), pch=16, cex=0.001, alpha=0.1) + 
  
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", 
                                 limits = c(-1.85, 1.85), oob = scales::squish ) + 
  
  labs(col = "G content on probe (scaled)") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))




### opposite to probe ----



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(sequence_post_48_g_content)) |> 
  dplyr::mutate(sequence_post_48_g_content = scale(sequence_post_48_g_content))



plt <- plt |> 
  dplyr::mutate(
    sequence_post_48_g_content_f = cut(sequence_post_48_g_content,
                                       breaks=quantile(plt$sequence_post_48_g_content,
                                                       probs = seq(0,1, by=0.2),
                                                       na.rm=TRUE, 
                                                       names=TRUE,
                                                       include.lowest=TRUE, 
                                                       right = TRUE)))



p4 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=sequence_post_48_g_content)) +
  #facet_grid(cols = vars(Infinium_Design_Type))  +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.1) + 
  
  geom_point(data=subset(plt, sequence_post_48_g_content_f == levels(plt$sequence_post_48_g_content_f)[5]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_post_48_g_content_f == levels(plt$sequence_post_48_g_content_f)[4]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_post_48_g_content_f == levels(plt$sequence_post_48_g_content_f)[3]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_post_48_g_content_f == levels(plt$sequence_post_48_g_content_f)[2]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_post_48_g_content_f == levels(plt$sequence_post_48_g_content_f)[1]), pch=16, cex=0.001, alpha=0.1) + 
  
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", 
                                 limits = c(-1.85, 1.85), oob = scales::squish ) + 
  
  labs(col = "G content after probe (scaled)") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))







## GC content ----



### on probe ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(sequence_pre_48_gc_content)) |> 
  dplyr::mutate(sequence_pre_48_gc_content = scale(sequence_pre_48_gc_content))


plt <- plt |> 
  dplyr::mutate(
    sequence_pre_48_gc_content_f = cut(sequence_pre_48_gc_content,
                                       breaks=quantile(plt$sequence_pre_48_gc_content,
                                                       probs = seq(0,1, by=0.2),
                                                       na.rm=TRUE, 
                                                       names=TRUE,
                                                       include.lowest=TRUE, 
                                                       right = TRUE)))



p5 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=sequence_pre_48_gc_content)) +
  #facet_grid(cols = vars(Infinium_Design_Type))  +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.1) + 
  
  geom_point(data=subset(plt, sequence_pre_48_gc_content_f == levels(plt$sequence_pre_48_gc_content_f)[5]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_pre_48_gc_content_f == levels(plt$sequence_pre_48_gc_content_f)[4]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_pre_48_gc_content_f == levels(plt$sequence_pre_48_gc_content_f)[3]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_pre_48_gc_content_f == levels(plt$sequence_pre_48_gc_content_f)[2]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_pre_48_gc_content_f == levels(plt$sequence_pre_48_gc_content_f)[1]), pch=16, cex=0.001, alpha=0.1) + 
  
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", 
                                 limits = c(-1.85, 1.85), oob = scales::squish ) + 
  
  labs(col = "GC content on probe (scaled)") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))




### opposite to probe ----



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(sequence_post_48_gc_content)) |> 
  dplyr::mutate(sequence_post_48_gc_content = scale(sequence_post_48_gc_content))



plt <- plt |> 
  dplyr::mutate(
    sequence_post_48_gc_content_f = cut(sequence_post_48_gc_content,
                                        breaks=quantile(plt$sequence_post_48_gc_content,
                                                        probs = seq(0,1, by=0.2),
                                                        na.rm=TRUE, 
                                                        names=TRUE,
                                                        include.lowest=TRUE, 
                                                        right = TRUE)))



p6 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=sequence_post_48_gc_content)) +
  #facet_grid(cols = vars(Infinium_Design_Type))  +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.1) + 
  
  geom_point(data=subset(plt, sequence_post_48_gc_content_f == levels(plt$sequence_post_48_gc_content_f)[5]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_post_48_gc_content_f == levels(plt$sequence_post_48_gc_content_f)[4]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_post_48_gc_content_f == levels(plt$sequence_post_48_gc_content_f)[3]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_post_48_gc_content_f == levels(plt$sequence_post_48_gc_content_f)[2]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, sequence_post_48_gc_content_f == levels(plt$sequence_post_48_gc_content_f)[1]), pch=16, cex=0.001, alpha=0.1) + 
  
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", 
                                 limits = c(-1.85, 1.85), oob = scales::squish ) + 
  
  labs(col = "GC content after probe (scaled)") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))





(p1 + p2) / (p3 + p4) / (p5 + p6)



## H: probeCpGcnt ----


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=probeCpGcnt)) +
  #facet_grid(cols = vars(Infinium_Design_Type))  +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, probeCpGcnt == 0), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, probeCpGcnt == 1), pch=16, cex=0.0025, alpha=0.15) + 
  geom_point(data=subset(plt, probeCpGcnt == 2), pch=16, cex=0.005, alpha=0.20) + 
  geom_point(data=subset(plt, probeCpGcnt == 3), pch=16, cex=0.05, alpha=0.25) + 
  geom_point(data=subset(plt, probeCpGcnt == 4), pch=16, cex=0.01, alpha=0.30) + 
  geom_point(data=subset(plt, probeCpGcnt == 5), pch=16, cex=0.01, alpha=0.35) + 
  geom_point(data=subset(plt, probeCpGcnt == 6), pch=16, cex=0.01, alpha=0.40) + 
  geom_point(data=subset(plt, probeCpGcnt == 7), pch=16, cex=0.025, alpha=0.45) + 
  geom_point(data=subset(plt, probeCpGcnt == 8), pch=16, cex=0.05, alpha=0.50) + 
  geom_point(data=subset(plt, probeCpGcnt == 9), pch=16, cex=0.1, alpha=0.55) + 
  geom_point(data=subset(plt, probeCpGcnt == 10), pch=16, cex=0.1, alpha=0.60) + 
  geom_point(data=subset(plt, probeCpGcnt == 11), pch=16, cex=0.1, alpha=0.65) + 
  
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", trans = "log1p",
                                 breaks = 0:11,
                                 labels = c(0,"",2,"","","",6,"","","","",11)
  ) + # , oob = scales::squish
  
  labs(col = "CpG's per probe") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))


ggsave("output/figures/vis_differential__g23_prim-rec__probe_CpG_count.png", width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))




## I: closest_CG  ----

plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::mutate(closest_CG_25 = dplyr::case_when(
    closest_CG > 25 ~ 26,
    closest_CG < -25 ~ -26,
    T ~ closest_CG
  ))




ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=abs(closest_CG_25))) +
  
  facet_wrap(~abs(closest_CG_25) > 8, scales="free",ncol=5) +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  scale_color_gradientn(colours = rev(col3(200))[1:175], na.value = "grey50", 
                                 #trans = "log1p",
                                 breaks = c(1,5,10,15,20,26),
                                 labels = c(1,"",10,"",20,"  25+")
  ) + # , oob = scales::squish
  
  labs(col = "distance to next CpG (bp)") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))



ggsave("output/figures/vis_differential__g23_prim-rec__closest_CG.png", width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))


### b: in quadrant ----


ggplot(plt, aes(y=DMP__g2_g3__pp_nc__t, x=as.factor(closest_CG))) +
  #geom_point(pch=16, cex=0.001, alpha=0.15, col="darkgray") +
  ggplot2::geom_violin(draw_quantiles = c(0.5), linewidth=theme_nature_lwd, col = "darkgray", adjust = 1.95) +
  theme_nature +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 




## J: solo-wCGw ----
# >= 35 & == wCGw [https://www.nature.com/articles/s41467-022-34268-8]


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=is_solo_WCGW)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, is_solo_WCGW==F) , pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=subset(plt, is_solo_WCGW==T) , pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=0) , pch=16) + # hack plotting empty data frame, but with right cex and alpha for guide
  

  labs(col = "is solo-WCGW", x = "Per probe t-score Grade 2 ~ Grade 3", y="Per probe t-score Primary ~ Recurrent") +
  scale_color_manual(values=c('TRUE'= col3(3)[1], 'FALSE'='gray80')) +

  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))

ggsave("output/figures/vis_differential__g23_prim-rec__solo-WCGW.png", width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))



## K: discrete motifs ----
### a: motifs alone ----

plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |>
  dplyr::filter(!is.na(gc_sequence_context_2)) |> 
  
  dplyr::mutate(cg_prefix_suffix = grepl("^CG", gc_sequence_context_2) | grepl("CG$", gc_sequence_context_2)) |> 
  dplyr::mutate(top_affinity = gc_sequence_context_2 %in% c(
    "GA[CG]AG", "GA[CG]TC","CT[CG]TT","GA[CG]TG", "CT[CG]TT","CA[CG]TG","CA[CG]TT","GA[CG]TT","AA[CG]TT", 
    #            palin                                        palin                            palin
    "CT[CG]TC", "GA[CG]TC","CA[CG]AG","CA[CG]TC", "AA[CG]AG","CA[CG]TG","AA[CG]TG","AA[CG]TC","AA[CG]TT" )) |> 
  
  dplyr::mutate(col_facet1 = ifelse(cg_prefix_suffix, "CG prefix/suffix", "other")) |> 
  dplyr::mutate(col_facet2 = ifelse(top_affinity, "top affinity", "other")) 
  

plt <- rbind(
  plt |> 
    dplyr::mutate(col = col_facet1) |> 
    dplyr::mutate(facet = "CG prefix/suffix")
  ,
  plt |> 
    dplyr::mutate(col = col_facet2) |> 
    dplyr::mutate(facet = "top affinity")
) |> 
  dplyr::mutate(col = factor(col, levels=c(   "top affinity"  ,"CG prefix/suffix" ,"other"   )))


levels(plt$col)


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=col)) +
  facet_wrap(~facet, scales="free") +

  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, col == "other"), pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=subset(plt, col != "other"), pch=16, cex=0.001, alpha=0.15, show.legend = F) + 

  geom_point(data=head(plt, n=0),pch=16) + 
  
  labs(col = "sequence context", x = "Per probe t-score Grade 2 ~ Grade 3", y="Per probe t-score Primary ~ Recurrent") +
  scale_color_manual(values=c('top affinity'= col3(3)[1],
                              'other'='gray80',
                              'CG prefix/suffix'= col3(3)[3])) +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white", colour=NA))



### b: motifs x min closest_CG ----



## G: 1P ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()

ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=is_1P)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, !is_1P), pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=subset(plt,  is_1P), pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white", colour=NA))




## G: 19Q ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()

ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=is_19Q)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, !is_19Q), pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=subset(plt,  is_19Q), pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white", colour=NA))







# ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=probeCpGcnt)) +
#   facet_grid(cols = vars(facet)) +
#   
#   geom_vline(xintercept=0, col="red") +
#   geom_hline(yintercept=0, col="red") +
#   
#   geom_point(data=subset(plt, probeCpGcnt == 0), pch=16, cex=0.001, alpha=0.1) + 
#   geom_point(data=subset(plt, probeCpGcnt == 1), pch=16, cex=0.0025, alpha=0.15) + 
#   geom_point(data=subset(plt, probeCpGcnt == 2), pch=16, cex=0.005, alpha=0.20) + 
#   geom_point(data=subset(plt, probeCpGcnt == 3), pch=16, cex=0.05, alpha=0.25) + 
#   geom_point(data=subset(plt, probeCpGcnt == 4), pch=16, cex=0.01, alpha=0.30) + 
#   geom_point(data=subset(plt, probeCpGcnt == 5), pch=16, cex=0.01, alpha=0.35) + 
#   geom_point(data=subset(plt, probeCpGcnt == 6), pch=16, cex=0.01, alpha=0.40) + 
#   geom_point(data=subset(plt, probeCpGcnt == 7), pch=16, cex=0.025, alpha=0.45) + 
#   geom_point(data=subset(plt, probeCpGcnt == 8), pch=16, cex=0.05, alpha=0.50) + 
#   geom_point(data=subset(plt, probeCpGcnt == 9), pch=16, cex=0.1, alpha=0.55) + 
#   geom_point(data=subset(plt, probeCpGcnt == 10), pch=16, cex=0.1, alpha=0.60) + 
#   geom_point(data=subset(plt, probeCpGcnt == 11), pch=16, cex=0.1, alpha=0.65) + 
#   geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
#   
#   scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", trans = "log1p",
#                         breaks = 0:11,
#                         labels = c(0,"",2,"","","",6,"","","","",11)
#   ) + # , oob = scales::squish
#   
#   labs(col = "Chromosome of CpG", 
#        subtitle=format_subtitle("chromosomal differences")) +
#   
#   
#   theme_nature +
#   theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
#   theme(plot.background = element_rect(fill="white", colour=NA))
# 



## L: PC1 ----

## M: PC2 ----

## N: PC3 ----

## N: decay methy unmethy + unmethylated combined ----
### FFPE decay ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id) & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  assertr::verify(!is.na(DPI__FFPE_decay_time__pp_nc__t)) |> 
  assertr::verify(!is.na(DMPI__FFPE_decay_time__pp_nc__t)) |> 
  assertr::verify(!is.na(DUPI__FFPE_decay_time__pp_nc__t))







ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t,
                col=DPI__FFPE_decay_time__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  labs(col="Total probe intensity ~ FFPE decay time") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-8, 8), oob = scales::squish)







ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t,
                y=DMP__primary_recurrence__pp_nc__t,
              col=DMPI__FFPE_decay_time__pp_nc__t)) +
  geom_point(pch=16, cex=0.001 , alpha=0.15) + 
  
  labs(x="Per probe t-score Grade 2 ~ Grade 3",
       y="Per probe t-score Primary ~ Recurrent",
       col="",
       caption=paste0("todo")) +
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))))+
  coord_cartesian(ylim=c(-max(abs(plt$DMP__primary_recurrence__pp_nc__t)), max(abs(plt$DMP__primary_recurrence__pp_nc__t))))+
  
  #theme(plot.background = element_rect(fill="white", colour = NA)) + # png export
  
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
  
  theme(axis.line =     element_line(linewidth = theme_nature_lwd)) # somehow w/ png export






ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t,
                     col=DUPI__FFPE_decay_time__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  labs(col="Unmethylated probe intensity ~ FFPE decay time") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-8, 8), oob = scales::squish)




ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t,
                     col=DMP__FFPE_decay_time__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  labs(col="Probe methylation ratio ~ FFPE decay time") +
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-8, 8), oob = scales::squish)


p1 + p2 + p3 + p4





## O: intensities ----


### FFPE decay ----
plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(DPI__FFPE_decay_time__pp_nc__t))


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t,
                col=DPI__FFPE_decay_time__pp_nc__t)) +
  #facet_wrap(~probe_type, scales="free",ncol=2) +
  #facet_wrap(~map, scales="free",ncol=2) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-8, 8), oob = scales::squish)




### primary - recurrence ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(DPI__primary_recurrence__pp_nc__t))


p2 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t,
                     col=DPI__primary_recurrence__pp_nc__t)) +
  #facet_wrap(~probe_type, scales="free",ncol=2) +
  #facet_wrap(~map, scales="free",ncol=2) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-8, 8), oob = scales::squish)





### g2 - g3 ----



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(DPI__g2_g3__pp_nc__t))


p3 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t,
                     col=DPI__g2_g3__pp_nc__t)) +
  #facet_wrap(~probe_type, scales="free",ncol=2) +
  #facet_wrap(~map, scales="free",ncol=2) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-8, 8), oob = scales::squish)




p1 + p2 + p3 



## X: Infinium_Design_Type ----
#' not need to be shown, effect is by CpG's per probe,
#' design type is defined by number of probes



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=probe_type)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, probe_type == "II (ligation Allele-A)"), pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, probe_type== "I (red & green)"), pch=16, cex=0.001, alpha=0.15) +   
  
  theme_nature +
  theme(legend.key.size = unit(0.6, 'lines')) +
  guides(color = guide_legend(override.aes = list(size = 1, alpha=1) ) )


ggsave(paste0("output/figures/vis_differential__probe_type_II_I.png"), width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))





ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=probe_type)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, probe_type== "I (red & green)"), pch=16, cex=0.001, alpha=0.15) +   
  geom_point(data = subset(plt, probe_type == "II (ligation Allele-A)"), pch=16, cex=0.001, alpha=0.15) + 
  
  theme_nature + theme(legend.key.size = unit(0.6, 'lines'))  +
  guides(color = guide_legend(override.aes = list(size = 1, alpha=1) ) )

ggsave(paste0("output/figures/vis_differential__probe_type_I_II.png"), width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))





## motif: xx[CG]xx violin(s) ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |>
  dplyr::filter(!is.na(gc_sequence_context_2))


plt.per.context <- plt |> 
  dplyr::select(gc_sequence_context_2) |> 
  dplyr::filter(!duplicated(gc_sequence_context_2)) |> 
  dplyr::mutate(gc_sequence_context_2_s = gsub("[CG]","CG",gc_sequence_context_2,fixed=T)) |> 
  dplyr::mutate(gc_sequence_context_2_rc = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(gc_sequence_context_2_s)))) |> 
  dplyr::mutate(facet_name = dplyr::case_when(
    gc_sequence_context_2_s < gc_sequence_context_2_rc ~ paste0(gc_sequence_context_2_s ,"/",gc_sequence_context_2_rc),
    gc_sequence_context_2_s > gc_sequence_context_2_rc ~ paste0(gc_sequence_context_2_rc,"/",gc_sequence_context_2_s ),
    gc_sequence_context_2_s == gc_sequence_context_2_rc ~ paste0("*",gc_sequence_context_2_s,"/",gc_sequence_context_2_s)
  ))
#dplyr::mutate(gc_sequence_context_2_s = NULL, gc_sequence_context_2_rc = NULL)


plt <- plt |> 
  dplyr::left_join(plt.per.context, by=c('gc_sequence_context_2'='gc_sequence_context_2'), suffix=c('','')) |> 
  dplyr::group_by(facet_name) |> 
  dplyr::mutate(facet_rank = median(DMP__g2_g3__pp_nc__t)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(facet_name = tolower(facet_name))


ggplot(plt, aes(x=reorder(facet_name, -facet_rank), y=DMP__g2_g3__pp_nc__t)) +
  #facet_wrap(~reorder(facet_name, -facet_rank), scales="free_x", ncol=length(unique(plt$facet_name))) +
  #ggbeeswarm::geom_quasirandom(size=0.01, alpha=0.65) +
  #ggplot2::geom_violin(draw_quantiles = c(), linewidth=theme_nature_lwd, col = "white", fill="darkgray", adjust = 1.95) +
  geom_boxplot(width=0.75, outlier.shape=NA, outlier.color=NA, col="darkgray",
              ,coef=0.5, fill=NA,linewidth=theme_nature_lwd) +
  
  stat_summary(fun.y = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5, col="red", width=0.85, linewidth=theme_nature_lwd) +
  
  labs(x = NULL, y="Per probe t-score Grade 2 ~ Grade 3") +
  coord_cartesian(ylim = c(-6.75, 3.5)) + # soft clip
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, family = "mono"))


ggsave("output/figures/vis_differential__xxCGxx_violin.pdf", width = 11 * 0.97, height=3.75)



### motif: gc_sequence_context_1 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_1)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) 


stat <- plt |> 
  dplyr::select(gc_sequence_context_1, DMP__g2_g3__pp_nc__t)

ggplot(stat, aes(x=gc_sequence_context_1, y=DMP__g2_g3__pp_nc__t)) +
  ggbeeswarm::geom_quasirandom(size=theme_nature_size/2) +
  theme_bw()


#### RC ----






### motif: gc_sequence_context_1_ed ----


plt <- plt |> 
  dplyr::mutate(gc_sequence_context_1_edit = gsub("[CG]",":", gc_sequence_context_1, fixed=T)) |> 
  dplyr::mutate(gc_sequence_context_1_edit = dplyr::case_when(
    gc_sequence_context_1_edit == "A:A" ~ "A:A | T:T | A:T | T:A",
    gc_sequence_context_1_edit == "T:T" ~ "A:A | T:T | A:T | T:A",
    
    gc_sequence_context_1_edit == "A:T" ~ "A:A | T:T | A:T | T:A",
    gc_sequence_context_1_edit == "T:A" ~ "A:A | T:T | A:T | T:A",
    
    gc_sequence_context_1_edit == "G:G" ~ "G:G | C:C | C:G | G:C",
    gc_sequence_context_1_edit == "C:C" ~ "G:G | C:C | C:G | G:C",
    gc_sequence_context_1_edit == "C:G" ~ "G:G | C:C | C:G | G:C",
    gc_sequence_context_1_edit == "G:C" ~ "G:G | C:C | C:G | G:C",
    
    
    T ~ "G:* | C:*"
  ))
  #dplyr::filter(gc_sequence_context_1_edit %in% c("A:A | T:T", "A:T | T:A"))


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_1_edit)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  + 
  scale_colour_Publication() 



### motif: gc_sequence_context_2 ----


stat <- plt |> 
  dplyr::select(gc_sequence_context_2, DMP__g2_g3__pp_nc__t)

ggplot(stat, aes(x=gc_sequence_context_2, y=DMP__g2_g3__pp_nc__t)) +
  ggbeeswarm::geom_quasirandom(size=0.1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))




### m1 m2 m3 ----


pplt1 <- plt |> 
  dplyr::mutate(label = dplyr::case_when(
    gc_sequence_context_2 %in% c("AA[CG]TT") ~ "AA[CG]TT (*palin)",
    gc_sequence_context_2 %in% c("AA[CG]TG", "CA[CG]TT") ~ "AA[CG]TG/CA[CG]TT",
    gc_sequence_context_2 %in% c("AA[CG]TC", "GA[CG]TT") ~ "AA[CG]TC/GA[CG]TT",
    T ~ "other"
  )) 


p1 = ggplot(pplt1 , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=label)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(pplt1, label == "other"), pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(pplt1, label != "other"), pch=16, cex=0.15, alpha=0.15) + 
  
  geom_point(data=head(pplt1, n=1),pch=16, cex=1, alpha=0.8) + # for label
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y = "Per probe t-score Primary ~ Recurrent") +
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) +
  
  scale_color_manual(
    
    values = c(
      "AA[CG]TC/GA[CG]TT" = "darkblue",
               "AA[CG]TG/CA[CG]TT" = "darkblue",
               "AA[CG]TT (*palin)" = "darkblue",
               "other" = "#F8766D")
  )



pplt2 <- plt |> 
  dplyr::mutate(label = dplyr::case_when(
    grepl("CG$",gc_sequence_context_2) ~ "xx[CG]CG suffix",
    grepl("^CG",gc_sequence_context_2) ~ "CG[CG]xx prefix",
    T ~ "other"
  )) 

p2 = ggplot(pplt2 , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=label)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(pplt2, label == "other"), pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(pplt2, label != "other"), pch=16, cex=0.15, alpha=0.15) + 
  
  geom_point(data=head(pplt2, n=1),pch=16, cex=1, alpha=0.8) + # for label
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y = "Per probe t-score Primary ~ Recurrent") +
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) +
  
  scale_color_manual(
    
    values = c(
      "xx[CG]CG suffix" = "darkblue",
      "CG[CG]xx prefix" = "darkblue",
      "other" = "#F8766D")
  )


p1 + p2




### motif: gc_sequence_context_2_prefix ----


stat <- plt |> 
  dplyr::select(gc_sequence_context_2, DMP__g2_g3__pp_nc__t) |> 
  dplyr::mutate(gc_sequence_context_2_prefix = gsub("..$","", gc_sequence_context_2))

ggplot(stat, aes(x=gc_sequence_context_2_prefix, y=DMP__g2_g3__pp_nc__t)) +
  ggbeeswarm::geom_quasirandom(size=0.1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))


### motif: gc_sequence_context_2_suffix ----


stat <- plt |> 
  dplyr::select(gc_sequence_context_2, DMP__g2_g3__pp_nc__t) |> 
  dplyr::mutate(gc_sequence_context_2_suffix = gsub("^..","", gc_sequence_context_2))

ggplot(stat, aes(x=gc_sequence_context_2_suffix, y=DMP__g2_g3__pp_nc__t)) +
  ggbeeswarm::geom_quasirandom(size=0.1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))



### motif: gc_sequence_context_l1 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l1 %in% c('A','T'))) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_l2 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l2 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_l3 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l3 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 




### motif: gc_sequence_context_l4 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l4 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_l5 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l5 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 


### motif: gc_sequence_context_l6 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l6 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 


### motif: gc_sequence_context_l7 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l7 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 


### motif: gc_sequence_context_l8 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l8 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_r1 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r1 %in% c('A','T'))) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
  


### motif: gc_sequence_context_r2 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r2 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_r3 ----




ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r3 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox
  



### motif: gc_sequence_context_r4 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r4 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox




### motif: gc_sequence_context_r5 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r5 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox




### motif: gc_sequence_context_r6 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r6 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox



### motif: gc_sequence_context_r7 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r7 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox




### motif: gc_sequence_context_r8 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r8 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_nature + 
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox




## x: all epigenetic clocks ----


clocks <- plt |> 
  dplyr::select(contains("dnaMethyAge") | contains("epiTOC2") | contains("DMP__RepliTali__up_nc__t")) |> 
  dplyr::select(!contains("__adj.P.Val")) |> 
  colnames()


for(clock in clocks) {
  
  clock_name <- gsub("__up_nc__t","",gsub("DMP__","", clock))
  
  plt.p <- plt |> 
    data.table::copy() |> # odd hack needed, because setnames also affects the former "glass_od.metadata.array_samples" object...
    data.table::setnames(old = c(clock), new = c("array_current_clock"))
    
  
  ggplot(plt.p, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=array_current_clock)) +
    geom_vline(xintercept=0, col="red") +
    geom_hline(yintercept=0, col="red") +
    
    geom_point(pch=16, cex=0.001, alpha=0.15) + 
    
    labs(col = gsub("_", " ",clock_name)) +
    
    theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
    
    theme(plot.background = element_rect(fill="white", colour=NA)) + # png export
    
    scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)
  
  
  ggsave(paste0("output/figures/vis_differential__epiGenetic_clock__",clock_name,".png"), width = (8.5*0.975/2), height = 4.3)

}




## early ~ late repli ----

plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(wgEncodeUwRepliSeqBjWaveSignalRep1)) |> 
  dplyr::filter(wgEncodeUwRepliSeqBjWaveSignalRep2 < 79) |>  # outliers that mess up linear statistics
  dplyr::filter(wgEncodeUwRepliSeqBg02esWaveSignalRep1 < 90) |> 
  dplyr::filter(wgEncodeUwRepliSeqBjWaveSignalRep1 < 95) |> 
  dplyr::filter(wgEncodeUwRepliSeqBjWaveSignalRep1 > 10) |> 
  dplyr::filter(wgEncodeUwRepliSeqNhekWaveSignalRep1 < 80)


ggplot(plt, aes(y= wgEncodeUwRepliSeqBjWaveSignalRep1, x=glass_nl_prim_rec__deep_significant)) +
  ggbeeswarm::geom_quasirandom(size=0.1) +
  ggpubr::stat_compare_means(aes(group=glass_nl_prim_rec__deep_significant), label.x.npc=0.1, method = "wilcox.test", show.lengend  = FALSE,  size=theme_nature_size)


wilcox.test(
  plt |>
    dplyr::filter(glass_nl_prim_rec__deep_significant == T) |> 
    dplyr::pull(wgEncodeUwRepliSeqBjWaveSignalRep1)
  ,
  plt |>
    dplyr::filter(glass_nl_prim_rec__deep_significant == F) |> 
    dplyr::pull(wgEncodeUwRepliSeqBjWaveSignalRep1)
)



ggplot(plt, aes(x=wgEncodeUwRepliSeqBg02esWaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqBg02esWaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqBjWaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqBjWaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqBjWaveSignalRep2, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA))  # png export
  #geom_vline(xintercept=79, col="red", lwd=0.5)
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqBjWaveSignalRep2.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm06990WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm06990WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm12801WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm12801WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm12812WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm12812WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm12813WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm12813WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm12878WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm12878WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqHelas3WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqHelas3WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqHepg2WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqHepg2WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqHuvecWaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqHuvecWaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqImr90WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqImr90WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqK562WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqK562WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqMcf7WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqMcf7WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqNhekWaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqNhekWaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqSknshWaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_nature_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_nature + theme(plot.background = element_rect(fill="white", colour=NA)) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqSknshWaveSignalRep1.png",width=8.5/2, height=8.5/2)



# Figure S2A: corr covars per sample ----
# corr qc ----

#p1 = ggcorrplot2::ggcorrplot(h$corr) +
#  theme(axis.title.x=element_blank(),
#        axis.text.x=element_blank(),
#        axis.ticks.x=element_blank())


## qc covars


plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::mutate(time_tissue_in_ffpe =  ifelse(isolation_material == "ffpe", time_between_resection_and_array, 0)) |> 
  dplyr::select(contains("array_qc_"), 
                array_percentage.detP.signi, 
                array_PC1, 
                time_tissue_in_ffpe,
                #array_GLASS_NL_g2_g3_sig,
                array_PC2, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit,
                
                array_PC3
  ) |> 
  dplyr::filter(!is.na(time_tissue_in_ffpe)) |> 
  
  dplyr::mutate(array_GLASS_OD_g2_g3_sig2 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig2 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_g2_g3_sig3 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig3 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_g2_g3_sig4 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig4 = NULL) |> 
  
  dplyr::mutate(array_PC3 = -1 * array_PC3) |> 
  
  #dplyr::mutate(AcCGAP = -1 * array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit = NULL) |> 
  #dplyr::mutate(`-1 * array_dnaMethyAge__ZhangY2017` = -1 * array_dnaMethyAge__ZhangY2017 , array_dnaMethyAge__ZhangY2017 = NULL) |> 
  #dplyr::mutate(`-1 * array_epiTOC2_hypoSC` = -1 * array_epiTOC2_hypoSC, array_epiTOC2_hypoSC = NULL) |>  # seems at inversed scale
  #dplyr::mutate(`-1 * array_dnaMethyAge__LuA2019` = -1 * array_dnaMethyAge__LuA2019, array_dnaMethyAge__LuA2019 = NULL) |>  # seems at inversed scale
  #dplyr::mutate(`-1 * QC: SPECIFICITY I GT MM 6` = -1 * `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`, `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`=NULL) |> 
  #dplyr::select(-array_dnaMethyAge__YangZ2016) |>  # epiTOC1, near identical to epiTOC2
  #dplyr::select(-array_dnaMethyAge__PCHorvathS2013) |>  # very similar to its 2018 equivalent
  dplyr::mutate(`array_qc_BISULFITE_CONVERSION_I_Beta_I-3_Beta_larger_0.12_0.18` = NULL) |> # contains N/A's
  dplyr::mutate(`array_qc_BISULFITE_CONVERSION_I_Beta_I-6_Beta_larger_0.2_0.3` = NULL)# contains N/A's


colnames(plt) <- gsub("array_","",colnames(plt))
colnames(plt) <- gsub("_"," ",colnames(plt))

corrplot::corrplot(abs(cor(plt, method="spearman")), order="hclust", tl.cex=0.75, tl.pos="l")





# AD alzheimer disease ----


tmp <- readRDS("cache/analysis_differential__ad_co__stats.Rds") |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::rename_with( ~ paste0("DMP__ad_brain__", .x)) |> 
  tibble::rownames_to_column('probe_id')




plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id')) |> 
  dplyr::filter(!is.na(DMP__ad_brain__t))



source("scripts/load_gene_annotations.R")



polycomb_tfs <- plt |> 
  dplyr::rename(gene = GencodeCompV12_NAME) |> 
  dplyr::select(probe_id, gene) |>
  tibble::tibble() |> 
  tidyr::separate_longer_delim(gene, delim = ";") |> 
  dplyr::filter(gene != "") |> 
  dplyr::distinct() |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 323685)
    return(.)
  })() |> 
  dplyr::mutate(polycomb_homeobox = 
                  gene %in% c(genes_polycomb_suz12_homeobox, genes_polycomb_eed_homeobox, genes_polycomb_h3k27_homeobox, genes_polycomb_prc2_homeobox)
  ) |> 
  dplyr::filter(polycomb_homeobox) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 3105)
    return(.)
  })()



plt <- plt |> 
  dplyr::mutate(polycomb_homeobox_tf = probe_id %in% polycomb_tfs$probe_id )



ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                y=DMP__ad_brain__t,
                col=polycomb_homeobox_tf,
                alpha=polycomb_homeobox_tf,
                size=polycomb_homeobox_tf
)) +
  
  geom_point(pch=16) +
  
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)), col="black", size=theme_nature_size, label.x=5) +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y = "Per probe t-score Alzheimer - Control",
       caption=paste0("Included samples per test: ")
  ) +
  
  theme_nature + 
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))))+
  coord_cartesian(ylim=c(-max(abs(plt$DMP__ad_brain__t)), max(abs(plt$DMP__ad_brain__t)))) +
  
  scale_color_manual(values=c("darkgray","red")) +
  scale_size_manual(values=c(0.001, 0.2)) +
  scale_alpha_manual(values=c(`TRUE`=1.0, `FALSE`=0.6)) + 
  labs(subtitle=format_subtitle("Oligodendroglioma vs. Alzheimer Disease")) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export





ggsave(paste0("output/figures/vis_differential__grade_PC1__alzheimer_disease.png"), width=(8.5 * 0.975 * (2/5)), height=3.48, dpi=300)




# AD alzheimer disease x aging ----


tmp <- readRDS("cache/analysis_differential__ad_co__stats.Rds") |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::rename_with( ~ paste0("DMP__ad_brain__", .x)) |> 
  tibble::rownames_to_column('probe_id')




plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id')) |> 
  dplyr::filter(!is.na(DMP__ad_brain__t)) 
#dplyr::arrange(abs(DMP__dnaMethyAge__PCHorvathS2018__up_nc__t))




p1 <- ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                      y=DMP__ad_brain__t,
                      #col=DMP__AcCGAP__pp_nc__t
                      #col=DMP__dnaMethyAge__PCHorvathS2018__up_nc__t
                      col=DMP__PC1_PC2_PC3_multivariate__t_PC2
)) +
  
  geom_point(pch=16, cex=0.01 , alpha=0.45) + 
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y = "Per probe t-score Alzheimer - Control",
       caption=paste0("Included samples per test: ")
  ) +
  
  theme_nature + 
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))),
                  ylim=c(-max(abs(plt$DMP__ad_brain__t)), max(abs(plt$DMP__ad_brain__t)))) +
  
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
  
  labs(subtitle=format_subtitle("Oligodendroglioma vs. Alzheimer Disease")) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export


p2 <- ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                      y=DMP__ad_brain__t,
                      #col=DMP__AcCGAP__pp_nc__t
                      #col=DMP__dnaMethyAge__PCHorvathS2018__up_nc__t
                      col=DMP__PC1_PC2_PC3_multivariate__t_PC3
)) +
  
  geom_point(pch=16, cex=0.01 , alpha=0.45) + 
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y = "Per probe t-score Alzheimer - Control",
       caption=paste0("Included samples per test: ")
  ) +
  
  theme_nature + 
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))),
                  ylim=c(-max(abs(plt$DMP__ad_brain__t)), max(abs(plt$DMP__ad_brain__t)))) +
  
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
  
  labs(subtitle=format_subtitle("Oligodendroglioma vs. Alzheimer Disease")) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export


p1 + p2


ggsave(paste0("output/figures/vis_differential__grade_PC1__alzheimer_disease__aging.png"), width=(8.5 * 0.975 * (3.5/5)), height=3.48, dpi=300)







# AD alzheimer disease x aging ----


tmp <- readRDS("cache/analysis_differential__ad_co__stats.Rds") |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::rename_with( ~ paste0("DMP__ad_brain__", .x)) |> 
  tibble::rownames_to_column('probe_id')




plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id')) |> 
  dplyr::filter(!is.na(DMP__ad_brain__t)) 
#dplyr::arrange(abs(DMP__dnaMethyAge__PCHorvathS2018__up_nc__t))




ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t,
                y=DMP__ad_brain__t,
                col=DMP__AcCGAP__pp_nc__t
                #col=DMP__dnaMethyAge__PCHorvathS2018__up_nc__t
)) +
  
  geom_point(pch=16, cex=0.01 , alpha=0.45) + 
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y = "Per probe t-score Alzheimer - Control",
       caption=paste0("Included samples per test: ")
  ) +
  
  theme_nature + 
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))),
                  ylim=c(-max(abs(plt$DMP__ad_brain__t)), max(abs(plt$DMP__ad_brain__t)))) +
  
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish) +
  
  labs(subtitle=format_subtitle("Oligodendroglioma vs. Alzheimer Disease")) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export




ggsave(paste0("output/figures/vis_differential__grade_PC1__alzheimer_disease__AcCGAP.png"), width=(8.5 * 0.975 * (2/5)), height=3.48, dpi=300)





# polycomb x aging ----




tmp <- readRDS("cache/analysis_differential__ad_co__stats.Rds") |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::rename_with( ~ paste0("DMP__ad_brain__", .x)) |> 
  tibble::rownames_to_column('probe_id')




plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id')) |> 
  dplyr::filter(!is.na(DMP__ad_brain__t))



#source("scripts/load_gene_annotations.R")



polycomb_tfs <- plt |> 
  dplyr::rename(gene = GencodeCompV12_NAME) |> 
  dplyr::select(probe_id, gene) |>
  tibble::tibble() |> 
  tidyr::separate_longer_delim(gene, delim = ";") |> 
  dplyr::filter(gene != "") |> 
  dplyr::distinct() |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 323685)
    return(.)
  })() |> 
  dplyr::mutate(polycomb_homeobox = 
                  gene %in% c(genes_polycomb_suz12_homeobox, genes_polycomb_eed_homeobox, genes_polycomb_h3k27_homeobox, genes_polycomb_prc2_homeobox)
  ) |> 
  dplyr::filter(polycomb_homeobox) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 3105)
    return(.)
  })()



plt <- plt |> 
  dplyr::mutate(polycomb_homeobox_tf = probe_id %in% polycomb_tfs$probe_id )



ggplot(plt, aes(x=polycomb_homeobox_tf, y=DMP__dnaMethyAge__PCHorvathS2018__up_nc__t)) +
  geom_violin() +
  #geom_boxplot(width=0.75, outlier.shape=NA, outlier.color=NA, col="darkgray") +
  theme_nature


ggplot(plt, aes(x=polycomb_homeobox_tf, y=DMP__AcCGAP__pp_nc__t)) +
  geom_boxplot(width=0.75, outlier.shape=NA, outlier.color=NA, col="darkgray")


wilcox.test(
  plt |> dplyr::filter(polycomb_homeobox_tf)      |> dplyr::pull(DMP__dnaMethyAge__PCHorvathS2018__up_nc__t),
  plt |> dplyr::filter(polycomb_homeobox_tf == F) |> dplyr::pull(DMP__dnaMethyAge__PCHorvathS2018__up_nc__t)
)





