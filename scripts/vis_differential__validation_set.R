#!/usr/bin/env R



# match both back to grade in main dataset ----


readRDS("cache/analysis_differential__g2_g3__partial_paired_nc__validationset__stats.Rds") |>
  dplyr::mutate(signi = adj.P.Val < 0.01) |>
  dplyr::pull(signi) |>
  table()


readRDS("cache/analysis_differential__g2_g3__PC1__partial_paired_nc__validationset__stats.Rds") |>
  dplyr::mutate(signi = adj.P.Val < 0.01) |>
  dplyr::pull(signi) |>
  table()


readRDS("cache/analysis_differential__p_r__partial_paired_nc__validationset__stats.Rds") |>
  dplyr::mutate(signi = adj.P.Val < 0.01) |>
  dplyr::pull(signi) |>
  table()


readRDS("cache/analysis_differential__p_r__PC1__partial_paired_nc__validationset__stats.Rds") |>
  dplyr::mutate(signi = adj.P.Val < 0.01) |>
  dplyr::pull(signi) |>
  table()


# plt ----


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

  
  dplyr::select(probe_id,
                DMP__primary_recurrence__pp_nc__t,
                DMP__primary_recurrence__pp_nc__adj.P.Val,
                DMP__primary_recurrence__pp_nc__logFC,
                
                DMP__primary_recurrence__pp_nc_PC1__t,
                DMP__primary_recurrence__pp_nc_PC1__adj.P.Val,
                DMP__primary_recurrence__pp_nc_PC1__logFC,
                
                
                DMP__g2_g3__pp_nc__t,
                DMP__g2_g3__pp_nc__adj.P.Val,
                DMP__g2_g3__pp_nc__logFC,
                
                DMP__g2_g3__pp_nc_PC1__t,
                DMP__g2_g3__pp_nc_PC1__adj.P.Val,
                DMP__g2_g3__pp_nc_PC1__logFC
                ) |> 
  dplyr::left_join(
    readRDS("cache/analysis_differential__g2_g3__partial_paired_nc__validationset__stats.Rds") |> 
      dplyr::select(probe_id, t, adj.P.Val, logFC) |> 
      dplyr::rename(t_validation_g2_g3 = t,
                    logFC_validation_g2_g3 = logFC,
                    adj.P.Val_validation_g2_g3 = adj.P.Val),
    by=c('probe_id'='probe_id')
  )|> 
  dplyr::left_join(
    readRDS("cache/analysis_differential__p_r__partial_paired_nc__validationset__stats.Rds") |> 
      dplyr::select(probe_id, t, adj.P.Val, logFC) |> 
      dplyr::rename(t_validation_p_r = t,
                    logFC_validation_p_r = logFC,
                    adj.P.Val_validation_p_r = adj.P.Val,
      ),
    by=c('probe_id'='probe_id')
  ) |> 
  dplyr::left_join(
    readRDS("cache/analysis_differential__g2_g3__PC1__partial_paired_nc__validationset__stats.Rds") |> 
      dplyr::select(probe_id, t, adj.P.Val, logFC) |> 
      dplyr::rename(t_validation_g2_g3__PC1 = t,
                    logFC_validation_g2_g3__PC1 = logFC,
                    adj.P.Val_validation_g2_g3__PC1 = adj.P.Val,
      ),
    by=c('probe_id'='probe_id')
  )|> 
  dplyr::left_join(
    readRDS("cache/analysis_differential__p_r__PC1__partial_paired_nc__validationset__stats.Rds") |> 
      dplyr::select(probe_id, t, adj.P.Val, logFC) |> 
      dplyr::rename(t_validation_p_r__PC1 = t,
                    adj.P.Val_validation_p_r__PC1 = adj.P.Val,
                    logFC_validation_p_r__PC1 = logFC),
    by=c('probe_id'='probe_id')
  ) |> 
  dplyr::filter(!is.na(DMP__primary_recurrence__pp_nc__t)) |> 
  dplyr::mutate(col=factor("black", levels=c("black","white")))




# coorr stuff ----

c = cor(plt |> tibble::column_to_rownames('probe_id'))

c = cor(plt |> dplyr::select(probe_id,
                             DMP__g2_g3__pp_nc__t, DMP__g2_g3__pp_nc_PC1__t,
                             t_validation_p_r,t_validation_p_r__PC1) |> tibble::column_to_rownames('probe_id'))


corrplot::corrplot(c, order="hclust", shade.lwd=0.5, tl.cex=0.4, cl.cex=0.4, tl.pos="l" )



# scatter GLASS-OD x validatie naive ----



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=t_validation_g2_g3)) +
  geom_vline(xintercept=0, col="red", lwd=theme_nature_lwd) +
  geom_hline(yintercept=0, col="red", lwd=theme_nature_lwd) +
  
  geom_point(pch=16, cex=0.001, alpha=0.01, col="black") +  
  
  geom_vline(xintercept=0, col="red", alpha=0.1, lwd=theme_nature_lwd) +
  geom_hline(yintercept=0, col="red", alpha=0.1, lwd=theme_nature_lwd) +
  
  #geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.1), lwd=theme_nature_lwd * 4) +
  #geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.1), lwd=theme_nature_lwd * 2) +
  #geom_smooth(method='lm', formula= y~x, se = F, lty=2, col="#6ba6e5", lwd=theme_nature_lwd) +
  
  ggpubr::stat_cor(
    label.x = -8,
    label.y =  5,
    method = "pearson", aes(label = after_stat(r.label)), col="black",
    size = theme_nature_size,
    family = theme_nature_font_family) +
  
  labs(x = "Per probe t-score GLASS-OD Grade 2 ~ Grade 3",
       y = "Per probe t-score GLASS-NL Grade 2 & 3 ~ Grade 4",
       subtitle=format_subtitle("Overlap outcome oligodendroglioma vs. astrocytoma"),
       caption="Overlap Grade associated differences oligodendroglioma & astrocytoma (patient and PC1 corrected)"
  ) +
  
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc__t)), max(abs(plt$DMP__g2_g3__pp_nc__t))))+
  coord_cartesian(ylim=c(-max(abs(plt$t_validation_g2_g3)), max(abs(plt$t_validation_g2_g3)))) +
  
  theme_nature +
  theme(legend.key.size = unit(theme_nature_lwd * 1.5, 'lines')) + # resize colbox
  #scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(0, 1), oob = scales::squish) +
  theme(plot.background = element_rect(fill="white", colour=NA))  # png export





ggsave(paste0("output/figures/vis_differential__GLASS-OD_x_validation__naive.png"),
       width=2.115, height=2.385, dpi=600)





# scatter GLASS-OD x validatie PC1 ----

# contingency table
tab <- plt |> 
  dplyr::filter(!is.na(DMP__g2_g3__pp_nc_PC1__adj.P.Val)) |> 
  dplyr::filter(!is.na(DMP__g2_g3__pp_nc_PC1__logFC)) |> 
  dplyr::filter(!is.na(adj.P.Val_validation_g2_g3__PC1)) |> 
  dplyr::filter(!is.na(logFC_validation_g2_g3__PC1)) |> 
  
  dplyr::mutate(significant_glass_od = 
                DMP__g2_g3__pp_nc_PC1__adj.P.Val < 0.01 & abs(DMP__g2_g3__pp_nc_PC1__logFC) > 0.5
                ) |> 
  dplyr::mutate(significant_validation = adj.P.Val_validation_g2_g3__PC1 < 0.01 & abs(logFC_validation_g2_g3__PC1) > 0.5) |> 
  dplyr::select(significant_glass_od, significant_validation)


ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t, y=t_validation_g2_g3__PC1, col=col)) +
  geom_vline(xintercept=0, col="red", lwd=theme_nature_lwd) +
  geom_hline(yintercept=0, col="red", lwd=theme_nature_lwd) +
  
  #geom_point(pch=16, cex=0.001, alpha=0.01, col="black") +  
  geom_point(pch=16, cex=0.001, alpha=0.0085, col="black") +
  
  geom_vline(xintercept=0, col="red", alpha=0.1, lwd=theme_nature_lwd) +
  geom_hline(yintercept=0, col="red", alpha=0.1, lwd=theme_nature_lwd) +
  
  #geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.1), lwd=theme_nature_lwd * 4) +
  #geom_smooth(method='lm', formula= y~x, se = F, lty=1, col=alpha("white",0.1), lwd=theme_nature_lwd * 2) +
  #geom_smooth(method='lm', formula= y~x, se = F, lty=2, col="#6ba6e5", lwd=theme_nature_lwd) +
  
  ggpubr::stat_cor(
    #label.x = -8,
    #label.y =  4,
    method = "pearson", aes(label = after_stat(r.label)), col="#444444",
    size = theme_nature_size,
    family = theme_nature_font_family) +
  
  labs(x = "Per probe t-score GLASS-OD Grade 2 ~ Grade 3 (quality corrected)",
       y = "Per probe t-score validation set Grade 2 ~ Grade 3 (quality corrected)"
       #subtitle=format_subtitle("Overlap outcome GLASS-OD vs. validation set"),
       #caption="caption"
  ) +
  
  coord_cartesian(xlim=c(-max(abs(plt$DMP__g2_g3__pp_nc_PC1__t)),
                          max(abs(plt$DMP__g2_g3__pp_nc_PC1__t))),
                  ylim=c(-max(abs(plt$t_validation_g2_g3__PC1)),
                          max(abs(plt$t_validation_g2_g3__PC1)))) +
  
  theme_nature +
  theme(legend.key.size = unit(theme_nature_lwd * 1.5, 'lines')) + # resize colbox
  
  #scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(0, 1), oob = scales::squish) +
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +
  theme(aspect.ratio=1)




ggsave(paste0("output/figures/vis_differential__GLASS-OD_x_validation__PC1.png"),
       width=2.1, height=3.9 , dpi=1200)




