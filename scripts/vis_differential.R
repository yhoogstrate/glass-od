#!/usr/bin/env R


# load ----


library(ggplot2)
library(patchwork)

source('scripts/load_palette.R')
source('scripts/load_themes.R')

if(!exixts('data.mvalues.probes')) {
  source('scripts/load_mvalues_hq_samples.R')
}




# plot ----

## A: overall density ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.0075, col="black") + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)





## B: AcCGAP in OD ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
    return(.)
  })()

ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__AcCGAP__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)




## C: GLASS-NL probes ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=-10*glass_nl_prim_rec__deep_significant)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == F),pch=19, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == T),pch=19, cex=0.1, alpha=0.35) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)





## D: FFPE | Tissue ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__FFPE_or_FF__pp_nc__t)) + # 
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)





## E: FFPE decay time ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__FFPE_decay_time__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)





## F: epiTOC2 / tsnc ----
#' normal epiTOC


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__epiTOC2_tnsc__up_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15) + 

  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)





## G: epiTOC2 / [hypoSC!!] ----

ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__epiTOC2_hypoSC__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == F),pch=19, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == T),pch=19, cex=0.1, alpha=0.35) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)



## H: probeCpGcnt ----

ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=probeCpGcnt)) +
  #facet_grid(cols = vars(Infinium_Design_Type))  +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, probeCpGcnt == 0), pch=19, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, probeCpGcnt == 1), pch=19, cex=0.0025, alpha=0.15) + 
  geom_point(data=subset(plt, probeCpGcnt == 2), pch=19, cex=0.005, alpha=0.20) + 
  geom_point(data=subset(plt, probeCpGcnt == 3), pch=19, cex=0.05, alpha=0.25) + 
  geom_point(data=subset(plt, probeCpGcnt == 4), pch=19, cex=0.01, alpha=0.30) + 
  geom_point(data=subset(plt, probeCpGcnt == 5), pch=19, cex=0.01, alpha=0.35) + 
  geom_point(data=subset(plt, probeCpGcnt == 6), pch=19, cex=0.01, alpha=0.40) + 
  geom_point(data=subset(plt, probeCpGcnt == 7), pch=19, cex=0.025, alpha=0.45) + 
  geom_point(data=subset(plt, probeCpGcnt == 8), pch=19, cex=0.05, alpha=0.50) + 
  geom_point(data=subset(plt, probeCpGcnt == 9), pch=19, cex=0.1, alpha=0.55) + 
  geom_point(data=subset(plt, probeCpGcnt == 10), pch=19, cex=0.1, alpha=0.60) + 
  geom_point(data=subset(plt, probeCpGcnt == 11), pch=19, cex=0.1, alpha=0.65) + 
  
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", trans = "log1p",
                                 breaks = 0:11,
                                 labels = c(0,"",2,"","","",6,"","","","",11)
  ) + # , oob = scales::squish
  
  labs(col = "CpG's per probe") +
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white")) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))

ggsave("output/figures/vis_differential__g23_prim-rec__probe_CpG_count.png", width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))


## I: closest_CG  ----

plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
    return(.)
  })() |> 
  dplyr::mutate(closest_CG_25 = dplyr::case_when(
    closest_CG > 25 ~ 26,
    closest_CG < -25 ~ -26,
    T ~ closest_CG
  ))

### a: closest_CG distance x effect  ----


ggplot(plt, aes(y=DMP__g2_g3__pp_nc__t, x=as.factor(closest_CG))) +
  #geom_point(pch=19, cex=0.001, alpha=0.15, col="darkgray") +
  ggplot2::geom_violin(draw_quantiles = c(0.5), linewidth=theme_cellpress_lwd, col = "darkgray", adjust = 1.95) +
  theme_cellpress +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


### b: in quadrant ----



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=abs(closest_CG_25))) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15) + 
  
  ggplot2::scale_color_gradientn(colours = rev(col3(200))[1:175], na.value = "grey50", 
                                 #trans = "log1p",
                                 breaks = c(1,5,10,15,20,26),
                                 labels = c(1,"",10,"",20,"  25+")
  ) + # , oob = scales::squish
  
  labs(col = "distance to next CpG (bp)") +
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white")) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))

ggsave("output/figures/vis_differential__g23_prim-rec__closest_CG.png", width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))



## J: solo-wCGw ----
# >= 35 & == wCGw [https://www.nature.com/articles/s41467-022-34268-8]


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=is_solo_WCGW)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, is_solo_WCGW==F) , pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=subset(plt, is_solo_WCGW==T) , pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=0) , pch=19) + # hack plotting empty data frame, but with right cex and alpha for guide
  

  labs(col = "is solo-WCGW", x = "Per probe t-score Grade2 ~ Grade 3", y="Per probe t-score Primary ~ Recurrence") +
  scale_color_manual(values=c('TRUE'= col3(3)[1], 'FALSE'='gray80')) +

  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white")) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))

ggsave("output/figures/vis_differential__g23_prim-rec__solo-WCGW.png", width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))



## K: discrete motifs ----
### a: motifs alone ----

plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
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



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=col)) +
  facet_wrap(~facet, scales="free") +

  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, col == "other"), pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=subset(plt, col != "other"), pch=19, cex=0.001, alpha=0.15, show.legend = F) + 

  geom_point(data=head(plt, n=0),pch=19) + 
  
  labs(col = "sequence context", x = "Per probe t-score Grade2 ~ Grade 3", y="Per probe t-score Primary ~ Recurrence") +
  scale_color_manual(values=c('top affinity'= col3(3)[1],
                              'other'='gray80',
                              'CG prefix/suffix'= col3(3)[3])) +
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white"))



### b: motifs x min closest_CG ----



## G: 1P ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
    return(.)
  })()

ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=is_1P)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, !is_1P), pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=subset(plt,  is_1P), pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white"))




## G: 19Q ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
    return(.)
  })()

ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=is_19Q)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, !is_19Q), pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=subset(plt,  is_19Q), pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white"))


## G: 1P | 19Q ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
    return(.)
  })() |> 
  dplyr::mutate(col = dplyr::case_when(
    CpG_chrm == "chr1" ~ CpG_chrm,
    CpG_chrm == "chr19" ~ CpG_chrm,
    T ~ "other"
  )) |> 
  dplyr::mutate(facet = ifelse(is_1P | is_19Q, "on 1p | 19q","other")) |> 
  dplyr::mutate( CpG_chrm = dplyr::case_when(
    CpG_chrm == "chr1" & is_1P ~ "chr1p",
    CpG_chrm == "chr1" & !is_1P ~ "chr1q",
  
    CpG_chrm == "chr19" & is_19Q ~ "chr19q",
    CpG_chrm == "chr19" & !is_19Q ~ "chr19p",
    T ~ CpG_chrm
  )) |> 
  dplyr::mutate(chr = factor(CpG_chrm, levels=gtools::mixedsort(unique(as.character(CpG_chrm))) ))



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=probeCpGcnt)) +
  facet_grid(cols = vars(facet)) +

  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, probeCpGcnt == 0), pch=19, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, probeCpGcnt == 1), pch=19, cex=0.0025, alpha=0.15) + 
  geom_point(data=subset(plt, probeCpGcnt == 2), pch=19, cex=0.005, alpha=0.20) + 
  geom_point(data=subset(plt, probeCpGcnt == 3), pch=19, cex=0.05, alpha=0.25) + 
  geom_point(data=subset(plt, probeCpGcnt == 4), pch=19, cex=0.01, alpha=0.30) + 
  geom_point(data=subset(plt, probeCpGcnt == 5), pch=19, cex=0.01, alpha=0.35) + 
  geom_point(data=subset(plt, probeCpGcnt == 6), pch=19, cex=0.01, alpha=0.40) + 
  geom_point(data=subset(plt, probeCpGcnt == 7), pch=19, cex=0.025, alpha=0.45) + 
  geom_point(data=subset(plt, probeCpGcnt == 8), pch=19, cex=0.05, alpha=0.50) + 
  geom_point(data=subset(plt, probeCpGcnt == 9), pch=19, cex=0.1, alpha=0.55) + 
  geom_point(data=subset(plt, probeCpGcnt == 10), pch=19, cex=0.1, alpha=0.60) + 
  geom_point(data=subset(plt, probeCpGcnt == 11), pch=19, cex=0.1, alpha=0.65) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", trans = "log1p",
                                 breaks = 0:11,
                                 labels = c(0,"",2,"","","",6,"","","","",11)
  ) + # , oob = scales::squish
  
  labs(col = "CpG's per probe") +
  
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white"))


ggplot(plt, aes(x = chr, y=DMP__g2_g3__pp_nc__t, fill=col)) +
  ggplot2::geom_violin(draw_quantiles = c(0.5), linewidth=theme_cellpress_lwd, col = "white", adjust = 1.95) +
  scale_fill_manual(values=c('chr1'='red','chr19'='blue','other'='darkgray')) +
  theme_cellpress +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 



## X: Infinium_Design_Type ----
#' not need to be shown, effect is by CpG's per probe,
#' design type is defined by number of probes



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=Infinium_Design_Type)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, Infinium_Design_Type == "II"), pch=19, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt,Infinium_Design_Type == "I"), pch=19, cex=0.001, alpha=0.15) +   
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) 




## motif: xx[CG]xx violin(s) ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
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
  ggplot2::geom_violin(draw_quantiles = c(0.5), linewidth=theme_cellpress_lwd, col = "white", fill = "darkgray", adjust = 1.95) +
  labs(x = NULL) +
  coord_cartesian(ylim = c(-6.75, 3.5)) + # soft clip
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, family = "mono"))


ggsave("output/figures/vis_differential__xxCGxx_violin.pdf", width = 11 * 0.97, height=3.75)





## motif: xx[CG]xx * TET1-2 NatCom ----



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
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



plt <- plt |> 
  dplyr::left_join(plt.per.context, by=c('gc_sequence_context_2'='gc_sequence_context_2'), suffix=c('','')) |> 
  dplyr::group_by(facet_name) |> 
  dplyr::mutate(context_median__DMP__g2_g3__pp_nc__t = median(DMP__g2_g3__pp_nc__t)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(facet_name = tolower(facet_name)) |> 
  dplyr::select(facet_name, gc_sequence_context_2, context_median__DMP__g2_g3__pp_nc__t) |> 
  dplyr::distinct() |> 
  dplyr::mutate(adjacent_CG = grepl("^cg",facet_name) | grepl("cg$",facet_name)) 





tmp <- readxl::read_xlsx('data/42003_2022_3033_MOESM4_ESM.xlsx', skip=3) |> 
  dplyr::rename(sequence = 1) |> 
  dplyr::rename(TET1_5mC_oxidation_kinetics  = 2) |> 
  dplyr::rename(TET1_5hmC_oxidation_kinetics = 3) |> 
  dplyr::rename(TET2_5mC_oxidation_kinetics  = 4) |> 
  dplyr::rename(TET2_5hmC_oxidation_kinetics = 5) |> 
  dplyr::mutate(sequence = gsub(" ", "" ,sequence))



plt <- plt |> 
  dplyr::mutate(sequence = gsub("[","",gsub("]","",gc_sequence_context_2),fixed=T),fixed=T) |> 
  dplyr::left_join(tmp, by=c('sequence'='sequence'), suffix=c('','')) 


plt.per.facet <- plt |> 
  dplyr::mutate(gc_sequence_context_2 = NULL) |> 
  dplyr::mutate(sequence = NULL) |> 
  dplyr::group_by(facet_name) |> 
  dplyr::mutate(TET1_5mC_oxidation_kinetics_per_facet = mean(TET1_5mC_oxidation_kinetics)) |> 
  dplyr::mutate(TET1_5hmC_oxidation_kinetics_per_facet = mean(TET1_5hmC_oxidation_kinetics)) |> 
  dplyr::mutate(TET2_5mC_oxidation_kinetics_per_facet = mean(TET2_5mC_oxidation_kinetics)) |> 
  dplyr::mutate(TET2_5hmC_oxidation_kinetics_per_facet = mean(TET2_5hmC_oxidation_kinetics)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(TET1_5mC_oxidation_kinetics = NULL) |> 
  dplyr::mutate(TET1_5hmC_oxidation_kinetics = NULL) |> 
  dplyr::mutate(TET2_5mC_oxidation_kinetics = NULL) |> 
  dplyr::mutate(TET2_5hmC_oxidation_kinetics = NULL) |> 
  dplyr::distinct() |> 
  tidyr::pivot_longer(cols=contains("kinetics"), names_to="oxidation_kinetics")

# TET1_5mC_oxidation_kinetics_per_facet

ggplot(plt.per.facet, aes(x=context_median__DMP__g2_g3__pp_nc__t, y=value, label=facet_name, col=adjacent_CG)) +
  facet_wrap(~oxidation_kinetics, scales="free",ncol=5) +
  geom_point() +
  ggrepel::geom_text_repel(data = subset(plt.per.facet, grepl("aacgtt|aacgtc|aacgtg", facet_name)), col="black", size=theme_cellpress_size, alpha=0.5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho", size=theme_cellpress_size, label.x.npc="right", hjust=1) +
  labs(col = "cg[CG]xx or xx[CG]cg", 
       x = "median t-statistic for probes with given sequence context",
       y = "oxidation kinetics value\n10.1038/s42003-022-03033-4\nAdams et al. NatCom 2021") +
  theme_cellpress

ggsave("output/figures/vis_differential__xxCGxx_x_TET1_2_NatCom.pdf", width = 11*0.975, height=3)


## motif: xx[CG]xx * TET1-2-3 SciAdv ----



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
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



plt <- plt |> 
  dplyr::left_join(plt.per.context, by=c('gc_sequence_context_2'='gc_sequence_context_2'), suffix=c('','')) |> 
  dplyr::group_by(facet_name) |> 
  dplyr::mutate(context_median__DMP__g2_g3__pp_nc__t = median(DMP__g2_g3__pp_nc__t)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(facet_name = tolower(facet_name)) |> 
  dplyr::select(facet_name, gc_sequence_context_2, context_median__DMP__g2_g3__pp_nc__t) |> 
  dplyr::distinct() |> 
  dplyr::mutate(adjacent_CG = grepl("^cg",facet_name) | grepl("cg$",facet_name)) 





tmp <- read.table("data/10.1126_sciadv.abm2427_Slope_summary.txt") |> 
  tibble::rownames_to_column('Motif2') |> 
  assertr::verify(Motif == Motif2) |> 
  dplyr::mutate(Motif2 = NULL) |> 
  dplyr::rename(sequence = Motif)



plt <- plt |> 
  dplyr::mutate(sequence = gsub("[","",gsub("]","",gc_sequence_context_2),fixed=T),fixed=T) |> 
  dplyr::left_join(tmp, by=c('sequence'='sequence'), suffix=c('','')) 


plt.per.facet <- plt |> 
  dplyr::mutate(gc_sequence_context_2 = NULL) |> 
  dplyr::mutate(sequence = NULL) |> 
  dplyr::group_by(facet_name) |> 
  #dplyr::mutate(ALL_72_per_facet = mean(ALL_72)) |> # 72H, not as powerful 
  dplyr::mutate(ALL_84_per_facet = mean(ALL_84)) |> 
  #dplyr::mutate(ALL_Dec_per_facet = mean(ALL_Dec)) |> 
  
  #dplyr::mutate(CGI_72_per_facet = mean(CGI_72)) |> 
  #dplyr::mutate(CGI_84_per_facet = mean(CGI_84)) |> 
  #dplyr::mutate(CGI_Dec_per_facet = mean(CGI_Dec)) |> 
  
  #dplyr::mutate(non_CGI_72_per_facet = mean(non_CGI_72)) |> 
  #dplyr::mutate(non_CGI_84_per_facet = mean(non_CGI_84)) |> 
  #dplyr::mutate(non_CGI_Dec_per_facet = mean(non_CGI_Dec)) |> 
  
  
  dplyr::ungroup() |> 
  dplyr::mutate(ALL_84 = NULL) |> 
  dplyr::mutate(ALL_Dec = NULL) |> 
  dplyr::mutate(CGI_Dec = NULL) |> 
  dplyr::mutate(non_CGI_Dec = NULL) |> 
  dplyr::mutate(ALL_72 = NULL) |> 
  dplyr::mutate(CGI_84 = NULL) |> 
  dplyr::mutate(non_CGI_84 = NULL) |> 
  dplyr::mutate(CGI_72 = NULL) |> 
  dplyr::mutate(non_CGI_72 = NULL) |> 
  dplyr::distinct() |> 
  tidyr::pivot_longer(cols=contains("_per_facet"), names_to="oxidation_kinetics")


ggplot(plt.per.facet, aes(x=context_median__DMP__g2_g3__pp_nc__t, y=value, label=facet_name, col=adjacent_CG)) +
  facet_wrap(~oxidation_kinetics, scales="free",ncol=5) +
  geom_point() +
  ggrepel::geom_text_repel(data = subset(plt.per.facet, grepl("aacgtt|aacgtc|aacgtg", facet_name)), col="black", size=theme_cellpress_size, alpha=0.5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho", size=theme_cellpress_size, label.x.npc="right", hjust=1) +
  labs(col = "cg[CG]xx or xx[CG]cg", 
       x = "median t-statistic for probes with given sequence context",
       y = "TET?3? de-meth velocity ALL 84h\nSciAdv 10.1126/sciadv.abm2427") +
  theme_cellpress

ggsave("output/figures/vis_differential__xxCGxx_x_TET1_2_3_SciAdv.pdf", width = 11*0.975/4, height=3)




plt <- plt |> 
  dplyr::left_join(tmp, by=c('sequence'='Motif'), suffix=c('','')) |> 
  tibble::column_to_rownames('sequence')


corrplot::corrplot(cor(plt, method="spearman"))





### motif: gc_sequence_context_1 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_1)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) 


stat <- plt |> 
  dplyr::select(gc_sequence_context_1, DMP__g2_g3__pp_nc__t)

ggplot(stat, aes(x=gc_sequence_context_1, y=DMP__g2_g3__pp_nc__t)) +
  ggbeeswarm::geom_quasirandom(size=theme_cellpress_size/2) +
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
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
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
  
  geom_point(data = subset(pplt1, label == "other"), pch=19, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(pplt1, label != "other"), pch=19, cex=0.15, alpha=0.15) + 
  
  geom_point(data=head(pplt1, n=1),pch=19, cex=1, alpha=0.8) + # for label
  
  labs(x = "Per probe t-score Grade2 ~ Grade 3",
       y = "Per probe t-score Primary ~ Recurrence") +
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
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
  
  geom_point(data = subset(pplt2, label == "other"), pch=19, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(pplt2, label != "other"), pch=19, cex=0.15, alpha=0.15) + 
  
  geom_point(data=head(pplt2, n=1),pch=19, cex=1, alpha=0.8) + # for label
  
  labs(x = "Per probe t-score Grade2 ~ Grade 3",
       y = "Per probe t-score Primary ~ Recurrence") +
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
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
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_l2 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l2 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_l3 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l3 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 




### motif: gc_sequence_context_l4 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l4 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_l5 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l5 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 


### motif: gc_sequence_context_l6 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l6 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 


### motif: gc_sequence_context_l7 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l7 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 


### motif: gc_sequence_context_l8 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l8 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_r1 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r1 %in% c('A','T'))) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
  


### motif: gc_sequence_context_r2 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r2 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_r3 ----




ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r3 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox
  



### motif: gc_sequence_context_r4 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r4 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox




### motif: gc_sequence_context_r5 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r5 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox




### motif: gc_sequence_context_r6 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r6 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox



### motif: gc_sequence_context_r7 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r7 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox




### motif: gc_sequence_context_r8 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r8 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox




## x: all epigenetic clocks ----


clocks <- plt |> 
  dplyr::select(contains("dnaMethyAge") | contains("epiTOC2")) |> 
  colnames()


for(clock in clocks) {
  
  clock_name <- gsub("__up_nc__t","",gsub("DMP__","", clock))
  
  plt.p <- plt |> 
    data.table::copy() |> # odd hack needed, because setnames also affects the former "glass_od.metadata.array_samples" object...
    data.table::setnames(old = c(clock), new = c("array_current_clock"))
    
  
  ggplot(plt.p, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=array_current_clock)) +
    geom_vline(xintercept=0, col="red") +
    geom_hline(yintercept=0, col="red") +
    
    geom_point(pch=19, cex=0.001, alpha=0.15) + 
    
    labs(col = gsub("_", " ",clock_name)) +
    
    theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
    
    theme(plot.background = element_rect(fill="white")) + # png export
    
    ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", oob = scales::squish)
  
  
  ggsave(paste0("output/figures/vis_differential__epiGenetic_clock__",clock_name,".png"), width = (8.5*0.975/2), height = 4.3)
  
  
  
  
}

## early ~ late repli ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (760405)) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693017)) 
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
  ggpubr::stat_compare_means(aes(group=glass_nl_prim_rec__deep_significant), label.x.npc=0.1, method = "wilcox.test", show.lengend  = FALSE,  size=theme_cellpress_size)


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
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqBg02esWaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqBjWaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqBjWaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqBjWaveSignalRep2, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white"))  # png export
  #geom_vline(xintercept=79, col="red", lwd=0.5)
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqBjWaveSignalRep2.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm06990WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm06990WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm12801WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm12801WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm12812WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm12812WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm12813WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm12813WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm12878WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm12878WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqHelas3WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqHelas3WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqHepg2WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqHepg2WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqHuvecWaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqHuvecWaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqImr90WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqImr90WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqK562WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqK562WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqMcf7WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqMcf7WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqNhekWaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqNhekWaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqSknshWaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=19,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqSknshWaveSignalRep1.png",width=8.5/2, height=8.5/2)




