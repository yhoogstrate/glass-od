#!/usr/bin/env R


# load ----


source('scripts/load_palette.R')

if(!exixts('data.mvalues.probes')) {
  source('scripts/load_mvalues_hq_samples.R')
}



# plot ----



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



p1 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.015, col="black") + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)
#p1




p2 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__FFPE_decay_time__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)
#p2




p3 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__AcCGAP__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)
#p3



p4 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=-10*glass_nl_prim_rec__deep_significant)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == F),pch=19, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == T),pch=19, cex=0.1, alpha=0.35) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)
#p4



p5 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__epiTOC2_hypoSC__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == F),pch=19, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == T),pch=19, cex=0.1, alpha=0.35) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)
#p5



(p1 + p3 + p4) / (p2 + p5 + p1)


## n seq diff ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=log(1+probeCpGcnt))) +
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
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits=c(log(2),log(13))) + # , oob = scales::squish
  theme(legend.key.size = unit(0.6, 'lines')) 



## sequence motifs ----


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


a = plt |> 
  dplyr::filter(!is.na(gc_sequence_context_1)) |> 
  dplyr::select(gc_sequence_context_1) |> 
  dplyr::filter(!duplicated(gc_sequence_context_1)) |> 
  dplyr::mutate(gc_sequence_context_1_s = gsub("[CG]","CG",gc_sequence_context_1,fixed=T)) |> 
  dplyr::mutate(gc_sequence_context_1_rc = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(gc_sequence_context_1_s)))) |> 
  dplyr::mutate(facet_name = ifelse(gc_sequence_context_1_s <= gc_sequence_context_1_rc,
                                    paste0(gc_sequence_context_1_s ,"_",gc_sequence_context_1_rc),
                                    paste0(gc_sequence_context_1_rc,"_",gc_sequence_context_1_s ))) |> 
  dplyr::arrange(facet_name) |> 
  dplyr::mutate(gc_sequence_context_1_s = NULL, gc_sequence_context_1_rc = NULL)


pplt <- plt |> 
  dplyr::left_join(a, by=c('gc_sequence_context_1'='gc_sequence_context_1')) 



ggplot(pplt, aes(x=gc_sequence_context_1, y=DMP__g2_g3__pp_nc__t)) +
  facet_wrap(~facet_name, scales="free_x", ncol=136) +
  ggbeeswarm::geom_quasirandom(size=0.1) +
  ggplot2::geom_violin(draw_quantiles = c(0.5), col = "white", fill = alpha("white", 0)) +
  theme_cellpress + 
  theme(axis.text.x = element_text(angle = 90))



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


#### RC ----


a = plt |> 
  dplyr::filter(!is.na(gc_sequence_context_2)) |> 
  dplyr::select(gc_sequence_context_2) |> 
  dplyr::filter(!duplicated(gc_sequence_context_2)) |> 
  dplyr::mutate(gc_sequence_context_2_s = gsub("[CG]","CG",gc_sequence_context_2,fixed=T)) |> 
  dplyr::mutate(gc_sequence_context_2_rc = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(gc_sequence_context_2_s)))) |> 
  dplyr::mutate(palindrom = gc_sequence_context_2_s == gc_sequence_context_2_rc) |> 
  dplyr::mutate(facet_name = dplyr::case_when(
    gc_sequence_context_2_s < gc_sequence_context_2_rc ~ paste0(gc_sequence_context_2_s ,"_",gc_sequence_context_2_rc),
    gc_sequence_context_2_s > gc_sequence_context_2_rc ~ paste0(gc_sequence_context_2_rc,"_",gc_sequence_context_2_s ),
    gc_sequence_context_2_s == gc_sequence_context_2_rc ~ paste0(gc_sequence_context_2_s, " *")
  )) |> 
  dplyr::arrange(facet_name) |> 
  dplyr::mutate(gc_sequence_context_2_s = NULL, gc_sequence_context_2_rc = NULL)


pplt <- plt |> 
  dplyr::left_join(a, by=c('gc_sequence_context_2'='gc_sequence_context_2')) 


ggplot(pplt, aes(x=facet_name, y=DMP__g2_g3__pp_nc__t)) +
  ggbeeswarm::geom_quasirandom(size=0.1) +
  ggplot2::geom_violin(draw_quantiles = c(0.5), col = "white", fill = alpha("white", 0)) +
  theme_cellpress + 
  theme(axis.text.x = element_text(angle = 90))




ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l4 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=19, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox



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



