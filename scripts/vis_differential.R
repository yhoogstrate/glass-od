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



