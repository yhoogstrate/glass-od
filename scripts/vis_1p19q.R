#!/usr/bin/env R

# load ----



# plt ----


tmp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG", "OLIGOSARC_IDH")) |> 
  dplyr::arrange(array_mnp_predictBrain_v12.8_cal_class, resection_id) 



plt.bins <- do.call(rbind, pbapply::pblapply(
  tmp$array_mnp_CNVP_v12.8_v5.2_CNVP_bins,
  function(x) {
  dat <- read.delim(x)
  
  dat <- dat |> 
    dplyr::mutate(array_sentrix_id = gsub("^X","", colnames(dat)[5])) |> 
    dplyr::rename(log2fc = 5) |> 
    dplyr::mutate(Feature = NULL) |> 
    dplyr::rename(chr = Chromosome, start = Start, end = End)

  return(dat)
}
))



plt <- plt.bins |> 
  dplyr::filter(chr %in% c("chr1", "chr19")) |> 
  dplyr::rename(col = chr) |> 
  dplyr::left_join(tmp |> dplyr::select(array_sentrix_id, array_mnp_predictBrain_v12.8_cal_class, resection_id), by=c('array_sentrix_id' = 'array_sentrix_id')) |> 
  dplyr::mutate(group = paste0("line-",1:dplyr::n()), type="bin") |>
  tidyr::pivot_longer(cols = c(start, end), values_to = "pos",names_to = "segment_point")  |> 
  
  dplyr::mutate(col = dplyr::case_when(
    col == "chr1" & pos <  132000000 ~ "chr1p",
    col == "chr1" & pos >= 132000000 ~ "chr1q",
    col == "chr19" & pos <  26000000 ~ "chr19p",
    col == "chr19" & pos >= 26000000 ~ "chr19q",
    T  ~ "?"
    )) |> 
  
  dplyr::mutate(pos = ifelse(col %in% c("chr19p", "chr19q"), pos +  260250000, pos)) |> 
  dplyr::mutate(facet = paste0(resection_id, " [", array_mnp_predictBrain_v12.8_cal_class,"]")) |> 
  dplyr::mutate(pos = pos/1000000)


plt.purity <- tmp |> 
  dplyr::mutate(col="purity") |> 
  dplyr::rename(log2fc = array_methylation_bins_1p19q_median.lfc) |> 
  dplyr::mutate(group = array_sentrix_id) |> 
  dplyr::mutate(type="purity estimate") |> 
  dplyr::mutate(facet = paste0(resection_id, " [", array_mnp_predictBrain_v12.8_cal_class,"]"))

plt.purity <- rbind(
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 1")) |> 
    dplyr::mutate(type = paste0(type," 1")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet) |> 
    dplyr::mutate(segment_point = "start") |> 
    dplyr::mutate(pos = 0)
  ,
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 1")) |> 
    dplyr::mutate(type = paste0(type," 1")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet) |> 
    dplyr::mutate(segment_point = "end") |> 
    dplyr::mutate(pos = 130),
  
  
  
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 2")) |> 
    dplyr::mutate(type = paste0(type," 2")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet) |> 
    dplyr::mutate(segment_point = "start") |> 
    dplyr::mutate(pos = 130)
  ,
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 2")) |> 
    dplyr::mutate(type = paste0(type," 2")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet) |> 
    dplyr::mutate(segment_point = "end") |> 
    dplyr::mutate(pos = 287),
  
  
  
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 3")) |> 
    dplyr::mutate(type = paste0(type," 1")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet) |> 
    dplyr::mutate(segment_point = "start") |> 
    dplyr::mutate(pos = 287)
  ,
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 3")) |> 
    dplyr::mutate(type = paste0(type," 1")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet) |> 
    dplyr::mutate(segment_point = "end") |> 
    dplyr::mutate(pos = 323)
)


plt <- rbind(plt, plt.purity)




plt.ybar <- tmp |> 
  dplyr::filter(resection_id %in% c("0003-R3", "0024-R2", "0047-R3", "0070-R3", "0086-R2") == F) |> 
  dplyr::mutate(col="y-bar") |> 
  dplyr::mutate(pos=0) |> 
  dplyr::mutate(group = paste0(array_sentrix_id,"-y-bar")) |> 
  dplyr::mutate(type="y-bar") |> 
  dplyr::mutate(facet = paste0(resection_id, " [", array_mnp_predictBrain_v12.8_cal_class,"]"))

plt.ybar <- rbind(
  plt.ybar |> 
    dplyr::select(col, pos, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet) |> 
    dplyr::mutate(segment_point = "start") |> 
    dplyr::mutate(log2fc = 2)
  ,
  plt.ybar |> 
    dplyr::select(col, pos, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet) |> 
    dplyr::mutate(segment_point = "end") |> 
    dplyr::mutate(log2fc = -2)
)



plt <- rbind(plt, plt.ybar)


ggplot(plt, aes(x=pos, y=log2fc, group=group, col=col)) +
  facet_wrap(~facet) + # , scales="free"
  coord_cartesian(ylim = c(-1.2, 1.7)) +
  
  
  geom_hline(yintercept = 0, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = -1, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = 1.5, col="gray80", lty=2, lwd=theme_nature_lwd) +
  
  
  geom_line(data=subset(plt, type=="bin"), lwd=theme_nature_lwd * 6) +
  geom_line(data=subset(plt, type=="purity estimate 1"), lwd=theme_nature_lwd, lty=1, col="red") +
  geom_line(data=subset(plt, type=="purity estimate 2"), lwd=theme_nature_lwd, lty=3, col="red") +
  
  geom_line(data=subset(plt, type=="y-bar"), lwd=theme_nature_lwd, lty=1, col="black") +
  
  labs(x = "million bases", col="") +
  
  scale_color_manual(values=c(`chr1p` = "darkblue",
                              `chr1q` = "lightblue",
                              `chr19p` = "lightgreen",
                              `chr19q` = "darkgreen")) +

  
  theme_nature +
  theme(plot.background = element_rect(fill="white", colour = NA)) # png export, white background, no border


ggsave("output/figures/vis_1p19q.png", width=8.5 * 0.975, height = 4, dpi=600)


