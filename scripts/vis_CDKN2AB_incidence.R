#!/usr/bin/env R

# load ----



# plt ----


tmp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::arrange(patient_id, resection_id)



plt.bins <- do.call(rbind, pbapply::pblapply(
  tmp$array_mnp_CNVP_v12.8_v5.2_CNVP_bins,
  function(x) {
  dat <- read.delim(x)
  
  dat <- dat |> 
    dplyr::mutate(array_sentrix_id = gsub("^X","", colnames(dat)[5])) |> 
    dplyr::rename(log2fc = 5) |> 
    dplyr::mutate(Feature = NULL) |> 
    dplyr::rename(chr = Chromosome, start_hg19 = Start, end_hg19 = End) # checked size of chr1 in this file is 249240621, identical to https://hgdownload.soe.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes

  return(dat)
}
)) |> 
  dplyr::left_join(tmp |> dplyr::select(array_sentrix_id, array_mnp_predictBrain_v12.8_cal_class, patient_id, resection_id), by=c('array_sentrix_id' = 'array_sentrix_id'))



pats <- plt.bins |> 
  dplyr::pull(patient_id) |> 
  unique()



## 1 - 19 ----


plt <- plt.bins |> 
  dplyr::filter(patient_id %in% pats[(1:19)]) |> 
  
  dplyr::filter(chr %in% c("chr9")) |> 
  dplyr::rename(col = chr) |> 
  dplyr::mutate(group = paste0("line-",1:dplyr::n()), type="bin") |> 
  dplyr::mutate(CDKN2AB = ifelse(
    #(start_hg19 >= 21967751 & start_hg19 <= 22009304) | (end_hg19 >= 21967751 & end_hg19 <= 22009304)
    (start_hg19 >= 21780758 & start_hg19 <= 22196297) | (end_hg19 >= 21780758 & end_hg19 <= 22196297)
    ,"CDKN2A/B","chr9")) |> 
  tidyr::pivot_longer(cols = c(start_hg19, end_hg19), values_to = "pos",names_to = "segment_point") |> 
  dplyr::mutate(facet = paste0(resection_id, "")) |> 
  dplyr::mutate(pos = pos/1000000)
#  dplyr::filter(pos < 50) # 9p


plt.purity <- tmp |> 
  dplyr::filter(resection_id %in% plt$resection_id) |> 
  dplyr::mutate(col="purity") |> 
  dplyr::rename(log2fc = array_methylation_bins_1p19q_median.lfc) |> 
  dplyr::mutate(group = array_sentrix_id) |> 
  dplyr::mutate(type="purity estimate") |> 
  dplyr::mutate(facet = paste0(resection_id, "")) |> 
  dplyr::mutate(CDKN2AB = F)

plt.purity <- rbind(
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 1")) |> 
    dplyr::mutate(type = paste0(type,"")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet, CDKN2AB, patient_id) |> 
    dplyr::mutate(segment_point = "start") |> 
    dplyr::mutate(pos = 0)
  ,
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 1")) |> 
    dplyr::mutate(type = paste0(type,"")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet, CDKN2AB, patient_id) |> 
    dplyr::mutate(segment_point = "end") |> 
    dplyr::mutate(pos = max(plt$pos))
)


plt <- rbind(plt, plt.purity)





plt <- plt |> 
  dplyr::mutate(resection = gsub("^.+\\-R","R", facet))


ggplot(plt, aes(x=pos, y=log2fc, group=group, col=CDKN2AB)) +
  facet_grid(rows = vars(patient_id), cols=vars(resection)) +
  #facet_wrap(~facet) + # , scales="free"
  coord_cartesian(ylim = c(-1.5, 1.3)) +
  
  geom_hline(yintercept = 0, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = -1, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = 1.5, col="gray80", lty=2, lwd=theme_nature_lwd) +
  
  #geom_vline(xintercept = (21.967751+21.994391)/2, col="lightgreen", lty=1, lwd=theme_nature_lwd) +
  #geom_vline(xintercept = (22.002902+22.009304)/2, col="lightgreen", lty=1, lwd=theme_nature_lwd) +
  
  
  geom_line(data=subset(plt, type=="bin" & CDKN2AB == "chr9"), lwd=theme_nature_lwd * 6) +
  geom_line(data=subset(plt, type=="bin" & CDKN2AB != "chr9"), lwd=theme_nature_lwd * 15) +
  geom_line(data=subset(plt, type=="purity estimate"), lwd=theme_nature_lwd, lty=1, col="darkblue") +
  
  #geom_line(data=subset(plt, type=="y-bar"), lwd=theme_nature_lwd, lty=1, col="black") +
  
  labs(x = "million bases", col="") +
  
  scale_color_manual(values=c(`CDKN2A/B` = "red",
                              `FALSE` = "darkblue",
                              `chr9` = "#DDDDDD")) +
  
  
  theme_nature +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=theme_nature_lwd)) +
  theme(plot.background = element_rect(fill="white", colour = NA)) # png export, white background, no border



ggsave("output/figures/vis_CDKN2AB_incidence__1_19.png", width=8.5 * 0.975, height = 11 * 0.975, dpi=600)






## 20 - 38 ----


plt <- plt.bins |> 
  dplyr::filter(patient_id %in% pats[(1:19) + 19]) |> 
  
  dplyr::filter(chr %in% c("chr9")) |> 
  dplyr::rename(col = chr) |> 
  dplyr::mutate(group = paste0("line-",1:dplyr::n()), type="bin") |> 
  dplyr::mutate(CDKN2AB = ifelse(
    #(start_hg19 >= 21967751 & start_hg19 <= 22009304) | (end_hg19 >= 21967751 & end_hg19 <= 22009304)
    (start_hg19 >= 21780758 & start_hg19 <= 22196297) | (end_hg19 >= 21780758 & end_hg19 <= 22196297)
    ,"CDKN2A/B","chr9")) |> 
  tidyr::pivot_longer(cols = c(start_hg19, end_hg19), values_to = "pos",names_to = "segment_point") |> 
  dplyr::mutate(facet = paste0(resection_id, "")) |> 
  dplyr::mutate(pos = pos/1000000)
#  dplyr::filter(pos < 50) # 9p


plt.purity <- tmp |> 
  dplyr::filter(resection_id %in% plt$resection_id) |> 
  dplyr::mutate(col="purity") |> 
  dplyr::rename(log2fc = array_methylation_bins_1p19q_median.lfc) |> 
  dplyr::mutate(group = array_sentrix_id) |> 
  dplyr::mutate(type="purity estimate") |> 
  dplyr::mutate(facet = paste0(resection_id, "")) |> 
  dplyr::mutate(CDKN2AB = F)


plt.purity <- rbind(
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 1")) |> 
    dplyr::mutate(type = paste0(type,"")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet, CDKN2AB, patient_id) |> 
    dplyr::mutate(segment_point = "start") |> 
    dplyr::mutate(pos = 0)
  ,
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 1")) |> 
    dplyr::mutate(type = paste0(type,"")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet, CDKN2AB, patient_id) |> 
    dplyr::mutate(segment_point = "end") |> 
    dplyr::mutate(pos = max(plt$pos))
)


plt <- rbind(plt, plt.purity)





plt <- plt |> 
  dplyr::mutate(resection = gsub("^.+\\-R","R", facet))


ggplot(plt, aes(x=pos, y=log2fc, group=group, col=CDKN2AB)) +
  facet_grid(rows = vars(patient_id), cols=vars(resection)) +
  #facet_wrap(~facet) + # , scales="free"
  coord_cartesian(ylim = c(-1.5, 1.3)) +
  
  geom_hline(yintercept = 0, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = -1, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = 1.5, col="gray80", lty=2, lwd=theme_nature_lwd) +
  
  #geom_vline(xintercept = (21.967751+21.994391)/2, col="lightgreen", lty=1, lwd=theme_nature_lwd) +
  #geom_vline(xintercept = (22.002902+22.009304)/2, col="lightgreen", lty=1, lwd=theme_nature_lwd) +
  
  
  geom_line(data=subset(plt, type=="bin" & CDKN2AB == "chr9"), lwd=theme_nature_lwd * 6) +
  geom_line(data=subset(plt, type=="bin" & CDKN2AB != "chr9"), lwd=theme_nature_lwd * 15) +
  geom_line(data=subset(plt, type=="purity estimate"), lwd=theme_nature_lwd, lty=1, col="darkblue") +
  
  #geom_line(data=subset(plt, type=="y-bar"), lwd=theme_nature_lwd, lty=1, col="black") +
  
  labs(x = "million bases", col="") +
  
  scale_color_manual(values=c(`CDKN2A/B` = "red",
                              `FALSE` = "darkblue",
                              `chr9` = "#DDDDDD")) +
  
  
  theme_nature +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=theme_nature_lwd)) +
  theme(plot.background = element_rect(fill="white", colour = NA)) # png export, white background, no border



ggsave("output/figures/vis_CDKN2AB_incidence__20_38.png", width=8.5 * 0.975, height = 11 * 0.975, dpi=600)


## 39 - 57 ----


plt <- plt.bins |> 
  dplyr::filter(patient_id %in% pats[(1:19) + 19 + 19]) |> 
  
  dplyr::filter(chr %in% c("chr9")) |> 
  dplyr::rename(col = chr) |> 
  dplyr::mutate(group = paste0("line-",1:dplyr::n()), type="bin") |> 
  dplyr::mutate(CDKN2AB = ifelse(
    #(start_hg19 >= 21967751 & start_hg19 <= 22009304) | (end_hg19 >= 21967751 & end_hg19 <= 22009304)
    (start_hg19 >= 21780758 & start_hg19 <= 22196297) | (end_hg19 >= 21780758 & end_hg19 <= 22196297)
    ,"CDKN2A/B","chr9")) |> 
  tidyr::pivot_longer(cols = c(start_hg19, end_hg19), values_to = "pos",names_to = "segment_point") |> 
  dplyr::mutate(facet = paste0(resection_id, "")) |> 
  dplyr::mutate(pos = pos/1000000)
#  dplyr::filter(pos < 50) # 9p


plt.purity <- tmp |> 
  dplyr::filter(resection_id %in% plt$resection_id) |> 
  dplyr::mutate(col="purity") |> 
  dplyr::rename(log2fc = array_methylation_bins_1p19q_median.lfc) |> 
  dplyr::mutate(group = array_sentrix_id) |> 
  dplyr::mutate(type="purity estimate") |> 
  dplyr::mutate(facet = paste0(resection_id, "")) |> 
  dplyr::mutate(CDKN2AB = F)


plt.purity <- rbind(
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 1")) |> 
    dplyr::mutate(type = paste0(type,"")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet, CDKN2AB, patient_id) |> 
    dplyr::mutate(segment_point = "start") |> 
    dplyr::mutate(pos = 0)
  ,
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 1")) |> 
    dplyr::mutate(type = paste0(type,"")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet, CDKN2AB, patient_id) |> 
    dplyr::mutate(segment_point = "end") |> 
    dplyr::mutate(pos = max(plt$pos))
)


plt <- rbind(plt, plt.purity)





plt <- plt |> 
  dplyr::mutate(resection = gsub("^.+\\-R","R", facet))


ggplot(plt, aes(x=pos, y=log2fc, group=group, col=CDKN2AB)) +
  facet_grid(rows = vars(patient_id), cols=vars(resection)) +
  #facet_wrap(~facet) + # , scales="free"
  coord_cartesian(ylim = c(-1.5, 1.3)) +
  
  geom_hline(yintercept = 0, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = -1, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = 1.5, col="gray80", lty=2, lwd=theme_nature_lwd) +
  
  #geom_vline(xintercept = (21.967751+21.994391)/2, col="lightgreen", lty=1, lwd=theme_nature_lwd) +
  #geom_vline(xintercept = (22.002902+22.009304)/2, col="lightgreen", lty=1, lwd=theme_nature_lwd) +
  
  
  geom_line(data=subset(plt, type=="bin" & CDKN2AB == "chr9"), lwd=theme_nature_lwd * 6) +
  geom_line(data=subset(plt, type=="bin" & CDKN2AB != "chr9"), lwd=theme_nature_lwd * 15) +
  geom_line(data=subset(plt, type=="purity estimate"), lwd=theme_nature_lwd, lty=1, col="darkblue") +
  
  #geom_line(data=subset(plt, type=="y-bar"), lwd=theme_nature_lwd, lty=1, col="black") +
  
  labs(x = "million bases", col="") +
  
  scale_color_manual(values=c(`CDKN2A/B` = "red",
                              `FALSE` = "darkblue",
                              `chr9` = "#DDDDDD")) +
  
  
  theme_nature +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=theme_nature_lwd)) +
  theme(plot.background = element_rect(fill="white", colour = NA)) # png export, white background, no border



ggsave("output/figures/vis_CDKN2AB_incidence__39_57.png", width=8.5 * 0.975, height = 11 * 0.975, dpi=600)


## 58 - 76 ----


plt <- plt.bins |> 
  dplyr::filter(patient_id %in% pats[(1:19) + 19 + 19 + 19]) |> 
  
  dplyr::filter(chr %in% c("chr9")) |> 
  dplyr::rename(col = chr) |> 
  dplyr::mutate(group = paste0("line-",1:dplyr::n()), type="bin") |> 
  dplyr::mutate(CDKN2AB = ifelse(
    #(start_hg19 >= 21967751 & start_hg19 <= 22009304) | (end_hg19 >= 21967751 & end_hg19 <= 22009304)
    (start_hg19 >= 21780758 & start_hg19 <= 22196297) | (end_hg19 >= 21780758 & end_hg19 <= 22196297)
    ,"CDKN2A/B","chr9")) |> 
  tidyr::pivot_longer(cols = c(start_hg19, end_hg19), values_to = "pos",names_to = "segment_point") |> 
  dplyr::mutate(facet = paste0(resection_id, "")) |> 
  dplyr::mutate(pos = pos/1000000)
#  dplyr::filter(pos < 50) # 9p


plt.purity <- tmp |> 
  dplyr::filter(resection_id %in% plt$resection_id) |> 
  dplyr::mutate(col="purity") |> 
  dplyr::rename(log2fc = array_methylation_bins_1p19q_median.lfc) |> 
  dplyr::mutate(group = array_sentrix_id) |> 
  dplyr::mutate(type="purity estimate") |> 
  dplyr::mutate(facet = paste0(resection_id, "")) |> 
  dplyr::mutate(CDKN2AB = F)


plt.purity <- rbind(
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 1")) |> 
    dplyr::mutate(type = paste0(type,"")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet, CDKN2AB, patient_id) |> 
    dplyr::mutate(segment_point = "start") |> 
    dplyr::mutate(pos = 0)
  ,
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 1")) |> 
    dplyr::mutate(type = paste0(type,"")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet, CDKN2AB, patient_id) |> 
    dplyr::mutate(segment_point = "end") |> 
    dplyr::mutate(pos = max(plt$pos))
)


plt <- rbind(plt, plt.purity)





plt <- plt |> 
  dplyr::mutate(resection = gsub("^.+\\-R","R", facet))


ggplot(plt, aes(x=pos, y=log2fc, group=group, col=CDKN2AB)) +
  facet_grid(rows = vars(patient_id), cols=vars(resection)) +
  #facet_wrap(~facet) + # , scales="free"
  coord_cartesian(ylim = c(-1.5, 1.3)) +
  
  geom_hline(yintercept = 0, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = -1, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = 1.5, col="gray80", lty=2, lwd=theme_nature_lwd) +
  
  #geom_vline(xintercept = (21.967751+21.994391)/2, col="lightgreen", lty=1, lwd=theme_nature_lwd) +
  #geom_vline(xintercept = (22.002902+22.009304)/2, col="lightgreen", lty=1, lwd=theme_nature_lwd) +
  
  
  geom_line(data=subset(plt, type=="bin" & CDKN2AB == "chr9"), lwd=theme_nature_lwd * 6) +
  geom_line(data=subset(plt, type=="bin" & CDKN2AB != "chr9"), lwd=theme_nature_lwd * 15) +
  geom_line(data=subset(plt, type=="purity estimate"), lwd=theme_nature_lwd, lty=1, col="darkblue") +
  
  #geom_line(data=subset(plt, type=="y-bar"), lwd=theme_nature_lwd, lty=1, col="black") +
  
  labs(x = "million bases", col="") +
  
  scale_color_manual(values=c(`CDKN2A/B` = "red",
                              `FALSE` = "darkblue",
                              `chr9` = "#DDDDDD")) +
  
  
  theme_nature +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=theme_nature_lwd)) +
  theme(plot.background = element_rect(fill="white", colour = NA)) # png export, white background, no border



ggsave("output/figures/vis_CDKN2AB_incidence__58_76.png", width=8.5 * 0.975, height = 11 * 0.975, dpi=600)




## 77 - 95 ----


plt <- plt.bins |> 
  dplyr::filter(patient_id %in% pats[(1:19) + 19 + 19 + 19 + 19]) |> 
  
  dplyr::filter(chr %in% c("chr9")) |> 
  dplyr::rename(col = chr) |> 
  dplyr::mutate(group = paste0("line-",1:dplyr::n()), type="bin") |> 
  dplyr::mutate(CDKN2AB = ifelse(
    #(start_hg19 >= 21967751 & start_hg19 <= 22009304) | (end_hg19 >= 21967751 & end_hg19 <= 22009304)
    (start_hg19 >= 21780758 & start_hg19 <= 22196297) | (end_hg19 >= 21780758 & end_hg19 <= 22196297)
    ,"CDKN2A/B","chr9")) |> 
  tidyr::pivot_longer(cols = c(start_hg19, end_hg19), values_to = "pos",names_to = "segment_point") |> 
  dplyr::mutate(facet = paste0(resection_id, "")) |> 
  dplyr::mutate(pos = pos/1000000)
#  dplyr::filter(pos < 50) # 9p


plt.purity <- tmp |> 
  dplyr::filter(resection_id %in% plt$resection_id) |> 
  dplyr::mutate(col="purity") |> 
  dplyr::rename(log2fc = array_methylation_bins_1p19q_median.lfc) |> 
  dplyr::mutate(group = array_sentrix_id) |> 
  dplyr::mutate(type="purity estimate") |> 
  dplyr::mutate(facet = paste0(resection_id, "")) |> 
  dplyr::mutate(CDKN2AB = F)


plt.purity <- rbind(
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 1")) |> 
    dplyr::mutate(type = paste0(type,"")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet, CDKN2AB, patient_id) |> 
    dplyr::mutate(segment_point = "start") |> 
    dplyr::mutate(pos = 0)
  ,
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 1")) |> 
    dplyr::mutate(type = paste0(type,"")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet, CDKN2AB, patient_id) |> 
    dplyr::mutate(segment_point = "end") |> 
    dplyr::mutate(pos = max(plt$pos))
)


plt <- rbind(plt, plt.purity)





plt <- plt |> 
  dplyr::mutate(resection = gsub("^.+\\-R","R", facet))


ggplot(plt, aes(x=pos, y=log2fc, group=group, col=CDKN2AB)) +
  facet_grid(rows = vars(patient_id), cols=vars(resection)) +
  #facet_wrap(~facet) + # , scales="free"
  coord_cartesian(ylim = c(-1.5, 1.3)) +
  
  geom_hline(yintercept = 0, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = -1, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = 1.5, col="gray80", lty=2, lwd=theme_nature_lwd) +
  
  #geom_vline(xintercept = (21.967751+21.994391)/2, col="lightgreen", lty=1, lwd=theme_nature_lwd) +
  #geom_vline(xintercept = (22.002902+22.009304)/2, col="lightgreen", lty=1, lwd=theme_nature_lwd) +
  
  
  geom_line(data=subset(plt, type=="bin" & CDKN2AB == "chr9"), lwd=theme_nature_lwd * 6) +
  geom_line(data=subset(plt, type=="bin" & CDKN2AB != "chr9"), lwd=theme_nature_lwd * 15) +
  geom_line(data=subset(plt, type=="purity estimate"), lwd=theme_nature_lwd, lty=1, col="darkblue") +
  
  #geom_line(data=subset(plt, type=="y-bar"), lwd=theme_nature_lwd, lty=1, col="black") +
  
  labs(x = "million bases", col="") +
  
  scale_color_manual(values=c(`CDKN2A/B` = "red",
                              `FALSE` = "darkblue",
                              `chr9` = "#DDDDDD")) +
  
  
  theme_nature +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=theme_nature_lwd)) +
  theme(plot.background = element_rect(fill="white", colour = NA)) # png export, white background, no border



ggsave("output/figures/vis_CDKN2AB_incidence__77_95.png", width=8.5 * 0.975, height = 11 * 0.975, dpi=600)



## 96 - 114 ----


plt <- plt.bins |> 
  dplyr::filter(patient_id %in% pats[(1:19) + 19 + 19 + 19 + 19 + 19]) |> 
  
  dplyr::filter(chr %in% c("chr9")) |> 
  dplyr::rename(col = chr) |> 
  dplyr::mutate(group = paste0("line-",1:dplyr::n()), type="bin") |> 
  dplyr::mutate(CDKN2AB = ifelse(
    #(start_hg19 >= 21967751 & start_hg19 <= 22009304) | (end_hg19 >= 21967751 & end_hg19 <= 22009304)
    (start_hg19 >= 21780758 & start_hg19 <= 22196297) | (end_hg19 >= 21780758 & end_hg19 <= 22196297)
    ,"CDKN2A/B","chr9")) |> 
  tidyr::pivot_longer(cols = c(start_hg19, end_hg19), values_to = "pos",names_to = "segment_point") |> 
  dplyr::mutate(facet = paste0(resection_id, "")) |> 
  dplyr::mutate(pos = pos/1000000)
#  dplyr::filter(pos < 50) # 9p


plt.purity <- tmp |> 
  dplyr::filter(resection_id %in% plt$resection_id) |> 
  dplyr::mutate(col="purity") |> 
  dplyr::rename(log2fc = array_methylation_bins_1p19q_median.lfc) |> 
  dplyr::mutate(group = array_sentrix_id) |> 
  dplyr::mutate(type="purity estimate") |> 
  dplyr::mutate(facet = paste0(resection_id, "")) |> 
  dplyr::mutate(CDKN2AB = F)


plt.purity <- rbind(
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 1")) |> 
    dplyr::mutate(type = paste0(type,"")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet, CDKN2AB, patient_id) |> 
    dplyr::mutate(segment_point = "start") |> 
    dplyr::mutate(pos = 0)
  ,
  plt.purity |> 
    dplyr::mutate(group = paste0(group," 1")) |> 
    dplyr::mutate(type = paste0(type,"")) |> 
    dplyr::select(col, log2fc, array_sentrix_id,  array_mnp_predictBrain_v12.8_cal_class, resection_id, group, type, facet, CDKN2AB, patient_id) |> 
    dplyr::mutate(segment_point = "end") |> 
    dplyr::mutate(pos = max(plt$pos))
)


plt <- rbind(plt, plt.purity)





plt <- plt |> 
  dplyr::mutate(resection = gsub("^.+\\-R","R", facet))


ggplot(plt, aes(x=pos, y=log2fc, group=group, col=CDKN2AB)) +
  facet_grid(rows = vars(patient_id), cols=vars(resection)) +
  #facet_wrap(~facet) + # , scales="free"
  coord_cartesian(ylim = c(-1.5, 1.3)) +
  
  geom_hline(yintercept = 0, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = -1, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = 1.5, col="gray80", lty=2, lwd=theme_nature_lwd) +
  
  #geom_vline(xintercept = (21.967751+21.994391)/2, col="lightgreen", lty=1, lwd=theme_nature_lwd) +
  #geom_vline(xintercept = (22.002902+22.009304)/2, col="lightgreen", lty=1, lwd=theme_nature_lwd) +
  
  
  geom_line(data=subset(plt, type=="bin" & CDKN2AB == "chr9"), lwd=theme_nature_lwd * 6) +
  geom_line(data=subset(plt, type=="bin" & CDKN2AB != "chr9"), lwd=theme_nature_lwd * 15) +
  geom_line(data=subset(plt, type=="purity estimate"), lwd=theme_nature_lwd, lty=1, col="darkblue") +
  
  #geom_line(data=subset(plt, type=="y-bar"), lwd=theme_nature_lwd, lty=1, col="black") +
  
  labs(x = "million bases", col="") +
  
  scale_color_manual(values=c(`CDKN2A/B` = "red",
                              `FALSE` = "darkblue",
                              `chr9` = "#DDDDDD")) +
  
  
  theme_nature +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=theme_nature_lwd)) +
  theme(plot.background = element_rect(fill="white", colour = NA)) # png export, white background, no border



ggsave("output/figures/vis_CDKN2AB_incidence__96_114.png", width=8.5 * 0.975, height = 11 * 0.975, dpi=600)



# 9p arm medians ----


sel <- c("0007-R2","0012-R3","0018-R3","0024-R2","0031-R2","0032-R1","0036-R2","0038-R2","0042-R3","0045-R2","0051-R2","0054-R3","0055-R2","0061-R2","0064-R4","0070-R3","0090-R3","0092-R3","0094-R2","0099-R2","0111-R2","0114-R2","0116-R2","0117-R3","0119-R2","0124-R2")


plt <- plt.bins |> 
  dplyr::filter(resection_id %in% sel)  |> 
  dplyr::filter(chr %in% c("chr9")) |> 
  dplyr::rename(col = chr) |> 
  dplyr::mutate(group = paste0("line-",1:dplyr::n()), type="bin") |> 
  dplyr::mutate(CDKN2AB = ifelse(
    #(start_hg19 >= 21967751 & start_hg19 <= 22009304) | (end_hg19 >= 21967751 & end_hg19 <= 22009304)
    (start_hg19 >= 21780758 & start_hg19 <= 22196297) | (end_hg19 >= 21780758 & end_hg19 <= 22196297)
    ,"CDKN2A/B","chr9")) |> 
  
  dplyr::mutate(CDKN2AB = 
    dplyr::case_when(
      (start_hg19 >= 21780758 & start_hg19 <= 22196297) | (end_hg19 >= 21780758 & end_hg19 <= 22196297) ~ "CDKN2A/B",
      end_hg19 < 50000000 ~ "chr9p",
      T ~ "chr9q"
    )) |> 
  
  
  
  dplyr::mutate(bin = paste0(start_hg19, "-", end_hg19)) |> 
  dplyr::select(bin, log2fc, CDKN2AB) |> 
  dplyr::group_by(bin) |> 
  dplyr::summarise(mean_log2fc = mean(log2fc),
                   median_log2fc = median(log2fc),
                   CDKN2AB = unique(CDKN2AB)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(start_hg19 = as.numeric(gsub("\\-.+$","",bin))) |> 
  dplyr::mutate(end_hg19 = as.numeric(gsub("^.+\\-","",bin))) |> 
  tidyr::pivot_longer(cols = c(start_hg19, end_hg19), values_to = "pos",names_to = "segment_point") |> 
  dplyr::mutate(pos = pos/1000000) |> 
  
  dplyr::mutate(log2fc = mean_log2fc)




ggplot(plt, aes(x=pos, y=log2fc, col=CDKN2AB, group=bin)) +
  coord_cartesian(ylim = c(-0.65, 0.125)) +
  
  geom_hline(yintercept = 0, col="gray80", lty=1, lwd=theme_nature_lwd) +
  geom_hline(yintercept = -0.5, col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept = 1.5, col="gray80", lty=2, lwd=theme_nature_lwd) +
  
  #geom_vline(xintercept = 29.8, lwd=theme_nature_lwd) +
  
  
  geom_line(data=subset(plt, CDKN2AB == "chr9"), lwd=theme_nature_lwd * 6) +
  geom_line(data=subset(plt, CDKN2AB != "chr9"), lwd=theme_nature_lwd * 15) +
  
  labs(x = "million bases hg19", y="mean log2fc", col="", caption="last resections with 9p arm losses (n=26, incl 4x HD)") +
  
  scale_color_manual(values=c(`CDKN2A/B` = mixcol("red","black",0.2),
                              `chr9q` = "#CB75A4",
                              `chr9p` = "#59B2E6")) +
  
  scale_x_continuous(breaks=(0:15)*10) +
  scale_y_continuous(breaks=c(0,-0.25,-0.5)) +
  
  theme_nature +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=theme_nature_lwd)) +
  theme(plot.background = element_rect(fill="white", colour = NA)) + # png export, white background, no border
  theme(legend.position = "right")


ggsave("output/figures/vis_CDKN2AB_incidence__9p_arm__mean.png", width=6.5, height=2.25, dpi=600)

