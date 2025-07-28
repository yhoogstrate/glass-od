#!/usr/bin/env


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


#source('scripts/load_chrom_sizes.R')
source('scripts/load_themes.R')



library(ggplot2)
library(patchwork)



# obtain data ----



metadata.all <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_first_G2_and_last_G3(154) |> 
  dplyr::mutate(gr.status = factor(ifelse(resection_tumor_grade == 2, "Grade2", "Grade3"), levels=c("Grade2", "Grade3")))



## get raw cnv profiles ----


f <- function(fn, purity) {
  #fn = metadata$heidelberg_cnvp_bins[1]
  
  d <- read.delim(fn) |> 
    dplyr::mutate(cnvp_bin = paste0("cnvp_bin_",Chromosome,"_",Start,"_",End)) |> 
    tibble::column_to_rownames('cnvp_bin') |> 
    dplyr::mutate(Chromosome = NULL, Start = NULL, End=NULL, Feature=NULL) 

  # moet wss iets anders...
  d[,1] <- d[,1] / purity
  
  d <- d |>  
    t() |> 
    as.data.frame()
  
  return (d)
}



data.cnv.profiles <- metadata.all |> 
  dplyr::select(array_sentrix_id, 
                array_mnp_CNVP_v12.8_v5.2_CNVP_bins,
                array_mnp_CNVP_v12.8_v5.2_CNVP_bins, array_methylation_bins_1p19q_purity) |>
  #head(n=30) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = f(array_mnp_CNVP_v12.8_v5.2_CNVP_bins, array_methylation_bins_1p19q_purity)) |> 
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(array_mnp_CNVP_v12.8_v5.2_CNVP_bins = NULL) |> 
  dplyr::mutate(array_methylation_bins_1p19q_purity = NULL) |> 
  tibble::column_to_rownames('array_sentrix_id')


stopifnot(is.na(data.cnv.profiles) == F)



metadata.cnv.profiles <- data.frame(cnvp_bin = colnames(data.cnv.profiles)) |> 
  dplyr::mutate(chr = gsub("^.+_(chr[^_]+)_.+$","\\1",cnvp_bin)) |>
  dplyr::filter(chr %in% c("chrX", "chrY") == F ) |>
  dplyr::mutate(chr = factor(chr, levels = gtools::mixedsort(unique(as.character(chr))))) |>
  dplyr::mutate(start = as.numeric(gsub("^.+_chr[^_]+_([0-9]+)_.+$","\\1",cnvp_bin))) |>
  dplyr::mutate(end = as.numeric(gsub("^.+_chr[^_]+_[0-9]+_([0-9]+)$","\\1",cnvp_bin))) |>
  dplyr::mutate(pos = (start + end)/2)


data.cnv.profiles <- data.cnv.profiles |> 
  dplyr::select(metadata.cnv.profiles$cnvp_bin)



# layer g3 ~ g2 ----


metadata.g3.g2 <- metadata.all |> 
  dplyr::filter(!is.na(resection_tumor_grade)) |> 
  filter_first_G2_and_last_G3(154) |>
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder"))))



data.cnv.profiles.g3.g2 <- data.cnv.profiles |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(metadata.g3.g2$array_sentrix_id) |> 
  t() |> 
  as.data.frame() 


## limma: resection_grade ----


design.g3.g2 <- model.matrix(~patient + array_methylation_bins_1p19q_purity + gr.status, data=metadata.g3.g2 |> tibble::column_to_rownames('array_sentrix_id'))


fit <- limma::lmFit(t(data.cnv.profiles.g3.g2), design.g3.g2)
fit <- limma::eBayes(fit)
stats.g3.g2 <- limma::topTable(fit, n=nrow(fit), sort="p", coef="gr.statusGrade3") |> # adjust="BH"
  tibble::rownames_to_column('cnvp_bin')   |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('','')) 




## per group quantiles ----


quantiles.g2 <- metadata.g3.g2 |> 
  dplyr::filter(resection_tumor_grade == 2) |> 
  dplyr::select(array_sentrix_id) |> 
  dplyr::left_join(data.cnv.profiles.g3.g2 |> tibble::rownames_to_column('array_sentrix_id'),
                   by=c('array_sentrix_id' = 'array_sentrix_id'), suffix=c('','')) |> 
  tibble::column_to_rownames('array_sentrix_id') |> 
  pbapply::pblapply(function(x) {return(quantile(x))}) |> 
  as.data.frame() |> 
  t() |> 
  as.data.frame() |> 
  dplyr::rename_with( ~ paste0("quantile.", .x)) |> 
  tibble::rownames_to_column('cnvp_bin')


quantiles.g3 <- metadata.g3.g2 |> 
  dplyr::filter(resection_tumor_grade == 3) |> 
  dplyr::select(array_sentrix_id) |> 
  dplyr::left_join(data.cnv.profiles.g3.g2 |> tibble::rownames_to_column('array_sentrix_id'),
                   by=c('array_sentrix_id' = 'array_sentrix_id'), suffix=c('','')) |> 
  tibble::column_to_rownames('array_sentrix_id') |> 
  pbapply::pblapply(function(x) {return(quantile(x))}) |> 
  as.data.frame() |> 
  t() |> 
  as.data.frame() |> 
  dplyr::rename_with( ~ paste0("quantile.", .x)) |> 
  tibble::rownames_to_column('cnvp_bin')



quantiles.g3.g2 <- rbind(
  quantiles.g2 |> dplyr::mutate(grade = 2),
  quantiles.g3 |> dplyr::mutate(grade = 3)
)





## plot g3 ~ g2 ----


plt <- quantiles.g3.g2 |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('','')) |> 
  dplyr::mutate(grade = as.factor(paste0("Grade ",grade))) |> 
  dplyr::arrange(`quantile.50%`)



p1 <- ggplot(plt, aes(x=pos/100000000,y=`quantile.50%`, col=grade)) +
  geom_hline(yintercept = 0, col="gray40",lwd=theme_cellpress_lwd,lty=1, show_guide = FALSE) +
  facet_wrap(~chr, scales="free_x", ncol=22) +
  geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
  geom_point(pch=20,cex=0.4, stroke=0, col="gray70") +
  #geom_point(pch=19,cex=0.2,alpha=0.15, col="gray70") +
  geom_smooth(se=F, method = "gam", lwd=theme_cellpress_lwd * 2.5, alpha=0.8) +
  coord_cartesian(ylim=c(-1.05, 1.05)) +
  labs(x=NULL, col=NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) + 
  #scale_color_manual(values=palette_g2_g3) +
  scale_color_manual(values=c("Grade 3" = mixcol("red","black",0.1), "Grade 2" = as.character(palette_g2_g3['Grade 3'])))+
  theme_nature
p1




# p2a <- ggplot(stats.g3.g2, aes(x=pos/100000000,y=logFC, col=adj.P.Val)) +
#   facet_wrap(~chr, scales="free_x", ncol=22) +
#   geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
#   #geom_point(pch=19,cex=0.2,alpha=0.05) +
#   geom_point(pch=19,cex=0.2,alpha=0.15) +
#   geom_hline(yintercept = 0, col="red",lwd=theme_cellpress_lwd) +
#   geom_smooth(se=F, lwd=theme_cellpress_lwd * 2, col="black") +
#   coord_cartesian(ylim=c(-0.3, 0.3)) +
#   labs(x=NULL) +
#   scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) +
#   scale_color_manual(values=c(`TRUE`= 'red', `FALSE`='gray60')) +
#   theme_nature
# p2a



p2b <- ggplot(stats.g3.g2, aes(x=pos/100000000,y=t, col=adj.P.Val < 0.01)) +
  facet_wrap(~chr, scales="free_x", ncol=22) +
  geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
  #geom_point(pch=19,cex=0.2,alpha=0.05) +
  #geom_point(pch=19,cex=0.2,alpha=0.15) +
  geom_point(pch=20,cex=0.4, stroke=0) +
  geom_hline(yintercept = 0, col="red",lwd=theme_cellpress_lwd) +
  geom_smooth(se=F, lwd=theme_cellpress_lwd * 2, col="black") +
  coord_cartesian(ylim=c(-6.5, 6.5)) +
  labs(x=NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) +
  scale_color_manual(values=c(`TRUE`= 'red', `FALSE`='gray60'), guide="none") +
  theme_nature
p2b

p2 <- p2b # otherwise logfc between categorial grade and cgc are incomparible



plt <- rbind(stats.g3.g2 |> dplyr::mutate(y = -log(adj.P.Val), type="data"),
             stats.g3.g2 |> dplyr::mutate(y= 0, type="border")
             )



p3 <- ggplot(plt, aes(x=pos/100000000,y=y, group=cnvp_bin)) +
  facet_wrap(~chr, scales="free_x", ncol=22) +
  geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
  geom_line(alpha=0.05,col=mixcol("darkgreen", "white", 0.15),lwd=theme_cellpress_lwd) +
  theme_bw() +
  geom_smooth(data=subset(plt, type == "data"), group=1, se=F, lwd=theme_cellpress_lwd * 2, col="black") +
  geom_hline(yintercept = 0, col="red",lwd=theme_cellpress_lwd,lty=1) +
  geom_hline(yintercept = -log(0.05), col="red",lwd=theme_cellpress_lwd,lty=2) +
  labs(y = "-log(Padj)")+
  coord_cartesian(ylim=c(0, 7.75)) +
  labs(x="Chromosomal position (/ 100 MB)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) +
  theme_nature
p3



p1 / p2 / p3



ggsave("output/figures/vis_cnv_profiles_GISTIC_style__g2_x_g3.png", width=11 * 0.975, height = 3.75) # removal of ticks needs to be incorp





# CGC[Ac] + PC2 + PC3 ----




metadata.LR.PC2 <- metadata.g3.g2 |> 
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  assertr::verify(!is.na(array_PC1)) |> 
  assertr::verify(!is.na(array_PC2)) |> 
  assertr::verify(!is.na(array_PC3))



data.cnv.profiles.LR.PC2 <- data.cnv.profiles |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(metadata.LR.PC2$array_sentrix_id) |> 
  t() |> 
  as.data.frame() 



design.LR <- model.matrix(~patient + array_methylation_bins_1p19q_purity + scale(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit), data=metadata.LR.PC2 |> tibble::column_to_rownames('array_sentrix_id'))
design.PC2 <- model.matrix(~patient + array_methylation_bins_1p19q_purity + scale(-array_PC2),                            data=metadata.LR.PC2 |> tibble::column_to_rownames('array_sentrix_id'))
design.PC3 <- model.matrix(~patient + array_methylation_bins_1p19q_purity + scale(-array_PC3),                            data=metadata.LR.PC2 |> tibble::column_to_rownames('array_sentrix_id'))

design.PC2_PC3 <- model.matrix(~patient + array_methylation_bins_1p19q_purity  + scale(-array_PC2) +  scale(-array_PC3),  data=metadata.LR.PC2 |> tibble::column_to_rownames('array_sentrix_id'))



fit.LR <- limma::lmFit(t(data.cnv.profiles.LR.PC2), design.LR)
fit.LR <- limma::eBayes(fit.LR)
stats.LR <- limma::topTable(fit.LR, n=nrow(fit.LR), sort="p", coef="scale(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)") |> # adjust="BH"
  tibble::rownames_to_column('cnvp_bin')

# 
# fit.PC2 <- limma::lmFit(t(data.cnv.profiles.LR.PC2), design.PC2)
# fit.PC2 <- limma::eBayes(fit.PC2)
# stats.PC2 <- limma::topTable(fit.PC2, n=nrow(fit.PC2), sort="p", coef="scale(-array_PC2)") |> # adjust="BH"
#   tibble::rownames_to_column('cnvp_bin')
# 
# 
# fit.PC3 <- limma::lmFit(t(data.cnv.profiles.LR.PC2), design.PC3)
# fit.PC3 <- limma::eBayes(fit.PC3)
# stats.PC3 <- limma::topTable(fit.PC3, n=nrow(fit.PC3), sort="p", coef="scale(-array_PC3)") |> # adjust="BH"
#   tibble::rownames_to_column('cnvp_bin')
# 


fit.PC2_PC3 <- limma::lmFit(t(data.cnv.profiles.LR.PC2), design.PC2_PC3)
fit.PC2_PC3 <- limma::eBayes(fit.PC2_PC3)

stats.PC2 <- limma::topTable(fit.PC2_PC3, n=nrow(fit.PC2_PC3), sort="p", coef="scale(-array_PC2)") |> # adjust="BH"
  tibble::rownames_to_column('cnvp_bin')
stats.PC3 <- limma::topTable(fit.PC2_PC3, n=nrow(fit.PC2_PC3), sort="p", coef="scale(-array_PC3)") |> # adjust="BH"
  tibble::rownames_to_column('cnvp_bin')



## plot CGC[Ac] ----


plt.LR <- stats.LR |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('',''))


plt.LR.padj <- stats.LR |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('','')) |> 
  dplyr::mutate(minLogPadj = -log(adj.P.Val))

plt.LR.padj <- rbind(
  plt.LR.padj |> dplyr::mutate(type = "data"),
  plt.LR.padj |> dplyr::mutate(type = "border", minLogPadj = 0))


p4 <- ggplot(plt.LR, aes(x=pos/100000000,y=t, col=adj.P.Val<0.01)) +
  facet_wrap(~chr, scales="free_x", ncol=22) +
  geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
  geom_point(pch=20,cex=0.4, stroke=0) +
  #geom_point(pch=19,cex=0.2,alpha=0.05) +
  geom_hline(yintercept = 0, col="red",lwd=theme_cellpress_lwd,lty=1) +
  geom_smooth(se=F, lwd=theme_cellpress_lwd * 2, col="black") +
  coord_cartesian(ylim=c(-6.5, 6.5)) +
  labs(x=NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) +
  scale_color_manual(values=c(`TRUE`= mixcol('red', "gray60", 0.4), `FALSE`='gray60'), guide="none") +
  theme_nature
p4

p5 <- ggplot(plt.LR.padj, aes(x=pos/100000000,y=minLogPadj, group=cnvp_bin)) +
  facet_wrap(~chr, scales="free_x", ncol=22) +
  geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
  geom_line(alpha=0.05,col=mixcol("darkgreen", "white", 0.15),lwd=theme_cellpress_lwd) +
  geom_smooth(data=subset(plt.LR.padj, type == "data"), group=1, se=F, lwd=theme_cellpress_lwd * 2, col="black") +
  geom_hline(yintercept = 0, col="red",lwd=theme_cellpress_lwd,lty=1) +
  geom_hline(yintercept = -log(0.05), col="red",lwd=theme_cellpress_lwd,lty=2) +
  labs(y = "-log(Padj)") +
  coord_cartesian(ylim=c(0, 7.75)) +
  labs(x="Chromosomal position (/ 100 MB)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) +
  theme_nature
p5



p4 / p5


ggsave("output/figures/vis_cnv_profiles_GISTIC_style__AcCGAP.png", width=11 * 0.975, height = 2.46) # removal of ticks needs to be incorp



# ggplot(data.frame(x= c(1), y=c(3)), aes(x=x,y=y)) +
#   #geom_point(cex=1, pch=19, col = "red", fill="white", alpha=0.3) +
#   geom_point(pch=20,cex=15,alpha=0.3, col = "red",stroke=0) +
#   xlim(c(0.9999999,1.0000001)) +
#   ylim(c(2.9999999,3.0000001)) +
#   theme_nature



## PC2 ----


plt.PC2 <- stats.PC2 |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('',''))


plt.PC2.padj <- stats.PC2 |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('','')) |> 
  dplyr::mutate(minLogPadj = -log(adj.P.Val))

plt.PC2.padj <- rbind(
  plt.PC2.padj |> dplyr::mutate(type = "data"),
  plt.PC2.padj |> dplyr::mutate(type = "border", minLogPadj = 0))



p6 <- ggplot(plt.PC2, aes(x=pos/100000000, y=logFC, col=adj.P.Val<0.01)) +
  facet_wrap(~chr, scales="free_x", ncol=22) +
  geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
  geom_point(pch=20, cex=0.4, stroke=0) +
  geom_hline(yintercept = 0, col="red",lwd=theme_cellpress_lwd,lty=1) +
  geom_smooth(se=F, lwd=theme_cellpress_lwd * 2, col="black") +
  coord_cartesian(ylim=c(-0.3, 0.3)) +
  labs(x=NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) +
  scale_color_manual(values=c(`TRUE`= mixcol('red', "gray60", 0.4), `FALSE`='gray60'), guide="none") +
  theme_nature
p6


p7 <- ggplot(plt.PC2.padj, aes(x=pos/100000000,y=minLogPadj, group=cnvp_bin)) +
  facet_wrap(~chr, scales="free_x", ncol=22) +
  geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
  geom_line(alpha=0.05,col=mixcol("darkgreen", "white", 0.15),lwd=theme_cellpress_lwd) +
  geom_smooth(data=subset(plt.PC2.padj, type == "data"), group=1, se=F, lwd=theme_cellpress_lwd * 2, col="black") +
  geom_hline(yintercept = 0, col="red",lwd=theme_cellpress_lwd,lty=1) +
  geom_hline(yintercept = -log(0.05), col="red",lwd=theme_cellpress_lwd,lty=2) +
  labs(y = "-log(Padj)") +
  coord_cartesian(ylim=c(0, 7.75)) +
  labs(x="Chromosomal position (/ 100 MB)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) +
  theme_nature
p7


p6 / p7


ggsave("output/figures/vis_cnv_profiles_GISTIC_style__PC2.png", width=11 * 0.975, height = 2.46) # removal of ticks needs to be incorp




## PC3 ----



plt.PC3 <- stats.PC3 |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('',''))


plt.PC3.padj <- stats.PC3 |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('','')) |> 
  dplyr::mutate(minLogPadj = -log(adj.P.Val))


plt.PC3.padj <- rbind(
  plt.PC3.padj |> dplyr::mutate(type = "data"),
  plt.PC3.padj |> dplyr::mutate(type = "border", minLogPadj = 0))


p8 <- ggplot(plt.PC3, aes(x=pos/100000000,y=logFC, col=adj.P.Val<0.01)) +
  facet_wrap(~chr, scales="free_x", ncol=22) +
  geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
  geom_point(pch=20, cex=0.4, stroke=0) +
  #geom_point(pch=19,cex=0.2,alpha=0.05) +
  geom_hline(yintercept = 0, col="red",lwd=theme_cellpress_lwd,lty=1) +
  geom_smooth(se=F, lwd=theme_cellpress_lwd * 2, col="black") +
  coord_cartesian(ylim=c(-0.3, 0.3)) +
  labs(x=NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) +
  scale_color_manual(values=c(`TRUE`= mixcol('red', "gray60", 0.4), `FALSE`='gray60'), guide="none") +
  theme_nature
p8


p9 <- ggplot(plt.PC3.padj, aes(x=pos/100000000,y=minLogPadj, group=cnvp_bin)) +
  facet_wrap(~chr, scales="free_x", ncol=22) +
  geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
  geom_line(alpha=0.05,col=mixcol("darkgreen", "white", 0.15),lwd=theme_cellpress_lwd) +
  geom_smooth(data=subset(plt.PC3.padj, type == "data"), group=1, se=F, lwd=theme_cellpress_lwd * 2, col="black") +
  geom_hline(yintercept = 0, col="red",lwd=theme_cellpress_lwd,lty=1) +
  geom_hline(yintercept = -log(0.05), col="red",lwd=theme_cellpress_lwd,lty=2) +
  labs(y = "-log(Padj)") +
  coord_cartesian(ylim=c(0, 7.75)) +
  labs(x="Chromosomal position (/ 100 MB)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) +
  theme_nature
p9


p8 / p9



ggsave("output/figures/vis_cnv_profiles_GISTIC_style__PC3.png", width=11 * 0.975, height = 2.46) # removal of ticks needs to be incorp






# layer parimary recurrence ----


metadata.p.r <- metadata.all |> 
  filter_primaries_and_last_recurrences(136) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(patient_factor = as.factor(ifelse(dplyr::n() == 2, paste0("p",patient_id), "other"))) |> 
  dplyr::ungroup()


data.cnv.profiles.p.r <- data.cnv.profiles |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(metadata.p.r$array_sentrix_id) |> 
  t() |> 
  as.data.frame() 



## limma: resection_grade ----


design.p.r <- model.matrix(~patient_factor + pr.status, data=metadata.p.r |> tibble::column_to_rownames('array_sentrix_id'))


fit <- limma::lmFit(t(data.cnv.profiles.p.r), design.p.r)
fit <- limma::eBayes(fit)
stats.p.r <- limma::topTable(fit, n=nrow(fit), sort="p", coef="pr.statusrecurrence") |> # adjust="BH"
  tibble::rownames_to_column('cnvp_bin')  |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('','')) 




## per group quantiles ----


quantiles.p <- metadata.p.r |> 
  dplyr::filter(pr.status == "primary") |> 
  dplyr::select(array_sentrix_id) |> 
  dplyr::left_join(data.cnv.profiles.p.r |> tibble::rownames_to_column('array_sentrix_id'),
                   by=c('array_sentrix_id' = 'array_sentrix_id'), suffix=c('','')) |> 
  tibble::column_to_rownames('array_sentrix_id') |> 
  pbapply::pblapply(function(x) {return(quantile(x))}) |> 
  as.data.frame() |> 
  t() |> 
  as.data.frame() |> 
  dplyr::rename_with( ~ paste0("quantile.", .x)) |> 
  tibble::rownames_to_column('cnvp_bin')


quantiles.r <- metadata.p.r |> 
  dplyr::filter(pr.status == "recurrence") |> 
  dplyr::select(array_sentrix_id) |> 
  dplyr::left_join(data.cnv.profiles.p.r |> tibble::rownames_to_column('array_sentrix_id'),
                   by=c('array_sentrix_id' = 'array_sentrix_id'), suffix=c('','')) |> 
  tibble::column_to_rownames('array_sentrix_id') |> 
  pbapply::pblapply(function(x) {return(quantile(x))}) |> 
  as.data.frame() |> 
  t() |> 
  as.data.frame() |> 
  dplyr::rename_with( ~ paste0("quantile.", .x)) |> 
  tibble::rownames_to_column('cnvp_bin')


quantiles.p.r <- rbind(
  quantiles.p |> dplyr::mutate(pr.status = "Primary"),
  quantiles.r |> dplyr::mutate(pr.status = "Recurrence")
)








## plot ----



plt <- quantiles.p.r |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('','')) |> 
  dplyr::mutate(pr.status = factor(pr.status, levels=c("Primary","Recurrence")))


p6 <- ggplot(plt, aes(x=pos/100000000,y=`quantile.50%`, col=pr.status)) +
  geom_hline(yintercept = 0, col="gray40",lwd=theme_cellpress_lwd,lty=1, show_guide = FALSE) +
  facet_wrap(~chr, scales="free_x", ncol=22) +
  geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
  geom_point(pch=19,cex=0.2,alpha=0.05) +
  geom_smooth(se=F, method = "gam", lwd=theme_cellpress_lwd * 2.5) +
  coord_cartesian(ylim=c(-1.15, 1.15)) +
  labs(x=NULL, col=NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) +
  theme_cellpress + 
  scale_color_manual(values=palette_p_r)
p6





p7 <- ggplot(stats.p.r, aes(x=pos/100000000,y=logFC)) +
  facet_wrap(~chr, scales="free_x", ncol=22) +
  geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
  geom_point(pch=19,cex=0.2,alpha=0.05) +
  geom_hline(yintercept = 0, col="red",lwd=theme_cellpress_lwd) +
  geom_smooth(se=F, lwd=theme_cellpress_lwd * 2, col="black") +
  coord_cartesian(ylim=c(-0.3, 0.3)) +
  labs(x=NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) +
  theme_cellpress
p7







plt <- rbind(stats.p.r |> dplyr::mutate(y = -log(adj.P.Val), type="data"),
             stats.p.r |> dplyr::mutate(y= 0, type="border")
)


p8 <- ggplot(plt, aes(x=pos/100000000,y=y, group=cnvp_bin)) +
  facet_wrap(~chr, scales="free_x", ncol=22) +
  geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
  geom_line(alpha=0.05,col=mixcol("darkgreen", "white", 0.15),lwd=theme_cellpress_lwd) +
  theme_bw() +
  geom_smooth(data=subset(plt, type == "data"), group=1, se=F, lwd=theme_cellpress_lwd * 2, col="black") +
  geom_hline(yintercept = 0, col="red",lwd=theme_cellpress_lwd,lty=1) +
  geom_hline(yintercept = -log(0.05), col="red",lwd=theme_cellpress_lwd,lty=2) +
  labs(y = "-log(Padj)")+
  coord_cartesian(ylim=c(0, 6.25)) +
  labs(x="Chromosomal position (/ 100 MB)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) +
  theme_cellpress
p8



p6 / p7 / p8


ggsave("output/figures/vis_cnv_profiles_GISTIC_style__p_x_r.pdf", width=11 * 0.975, height = 3.75) # removal of ticks needs to be incorp



# CoxPH ----

## post rec svvl ----


metadata <- metadata.all |> 
  dplyr::filter(resection_number > 1) |> 
  dplyr::group_by(patient_id) |>
  dplyr::slice_max(resection_number, with_ties = FALSE) |>
  dplyr::ungroup() |> 
  assertr::verify(!duplicated(patient_id)) |> 
  dplyr::filter(!is.na(patient_last_follow_up_event)& !is.na(time_between_resection_and_last_event))



data.cnv <- data.cnv.profiles |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  t() |> 
  as.data.frame() 



stopifnot(metadata$array_sentrix_id == rownames(data.cnv))

svvl <- survival::Surv(metadata$time_between_resection_and_last_event, metadata$patient_last_follow_up_event)


out <- data.frame()
for(bin in 1:ncol(data.cnv)) {
  
  name <- colnames(data.cnv)[bin]
  dat <- data.frame(data = data.cnv[,bin],
                    array_methylation_bins_1p19q_purity = metadata$array_methylation_bins_1p19q_purity
                    )
  

  cox_model <- survival::coxph(svvl ~ data +array_methylation_bins_1p19q_purity, data = dat)
  cox_model_s <- summary(cox_model)
  
  append <- cox_model_s$coefficients |> 
    as.data.frame() |>
    tibble::rownames_to_column('factor') |>
    dplyr::filter(factor == "data") |> 
    dplyr::mutate(cnvp_bin=name)
  
  out <- rbind(out, append)
}



plt <- out |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('','')) |> 
  dplyr::mutate(adj.P.Val = p.adjust(`Pr(>|z|)`, method = "BH")) |> 
  dplyr::mutate(minLogPadj = -log(adj.P.Val))



plt.padj <- rbind(
  plt |> dplyr::mutate(type = "data"),
  plt |> dplyr::mutate(type = "border", minLogPadj = 0)) |> 
  dplyr::mutate(col = dplyr::case_when(
    adj.P.Val < 0.01 ~ "significant",
    adj.P.Val < 0.05 & adj.P.Val >= 0.01  ~ "trend",
    T ~ "other"
  ))



p10 <- ggplot(plt.padj, aes(x=pos/100000000,y=z, col=col)) +
  facet_wrap(~chr, scales="free_x", ncol=22) +
  geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
  #geom_point(pch=19,cex=0.2,alpha=0.15) +
  geom_point(pch=20,cex=0.4, stroke=0) +
  geom_hline(yintercept = 0, col="red",lwd=theme_cellpress_lwd,lty=1) +
  geom_smooth(se=F, lwd=theme_cellpress_lwd * 2, col="black") +
  coord_cartesian(ylim=c(-5, 5)) +
  labs(x=NULL,y="z-score CoxPH") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) +
  scale_color_manual(values=c(`significant`= 'red', 
                              `trend` = 'deeppink',
                              `other`='gray60'), guide="none") +
  theme_nature
p10


p11 <- ggplot(plt.padj, aes(x=pos/100000000,y=minLogPadj, group=cnvp_bin)) +
  facet_wrap(~chr, scales="free_x", ncol=22) +
  geom_vline(xintercept = 0, col="black",lwd=theme_cellpress_lwd,lty=1) +
  geom_line(alpha=0.05,col=mixcol("darkgreen", "white", 0.15),lwd=theme_cellpress_lwd) +
  geom_smooth(data=subset(plt.padj, type == "data"), group=1, se=F, lwd=theme_cellpress_lwd * 2, col="black") +
  geom_hline(yintercept = 0, col="red",lwd=theme_cellpress_lwd,lty=1) +
  geom_hline(yintercept = -log(0.05), col="red",lwd=theme_cellpress_lwd,lty=2) +
  labs(y = "-log(Padj)") +
  coord_cartesian(ylim=c(0, 3.75)) +
  labs(x="Chromosomal position (/ 100 MB)") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c(0,"",1,"",2,"",3), limits=c(0, NA)) +
  theme_nature
p11


p10 / p11


ggsave("output/figures/vis_cnv_profiles_GISTIC_style__post_rec_svvl.png", width=11 * 0.975, height = 2.46) # removal of ticks needs to be incorp





# stack ----


p1 / p2 / p3 / p4 / p5 / p6 / p7 / p8


# sandbox ----

#ggplot(plt.quantiles.g2[1:160,], aes(x=pos/1000000, y=bin_lfc_normalised)) +
#  geom_rect(aes(NULL,NULL,xmin=start,xmax=end,ymin=`quantile.25%`,ymax=`quantile.75%`), size=0.5, alpha=0.2)




# plt.quantiles.g3 <-  data.cnv.profiles |> 
#   tibble::rownames_to_column('array_sentrix_id') |> 
#   dplyr::filter(sentrix_id %in% (metadata |> dplyr::filter(resection_tumor_grade == 2) |> dplyr::pull(array_sentrix_id))) |> 
#   tibble::column_to_rownames('array_sentrix_id') |> 
#   pbapply::pblapply(function(x) {return(quantile(x))}) |> 
#   as.data.frame() |> 
#   t() |> 
#   as.data.frame() |> 
#   dplyr::rename_with( ~ paste0("quantile.", .x)) |> 
#   tibble::rownames_to_column('cnvp_bin') |> 
#   tidyr::pivot_longer(cols = -c('cnvp_bin'), names_to = 'type', values_to='bin_lfc_normalised') |> 
#   dplyr::filter(type %in% c('quantile.25%','quantile.25%'))







