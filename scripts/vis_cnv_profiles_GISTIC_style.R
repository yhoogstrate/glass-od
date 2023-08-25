#!/usr/bin/env


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}

source('scripts/load_chrom_sizes.R')


library(ggplot2)
library(patchwork)



# obtain data ----


metadata.all <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(163) |> 
  filter_primaries_and_last_recurrences(136) |> 
  dplyr::mutate(resection_tumor_hg = ifelse(resection_tumor_grade == 3, 1, 0)) # 1 or 0 for regression multiplication factor



#'@todo remove interim samples, second g2's and non-last g3's etc.
#'@todo first grading needs to be complete


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
  dplyr::select(sentrix_id, heidelberg_cnvp_bins, methylation_bins_1p19q_purity) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = f(heidelberg_cnvp_bins, methylation_bins_1p19q_purity)) |> 
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(heidelberg_cnvp_bins = NULL) |> 
  dplyr::mutate(methylation_bins_1p19q_purity = NULL) |> 
  tibble::column_to_rownames('sentrix_id')


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
  dplyr::filter(!is.na(resection_tumor_grade))


data.cnv.profiles.g3.g2 <- data.cnv.profiles |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(metadata.g3.g2$sentrix_id) |> 
  t() |> 
  as.data.frame() 


## limma: resection_grade ----


design.g3.g2 <- model.matrix(~resection_tumor_hg, data=metadata.g3.g2 |> tibble::column_to_rownames('sentrix_id'))


fit <- limma::lmFit(t(data.cnv.profiles.g3.g2), design.g3.g2)
fit <- limma::eBayes(fit)
stats.g3.g2 <- limma::topTable(fit, n=nrow(fit), sort="p") |> # adjust="BH"
  tibble::rownames_to_column('cnvp_bin')   |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('','')) 



## per group quantiles ----


quantiles.g2 <- metadata.g3.g2 |> 
  dplyr::filter(resection_tumor_grade == 2) |> 
  dplyr::select(sentrix_id) |> 
  dplyr::left_join(data.cnv.profiles.g3.g2 |> tibble::rownames_to_column('sentrix_id'),
                   by=c('sentrix_id' = 'sentrix_id'), suffix=c('','')) |> 
  tibble::column_to_rownames('sentrix_id') |> 
  pbapply::pblapply(function(x) {return(quantile(x))}) |> 
  as.data.frame() |> 
  t() |> 
  as.data.frame() |> 
  dplyr::rename_with( ~ paste0("quantile.", .x)) |> 
  tibble::rownames_to_column('cnvp_bin')


quantiles.g3 <- metadata.g3.g2 |> 
  dplyr::filter(resection_tumor_grade == 3) |> 
  dplyr::select(sentrix_id) |> 
  dplyr::left_join(data.cnv.profiles.g3.g2 |> tibble::rownames_to_column('sentrix_id'),
                   by=c('sentrix_id' = 'sentrix_id'), suffix=c('','')) |> 
  tibble::column_to_rownames('sentrix_id') |> 
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
  dplyr::mutate(grade = as.factor(paste0("Grade ",grade))) 



p1 <- ggplot(plt, aes(x=pos/1000000,y=`quantile.50%`, col=grade)) +
  facet_grid(cols = vars(chr), scales = "free", space="free")  +
  geom_point(pch=19,cex=0.2,alpha=0.05) +
  theme_bw() +
  geom_smooth() +
  coord_cartesian(ylim=c(-1.25, 1.25)) +
  labs(x=NULL) +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300)) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.1, "lines"))




p2 <- ggplot(stats.g3.g2, aes(x=pos/1000000,y=logFC)) +
  facet_grid(cols = vars(chr), scales = "free", space="free")  +
  geom_point(pch=19,cex=0.2,alpha=0.05) +
  theme_bw() +
  geom_hline(yintercept = 0, col="red",lwd=0.25,lty=1) +
  geom_smooth(se=F, lwd=0.5, col="black") +
  coord_cartesian(ylim=c(-0.3, 0.3)) +
  labs(x=NULL) +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300)) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(panel.spacing = unit(0.1, "lines"))
p2



plt <- rbind(stats.g3.g2 |> dplyr::mutate(y = -log(adj.P.Val), type="data"),
             stats.g3.g2 |> dplyr::mutate(y= 0, type="border")
             )



p3 <- ggplot(plt, aes(x=pos/1000000,y=y, group=cnvp_bin)) +
  facet_grid(cols = vars(chr), scales = "free", space="free")  +
  geom_line(alpha=0.05,col="darkgreen",lwd=0.25) +
  theme_bw() +
  geom_smooth(data=subset(plt, type == "data"), group=1, se=F, lwd=0.5, col="black") +
  geom_hline(yintercept = 0, col="red",lwd=0.25,lty=1) +
  geom_hline(yintercept = -log(0.05), col="red",lwd=0.25,lty=2) +
  labs(y = "-log(Padj)")+
  coord_cartesian(ylim=c(0, 5)) +
  labs(x=NULL) +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300)) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(panel.spacing = unit(0.1, "lines"))
p3


p1 / p2 / p3



# layer A_IDH_HG__A_IDH_LG_lr__lasso_fit + PC2 ----


metadata.LR.PC2 <- metadata.all |> 
  dplyr::filter(!is.na(A_IDH_HG__A_IDH_LG_lr__lasso_fit) & !is.na(PC2))

data.cnv.profiles.LR.PC2 <- data.cnv.profiles |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(metadata.LR.PC2$sentrix_id) |> 
  t() |> 
  as.data.frame() 



design.LR <- model.matrix(~scale(A_IDH_HG__A_IDH_LG_lr__lasso_fit), data=metadata.LR.PC2 |> tibble::column_to_rownames('sentrix_id'))
design.PC2 <- model.matrix(~ scale(-PC2), data=metadata.LR.PC2 |> tibble::column_to_rownames('sentrix_id'))


fit.LR <- limma::lmFit(t(data.cnv.profiles.LR.PC2), design.LR)
fit.LR <- limma::eBayes(fit.LR)
stats.LR <- limma::topTable(fit.LR, n=nrow(fit.LR), sort="p") |> # adjust="BH"
  tibble::rownames_to_column('cnvp_bin')

fit.PC2 <- limma::lmFit(t(data.cnv.profiles.LR.PC2), design.PC2)
fit.PC2 <- limma::eBayes(fit.PC2)
stats.PC2 <- limma::topTable(fit.PC2, n=nrow(fit.PC2), sort="p") |> # adjust="BH"
  tibble::rownames_to_column('cnvp_bin')


## plot ----


plt.LR <- stats.LR |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('',''))


plt.LR.padj <- stats.LR |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('','')) |> 
  dplyr::mutate(minLogPadj = -log(adj.P.Val))

plt.LR.padj <- rbind(
  plt.LR.padj |> dplyr::mutate(type = "data"),
  plt.LR.padj |> dplyr::mutate(type = "border", minLogPadj = 0))

# plt.PC2 <- stats.PC2 |> 
#   dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('','')) |> 
#   dplyr::mutate(minLogPadj = -log(adj.P.Val))


p4 <- ggplot(plt.LR, aes(x=pos/1000000,y=logFC)) +
  facet_grid(cols = vars(chr), scales = "free", space="free")  +
  geom_point(pch=19,cex=0.2,alpha=0.05) +
  theme_bw() +
  geom_hline(yintercept = 0, col="red",lwd=0.25,lty=1) +
  geom_smooth(se=F, lwd=0.5, col="black") +
  coord_cartesian(ylim=c(-0.3, 0.3)) +
  labs(x=NULL) +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300)) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(panel.spacing = unit(0.1, "lines"))


p5 <- ggplot(plt.LR.padj, aes(x=pos/1000000,y=minLogPadj, group=cnvp_bin)) +
  facet_grid(cols = vars(chr), scales = "free", space="free")  +
  geom_line(alpha=0.05, col="darkgreen",lwd=0.25) +
  theme_bw() +
  geom_smooth(data=subset(plt.LR.padj, type == "data"), group=1, se=F, lwd=0.5, col="black") +
  geom_hline(yintercept = 0, col="red",lwd=0.25,lty=1) +
  geom_hline(yintercept = -log(0.05), col="red",lwd=0.25,lty=2) +
  labs(y = "-log(Padj)") +
  coord_cartesian(ylim=c(-0, 12)) +
  labs(x=NULL) +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300)) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(panel.spacing = unit(0.1, "lines"))


p4 / p5




# p2 <- ggplot(plt.PC2, aes(x=pos/1000000,y=logFC)) +
#   facet_grid(cols = vars(chr), scales = "free", space="free")  +
#   geom_point(pch=19,cex=0.2,alpha=0.05) +
#   theme_bw() +
#   geom_hline(yintercept = 0, col="red",lwd=0.25,lty=1) +
#   geom_smooth(se=F, lwd=0.5, col="black") +
#   coord_cartesian(ylim=c(-0.25, 0.25))



# p4 <- ggplot(plt.PC2, aes(x=pos/1000000,y=minLogPadj)) +
#   facet_grid(cols = vars(chr), scales = "free", space="free")  +
#   geom_point(pch=19,cex=0.2,alpha=0.05) +
#   theme_bw() +
#   geom_smooth(se=F, lwd=0.5, col="black") +
#   geom_hline(yintercept = 0, col="red",lwd=0.25,lty=1) +
#   geom_hline(yintercept = -log(0.01), col="red",lwd=0.25,lty=2) +
#   coord_cartesian(ylim=c(-0, 12))


# layer parimary recurrence ----

metadata.p.r <- metadata.all |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence")))

data.cnv.profiles.p.r <- data.cnv.profiles |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(metadata.p.r$sentrix_id) |> 
  t() |> 
  as.data.frame() 



## limma: resection_grade ----


design.p.r <- model.matrix(~pr.status, data=metadata.p.r |> tibble::column_to_rownames('sentrix_id'))


fit <- limma::lmFit(t(data.cnv.profiles.p.r), design.p.r)
fit <- limma::eBayes(fit)
stats.p.r <- limma::topTable(fit, n=nrow(fit), sort="p") |> # adjust="BH"
  tibble::rownames_to_column('cnvp_bin')  |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('','')) 




## per group quantiles ----


quantiles.p <- metadata.p.r |> 
  dplyr::filter(pr.status == "primary") |> 
  dplyr::select(sentrix_id) |> 
  dplyr::left_join(data.cnv.profiles.p.r |> tibble::rownames_to_column('sentrix_id'),
                   by=c('sentrix_id' = 'sentrix_id'), suffix=c('','')) |> 
  tibble::column_to_rownames('sentrix_id') |> 
  pbapply::pblapply(function(x) {return(quantile(x))}) |> 
  as.data.frame() |> 
  t() |> 
  as.data.frame() |> 
  dplyr::rename_with( ~ paste0("quantile.", .x)) |> 
  tibble::rownames_to_column('cnvp_bin')


quantiles.r <- metadata.p.r |> 
  dplyr::filter(pr.status == "recurrence") |> 
  dplyr::select(sentrix_id) |> 
  dplyr::left_join(data.cnv.profiles.p.r |> tibble::rownames_to_column('sentrix_id'),
                   by=c('sentrix_id' = 'sentrix_id'), suffix=c('','')) |> 
  tibble::column_to_rownames('sentrix_id') |> 
  pbapply::pblapply(function(x) {return(quantile(x))}) |> 
  as.data.frame() |> 
  t() |> 
  as.data.frame() |> 
  dplyr::rename_with( ~ paste0("quantile.", .x)) |> 
  tibble::rownames_to_column('cnvp_bin')


quantiles.p.r <- rbind(
  quantiles.p |> dplyr::mutate(pr.status = "primary"),
  quantiles.r |> dplyr::mutate(pr.status = "recurrence")
)








## plot ----



plt <- quantiles.p.r |> 
  dplyr::left_join(metadata.cnv.profiles, by=('cnvp_bin'='cnvp_bin'), suffix=c('','')) |> 
  dplyr::mutate(pr.status = factor(pr.status, levels=c("primary","recurrence")))


p6 <- ggplot(plt, aes(x=pos/1000000,y=`quantile.50%`, col=pr.status)) +
  facet_grid(cols = vars(chr), scales = "free", space="free")  +
  geom_point(pch=19,cex=0.2,alpha=0.05) +
  theme_bw() +
  geom_smooth() +
  coord_cartesian(ylim=c(-1.25, 1.25)) +
  labs(x=NULL) +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300)) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(panel.spacing = unit(0.1, "lines"))
p6




p7 <- ggplot(stats.p.r, aes(x=pos/1000000,y=logFC)) +
  facet_grid(cols = vars(chr), scales = "free", space="free")  +
  geom_point(pch=19,cex=0.2,alpha=0.05) +
  theme_bw() +
  geom_hline(yintercept = 0, col="red",lwd=0.25,lty=1) +
  geom_smooth(se=F, lwd=0.5, col="black") +
  coord_cartesian(ylim=c(-0.3, 0.3)) +
  labs(x=NULL) +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300)) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(panel.spacing = unit(0.1, "lines"))
p7




plt <- rbind(stats.p.r |> dplyr::mutate(y = -log(adj.P.Val), type="data"),
             stats.p.r |> dplyr::mutate(y= 0, type="border")
)


p8 <- ggplot(plt, aes(x=pos/1000000,y=y, group=cnvp_bin)) +
  facet_grid(cols = vars(chr), scales = "free", space="free")  +
  geom_line(alpha=0.033,col="darkgreen",lwd=0.25) +
  theme_bw() +
  geom_smooth(data=subset(plt, type == "data"), group=1, se=F, lwd=0.5, col="black") +
  geom_hline(yintercept = 0, col="red",lwd=0.25,lty=1) +
  geom_hline(yintercept = -log(0.01), col="red",lwd=0.25,lty=2) +
  labs(y = "-log(Padj)", x=NULL) +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300)) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(panel.spacing = unit(0.1, "lines"))
p8


p6 / p7 / p8



# stack ----


p1 / p2 / p3 / p4 / p5 / p6 / p7 / p8


# sandbox ----

#ggplot(plt.quantiles.g2[1:160,], aes(x=pos/1000000, y=bin_lfc_normalised)) +
#  geom_rect(aes(NULL,NULL,xmin=start,xmax=end,ymin=`quantile.25%`,ymax=`quantile.75%`), size=0.5, alpha=0.2)




# plt.quantiles.g3 <-  data.cnv.profiles |> 
#   tibble::rownames_to_column('sentrix_id') |> 
#   dplyr::filter(sentrix_id %in% (metadata |> dplyr::filter(resection_tumor_grade == 2) |> dplyr::pull(sentrix_id))) |> 
#   tibble::column_to_rownames('sentrix_id') |> 
#   pbapply::pblapply(function(x) {return(quantile(x))}) |> 
#   as.data.frame() |> 
#   t() |> 
#   as.data.frame() |> 
#   dplyr::rename_with( ~ paste0("quantile.", .x)) |> 
#   tibble::rownames_to_column('cnvp_bin') |> 
#   tidyr::pivot_longer(cols = -c('cnvp_bin'), names_to = 'type', values_to='bin_lfc_normalised') |> 
#   dplyr::filter(type %in% c('quantile.25%','quantile.25%'))







