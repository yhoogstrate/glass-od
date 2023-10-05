#!/usr/bin/env


library(ggplot2)
library(patchwork)


source("scripts/load_functions.R")
source("scripts/load_themes.R")
source("scripts/load_palette.R")


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# make as-much-as-possible covering CNV plot per sample ----

#' @todo: detP x usv-PC2 + _reason_excl
#' @todo: PC2 x AcCGAP
#' @todo: thumbnail of H&E

cnv_plot <- function(cur_sentrix_id) {
  cur_sentrix_id <- "207331540060_R01C01"
  print(cur_sentrix_id)
  #cur_sentrix_id <- "201496850071_R02C01"
  #cur_sentrix_id <- "201496850071_R02C01"
  #a = c("201496850071_R02C01", "203293640061_R08C01", "205828590003_R02C01", "205832320037_R07C01",
  #  "206238130171_R03C01", "207331540060_R02C01", "207331540060_R03C01",  "207331540060_R04C01")
  #cur_sentrix_id = a[6]
  
  cur_idat <- glass_od.metadata.array_samples |> 
    dplyr::filter(array_sentrix_id == cur_sentrix_id)
  

  ## p1: CNV plot ----
  bins <- read.delim(cur_idat$array_mnp_CNVP_v12.8_v5.2_CNVP_bins)
  segments <- read.delim(cur_idat$array_mnp_CNVP_v12.8_v5.2_CNVP_segments)
  
  plt.bins <- bins |> 
    dplyr::rename(log2fc = 5) |> 
    dplyr::mutate(Feature = NULL) |> 
    dplyr::rename(chr = Chromosome, start = Start, end = End) |> 
    dplyr::filter(end - start > 35) |> 
    dplyr::mutate(type="bin")
  
  plt.segments <- segments |> 
    dplyr::select(chrom, loc.start, loc.end, seg.median) |> 
    dplyr::rename(chr = chrom, start = loc.start, end = loc.end, log2fc=seg.median) |> 
    dplyr::mutate(type="segment")
  
  plt <- rbind(plt.bins, plt.segments) |>
    dplyr::mutate(group = paste0("line-",1:dplyr::n())) |> 
    dplyr::filter(log2fc >= -1.2 & log2fc <= 1.7) |> 
    tidyr::pivot_longer(cols = c(start, end), values_to = "pos",names_to = "segment_point") |> 
    dplyr::filter(chr %in% c("chrX","chrY","chrM") == F ) |> 
    dplyr::mutate(chr = factor(chr, levels = gtools::mixedsort(unique(as.character(chr))))) |> 
    dplyr::mutate(pos = round(pos/1000000))
  
  sel <- unique(c(cur_idat$array_mnp_predictBrain_v12.8_cal_class, "O_IDH", "A_IDH_LG", "A_IDH_HG", "OLIGOSARC_IDH"))
  pct <- cur_idat |> dplyr::select(paste0("array_mnp_predictBrain_v12.8_cal_",sel))
  class_txt <- data.frame(sel=sel, pct) |> 
    dplyr::mutate(txt = paste0(sel, ": ",round(pct,2),"")) |> 
    dplyr::pull(txt) |> 
    stringi::stri_paste(collapse='  -  ')
  
  p1 <- ggplot(plt, aes(x = pos, y=log2fc, group=group, col=log2fc)) +
    facet_grid(cols = vars(chr), scales = "free", space="free") +
    
    geom_vline(xintercept = 0, col="black", lwd=theme_cellpress_lwd) +
    
    geom_hline(yintercept = 0, col="gray80", lty=2, lwd=theme_cellpress_lwd) +
    geom_hline(yintercept = -1, col="gray80", lty=2, lwd=theme_cellpress_lwd) +
    geom_hline(yintercept = 1.5, col="gray80", lty=2, lwd=theme_cellpress_lwd) +

    geom_hline(yintercept = cur_idat$array_methylation_bins_1p19q_median.lfc, col="blue", lty=3, lwd=theme_cellpress_lwd) +
    
    geom_line(data = plt |> dplyr::filter(type =="bin"), lwd=1.25) +
    geom_line(data = plt |> dplyr::filter(type =="segment"), lwd=theme_cellpress_lwd * 2, col="red") +

    theme(
      axis.title = element_text(face = "bold", size = rel(1)),
      # axis.text.x = element_blank(),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.45, hjust = 1),
      
      panel.spacing = unit(0.1, "lines")
    ) +
    labs(
      x = NULL,
      y = paste0("CNVP v12.8 log2 ratio"),
      subtitle = paste0(cur_idat$resection_id, " / ", 
                       cur_idat$array_sentrix_id, "  -  ",
                       class_txt , 
                       "  -  purity: ", round(100 * cur_idat$array_methylation_bins_1p19q_purity,1), "%", 
                       "  -  grade: ", as.character(cur_idat$resection_tumor_grade))
    ) +
    ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", 
                                   limits = c(-abs(cur_idat$array_methylation_bins_1p19q_median.lfc), abs(cur_idat$array_methylation_bins_1p19q_median.lfc)),
                                   breaks= c( -abs(cur_idat$array_methylation_bins_1p19q_median.lfc) , abs(cur_idat$array_methylation_bins_1p19q_median.lfc)), oob = scales::squish) +
    scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300)) +
    coord_cartesian(ylim = c(-1.2, 1.7)) +
    theme_cellpress
  p1
  
  
  ## p2: purity x MNP class plot ----
  plt <- glass_od.metadata.array_samples |> 
    dplyr::filter(!is.na(array_methylation_bins_1p19q_median.lfc))|> 
    dplyr::filter(!is.na(array_methylation_bins_1p19q_sd))
  
  p2 <- ggplot(plt, aes(x=array_methylation_bins_1p19q_purity,
                  y=array_methylation_bins_1p19q_sd,
                  col=array_mnp_predictBrain_v12.8_cal_class,
                  label=isolation_id)
  ) + 
    geom_vline(xintercept = 0.1, col="red",lwd=theme_cellpress_lwd,lty=2)  + 
    geom_point(size=1) +
    geom_point(data = subset(plt, array_sentrix_id == cur_sentrix_id), size=2.2, col="black", pch=21) +
    ggrepel::geom_text_repel(data = cur_idat, col="black",
                             nudge_x = 0.1*0.4,
                             nudge_y = 0.05*1.2,
                             size=theme_cellpress_size) +
    labs(x = "Tumor purity [1P/19Q log2Fc based]",
         y="Standard deviation purity estimator") +
    coord_cartesian(ylim = c(0, 0.5)) +
    labs(col="v12.8 class") +
    theme_cellpress
  p2

  
  
  ## p3: QC plot ----
  plt <- glass_od.metadata.array_samples |> 
    dplyr::mutate(col = dplyr::case_when(
      array_sentrix_id == cur_sentrix_id ~ isolation_id,
      array_qc.pca.detP.outlier == T ~ "outlier: excluded",
      array_qc.pca.detP.outlier == F ~ "no outlier: good / sufficient quality"
    ))
  
  p3 <- ggplot(plt, aes(x=array_qc.pca.comp1, y=array_percentage.detP.signi, col = col, label=isolation_id)) +
    scale_y_continuous(trans = "log1p", limits=c(0,100), breaks=c(0,2.5,5,10,25,50,100)) +
    scale_x_continuous(breaks=c(-400,0,400,600,800,1200)) +
    geom_hline(yintercept=2.5, col="red", lty=2, lwd=theme_cellpress_lwd) +
    geom_vline(xintercept=600, col="red", lty=2, lwd=theme_cellpress_lwd) +
    geom_point(size=1) +
    geom_point(data = subset(plt, array_sentrix_id == cur_sentrix_id), size=2.2, col="black", pch=21) +
    ggrepel::geom_text_repel(data=subset(plt, array_sentrix_id == cur_sentrix_id), show_guide=F, 
                             col="black", size=theme_cellpress_size,nudge_y = 0.5) +
    labs(x = "PC1 all samples", y="percentage detP failed", col=NULL) +
    theme_cellpress
  p3
  
  
  
  ## p4: A/O/IDH/LG/HG ----
  
  sel <- unique(c("O_IDH", "A_IDH_LG", "A_IDH_HG", "OLIGOSARC_IDH"))
  plt <- cur_idat |> dplyr::select(paste0("array_mnp_predictBrain_v12.8_cal_",sel)) |> 
    t() |> 
    as.data.frame() |> 
    tibble::rownames_to_column('v12.8_predicted_class') |> 
    dplyr::mutate(v12.8_predicted_class = gsub("array_mnp_predictBrain_v12.8_cal_", "", v12.8_predicted_class)) |> 
    dplyr::mutate(tumor_type = ifelse(grepl("A_IDH",v12.8_predicted_class), "Astrocytoma", "Oligodendroglioma")) |> 
    dplyr::mutate(lg_hg = ifelse(v12.8_predicted_class %in% c("A_IDH","A_IDH_LG", "O_IDH"), "LG", "HG")) |> 
    dplyr::mutate(y_ref = ifelse(lg_hg == "LG", 0 , 1)) |> 
    dplyr::mutate(y_dat = ifelse(lg_hg == "HG", 1 - V1, V1)) |> 
    tidyr::pivot_longer(cols=c(y_ref,y_dat), values_to = 'y')
  
  p4 <- ggplot(plt, aes(x=tumor_type, y=y, group=v12.8_predicted_class, col=v12.8_predicted_class)) +
    geom_hline(yintercept = 1 , col="black", alpha=0.5, lwd=theme_cellpress_lwd) +
    geom_hline(yintercept = 0 , col="black", alpha=0.5, lwd=theme_cellpress_lwd) +
    geom_line(lwd=10, show_guide=F) +
    geom_point(data=subset(plt, y==9999)) +
    labs(x = NULL, y="Prediction probability MNP v12.8", col=NULL) +
    scale_color_manual(values = palette_mnp_12.8_6) +
    theme_cellpress
  p4
  
  
  ## p5: Reasons excluded plot ----
  plt <- cur_idat |> 
    dplyr::select(contains("excluded")) |> 
    t() |>
    as.data.frame() |> 
    tibble::rownames_to_column("reason") |> 
    dplyr::mutate(reason = gsub("_reason_excluded","",reason)) |> 
    dplyr::mutate(excluded = !is.na(V1)) |> 
    dplyr::mutate(V1 = ifelse(excluded, V1, "no reason for exclusion"))
  
  
  p5 <- ggplot(plt, aes(y = reason, label=V1, fill=excluded)) +
    geom_label(x = 0, hjust=0, size=theme_cellpress_size, show_guide=F) +
    scale_fill_manual(values=c('FALSE'=mixcol('green','darkgreen',0.4), 'TRUE'='red')) +
    labs(y=NULL, subtitle = paste0(cur_idat$isolation_id, ": reasons excluded"), fill=NULL, col=NULL, x=NULL) +
    theme_cellpress
  p5
  
  
  ## p6: H&E ------
  
  if(is.na(cur_idat$HE_coupe_thumbnail_roi)) {
    p6 <- grid::rasterGrob(width = 1.25, png::readPNG('data/gray-pix.png'), interpolate=TRUE)
    
  } else {
    p6 <- grid::rasterGrob(jpeg::readJPEG(cur_idat$HE_coupe_thumbnail_roi), interpolate=TRUE)
  }
  
  
  if(is.na(cur_idat$HE_coupe_thumbnail_full)) {
    p7 <- grid::rasterGrob(width = 1.25, png::readPNG('data/gray-pix.png'), interpolate=TRUE)
    
  } else {
    p7 <- grid::rasterGrob(png::readPNG(cur_idat$HE_coupe_thumbnail_full), interpolate=TRUE)
  }
  
  
  
  ## merge ----
  
  p1 / (p2 + p3 + p4) / (p5 + p6 + p7)
  
  
  ggsave(paste0("output/figures/cnv_profiles/",
                cur_idat$isolation_id,"_badge.png"), width = 8.5 * 0.975, height = 5.25, scale=1.5)
}


pbapply::pblapply(glass_od.metadata.array_samples |> 
                    dplyr::filter(!is.na(array_heidelberg_cnvp_bins)) |> 
                    dplyr::filter(!is.na(array_mnp_predictBrain_v12.8_cal_class)) |> 
                    
                    #dplyr::filter(is.na(reason_excluded_patient)) |> 
                    #dplyr::filter(is.na(reason_excluded_resection_isolation)) |> 
                    #dplyr::filter(is.na(reason_excluded_resection)) |> 
                    #dplyr::filter(is.na(reason_excluded_array_sample)) |> 
                    dplyr::pull(array_sentrix_id), cnv_plot)





tpc.estimate = data.frame()
for(bc in  glass_od.metadata.idat |> dplyr::filter(!is.na(heidelberg_CNV_segments)) |> dplyr::pull(Sentrix_ID) ) {
  
  seg <- glass_od.metadata.array_samples |> 
    dplyr::filter(`sentrix_id` == sentrix_id)
  
  
  segs <- read.delim(seg |> 
                        dplyr::pull(heidelberg_CNV_segments)) |> 
    dplyr::rename(seg.median.l2fc = seg.median) |> 
    dplyr::rename(seg.mean.l2fc = seg.mean) |> 
    dplyr::filter(num.mark >= 35)
  
  center <-  segs |> 
    dplyr::filter(is.na(pval))
  center <- mean(rep(center$seg.mean.l2fc, center$num.mark))
  
  a = segs |> 
    dplyr::filter(chrom %in% c('chrX','chrY') == F) |> 
    dplyr::mutate(chrom.offset = chrs_hg19_s[chrom] ) |> 
    dplyr::mutate(x1 = loc.start + chrom.offset, x2 = loc.end + chrom.offset) |> 
    
    #dplyr::mutate(outlier = num_points <= 60 ) |> 
    #dplyr::filter(outlier == F) |> 
    
    #dplyr::mutate(seg.median.l2fc = seg.median.l2fc - center) |> 
    dplyr::filter(seg.median.l2fc < 1.1) |> 
    dplyr::filter(seg.median.l2fc > -2.0)
  
  
  
 
  
  
  # chr17 84622 190.863.195
  if(T) {
    a <- a |> dplyr::filter(chrom %in% c('chr1', 'chr19'))
  }
  # } else if(bc == "206119350033_R08C01") { a <- a |> dplyr::filter(chrom %in% c('chr1', 'chr19'))
  # } else if(bc == "206137490109_R08C01") { a <- a |> dplyr::filter(chrom %in% c('chr1', 'chr19'))
  # } else if(bc == "206137490109_R04C01") { a <- a |> dplyr::filter(chrom %in% c('chr1', 'chr19'))
  # } else if(bc == "206137490109_R02C01") { a <- a |> dplyr::filter(chrom %in% c('chr1', 'chr19'))
  # } else if(bc == "206137490057_R08C01") { a <- a |> dplyr::filter(chrom %in% c('chr1', 'chr19'))

  
  
  
  #print(paste0(bc, " -> ", str_c(unique(a$chrom), collapse = ",")))
  
  
  # code
  out <- data.frame()
  for(frac in  1:100 / 100) {
    fc.p.4 <- ((1 - frac) * 2 + frac * 4) / 2
    fc.p.3 <- ((1 - frac) * 2 + frac * 3) / 2
    fc.n.1 <- ((1 - frac) * 2 + frac * 1) / 2
    
    lfc.p.4 <- log2(fc.p.4)
    lfc.p.3 <- log2(fc.p.3)
    lfc.n.1 <- log2(fc.n.1)
    
    dists <- c()
    dist <- 0
    for(i in 1:nrow(a)) {
      
      e <- a[i,]
      
      d <- c(
        (((e$seg.median.l2fc - 0) * e$num.mark)^2),
        (((e$seg.median.l2fc - lfc.p.4) * e$num.mark)^2) * 1.1 * 999, # penalize, do never prefer
        (((e$seg.median.l2fc - lfc.p.3) * e$num.mark)^2),
        (((e$seg.median.l2fc - lfc.n.1) * e$num.mark)^2)
      ) 
      d <- min(d)
      dists <- c(dists, d)
      
      dist <- dist + d
    }
    
    out <- rbind(out, data.frame(pct = frac * 100,
                                 lfc.3p = lfc.p.3,
                                 lfc.4p = lfc.p.4,
                                 lfc.n = lfc.n.1,
                                 dist = dist) )
  }
  
  
  r = out |> 
    dplyr::arrange(dist) |> 
    dplyr::slice(1) |> 
    dplyr::mutate(
                  Sentrix_ID = bc
                  #sample.short = gsub("^([^-]+-[^-]+-[^-]+-[^-]+).+$","\\1",bc),
                  #portion_barcode = gsub("^([^\\-]+)-([^\\-]+)-([^\\-]+)-([^\\-]+)-([0-9]+).+$$","\\1-\\2-\\3-\\4-\\5",bc),
                  #estimate.purity = p.tpc.estimate,
                  #hoogstrate.rf.purity.2021 = p.tpc.hoogstrate.2021
    )
  tpc.estimate <- rbind(tpc.estimate, r)
  

  
  
  plt <- segs |> 
    dplyr::mutate(chrom = as.character(chrom)) |> 
    dplyr::filter(chrom %in% c('chrX','chrY') == F) |> 
    dplyr::mutate(id = paste0("id.", 1:dplyr::n())) |>
    tidyr::pivot_longer(cols = c(`loc.start`, `loc.end`), names_to = "type", values_to = "pos") |>
    dplyr::mutate(chrom = factor(chrom, levels = gtools::mixedsort(unique(as.character(chrom))))) |>
    dplyr::arrange(chrom) |> 
    dplyr::mutate(pos = pos / 1000000) |> 
    dplyr::filter(num.mark >= 35) |> 
    dplyr::mutate(pos2=NA, seg.median.l2fc2=NA)

  
  plt <- rbind(plt,
               data.frame(
                 pos = c(125038530, 125038530,   26546521,26546521)/1000000,
                 seg.median.l2fc = c(-2, 2,   -2,2),
                 id=c("1p/q marker","1p/q marker", "19p/q marker","19p/q marker"),
                 chrom=c("chr1","chr1","chr19","chr19"),
                 ID=NA, num.mark=NA, bstat=NA,pval=NA, seg.mean.l2fc=NA,pos2=NA,seg.median.l2fc2=NA,type="vline")
          ,
          data.frame(
                pos =  c(0, 0)/1000000,
                pos2 = c(chrs_hg19["chr1"], chrs_hg19["chr19"])/1000000,
                
                seg.median.l2fc = c(-2, -2),
                seg.median.l2fc2 = c(2, 2),
                
                id=c("1p/q box","19p/q box"),
                chrom=c("chr1","chr19"),
                ID=NA, num.mark=NA, bstat=NA,pval=NA, seg.mean.l2fc=NA,type="box")
               )
  
  
  sel <- unique(c(seg$predictBrain_11_scores_cal_class, "O_IDH", "A_IDH", "A_IDH_HG"))
  pct <- seg |> dplyr::select(paste0("predictBrain_11_scores_cal_",sel))
  class_txt <- data.frame(sel=sel, pct) |> 
    dplyr::mutate(txt = paste0(sel, ": ",pct,"%")) |> 
    dplyr::pull(txt) |> 
    stringi::stri_paste(collapse='  -  ')
  
  
  
  purity <- r$pct / 100
  ggplot(plt |> dplyr::filter(type %in% c("loc.start", "loc.end" )),
         aes(x = pos, y = seg.median.l2fc, group = id, col = chrom)) +
    facet_grid(cols=vars(chrom), scales = "free", space = "free") +
    geom_rect(aes(xmin=pos, xmax=pos2, ymin=seg.median.l2fc, ymax=seg.median.l2fc2), 
              data =  plt |> dplyr::filter(type == "box"), col=NA, alpha=0.04,fill="green") +
    
    geom_hline(yintercept = 0, lwd = 0.5, lty = 3, col = "black", alpha = 0.5) +
    geom_line(lwd = 2) +
    geom_hline(yintercept = log2(((1 - purity) * 2 + purity * 4) / 2), lwd = 0.7, lty = 2, col = "black", alpha = 0.35) +
    geom_hline(yintercept = log2(((1 - purity) * 2 + purity * 3) / 2), lwd = 0.7, lty = 2, col = "black", alpha = 0.35) +
    geom_hline(yintercept = log2(((1 - purity) * 2 + purity * 1) / 2), lwd = 0.7, lty = 2, col = "black", alpha = 0.35) +
    geom_line(data =  plt |> dplyr::filter(type == "vline"),lty=2,col="black") +
    
    theme_bw() +
    theme(
      axis.title = element_text(face = "bold", size = rel(1)),
      # axis.text.x = element_blank(),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.45, hjust = 1),
      
      panel.spacing = unit(0.1, "lines")
    ) +
    scale_color_discrete(guide = "none") +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    labs(
      x = NULL,
      y = paste0("CNVP ",seg$heidelberg_cnvp_version," log2 ratio"),
      caption = paste0(seg$GLASS_OD_patient, " / ", bc, "  -  ", class_txt , "  -  purity estimate: ", purity * 100, "%")
    ) +
    scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300))
  
  
  
  ggsave(paste0("output/figures/purity_cnv-epic/",seg$GLASS_OD_patient,"_", bc ,".tpc.estimate.png"), width = 8.3, height = 3.5, scale = 1.75)
  
}


# per grade CNV plot normalised ----



## test with simple small test set [3x g2, 3x g3] ----





