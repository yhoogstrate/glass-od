#!/usr/bin/env


source('scripts/load_chrom_sizes.R')


# CNV based purity estimation

# type 1: segment based [excl] ----

estimate_1p19q_purity_segment_based <- function(fn) {
  #fn = glass_od.metadata.array_samples$heidelberg_cnvp_segments[1]
  
  segs <- read.delim(fn) |> 
    dplyr::rename(seg.median.l2fc = seg.median) |> 
    dplyr::rename(seg.mean.l2fc = seg.mean) |> 
    dplyr::filter(num.mark >= 35)
  
  
  center <-  segs |> 
    dplyr::filter(is.na(pval))
  center <- mean(rep(center$seg.mean.l2fc, center$num.mark))
}








tpc.estimate = data.frame()
for(bc in (glass_od.metadata.array_samples |> 
           dplyr::filter(arraychip_version != "450k" ) |> 
           dplyr::filter(!is.na(array_mnp_CNVP_v12.8_v5.2_CNVP_segments)) |> 
           dplyr::pull(array_array_sentrix_id))) {
  bc = "204808700074_R05C01"
  print(bc)
  seg <- glass_od.metadata.idat |> 
    dplyr::filter(`array_sentrix_id` == bc)
  
  
  segs <- read.delim(seg |> 
                        dplyr::pull(array_mnp_CNVP_v12.8_v5.2_CNVP_segments)) |> 
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
                  array_sentrix_id = bc
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



# type 2: bin based ----

## 1. preliminary analysis finding important bins

fun <- function(fn) {
  df <- read.delim(fn) |> 
    dplyr::filter(Chromosome %in% c("chr1", "chr19")) |> 
    tibble::column_to_rownames('Feature') |> 
    dplyr::select(c(4)) |> 
    t() |> 
    as.data.frame(stringsAsFactors=F)

  return(df)
}



df <- do.call(
  rbind,
  pbapply::pblapply(glass_od.metadata.array_samples |> 
                      dplyr::filter(!is.na(array_mnp_CNVP_v12.8_v5.2_CNVP_bins)) |> 
                      #dplyr::filter(is.na(reason_excluded_array_sample)) |> 
                      #dplyr::filter(is.na(reason_excluded_patient)) |> 
                      #dplyr::filter(is.na(reason_excluded_resection)) |> 
                      #dplyr::filter(is.na(reason_excluded_resection_isolation)) |> 
                      dplyr::pull(array_mnp_CNVP_v12.8_v5.2_CNVP_bins)
                    , FUN = fun)
) 



df[1:5,1:5]

pca.raw <- prcomp(df)


plt <- pca.raw$rotation |> 
  as.data.frame(stringsAsFactors = F) |> 
  dplyr::select(paste0("PC", 1:5)) |> 
  tibble::rownames_to_column('bin') |> 
  dplyr::mutate(chr = gsub("^([^\\-]+).+$","\\1",bin)) |> 
  tidyr::pivot_longer(cols = -c(bin, chr)) |> 
  dplyr::rename(pca_rotation_val = value) |> 
  dplyr::filter(name == "PC1")


ggplot(plt, aes(x = bin, y = pca_rotation_val, col=chr)) +
  facet_grid(rows = vars(name), scales = "free", space="free") +
  geom_hline(yintercept = -0.021,col="gray", lty=2) +
  geom_point(pch=19,cex=0.5) +
  geom_hline(yintercept = -0.0345,col="gray", lty=2)



features <- plt |> 
  dplyr::filter(pca_rotation_val >= -0.0345 & pca_rotation_val <= -0.021) |> 
  dplyr::pull(bin)

stopifnot(length(features) == 1013)

# saveRDS(features, "cache/1p_19q_purity_most_informative_bins.Rds")


## 2. calc


glass_od.data.1p_19q_purity_most_informative_bins <- readRDS("cache/1p_19q_purity_most_informative_bins.Rds")


estimate_1p19q_purity_bin_based <- function(fn) {
  #fn <- glass_od.metadata.array_samples$heidelberg_cnvp_bins[1]
  
  sentrix <- gsub("^.+\\/([0-9]+_R[0-9]+C[0-9]+)\\.bins\\.igv$","\\1", fn)
  
  bins <- read.delim(fn) |> 
    dplyr::filter(Feature %in% glass_od.data.1p_19q_purity_most_informative_bins) |> 
    (function(.) {
      assertthat::assert_that(nrow(.) == 1013)
      return(.)
    })()
    
  median.lfc <- median(bins |> dplyr::pull(c(5)))
  purity <- -2 * (2^{median.lfc} - 1)
  sd <- sd(bins |> dplyr::pull(c(5)))
  
  out <- data.frame(array_sentrix_id = sentrix,
                    array_methylation_bins_1p19q_median.lfc = median.lfc,
                    array_methylation_bins_1p19q_purity = purity,
                    array_methylation_bins_1p19q_sd = sd
                    )
  
  return (out)
}



purities.bin <- do.call(
rbind,
pbapply::pblapply(glass_od.metadata.array_samples |> 
         dplyr::filter(!is.na(array_mnp_CNVP_v12.8_v5.2_CNVP_bins)) |> 
         dplyr::pull(array_mnp_CNVP_v12.8_v5.2_CNVP_bins)
       , estimate_1p19q_purity_bin_based)
) |> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_median.lfc)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_purity)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_sd)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 373)
    return(.)
  })()


saveRDS(purities.bin, file="cache/analysis_tumor_purity_EPIC_bin-based.Rds")





# mixed-in cortex: bin based ----






glass_od.data.1p_19q_purity_most_informative_bins <- readRDS("cache/1p_19q_purity_most_informative_bins.Rds")







purities.bin <- do.call(
  rbind,
  pbapply::pblapply(rbind(
    normal_cortex.metadata.array_samples |> dplyr::select(array_mnp_CNVP_v12.8_v5.2_CNVP_bins),
    cortex_dilution.metadata.array_samples |> dplyr::select(array_mnp_CNVP_v12.8_v5.2_CNVP_bins)
  ) |> 
                      dplyr::pull(array_mnp_CNVP_v12.8_v5.2_CNVP_bins)
                    , estimate_1p19q_purity_bin_based)
) |> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_median.lfc)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_purity)) |> 
  assertr::verify(is.numeric(array_methylation_bins_1p19q_sd)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 14+148)
    return(.)
  })()


saveRDS(purities.bin, file="cache/analysis_tumor_purity_EPIC_bin-based__normal_cortex.Rds")








