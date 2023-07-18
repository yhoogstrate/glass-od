#!/usr/bin/env


if(!exists('glass_od.metadata.idats')) {
  source('scripts/load_metadata.R')
}

source('scripts/load_chrom_sizes.R')


# obtain data ----


metadata <- glass_od.metadata.idats |> 
  dplyr::filter(is.na(reason_excluded_patient)) |> 
  dplyr::filter(is.na(reason_excluded_resection)) |> 
  dplyr::filter(is.na(reason_excluded_resection_isolation)) |> 
  dplyr::filter(is.na(reason_excluded_array_sample)) |> 
  dplyr::filter(study_name == "GLASS-OD") |> 
  assertr::verify(!is.na(qc.pca.outlier)) |> 
  dplyr::filter(qc.pca.outlier == F) |> 
  
  dplyr::filter(!is.na(resection_tumor_grade)) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  
  dplyr::mutate(resection_tumor_hg = ifelse(resection_tumor_grade == 3, 1, 0)) # 1 or 0 for regression multiplication factor


#'@todo remove interim samples, second g2's and non-last g3's etc.
#'@todo first grading needs to be complete


# get cnv profiles ----

f <- function(fn, purity) {
  #fn = metadata$heidelberg_cnvp_bins[1]
  
  d <- read.delim(fn) |> 
    dplyr::mutate(bid = paste0("cnvp_bin_",Chromosome,"_",Start,"_",End)) |> 
    tibble::column_to_rownames('bid') |> 
    dplyr::mutate(Chromosome = NULL, Start = NULL, End=NULL, Feature=NULL) 

  # moet wss iets anders...
  d[,1] <- d[,1] / purity
  
  d <- d |>  
    t() |> 
    as.data.frame()
  
  return (d)
}


data.cnv.profiles <- metadata |> 
  dplyr::select(sentrix_id, heidelberg_cnvp_bins, methylation_bins_1p19q_purity) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = f(heidelberg_cnvp_bins, methylation_bins_1p19q_purity)) |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(heidelberg_cnvp_bins = NULL) |> 
  dplyr::mutate(methylation_bins_1p19q_purity = NULL) |> 
  tibble::column_to_rownames('sentrix_id')


# limma - g2~g3 as condition ----


design <- model.matrix(~0 + resection_tumor_hg, data=metadata)
fit <- limma::lmFit(t(data.cnv.profiles), design)
fit <- limma::eBayes(fit)
stats.conditional <- limma::topTable(fit, n=nrow(fit), sort="p") |> # adjust="BH"
  tibble::rownames_to_column('cnvp_bin')


# per group quantiles [only for as-condition] ----


plt.quantiles.g2 <-  data.cnv.profiles |> 
  tibble::rownames_to_column('sentrix_id') |> 
  dplyr::filter(sentrix_id %in% (metadata |> dplyr::filter(resection_tumor_grade == 2) |> dplyr::pull(sentrix_id))) |> 
  tibble::column_to_rownames('sentrix_id') |> 
  pbapply::pblapply(function(x) {return(quantile(x))}) |> 
  as.data.frame() |> 
  t() |> 
  as.data.frame() |> 
  dplyr::rename_with( ~ paste0("quantile.", .x)) |> 
  tibble::rownames_to_column('bid') |> 
  #tidyr::pivot_longer(cols = -c('bid'), names_to = 'type', values_to='bin_lfc_normalised') |> 
  #dplyr::filter(type %in% c('quantile.25%','quantile.25%')) |> 
  
  dplyr::mutate(chr = gsub("^.+_(chr[^_]+)_.+$","\\1",bid)) |> 
  dplyr::filter(chr %in% c("chrX", "chrY") == F ) |> 
  dplyr::mutate(chr = factor(chr, levels = gtools::mixedsort(unique(as.character(chr))))) |> 
  dplyr::mutate(start = as.numeric(gsub("^.+_chr[^_]+_([0-9]+)_.+$","\\1",bid))) |> 
  dplyr::mutate(end = as.numeric(gsub("^.+_chr[^_]+_[0-9]+_([0-9]+)$","\\1",bid))) |> 
  dplyr::mutate(pos = (start + end)/2) |> 
  dplyr::mutate(bin_lfc_normalised = NULL)


ggplot(plt.quantiles.g2[1:160,], aes(x=pos/1000000, y=bin_lfc_normalised)) +
  geom_rect(aes(NULL,NULL,xmin=start,xmax=end,ymin=`quantile.25%`,ymax=`quantile.75%`), size=0.5, alpha=0.2)

  
# plt.quantiles.g3 <-  data.cnv.profiles |> 
#   tibble::rownames_to_column('sentrix_id') |> 
#   dplyr::filter(sentrix_id %in% (metadata |> dplyr::filter(resection_tumor_grade == 2) |> dplyr::pull(sentrix_id))) |> 
#   tibble::column_to_rownames('sentrix_id') |> 
#   pbapply::pblapply(function(x) {return(quantile(x))}) |> 
#   as.data.frame() |> 
#   t() |> 
#   as.data.frame() |> 
#   dplyr::rename_with( ~ paste0("quantile.", .x)) |> 
#   tibble::rownames_to_column('bid') |> 
#   tidyr::pivot_longer(cols = -c('bid'), names_to = 'type', values_to='bin_lfc_normalised') |> 
#   dplyr::filter(type %in% c('quantile.25%','quantile.25%'))





# plot

plt <- data.cnv.profiles |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('bid') |> 
  tidyr::pivot_longer(cols = -c(bid),names_to = "sentrix_id",values_to = 'bin_lfc_normalised') |> 
  dplyr::mutate(chr = gsub("^.+_(chr[^_]+)_.+$","\\1",bid)) |> 
  dplyr::filter(chr %in% c("chrX", "chrY") == F ) |> 
  dplyr::mutate(chr = factor(chr, levels = gtools::mixedsort(unique(as.character(chr))))) |> 
  dplyr::mutate(start = as.numeric(gsub("^.+_chr[^_]+_([0-9]+)_.+$","\\1",bid))) |> 
  dplyr::mutate(end = as.numeric(gsub("^.+_chr[^_]+_[0-9]+_([0-9]+)$","\\1",bid))) |> 
  dplyr::mutate(pos = (start + end)/2) |> 
  dplyr::left_join(metadata |> dplyr::select(sentrix_id, resection_tumor_grade), by=c('sentrix_id'='sentrix_id'), suffix = c('','')) 




ggplot(plt, aes(x=pos/1000000, y=bin_lfc_normalised, col=as.factor(resection_tumor_grade))) +
  facet_grid(cols = vars(chr), scales = "free", space="free")  +
  geom_point(pch=19,cex=0.2,alpha=0.05) +
  theme_bw() +
  geom_smooth() +
  coord_cartesian(ylim=c(-2.5, 2.5))




  




