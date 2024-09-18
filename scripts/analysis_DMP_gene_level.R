#!/usr/bin/env

# data ----


if(!exists("data.mvalues.probes")) {
  source("scripts/load_mvalues_hq_samples.R")
}


# gene annotations ----


source("scripts/load_gene_annotations.R")


# libs ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')

library(ggplot2)



# funcs ----



gene_enrich_ref <- function(gene, stats, idx) {
  target_probes <- idx |> 
    dplyr::filter(gene == gene) |> 
    dplyr::pull(probe_id)
  
  
  split1 <- stat |> 
    dplyr::filter(probe_id %in% target_probes) |> 
    dplyr::pull(DMP__g2_g3__pp_nc_PC1__t)
  
  split2 <- stat |> 
    dplyr::filter(probe_id %in% target_probes == F) |> 
    dplyr::pull(DMP__g2_g3__pp_nc_PC1__t)
  
  #print(paste0(gene,": ", length(split1)))
  
  a = t.test(split1, split2)
  
  return( data.frame(gene = gene,
                     t= a$statistic,
                     p.value = a$p.value,
                     n.target = length(split1),
                     n.offtarget = length(split2)) |> 
            tibble::remove_rownames())
}


gene_enrich_mu <- function(gene_t, mu, stats, idx) {
  target_probes <- idx |> 
    dplyr::filter(gene == gene_t) |> 
    dplyr::pull(probe_id)
  
  split1 <- stat |> 
    dplyr::filter(probe_id %in% target_probes) |> 
    dplyr::pull(DMP__g2_g3__pp_nc_PC1__t)
  
  
  a = t.test(split1, mu=mu)
  
  
  return( data.frame(gene     = gene_t,
                     mean     = a$estimate,
                     t        = a$statistic,
                     p.value  = a$p.value,
                     n.target = length(split1)) |> 
            tibble::remove_rownames())
}



# stat ----




stat <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



# gene annotations: GencodeCompV12_NAME ----



idx.GencodeCompV12_NAME <- stat |>
  dplyr::rename(gene = GencodeCompV12_NAME) |> 
  dplyr::select(probe_id, gene) |>
  tibble::tibble() |> 
  tidyr::separate_longer_delim(gene, delim = ";") |> 
  dplyr::filter(gene != "") |> 
  dplyr::distinct() |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 604089)
    return(.)
  })() |> 
  dplyr::filter(!grepl("^RP[0-9]+-", gene)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 531913)
    return(.)
  })() 
  




## t-test per median t-value ----


data.per.gene.GencodeCompV12_NAME <- idx.GencodeCompV12_NAME |> 
  dplyr::full_join(
    stat |> dplyr::select(probe_id, DMP__g2_g3__pp_nc_PC1__t),
    by=c('probe_id'='probe_id')
  ) |> 
  dplyr::group_by(gene) |> 
  dplyr::summarise(
    DMP__g2_g3__pp_nc_PC1__t = mean(DMP__g2_g3__pp_nc_PC1__t)
  ) |> 
  dplyr::ungroup()


out.per.gene.GencodeCompV12_NAME <- data.per.gene.GencodeCompV12_NAME |> 
  dplyr::mutate(p.value = dnorm(DMP__g2_g3__pp_nc_PC1__t, mean=0, sd=sd(data.per.gene.GencodeCompV12_NAME$DMP__g2_g3__pp_nc_PC1__t))) |> 
  dplyr::mutate(p.adj = p.adjust(p.value, method="fdr"))


rm(data.per.gene.GencodeCompV12_NAME)


ggplot(out.per.gene.GencodeCompV12_NAME, aes(x = DMP__g2_g3__pp_nc_PC1__t, y=-log(p.value), col=grepl("^HOX", gene))) + 
  geom_point(data=subset(out.per.gene.GencodeCompV12_NAME, !grepl("^HOX", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(out.per.gene.GencodeCompV12_NAME, grepl("^HOX", gene)), pch=19,cex=1) +
  theme_nature


ggplot(out.per.gene.GencodeCompV12_NAME, aes(x = DMP__g2_g3__pp_nc_PC1__t, y=-log(p.value), col=grepl("^RP[0-9]+-", gene))) + 
  geom_point(data=subset(out.per.gene.GencodeCompV12_NAME, !grepl("^RP[0-9]+-", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(out.per.gene.GencodeCompV12_NAME, grepl("^RP[0-9]+-", gene)), pch=19,cex=1) +
  theme_nature



ggplot(out.per.gene.GencodeCompV12_NAME, aes(x = DMP__g2_g3__pp_nc_PC1__t, y=-log(p.adj), col=grepl("^COL", gene))) + 
  geom_point(data=subset(out.per.gene.GencodeCompV12_NAME, !grepl("^COL", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(out.per.gene.GencodeCompV12_NAME, grepl("^COL", gene)), pch=19,cex=1) +
  theme_nature


ggplot(out.per.gene.GencodeCompV12_NAME, aes(x = DMP__g2_g3__pp_nc_PC1__t, y=-log(p.adj), col=grepl("^RP", gene))) + 
  geom_point(data=subset(out.per.gene.GencodeCompV12_NAME, !grepl("^RP", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(out.per.gene.GencodeCompV12_NAME, grepl("^RP", gene)), pch=19,cex=1) +
  theme_nature


ggplot(out.per.gene.GencodeCompV12_NAME, aes(x = DMP__g2_g3__pp_nc_PC1__t, y=-log(p.adj), col=grepl("^H[0-9]+", gene))) + 
  geom_point(data=subset(out.per.gene.GencodeCompV12_NAME, !grepl("^H[0-9]+", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(out.per.gene.GencodeCompV12_NAME, grepl("^H[0-9]+", gene)), pch=19,cex=1) +
  theme_nature


ggplot(out.per.gene.GencodeCompV12_NAME, aes(x = DMP__g2_g3__pp_nc_PC1__t, y=-log(p.adj), col=grepl("^HIST[0-9]+", gene))) + 
  geom_point(data=subset(out.per.gene.GencodeCompV12_NAME, !grepl("^HIST", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(out.per.gene.GencodeCompV12_NAME, grepl("^HIST", gene)), pch=19,cex=1) +
  theme_nature






## t-test per all t-values ----
# more confident for higher number of probes per gene


genes.GencodeCompV12_NAME = idx.GencodeCompV12_NAME |> 
  dplyr::group_by(gene) |> 
  dplyr::filter(dplyr::n() > 1) |> 
  dplyr::ungroup() |> 
  #head(n=10) |> 
  dplyr::pull(gene) |> 
  unique()



if(!exists('gene_enrichment_0.GencodeCompV12_NAME')) {
  if(file.exists('cache/gene_enrichment_0.GencodeCompV12_NAME.Rds')) {
    gene_enrichment_0.GencodeCompV12_NAME <- readRDS("cache/gene_enrichment_0.GencodeCompV12_NAME.Rds")
  } else {
    gene_enrichment_0.GencodeCompV12_NAME = do.call(rbind, pbapply::pblapply(genes.GencodeCompV12_NAME, gene_enrich_mu, mu=0, stats=stat, idx=idx.GencodeCompV12_NAME))
    saveRDS(gene_enrichment_0.GencodeCompV12_NAME, "cache/gene_enrichment_0.GencodeCompV12_NAME.Rds")
  }
}



### plt HOX ----


panels <- c(
  "ALX",
  "BARX",
  "DBX",
  "DLX",
   "EMX[0-9]$",
   "EN[12]$",
   "EVX",
  # "FOXA", # no homeobox
  # "FOXC",
  # "FOXD",
  # "FOXG",
  # "FOXL",
  # "FOXO",
  # "FOXP",
  "HMX",
   "HOXA|HOTTIP",
   "HOXB",
   "HOXC|HOTAIR",
   "HOXD",
   "GBX",
   "IRX",
   "ISL",
   "LHX",
   "LMX",
   "MKX",
   "MNX",
   "MSX",
   "NKX",
   "OTX",
   "PAX[0-9]+",
  "PHOX",
  "PITX",
   #"SOX",
   "SHOX",
  "SIX",
   "TBX[0-9]$",
   "TLX",
  "VSX"
)

plt <- gene_enrichment_0.GencodeCompV12_NAME
plts <- data.frame()

for(tf in panels) {
  
  plts <- rbind(plts,
                
                plt |> 
                  #dplyr::mutate(color = ifelse(grepl(paste0("^", tf) & !grepl("-AS[0-9]$", tf), gene), tf, "other")) |> 
                  dplyr::mutate(is.tf = ifelse(grepl(paste0("^", tf), gene) & !grepl("-AS[0-9]$", gene)  & !grepl("-AS$", gene) & !grepl("-IT[0-9]$", gene) & !grepl("-OT$", gene), "Yes", "No")) |> 
                  dplyr::mutate(facet = gsub("^([A-Z]+)[^A-Z]*?$","\\1", tf))
                ) 
}

plts <- plts |> 
  dplyr::group_by(facet) |> 
  dplyr::mutate(facet_i = dplyr::cur_group_id()) |> 
  dplyr::mutate(n.tfs = sum(is.tf == "Yes")) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(facet_x = facet_i %/% 8) |> 
  dplyr::mutate(facet_y = facet_i %% 4)

plts <- plts |> 
  dplyr::mutate(
    facet = factor(facet, levels = plts |> 
                     dplyr::filter(!duplicated(facet)) |> 
                     dplyr::arrange(-n.tfs, facet) |> 
                     dplyr::pull(facet))
  ) |> 
  dplyr::mutate(`-log(p.value)` = -log(p.value)) |> 
  dplyr::mutate(`log_pval_truncated` = `-log(p.value)` > 125) |> 
  dplyr::mutate(`-log(p.value)` = ifelse(`log_pval_truncated`, 125, `-log(p.value)`))



ggplot(plts, aes(x=mean, y=`-log(p.value)`, col=is.tf, label=gene, shape=log_pval_truncated)) + 
  facet_wrap(~facet, scale="free",ncol=6,nrow=5) +
  geom_point(data=subset(plts, is.tf == "No" & !log_pval_truncated), cex=0.01, alpha=0.5) +
  geom_point(data=subset(plts, is.tf == "No" & log_pval_truncated), cex=theme_nature_size/6) +
  ggrepel::geom_text_repel(data=subset(plts, gene %in% c("HOXD13", "TMPRSS3", "DAXX")),
                           box.padding = 0.65,
                           col="black", 
                           size=theme_nature_size,
                           family = theme_nature_font_family,
                           segment.size=theme_nature_lwd) +
  geom_point(data=subset(plts, gene %in% c("HOXD13", "TMPRSS3", "DAXX")), col="black", cex=theme_nature_size/3, pch=1) +
  geom_point(data=subset(plts, is.tf == "Yes"), size=theme_nature_size/3) +
  coord_cartesian(xlim = c(-7.5, 7.5)) + 
  coord_cartesian(ylim = c(0, 125)) + 
  scale_shape_manual(values=c(T=8,F=19,'TRUE'=8,'FALSE'=16)) +
  scale_color_manual(values=c('Yes'='red', 'No'=mixcol(col3(11)[10], "white", 0.25))) + 
  labs(x = "Mean t-score for all probes per gene") +
  theme_nature +
  theme(plot.background = element_rect(fill="white", colour=NA))   # png export


ggsave("output/figures/vis_analysis_DMP_gene_level.png", width=(8.5 * 0.975), height=4.5, dpi=600)






plts |> dplyr::filter(facet == "TBX") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)
plts |> dplyr::filter(facet == "NKX") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)
plts |> dplyr::filter(facet == "SOX") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)
plts |> dplyr::filter(facet == "HOXA") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)
plts |> dplyr::filter(facet == "HOXB") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)
plts |> dplyr::filter(facet == "HOXC") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)
plts |> dplyr::filter(facet == "HOXD") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)
plts |> dplyr::filter(facet == "PAX") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)

plts |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)




plt <- plt |> 
  dplyr::mutate(col = dplyr::case_when(
    grepl("^HOX|HOTAIR|HOTTIP", gene) ~ "HOX",
    gene %in% c(genes_polycomb_eed_homeobox, 
                genes_polycomb_h3k27_homeobox,
                genes_polycomb_prc2_homeobox,
                genes_polycomb_suz12_homeobox
                ) ~ "polycomb TF",
    T ~ "other"
  )) |> 
  dplyr::mutate(`-log(p.value)` = -log(p.value)) |> 
  dplyr::mutate(`log_pval_truncated` = `-log(p.value)` > 125) |> 
  dplyr::mutate(`-log(p.value)` = ifelse(`log_pval_truncated`, 125, `-log(p.value)`))




ggplot(plt, aes(x=mean, y=`-log(p.value)`, col=col, label=gene, shape=log_pval_truncated)) + 
  geom_point(data=subset(plt, col == "other" & !log_pval_truncated), cex=0.01, alpha=0.7) +
  geom_point(data=subset(plt, log_pval_truncated), cex=theme_nature_size/2) +
  geom_point(data=subset(plt, col != "other"), cex=theme_nature_size/6) +
  ggrepel::geom_text_repel(data=subset(plt, gene %in% c("HOXD13", "TMPRSS3", "DAXX", "TERT")),
                           box.padding = 0.95,
                           col="black",
                           size=theme_nature_size,
                           family = theme_nature_font_family,
                           segment.size=theme_nature_lwd) +
  geom_point(data=subset(plt, gene %in% c("HOXD13", "TMPRSS3", "DAXX", "TERT")), col="black", cex=theme_nature_size/2, pch=1) +
  scale_shape_manual(values=c(T=8,F=19,'TRUE'=8,'FALSE'=19)) +
  scale_color_manual(values=c( col3(11)[10],'gray', "red")) +
  xlim(c(-8,8)) +
  labs(x = "mean t-score probes per gene") +
  theme_nature +
  theme(plot.background = element_rect(fill="white", colour=NA))   # png export



ggsave("output/figures/vis_DMP_gene_level__polycomb_hox.png", width=(8.5 * 0.975)/3, height=2.25, dpi=600)




ggplot(gene_enrichment_0, aes(x=mean, y=-log(p.value), col=grepl("^COL", gene))) + 
  geom_point(data=subset(gene_enrichment_0, !grepl("^COL", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(gene_enrichment_0, grepl("^COL", gene)), pch=19,cex=1) +
  xlim(c(-10,10)) +
  theme_nature


ggplot(gene_enrichment_0, aes(x=mean, y=-log(p.value), col=grepl("^HIST[0-9]+", gene))) + 
  geom_point(data=subset(gene_enrichment_0, !grepl("^HIST[0-9]+", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(gene_enrichment_0, grepl("^HIST[0-9]+", gene)), pch=19,cex=1) +
  xlim(c(-25,25)) +
  theme_nature


ggplot(gene_enrichment_0, aes(x=t, y=-log(p.value), col=grepl("^H[0-9]+", gene))) + 
  geom_point(data=subset(gene_enrichment_0, !grepl("^H[0-9]+", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(gene_enrichment_0, grepl("^H[0-9]+", gene)), pch=19,cex=1) +
  xlim(c(-25,25)) +
  theme_nature


ggplot(gene_enrichment_0, aes(x=t, y=-log(p.value), col=grepl("^DAXX$", gene))) + 
  geom_point(data=subset(gene_enrichment_0, !grepl("^DAXX$", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(gene_enrichment_0, grepl("^DAXX$", gene)), pch=19,cex=1) +
  xlim(c(-25,25)) +
  theme_nature


ggplot(gene_enrichment_0, aes(x=t, y=-log(p.value), col=grepl("^ATRX$", gene))) + 
  geom_point(data=subset(gene_enrichment_0, !grepl("^ATRX$", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(gene_enrichment_0, grepl("^ATRX$", gene)), pch=19,cex=1) +
  xlim(c(-25,25)) +
  theme_nature


ggplot(gene_enrichment_0, aes(x=t, y=-log(p.value), col=grepl("^MKI67$", gene))) + 
  geom_point(data=subset(gene_enrichment_0, !grepl("^MKI67$", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(gene_enrichment_0, grepl("^MKI67$", gene)), pch=19,cex=1) +
  xlim(c(-25,25)) +
  theme_nature


gene_enrichment_0 |>
  dplyr::filter(t > 0) |> 
  dplyr::filter(!grepl("^RP[0-9]", gene)) |> 
  dplyr::filter(!grepl("^AC[0-9]", gene)) |> 
  View()



gene_enrichment_0 |>
  dplyr::filter(t < 0) |> 
  dplyr::filter(!grepl("^RP[0-9]", gene)) |> 
  dplyr::filter(!grepl("^AC[0-9]", gene)) |> 
  View()




# UCSC_RefGene_Name ----


idx.UCSC_RefGene_Name <- stat |> 
  dplyr::rename(gene = UCSC_RefGene_Name) |> 
  dplyr::select(probe_id, gene) |>
  tibble::tibble() |> 
  tidyr::separate_longer_delim(gene, delim = ";") |> 
  dplyr::filter(gene != "") |> 
  dplyr::distinct() |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 555727)
    return(.)
  })()





## t-test per median t-value ----




data.per.gene.UCSC_RefGene_Name <- idx.UCSC_RefGene_Name |> 
  dplyr::full_join(
    stat |> dplyr::select(probe_id, DMP__g2_g3__pp_nc_PC1__t),
    by=c('probe_id'='probe_id')
  ) |> 
  dplyr::group_by(gene) |> 
  dplyr::summarise(
    DMP__g2_g3__pp_nc_PC1__t = mean(DMP__g2_g3__pp_nc_PC1__t)
  ) |> 
  dplyr::ungroup()



#m = median(data.per.gene$DMP__g2_g3__pp_nc_PC1__t)
sd = sd(data.per.gene$DMP__g2_g3__pp_nc_PC1__t)

dat <- data.per.gene |> 
  dplyr::mutate(p.value = dnorm(DMP__g2_g3__pp_nc_PC1__t, mean=0, sd=sd)) |> 
  dplyr::mutate(p.adj = p.adjust(p.value, method="fdr"))





ggplot(dat, aes(x = DMP__g2_g3__pp_nc_PC1__t, y=-log(p.adj), col=grepl("^HOX", gene))) + 
  geom_point(data=subset(dat, !grepl("^HOX", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(dat, grepl("^HOX", gene)), pch=19,cex=1) +
  theme_nature


ggplot(dat, aes(x = median.t, y=-log(p.adj), col=grepl("^COL", gene))) + 
  geom_point(data=subset(dat, !grepl("^COL", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(dat, grepl("^COL", gene)), pch=19,cex=1) +
  theme_nature


ggplot(dat, aes(x = median.t, y=-log(p.adj), col=grepl("^RP", gene))) + 
  geom_point(data=subset(dat, !grepl("^RP", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(dat, grepl("^RP", gene)), pch=19,cex=1) +
  theme_nature


ggplot(dat, aes(x = median.t, y=-log(p.adj), col=grepl("^H[0-9]+", gene))) + 
  geom_point(data=subset(dat, !grepl("^H[0-9]+", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(dat, grepl("^H[0-9]+", gene)), pch=19,cex=1) +
  theme_nature


ggplot(dat, aes(x = median.t, y=-log(p.adj), col=grepl("^HIST[0-9]+", gene))) + 
  geom_point(data=subset(dat, !grepl("^HIST", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(dat, grepl("^HIST", gene)), pch=19,cex=1) +
  theme_nature



## t-test per all t-values ----
# more confident for higher number of probes per gene


genes.UCSC_RefGene_Name = idx.UCSC_RefGene_Name |> 
  dplyr::group_by(gene) |> 
  dplyr::filter(dplyr::n() > 1) |> 
  dplyr::ungroup() |> 
  dplyr::pull(gene) |> 
  unique()



if(!exists('gene_enrichment_0.UCSC_RefGene_Name')) {
  if(file.exists('cache/gene_enrichment_0.UCSC_RefGene_Name.Rds')) {
    gene_enrichment_0.UCSC_RefGene_Name <- readRDS("cache/gene_enrichment_0.UCSC_RefGene_Name.Rds")
  } else {
    gene_enrichment_0.UCSC_RefGene_Name = do.call(rbind, pbapply::pblapply(genes.UCSC_RefGene_Name, gene_enrich_mu, mu=0, stats=stat, idx=idx.UCSC_RefGene_Name))
    saveRDS(gene_enrichment_0.UCSC_RefGene_Name, "cache/gene_enrichment_0.UCSC_RefGene_Name.Rds")
  }
}

### plt ----




# hg38.manifest.gencode.v36.genesUniq ----



idx.hg38.manifest.gencode.v36.genesUniq <- stat |> 
  dplyr::rename(gene = hg38.manifest.gencode.v36.genesUniq) |> 
  dplyr::select(probe_id, gene) |>
  tibble::tibble() |> 
  tidyr::separate_longer_delim(gene, delim = ";") |> 
  dplyr::filter(gene != "") |> 
  dplyr::distinct() |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 784207)
    return(.)
  })()






## t-test per median t-value ----



data.per.gene <- idx |> 
  dplyr::full_join(
    stat |> dplyr::select(probe_id, DMP__g2_g3__pp_nc_PC1__t),
    by=c('probe_id'='probe_id')
  ) |> 
  dplyr::group_by(gene) |> 
  dplyr::summarise(
    DMP__g2_g3__pp_nc_PC1__t = median(DMP__g2_g3__pp_nc_PC1__t)
  ) |> 
  dplyr::ungroup()



#m = median(data.per.gene$DMP__g2_g3__pp_nc_PC1__t)
sd = sd(data.per.gene$DMP__g2_g3__pp_nc_PC1__t)

dat <- data.per.gene |> 
  dplyr::mutate(p.value = dnorm(DMP__g2_g3__pp_nc_PC1__t, mean=0, sd=sd)) |> 
  dplyr::mutate(p.adj = p.adjust(p.value, method="fdr"))



ggplot(dat, aes(x = DMP__g2_g3__pp_nc_PC1__t, y=-log(p.adj), col=grepl("^HOX", gene))) + 
  geom_point(data=subset(dat, !grepl("^HOX", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(dat, grepl("^HOX", gene)), pch=19,cex=1) +
  theme_nature


ggplot(dat, aes(x = median.t, y=-log(p.adj), col=grepl("^COL", gene))) + 
  geom_point(data=subset(dat, !grepl("^COL", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(dat, grepl("^COL", gene)), pch=19,cex=1) +
  theme_nature


ggplot(dat, aes(x = median.t, y=-log(p.adj), col=grepl("^RP", gene))) + 
  geom_point(data=subset(dat, !grepl("^RP", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(dat, grepl("^RP", gene)), pch=19,cex=1) +
  theme_nature


ggplot(dat, aes(x = median.t, y=-log(p.adj), col=grepl("^H[0-9]+", gene))) + 
  geom_point(data=subset(dat, !grepl("^H[0-9]+", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(dat, grepl("^H[0-9]+", gene)), pch=19,cex=1) +
  theme_nature


ggplot(dat, aes(x = median.t, y=-log(p.adj), col=grepl("^HIST[0-9]+", gene))) + 
  geom_point(data=subset(dat, !grepl("^HIST", gene)), pch=19,cex=0.01) +
  geom_point(data=subset(dat, grepl("^HIST", gene)), pch=19,cex=1) +
  theme_nature


## t-test per all t-values ----
# more confident for higher number of probes per gene


genes.hg38.manifest.gencode.v36.genesUniq = idx.hg38.manifest.gencode.v36.genesUniq |> 
  dplyr::group_by(gene) |> 
  dplyr::filter(dplyr::n() > 1) |> 
  dplyr::ungroup() |> 
  dplyr::pull(gene) |> 
  unique()




if(!exists('gene_enrichment_0.hg38.manifest.gencode.v36.genesUniq')) {
  if(file.exists('cache/gene_enrichment_0.hg38.manifest.gencode.v36.genesUniq.Rds')) {
    gene_enrichment_0.hg38.manifest.gencode.v36.genesUniq <- readRDS("cache/gene_enrichment_0.hg38.manifest.gencode.v36.genesUniq.Rds")
  } else {
    gene_enrichment_0.hg38.manifest.gencode.v36.genesUniq = do.call(rbind, pbapply::pblapply(genes.hg38.manifest.gencode.v36.genesUniq, gene_enrich_mu, mu=0, stats=stat, idx=idx.hg38.manifest.gencode.v36.genesUniq))
    saveRDS(gene_enrichment_0.hg38.manifest.gencode.v36.genesUniq, "cache/gene_enrichment_0.hg38.manifest.gencode.v36.genesUniq.Rds")
  }
}



### plt ----

