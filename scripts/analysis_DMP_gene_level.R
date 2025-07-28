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
    assertthat::assert_that(nrow(.) == 604185)
    return(.)
  })() |> 
  dplyr::filter(!grepl("^RP[0-9]+-", gene)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 532000)
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


fn <- 'cache/gene_enrichment_0.GencodeCompV12_NAME.Rds'
if(!exists(fn)) {
  if(file.exists(fn)) {
    gene_enrichment_0.GencodeCompV12_NAME <- readRDS(fn)
  } else {
    gene_enrichment_0.GencodeCompV12_NAME = do.call(rbind, pbapply::pblapply(genes.GencodeCompV12_NAME, gene_enrich_mu, mu=0, stats=stat, idx=idx.GencodeCompV12_NAME))
    saveRDS(gene_enrichment_0.GencodeCompV12_NAME, fn)
  }
} else {
  gene_enrichment_0.GencodeCompV12_NAME <- readRDS(fn)
}
rm(fn)


gene_enrichment_0.GencodeCompV12_NAME <- gene_enrichment_0.GencodeCompV12_NAME |> 
  dplyr::mutate(q.value = p.adjust(p.value))



gene_enrichment_0.GencodeCompV12_NAME |>
  dplyr::filter(gene %in% c('TERT', 'DAXX'))



### plt HOX ----



panels <- c( # CATNON
  "HAND",
  #"HOXA10",
  #"HOXA11",
  #"HOXA6",
  #"HOXA7",
  "HOXA",
  #"HOXC4",
  "HOXC",
  #"HOXD10",
  #"HOXD3",
  #"HOXD4",
  "HOXD[0-9]$",
  "IRX[0-9]$",
  "ISL",
  "MNX",
  "OCIAD",
  #"OSR1",
  "OSR|OTP",
  #"PAX3",
  "PAX[0-9]$",
  "SHOX[0-9]$",
  "SIM",
  "SIX",
  "TBX[0-9]$",
  "TWIST[0-9]$",
  "VAX[0-9]$"
)

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
  "OTP|OTX",
  "OSR[0-9]+",
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




names <- plts |> 
  dplyr::filter(is.tf == "Yes") |> 
  dplyr::select(facet, gene) |> 
  dplyr::arrange(facet, gene) |> 
  dplyr::group_by(facet) |> 
  dplyr::summarise(genes = paste0(gene, collapse=", "))


plts <- plts |> 
  dplyr::left_join(names, by=c('facet'='facet'), suffix=c('',''))


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
  facet_wrap(~facet, scale="free") + #,ncol=6,nrow=5) +
  geom_point(data=subset(plts, is.tf == "No" & !log_pval_truncated), cex=0.01, alpha=0.5) +
  geom_point(data=subset(plts, is.tf == "No" & log_pval_truncated), cex=theme_nature_size/6) +
  ggrepel::geom_text_repel(data=subset(plts, gene %in% c("HOXD13", "TMPRSS3", "DAXX", "TERT")),
                           box.padding = 0.65,
                           col="black", 
                           size=theme_nature_size,
                           family = theme_nature_font_family,
                           segment.size=theme_nature_lwd) +
  geom_point(data=subset(plts, gene %in% c("HOXD13", "TMPRSS3", "DAXX", "TERT")), col="black", cex=theme_nature_size/3, pch=1) +
  geom_point(data=subset(plts, is.tf == "Yes"), size=theme_nature_size/3) +
  coord_cartesian(xlim = c(-7.5, 7.5)) + 
  coord_cartesian(ylim = c(0, 125)) + 
  scale_shape_manual(values=c(T=8,F=19,'TRUE'=8,'FALSE'=16)) +
  scale_color_manual(values=c('Yes'='red', 'No'=mixcol(col3(11)[10], "white", 0.25))) + 
  labs(x = "Mean t-score per gene") +
  theme_nature +
  theme(plot.background = element_rect(fill="white", colour=NA))   # png export



ggsave("output/figures/vis_analysis_DMP_gene_level.png", width=(8.5 * 0.975), height=4.5, dpi=600)



ggplot(plts, aes(x=genes,
                 y=`mean`
                 )) + 
  #ggplot2::geom_violin(draw_quantiles = c(0.5), linewidth=theme_nature_lwd, col = "darkgray", adjust = 1.95) +
  #see::geom_violinhalf() +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_nature_lwd, lty=2) +
  gghalves::geom_half_violin(fill=mixcol("#FFFFFF", "#EEEEEE"),side = "r", draw_quantiles = c(0.5), linewidth=theme_nature_lwd, col = "darkgray", width = 1.6) +
  #geom_point(data=subset(plts, is.tf == "Yes"), size=theme_nature_size/3, col="red") +
  #ggbeeswarm::geom_beeswarm(data=subset(plts, is.tf == "Yes"), size=theme_nature_size/3, col="red", side=1, method="compactswarm") +
  geom_point(data=subset(plts, is.tf == "Yes"), size=theme_nature_size, col="red", pch="|") +
  coord_flip() +
  labs(y= "Mean t-score for all CpGs per genes", x=NULL) +
  theme_nature

ggsave("output/figures/vis_analysis_DMP_gene_level__TFs_gghalve.pdf", width=(8.5 * 0.975), height=4.0, dpi=600)




ggplot(plts, aes(x=median, y=`-log(p.value)`, col=is.tf, label=gene, shape=log_pval_truncated)) + 
  facet_wrap(~facet, scale="free") + #,ncol=6,nrow=5) +
  geom_point(data=subset(plts, is.tf == "No" & !log_pval_truncated), cex=0.01, alpha=0.5) +
  geom_point(data=subset(plts, is.tf == "No" & log_pval_truncated), cex=theme_nature_size/6) +
  ggrepel::geom_text_repel(data=subset(plts, gene %in% c("HOXD13", "TMPRSS3", "DAXX", "TERT")),
                           box.padding = 0.65,
                           col="black", 
                           size=theme_nature_size,
                           family = theme_nature_font_family,
                           segment.size=theme_nature_lwd) +
  geom_point(data=subset(plts, gene %in% c("HOXD13", "TMPRSS3", "DAXX", "TERT")), col="black", cex=theme_nature_size/3, pch=1) +
  geom_point(data=subset(plts, is.tf == "Yes"), size=theme_nature_size/3) +
  coord_cartesian(xlim = c(-7.5, 7.5)) + 
  coord_cartesian(ylim = c(0, 125)) + 
  scale_shape_manual(values=c(T=8,F=19,'TRUE'=8,'FALSE'=16)) +
  scale_color_manual(values=c('Yes'='red', 'No'=mixcol(col3(11)[10], "white", 0.25))) + 
  labs(x = "Mean t-score for all probes per gene") +
  theme_nature +
  theme(plot.background = element_rect(fill="white", colour=NA))   # png export





plts |> dplyr::filter(facet == "TBX") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)
plts |> dplyr::filter(facet == "NKX") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)
plts |> dplyr::filter(facet == "SOX") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)
plts |> dplyr::filter(facet == "HOXA") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)
plts |> dplyr::filter(facet == "HOXB") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)
plts |> dplyr::filter(facet == "HOXC") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)
plts |> dplyr::filter(facet == "HOXD") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)
plts |> dplyr::filter(facet == "PAX") |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)

plts |> dplyr::filter(is.tf == "Yes") |> dplyr::pull(gene)




### combined ----


panels <- c(
  "ALX",
  "BARX",
  "DBX",
  "DLX",
  "EMX[0-9]$",
  "EN[12]$",
  "EVX",
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
  "OTP|OTX",
  "OSR[0-9]+",
  "PAX[0-9]+",
  "PHOX",
  "PITX",
  #"SOX",
  "SHOX",
  "SIX",
  "TBX[0-9]$",
  "TLX",
  "VSX",
  
  "^OR",
  
  "^IDH1$",
  "^IDH2$",
  #"^IDH3$",
  "DAXX",
  "ATRX",
  "TERT",
  "H3F3A",
  "CIC",
  "FUBP1",
  "TCF12"
  
  #"GATA6",
  #"GATA4",
  #"^KRT"
  
)


facettes <- rep("?", length(panels))
facettes[1:which(panels == "VSX")] <- "Homeobox TFs"
facettes[which(panels == "^OR")] <- "Olfactory Receptors"
facettes[which(panels == "^IDH1$"):length(panels)] <- "Disease implicated"





plt <- gene_enrichment_0.GencodeCompV12_NAME
plts <- data.frame()

for(i in 1:1:length(panels)) {
  tf <- panels[i]
  screen <- facettes[i]

  plts <- rbind(plts,
                plt |> 
                  dplyr::mutate(is.tf = ifelse(grepl(paste0("^", tf), gene) & !grepl("-AS[0-9]$", gene)  & !grepl("-AS$", gene) & !grepl("-IT[0-9]$", gene) & !grepl("-OT$", gene), "Yes", "No")) |> 
                  dplyr::mutate(facet = gsub("^([A-Z]+)[^A-Z]*?$","\\1", tf)) |> 
                  dplyr::mutate(screen = screen)
  ) 
}



a =               plt |> 
  dplyr::mutate(HOX = grepl("^HOX|HOTAIR|HOTTIP", gene) & !grepl("-AS[0-9]$", gene) & !grepl("-AS$", gene)) |> 
  dplyr::mutate(polycomb =  gene %in% c(genes_polycomb_eed_homeobox, 
                                        genes_polycomb_h3k27_homeobox,
                                        genes_polycomb_prc2_homeobox,
                                        genes_polycomb_suz12_homeobox
  )) |> 
  dplyr::mutate(is.tf = ifelse(HOX | polycomb, "Yes", "No")) |> 
  dplyr::mutate(col = dplyr::case_when(
    HOX & polycomb ~ "HOX & polycomb",
    HOX ~ "HOX",
    polycomb ~ "polycomb",
    T ~ "not selected"
  )) |> 
  dplyr::mutate(facet = "Polycomb TF & HOX") |> 
  dplyr::mutate(screen = "Polycomb TF & HOX") |> 
  dplyr::mutate(HOX = NULL) |> 
  dplyr::mutate(polycomb = NULL)


plts <- rbind(plts |> dplyr::mutate(col = "see y-axis" ), 

a
              
              )






plts <- plts |> 
  dplyr::arrange(facet, screen) |> 
  dplyr::mutate(screen = factor(screen, levels=c(
    "Homeobox TFs", "Disease implicated" , "Olfactory Receptors", "Polycomb TF & HOX"
  )))


names <- plts |> 
  dplyr::filter(is.tf == "Yes") |> 
  dplyr::select(facet, gene) |> 
  dplyr::arrange(facet, gene) |> 
  dplyr::group_by(facet) |> 
  dplyr::summarise(genes = paste0(gene, collapse=", ")) |> 
  dplyr::mutate(genes = ifelse(grepl("^OR", genes), "Olfactory Receptors", genes)) |> 
  dplyr::mutate(genes = ifelse(facet == "Polycomb TF & HOX", facet, genes))


plts <- plts |> 
  dplyr::left_join(names, by=c('facet'='facet'), suffix=c('',''))


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





plts |>
  dplyr::filter(gene %in% c("DAXX", "TERT", "CIC", "IDH2", "H3F3A"))  |> 
  dplyr::filter(facet == plts$facet[1]) |> 
  dplyr::select(gene, p.value, q.value) |> 
  dplyr::mutate(signi = q.value < 0.01) |> 
  dplyr::mutate(q.str = format.pval(q.value, digits=3)) |> 
  as.data.frame()




ggplot(plts, aes(x=reorder(genes, dplyr::desc(genes)),
                 y=`mean`,
                 color=col,
                 label=gene
)) + 
  facet_grid(rows = vars(screen), scales = "free", space="free") +
  #ggplot2::geom_violin(draw_quantiles = c(0.5), linewidth=theme_nature_lwd, col = "darkgray", adjust = 1.95) +
  #see::geom_violinhalf() +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_nature_lwd, lty=2) +
  gghalves::geom_half_violin(fill=mixcol("#FFFFFF", "#EEEEEE"),side = "r", draw_quantiles = c(0.5), linewidth=theme_nature_lwd, col = "darkgray", width = 1.6) +
  #geom_point(data=subset(plts, is.tf == "Yes"), size=theme_nature_size/3, col="red") +
  #ggbeeswarm::geom_beeswarm(data=subset(plts, is.tf == "Yes"), size=theme_nature_size/3, col="red", side=1, method="compactswarm") +
  geom_point(data=subset(plts, is.tf == "Yes"), size=theme_nature_size, pch="|") +
  ggrepel::geom_text_repel(segment.size=theme_nature_lwd,
                           data = plts |> dplyr::filter(is.tf == "Yes" &
                                                          gene %in% c("HOXD12", "HOXD13") &
                                                          screen == "Homeobox TFs"),
                           size=theme_nature_size,
                           family=theme_nature_font_family,
                           col="#444444",
                           nudge_x = -1.5,
                           nudge_y = 0.5,
                           min.segment.length = 0
                           ) +
  coord_flip() +
  labs(y= "Mean t-score (methylation change Grade 2 - Grade 3) for all CpGs per genes", x=NULL) +
  theme_nature +
  scale_color_manual(values=c(
    `see y-axis` = 'red',
    `polycomb` = col3(11)[10],
    `HOX` = 'red',
    `HOX & polycomb` = 'purple'
  ))



ggsave("output/figures/vis_analysis_DMP_gene_level__TFs_gghalve_l__combi.pdf", width=(8.5 * 0.975), height=3.50, dpi=600)







### polycomb ----



plt <- gene_enrichment_0.GencodeCompV12_NAME



plt <- plt |> 
  dplyr::mutate(HOX = grepl("^HOX|HOTAIR|HOTTIP", gene) & !grepl("-AS[0-9]$", gene) & !grepl("-AS$", gene)) |> 
  dplyr::mutate(polycomb =  gene %in% c(genes_polycomb_eed_homeobox, 
                                        genes_polycomb_h3k27_homeobox,
                                        genes_polycomb_prc2_homeobox,
                                        genes_polycomb_suz12_homeobox
  )) |> 
  dplyr::mutate(col = dplyr::case_when(
    HOX & polycomb ~ "HOX & polycomb",
    HOX ~ "HOX",
    polycomb ~ "polycomb",
    T ~ "other"
  )) |> 
  dplyr::mutate(`-log(p.value)` = -log(p.value)) |> 
  dplyr::mutate(`log_pval_truncated` = `-log(p.value)` > 60) |> 
  dplyr::mutate(`-log(p.value)` = ifelse(`log_pval_truncated`, 60, `-log(p.value)`))



ggplot(plt, aes(x=mean, y=`-log(p.value)`, col=col, label=gene, shape=log_pval_truncated)) + 
  geom_point(data=subset(plt, col == "other" & !log_pval_truncated), cex=0.01, alpha=0.7) +
  geom_point(data=subset(plt, log_pval_truncated), cex=theme_nature_size/4) +
  geom_point(data=subset(plt, col != "other"), cex=theme_nature_size/6) +
  ggrepel::geom_text_repel(data=subset(plt, gene %in% c("HOXD12", "HOXD13", "DAXX", "TERT")),
                           box.padding = 0.95,
                           col="black",
                           size=theme_nature_size,
                           family = theme_nature_font_family,
                           segment.size=theme_nature_lwd) +
  geom_point(data=subset(plt, gene %in% c("HOXD12", "HOXD13", "DAXX", "TERT")), col="black", cex=theme_nature_size/2, pch=1) +
  scale_shape_manual(values=c(T=8,F=19,'TRUE'=8,'FALSE'=19)) +
  scale_color_manual(values=c(
    `polycomb` = col3(11)[10],
    `HOX` = 'red',
    `HOX & polycomb` = 'purple',
    `other`='gray'
    )) +
  xlim(c(-8,8)) +
  labs(x = "mean t-score probes per gene") +
  theme_nature +
  theme(plot.background = element_rect(fill="white", colour=NA))   # png export



ggsave("output/figures/vis_DMP_gene_level__polycomb_hox.png", width=2.0925, height=2.28, dpi=600)






plt <- plt |> 
  dplyr::mutate(HOX = grepl("^HOX|HOTAIR|HOTTIP", gene) & !grepl("-AS[0-9]$", gene) & !grepl("-AS$", gene)) |> 
  dplyr::mutate(polycomb =  gene %in% c(genes_polycomb_eed_homeobox, 
                                        genes_polycomb_h3k27_homeobox,
                                        genes_polycomb_prc2_homeobox,
                                        genes_polycomb_suz12_homeobox
  )) |> 
  dplyr::mutate(col = grepl("^OR", gene)) |> 
  dplyr::mutate(`-log(p.value)` = -log(p.value)) |> 
  dplyr::mutate(`log_pval_truncated` = `-log(p.value)` > 60) |> 
  dplyr::mutate(`-log(p.value)` = ifelse(`log_pval_truncated`, 60, `-log(p.value)`))










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



### examples 09746sitable1.pdf ----
# https://www.pnas.org/doi/full/10.1073/pnas.0609746104#glossary



plt <- plt |> 
  dplyr::mutate(
    col = gene %in% c(
      "TBX15","PRRX1","LHX8","PBX1","LMO4","Sp5","HOXD13","MEIS1","MEIS1","PAX3","HOXD3","HOXD4","TBR1","DLX1",
      "DLX2","OTX1","SIX3","ID2","NEUROD1","EN1","GBX2","SOX2","SOX14","TBR2","SHOX2","ZIC4",
      "ZIC1","FOXP1","BAPX1","PITX2","HAND2","NEUROG2","IRX2","IRX1","ISL1","NKX2-5","MEF2C","OTP1","SIM1","TFAP2A","TFAP2D","TFAP2A","TWIST1","Sp8","MEOX2","DLX5,6","HOXA3","HOXA4","HOXA5","HOXA6","HOXA7","HOXA9","HOXA10","HOXA11","CBX3",
      "SOX7","SOX17","GATA4","OSR2","DMRT1","DMRT3","NR6A1","BNC2","GATA3","PCGF4","ARID5B","KLF6",
      "EMX2","PAX6","FLI1","TBX5","CART1","TBX3","LHX5","HOXC","ZIC2","ZIC5","CDX2",
      "POU4F1","FOXG1B","SIX1","OTX2","FOXA1","TITF1","SIX6","GSC","MEIS2","ONECUT1","LBXCOR1","MEIS2","FOXB1","FOXF1","IRX3","IRX5","CBX2,CBX6","HOXB","HOXB","SALL3","ONECUT2","NKX2-2","PAX1","MAFB","SIM2"
    )) |> 
  dplyr::mutate(`-log(p.value)` = -log(p.value)) |> 
  dplyr::mutate(`log_pval_truncated` = `-log(p.value)` > 60) |> 
  dplyr::mutate(`-log(p.value)` = ifelse(`log_pval_truncated`, 60, `-log(p.value)`)) |> 
  dplyr::mutate(screen=1)



panels <- c(
  "ALX",
  "BARX",
  "DBX",
  "DLX",
  "EMX[0-9]$",
  "EN[12]$",
  "EVX",
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
  "OTP|OTX",
  "OSR[0-9]+",
  "PAX[0-9]+",
  "PHOX",
  "PITX",
  #"SOX",
  "SHOX",
  "SIX",
  "TBX[0-9]$",
  "TLX",
  "VSX",
  
  "^OR",
  
  "^IDH1$",
  "^IDH2$",
  #"^IDH3$",
  "DAXX",
  "ATRX",
  "TERT",
  "H3F3A",
  "CIC",
  "FUBP1",
  "TCF12",
  
  "GATA6",
  "GATA4",
  "^KRT"
  
)


facettes <- rep("?", length(panels))
facettes[1:which(panels == "VSX")] <- "Homeobox TFs"
facettes[which(panels == "^OR")] <- "Olfactory Receptors"
facettes[which(panels == "^IDH1$"):length(panels)] <- "Disease implicated"





plt <- gene_enrichment_0.GencodeCompV12_NAME
plts <- data.frame()

for(i in 1:1:length(panels)) {
  tf <- panels[i]
  screen <- facettes[i]
  
  plts <- rbind(plts,
                plt |> 
                  dplyr::mutate(is.tf = ifelse(grepl(paste0("^", tf), gene) & !grepl("-AS[0-9]$", gene)  & !grepl("-AS$", gene) & !grepl("-IT[0-9]$", gene) & !grepl("-OT$", gene), "Yes", "No")) |> 
                  dplyr::mutate(facet = gsub("^([A-Z]+)[^A-Z]*?$","\\1", tf)) |> 
                  dplyr::mutate(screen = screen)
  ) 
}



a =               plt |> 
  dplyr::mutate(HOX = grepl("^HOX|HOTAIR|HOTTIP", gene) & !grepl("-AS[0-9]$", gene) & !grepl("-AS$", gene)) |> 
  dplyr::mutate(polycomb =  gene %in% c(genes_polycomb_eed_homeobox, 
                                        genes_polycomb_h3k27_homeobox,
                                        genes_polycomb_prc2_homeobox,
                                        genes_polycomb_suz12_homeobox
  )) |> 
  dplyr::mutate(is.tf = ifelse(HOX | polycomb, "Yes", "No")) |> 
  dplyr::mutate(col = dplyr::case_when(
    HOX & polycomb ~ "HOX & polycomb",
    HOX ~ "HOX",
    polycomb ~ "polycomb",
    T ~ "not selected"
  )) |> 
  dplyr::mutate(facet = "Polycomb TF & HOX") |> 
  dplyr::mutate(screen = "Polycomb TF & HOX") |> 
  dplyr::mutate(HOX = NULL) |> 
  dplyr::mutate(polycomb = NULL)


plts <- rbind(plts |> dplyr::mutate(col = "see y-axis" ), 
              
              a
              
)






plts <- plts |> 
  dplyr::arrange(facet, screen) |> 
  dplyr::mutate(screen = factor(screen, levels=c(
    "Homeobox TFs", "Disease implicated" , "Olfactory Receptors", "Polycomb TF & HOX"
  )))


names <- plts |> 
  dplyr::filter(is.tf == "Yes") |> 
  dplyr::select(facet, gene) |> 
  dplyr::arrange(facet, gene) |> 
  dplyr::group_by(facet) |> 
  dplyr::summarise(genes = paste0(gene, collapse=", ")) |> 
  dplyr::mutate(genes = ifelse(grepl("^OR", genes), "Olfactory Receptors", genes)) |> 
  dplyr::mutate(genes = ifelse(facet == "Polycomb TF & HOX", facet, genes))


plts <- plts |> 
  dplyr::left_join(names, by=c('facet'='facet'), suffix=c('',''))


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





plts |>
  dplyr::filter(gene %in% c("DAXX", "TERT", "CIC", "IDH2", "H3F3A"))  |> 
  dplyr::filter(facet == plts$facet[1]) |> 
  dplyr::select(gene, p.value, q.value) |> 
  dplyr::mutate(signi = q.value < 0.01) |> 
  dplyr::mutate(q.str = format.pval(q.value, digits=3)) |> 
  as.data.frame()




ggplot(plts, aes(x=reorder(genes, dplyr::desc(genes)),
                 y=`mean`,
                 color=col,
                 label=gene
)) + 
  facet_grid(rows = vars(screen), scales = "free", space="free") +
  #ggplot2::geom_violin(draw_quantiles = c(0.5), linewidth=theme_nature_lwd, col = "darkgray", adjust = 1.95) +
  #see::geom_violinhalf() +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_nature_lwd, lty=2) +
  gghalves::geom_half_violin(fill=mixcol("#FFFFFF", "#EEEEEE"),side = "r", draw_quantiles = c(0.5), linewidth=theme_nature_lwd, col = "darkgray", width = 1.6) +
  #geom_point(data=subset(plts, is.tf == "Yes"), size=theme_nature_size/3, col="red") +
  #ggbeeswarm::geom_beeswarm(data=subset(plts, is.tf == "Yes"), size=theme_nature_size/3, col="red", side=1, method="compactswarm") +
  geom_point(data=subset(plts, is.tf == "Yes"), size=theme_nature_size, pch="|") +
  ggrepel::geom_text_repel(segment.size=theme_nature_lwd,
                           data = plts |> dplyr::filter(is.tf == "Yes" & 
                                                          gene %in% c("HOXD12", "HOXD13") & 
                                                          screen == "Homeobox TFs"),
                           size=theme_nature_size,
                           family=theme_nature_font_family,
                           col="#444444",
                           nudge_x = -1.5,
                           nudge_y = 0.5,
                           min.segment.length = 0
  ) +
  coord_flip() +
  labs(y= "Mean t-score (methylation change Grade 2 - Grade 3) for all CpGs per genes", x=NULL) +
  theme_nature +
  scale_color_manual(values=c(
    `see y-axis` = 'red',
    `polycomb` = col3(11)[10],
    `HOX` = 'red',
    `HOX & polycomb` = 'purple'
  ))






### Olfactory receptors ----

plt <- plt |> 
  dplyr::mutate(HOX = grepl("^HOX|HOTAIR|HOTTIP", gene) & !grepl("-AS[0-9]$", gene) & !grepl("-AS$", gene)) |> 
  dplyr::mutate(polycomb =  gene %in% c(genes_polycomb_eed_homeobox, 
                                        genes_polycomb_h3k27_homeobox,
                                        genes_polycomb_prc2_homeobox,
                                        genes_polycomb_suz12_homeobox
  )) |> 
  dplyr::mutate(col = grepl("^OR", gene)) |> 
  dplyr::mutate(`-log(p.value)` = -log(p.value)) |> 
  dplyr::mutate(`log_pval_truncated` = `-log(p.value)` > 60) |> 
  dplyr::mutate(`-log(p.value)` = ifelse(`log_pval_truncated`, 60, `-log(p.value)`))


ggplot(plt, aes(x=mean, y=`-log(p.value)`, col=col, label=gene, shape=log_pval_truncated)) + 
  geom_point(data=subset(plt, col == F ), cex=0.01, alpha=0.7) +
  #geom_point(data=subset(plt, log_pval_truncated), cex=theme_nature_size/4) +
  geom_point(data=subset(plt, col == T), cex=theme_nature_size/6) +
  scale_shape_manual(values=c(T=8,F=19,'TRUE'=8,'FALSE'=19)) +
  scale_color_manual(values=c(
    'TRUE' = 'red',
    'FALSE' ='gray'
  )) +
  xlim(c(-8,8)) +
  labs(x = "mean t-score probes per gene") +
  theme_nature +
  theme(plot.background = element_rect(fill="white", colour=NA))   # png export


ggsave("output/figures/vis_DMP_gene_level__olfactory.png", width=2.0925, height=2.28, dpi=600)



### test method ----



tmp <- idx.GencodeCompV12_NAME |> 
  dplyr::full_join(
    stat |> dplyr::select(probe_id, DMP__g2_g3__pp_nc_PC1__t),
    by=c('probe_id'='probe_id')
  ) |> 
  dplyr::group_by(gene) |> 
  dplyr::summarise(
    DMP__g2_g3__pp_nc_PC1__t = mean(DMP__g2_g3__pp_nc_PC1__t)
  ) |> 
  dplyr::ungroup()





## GSEA analysis ----

tmp <- gene_enrichment_0.GencodeCompV12_NAME |> 
  dplyr::left_join(out.per.gene.GencodeCompV12_NAME, by=c('gene'='gene')) |> 
  dplyr::filter(!is.na(DMP__g2_g3__pp_nc_PC1__t)) |> 
  dplyr::filter(!is.na(t))


ggplot(tmp, aes(x=t, y=DMP__g2_g3__pp_nc_PC1__t)) +
  geom_point(size=theme_nature_size) + 
  xlim(-60,60)


# https://rdrr.io/bioc/missMethyl/src/R/gometh.R
## out per gene seems better statistic
out.per.gene.GencodeCompV12_NAME |> 
  dplyr::filter(p.adj < 0.01) |> 
  View()





sel_signi <- data.mvalues.probes |> dplyr::filter(!is.na(DMP__g2_g3__pp_nc_PC1__adj.P.Val)) |> dplyr::filter(DMP__g2_g3__pp_nc_PC1__adj.P.Val < 0.01) |> dplyr::pull(probe_id)
all_cpgs = data.mvalues.probes |> dplyr::filter(!is.na(DMP__g2_g3__pp_nc_PC1__adj.P.Val)) |> dplyr::pull(probe_id)

# ? GOmeth - https://rdrr.io/bioc/missMethyl/man/gometh.html
gometh <- missMethyl::gometh(sig.cpg = sel_signi,
                             all.cpg = all_cpgs,
                             collection = "GO", # KEGG
                             array.type = "EPIC", 
                             prior.prob=FALSE,
                             anno = minfi::getAnnotation('IlluminaHumanMethylationEPICanno.ilm10b4.hg19'))


sel_signi_up <- data.mvalues.probes |> dplyr::filter(!is.na(DMP__g2_g3__pp_nc_PC1__adj.P.Val)) |> dplyr::filter(DMP__g2_g3__pp_nc_PC1__adj.P.Val < 0.0001 & DMP__g2_g3__pp_nc_PC1__t > 0) |> dplyr::pull(probe_id)
gometh_up <- missMethyl::gometh(sig.cpg = sel_signi_up,
                                all.cpg = all_cpgs,
                                collection = "GO", # KEGG
                             array.type = "EPIC", 
                             prior.prob=FALSE,
                             sig.genes=T,
                             anno = minfi::getAnnotation('IlluminaHumanMethylationEPICanno.ilm10b4.hg19'))


#sel_signi_down <- data.mvalues.probes |> dplyr::filter(runif(nrow(data.mvalues.probes),0,1) < (20*0.01)) |> dplyr::filter(!is.na(DMP__g2_g3__pp_nc_PC1__adj.P.Val)) |> dplyr::pull(probe_id)
sel_signi_down <- data.mvalues.probes |> dplyr::filter(!is.na(DMP__g2_g3__pp_nc_PC1__adj.P.Val)) |> dplyr::filter(DMP__g2_g3__pp_nc_PC1__adj.P.Val < 0.0001 & DMP__g2_g3__pp_nc_PC1__t < 0) |> dplyr::pull(probe_id)
gometh_down <- missMethyl::gometh(sig.cpg = sel_signi_down,
                             all.cpg = all_cpgs,
                             collection = "KEGG",  # GO
                             array.type = "EPIC", 
                             prior.prob=FALSE,
                             sig.genes=T,
                             anno = minfi::getAnnotation('IlluminaHumanMethylationEPICanno.ilm10b4.hg19'))


View(gometh_up)
View(gometh_down)


plt <- gometh_up |> 
  dplyr::rename_with(~paste0(.x,"_up")) |> 
  tibble::rownames_to_column('hsa_id') |> 
  dplyr::left_join(
    gometh_down |> 
      dplyr::rename_with(~paste0(.x,"_down")) |> 
      tibble::rownames_to_column('hsa_id'),
    by=c('hsa_id'='hsa_id')
  )


ggplot(plt, aes(x=-log10(P.DE_up),
                y=-log10(P.DE_down)
                )
       )  + 
  geom_point(pch=16)



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

