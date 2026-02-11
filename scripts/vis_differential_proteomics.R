#!/usr/bin/env

# libs ----


library(limma)
library(ggplot2)
library(patchwork)

source('scripts/load_palette.R')
source('scripts/load_themes.R')
source('scripts/load_gene_annotations.R')



# load data ----


if(!exists('metadata.proteomics.glass_od')) {
  source('scripts/load_GLASS-OD_proteomics.R')
}


if(!exists('glass_od.metadata.proteomics')) {
  source('scripts/load_GLASS-OD_proteomics.R')
}


if(!exists('glass_nl.proteomics')) {
  source('scripts/load_GLASS-NL_proteomics.R')
}



## plt tmp ----
### cell cycling ----


plt <- stats.time |> 
  dplyr::rename_with( ~ paste0(.x, "_time"), .cols=!matches("^protein_id$",perl = T)) |> 
  dplyr::full_join(stats.grade |> dplyr::rename_with( ~ paste0(.x, "_grade"), .cols=!matches("^protein_id$",perl = T)), by=c('protein_id'='protein_id')) |> 
  dplyr::full_join(stats.cgc |> dplyr::rename_with( ~ paste0(.x, "_cgc"), .cols=!matches("^protein_id$",perl = T)), by=c('protein_id'='protein_id')) |> 
  
  dplyr::full_join(stats.time.pat.corrected |> tibble::rownames_to_column("protein_id") |> dplyr::rename_with( ~ paste0(.x, "_time.cor"), .cols=!matches("^protein_id$",perl = T)), by=c('protein_id'='protein_id')) |> 
  dplyr::full_join(stats.grade.pat.corrected |> tibble::rownames_to_column("protein_id") |> dplyr::rename_with( ~ paste0(.x, "_grade.cor"), .cols=!matches("^protein_id$",perl = T)), by=c('protein_id'='protein_id')) |>  
  
  dplyr::mutate(cellcycling_go_ = protein_id %in% clng) |> 
  #dplyr::mutate(fn1 = grepl("FN1|FN|CIG|FINC|FIBRO", protein_id))
  dplyr::mutate(fn1 = grepl("FINC", protein_id))


#plt |> dplyr::filter(fn1) |> head()


#saveRDS(plt, file="cache/analysis_differential_proteomics__plt.Rds")


plt <- readRDS(file="cache/analysis_differential_proteomics__plt.Rds")





plt <- plt |> dplyr::mutate(label = cellcycling_go_)
ggplot(plt, aes(x=t_time, y=t_grade, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  theme_cellpress



plt <- plt |> dplyr::mutate(label = cellcycling_go_)
ggplot(plt, aes(x=t_time, y=t_time.cor, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  theme_cellpress


plt <- plt |> dplyr::mutate(label = cellcycling_go_)
ggplot(plt, aes(x=t_grade, y=t_grade.cor, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  theme_cellpress



pltcor <- plt |> 
  dplyr::select(starts_with("t")) |> 
  dplyr::filter(!is.na(t_time)) |> 
  dplyr::filter(!is.na(t_grade)) |> 
  dplyr::filter(!is.na(t_cgc))|> 
  dplyr::filter(!is.na(t_time.cor)) |> 
  dplyr::filter(!is.na(t_grade.cor))
#dplyr::filter(!is.na(t_cgc.cor))

corrplot::corrplot(cor(pltcor), order="hclust")


sum((plt |> dplyr::filter(!is.na(adj.P.Val_time)) |> dplyr::pull(adj.P.Val_time)) < 0.01)
sum((plt |> dplyr::filter(!is.na(adj.P.Val_time.cor)) |> dplyr::pull(adj.P.Val_time.cor)) < 0.01)
sum((plt |> dplyr::filter(!is.na(adj.P.Val_grade)) |> dplyr::pull(adj.P.Val_grade)) < 0.01)
sum((plt |> dplyr::filter(!is.na(adj.P.Val_grade.cor)) |> dplyr::pull(adj.P.Val_grade.cor)) < 0.01)
sum((plt |> dplyr::filter(!is.na(adj.P.Val_cgc)) |> dplyr::pull(adj.P.Val_cgc)) < 0.01)
sum((plt |> dplyr::filter(!is.na(adj.P.Val_cgc.cor)) |> dplyr::pull(adj.P.Val_cgc.cor)) < 0.01)




plt <- plt |> dplyr::mutate(label = protein_id %in% c("IDH1"))
p1 = ggplot(plt, aes(x=t, y=t.grade, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  ggrepel::geom_text_repel(data=subset(plt,  label == T), size=3, col="red") + 
  labs(x = "t-score CGC (continuous)", y="t-score WHO grade 2 - 3") + 
  theme_cellpress

plt <- plt |> dplyr::mutate(label = protein_id %in% c("IDH2"))
p2 = ggplot(plt, aes(x=t, y=t.grade, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  ggrepel::geom_text_repel(data=subset(plt,  label == T), size=3, col="red")  + 
  labs(x = "t-score CGC (continuous)", y="t-score WHO grade 2 - 3") + 
  theme_cellpress

plt <- plt |> dplyr::mutate(label = protein_id %in% c("GLUD1"))
p3 = ggplot(plt, aes(x=t, y=t.grade, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  ggrepel::geom_text_repel(data=subset(plt,  label == T), size=3, col="red") + 
  labs(x = "t-score CGC (continuous)", y="t-score WHO grade 2 - 3") +
  theme_cellpress


plt <- plt |> dplyr::mutate(label = protein_id %in% c("FINC"))
p4 = ggplot(plt, aes(x=t, y=t.grade, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  ggrepel::geom_text_repel(data=subset(plt,  label == T), size=3, col="red") + 
  labs(x = "t-score CGC (continuous)", y="t-score WHO grade 2 - 3") +
  theme_cellpress




p1 + p2 + p3 + p4

### corr t's ----


plt.c <- plt |> 
  dplyr::select(starts_with("t_")) |> 
  dplyr::filter(!is.na(t_time)) |> 
  dplyr::mutate(t_time = NULL ) |> 
  dplyr::mutate(t_grade = NULL) |> 
  cor()


corrplot::corrplot(plt.c, order="hclust") 





### FN1 ----



ggplot(plt, aes(x=t_cgc, y=t_time.cor, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_point(data=subset(plt, fn1 == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, fn1 == T), pch=19, cex=0.6, col="red") +
  ggrepel::geom_text_repel(data=subset(plt,  fn1 == T), size=3, col="red") + 
  labs(x = "t-score CGC (continuous)", y="t-score WHO grade 2 - 3") +
  theme_cellpress



#### top100 -----


top100 <- plt |> 
  dplyr::filter(!is.na(protein_id)) |> 
  dplyr::filter(adj.P.Val_cgc.cor < 0.01) |> 
  dplyr::arrange(abs(adj.P.Val_cgc.cor)) |> 
  dplyr::select(protein_id, adj.P.Val_cgc.cor) |> 
  head(n=100) |> 
  dplyr::pull(protein_id) |> 
  unique()


plt <- plt |> dplyr::mutate(label = protein_id %in% top100)
ggplot(plt, aes(x=t_time.cor, y=t_cgc.cor, label= protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_cellpress_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.15) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  ggrepel::geom_text_repel(data=subset(plt,  label == T), size=3, col="red") + 
  labs(x = "t-score CGC (continuous)", y="t-score WHO grade 2 - 3") +
  theme_cellpress


#### enhanced volcano ----



EnhancedVolcano(plt,
                lab = plt$protein_id,
                FCcutoff = 0.5,
                x = 'logFC_cgc.cor',
                y = 'adj.P.Val_cgc.cor')




plt <- plt |> dplyr::mutate(significant = adj.P.Val_cgc.cor < 0.01 & abs(logFC_cgc.cor) > lfc_cut)
ggplot(plt, aes(x= `logFC_cgc.cor`, y= -log(`adj.P.Val_cgc.cor`), col=significant, label=protein_id)) +
  
  geom_hline(yintercept=0,          col="black", lwd=theme_cellpress_lwd, lty=1) +
  geom_hline(yintercept=-log(0.01), col="black", lwd=theme_cellpress_lwd, lty=2) +
  #geom_hline(yintercept=-log(0.05), col="gray", lwd=theme_cellpress_lwd) +
  
  geom_vline(xintercept=0,        col="black", lwd=theme_cellpress_lwd, lty=2) +
  geom_vline(xintercept=-lfc_cut, col="black", lwd=theme_cellpress_lwd, lty=2) +
  geom_vline(xintercept= lfc_cut, col="black", lwd=theme_cellpress_lwd, lty=2) +
  
  geom_point(pch=19, cex=0.5) +
  ggrepel::geom_text_repel(data=subset(plt, significant==T), size=theme_cellpress_size, col="black") +
  scale_color_manual(values=c(`TRUE`='red',`FALSE`='darkgray')) +
  theme_cellpress




lfc_cut <- log2(1.25)
plt <- plt |> dplyr::mutate(significant = adj.P.Val_cgc.cor < 0.01 & abs(logFC_cgc.cor) > lfc_cut)
ggplot(plt, aes(x= `logFC_cgc.cor`, y= -log(`adj.P.Val_cgc.cor`), col=significant, label=protein_id)) +
  
  geom_hline(yintercept=0,          col="black", lwd=theme_cellpress_lwd, lty=1) +
  geom_hline(yintercept=-log(0.01), col="black", lwd=theme_cellpress_lwd, lty=2) +
  #geom_hline(yintercept=-log(0.05), col="gray", lwd=theme_cellpress_lwd) +
  
  geom_vline(xintercept=0,        col="black", lwd=theme_cellpress_lwd, lty=2) +
  geom_vline(xintercept=-lfc_cut, col="black", lwd=theme_cellpress_lwd, lty=2) +
  geom_vline(xintercept= lfc_cut, col="black", lwd=theme_cellpress_lwd, lty=2) +
  
  geom_point(pch=19, cex=0.5) +
  ggrepel::geom_text_repel(data=subset(plt, protein_id %in% c("GLUD2", "IDH1","IDH2", "GLUL", "GLUD1", "ALDH18A1", "GLS", "PC", "GLS2", "GOT1L1", "GOT1", "CPS1", "GOT2",
                                                              "IDH3B",  "IDH3A", "IDH3G", "LDHA", "LDHB", "OGDH")), size=theme_cellpress_size, col="black") +
  scale_color_manual(values=c(`TRUE`='red',`FALSE`='darkgray')) +
  theme_cellpress





lfc_cut <- log2(1.25)
plt <- plt |> dplyr::mutate(significant = adj.P.Val_cgc.cor < 0.01 & abs(logFC_cgc.cor) > lfc_cut)
ggplot(plt, aes(x= `logFC_cgc.cor`, y= -log(`adj.P.Val_cgc.cor`), col=significant, label=protein_id)) +
  
  geom_hline(yintercept=0,          col="black", lwd=theme_cellpress_lwd, lty=1) +
  geom_hline(yintercept=-log(0.01), col="black", lwd=theme_cellpress_lwd, lty=2) +
  #geom_hline(yintercept=-log(0.05), col="gray", lwd=theme_cellpress_lwd) +
  
  geom_vline(xintercept=0,        col="black", lwd=theme_cellpress_lwd, lty=2) +
  geom_vline(xintercept=-lfc_cut, col="black", lwd=theme_cellpress_lwd, lty=2) +
  geom_vline(xintercept= lfc_cut, col="black", lwd=theme_cellpress_lwd, lty=2) +
  
  geom_point(pch=19, cex=0.5) +
  ggrepel::geom_text_repel(data=subset(plt, protein_id %in% c("DNMT1",  "DNMT3A", "APOBEC3C", "APOBEC3B", "APOBEC2"
                                                              # no TET's
                                                              # no TDG's
  )), size=theme_cellpress_size, col="black") +
  scale_color_manual(values=c(`TRUE`='red',`FALSE`='darkgray')) +
  theme_cellpress




plt |> dplyr::filter(grepl("BER", protein_id)) |> dplyr::pull(protein_id)



plt |> dplyr::filter(significant) |> dplyr::arrange(adj.P.Val_cgc.cor) |>  dplyr::pull(protein_id)



### HOX + polycomb ----


selected_genes <- c(genes_homeobox_nonhox,
                    genes_polycomb_eed_homeobox,
                    genes_polycomb_h3k27_homeobox,
                    genes_polycomb_prc2_homeobox,
                    genes_polycomb_suz12_homeobox)


metadata.proteomics.glass_od |>
  dplyr::rowwise() |> 
  dplyr::mutate(
    is_selected1 = stringr::str_split(Genes, ";") |>
      unlist() |>
      (\(x) any(x %in% selected_genes))(),
    is_selected2 = stringr::str_split(Protein.Names, ";") |>
      unlist() |>
      (\(x) any(gsub("_HUMAN", "", x) %in% selected_genes))()
                ) |> 
  dplyr::mutate(is_selected = is_selected1 | is_selected2) |> 
  dplyr::ungroup() |> 
  dplyr::select(Genes, is_selected) |>
  dplyr::pull(is_selected) |> 
  table()



lfc_cut <- 0.5

plt <- metadata.proteomics.glass_od |>
  dplyr::rowwise() |> 
  dplyr::mutate(
    is_selected1 = stringr::str_split(Genes, ";") |>
      unlist() |>
      (\(x) any(x %in% selected_genes))(),
    is_selected2 = stringr::str_split(Protein.Names, ";") |>
      unlist() |>
      (\(x) any(gsub("_HUMAN", "", x) %in% selected_genes))()
  ) |> 
  dplyr::mutate(label = is_selected1 | is_selected2) |> 
  dplyr::ungroup() |> 
  
  dplyr::filter(!is.na(DPA__GLASS_OD__CGC__t)) |> 
  dplyr::mutate(significant = DPA__GLASS_OD__CGC__adj.P.Val < 0.01 & abs(DPA__GLASS_OD__CGC__logFC) > lfc_cut)




ggplot(plt, aes(x= `DPA__GLASS_OD__CGC__logFC`, y= -log(`DPA__GLASS_OD__CGC__adj.P.Val`), col=label, label=Genes)) +
  
  geom_point(data=subset(plt, label==F), size=theme_cellpress_lwd/3) +
  
  geom_hline(yintercept=0,          col="black", lwd=theme_cellpress_lwd, lty=1) +
  geom_hline(yintercept=-log(0.01), col="black", lwd=theme_cellpress_lwd, lty=2) +

  geom_vline(xintercept=0,        col="black", lwd=theme_cellpress_lwd, lty=2) +
  geom_vline(xintercept=-lfc_cut, col="black", lwd=theme_cellpress_lwd, lty=2) +
  geom_vline(xintercept= lfc_cut, col="black", lwd=theme_cellpress_lwd, lty=2) +
  
  ggrepel::geom_text_repel(data=subset(plt, label == T),
                           max.overlaps = Inf,
                           #box.padding = 0.25,
                           #nudge_x = 0.135,
                           nudge_y = 7.5,
                           segment.color = "darkgreen",
                           col="black",
                           size=theme_cellpress_size, segment.size=theme_cellpress_lwd, family = theme_cellpress_font_family
                           ) +
  
  geom_point(data=subset(plt, label==T), size=theme_cellpress_lwd/1.5) +
  
  scale_color_manual(values=c(`TRUE`='red',`FALSE`='darkgray')) +
  labs(subtitle = "Differential protein analysis: GLASS-OD ~ CGC, colored by Homeobox TF's", col="Homeobox TF", x="log2FC", y="-log(adj. P-value)") +
  theme_cellpress +
  theme_cellpress + theme(plot.background = element_rect(fill="white", colour=NA)) 

ggsave("output/figures/analysis_differential_proteomics__homeobox_tfs_2.pdf", width = 3.5, height=2)
ggsave("output/figures/analysis_differential_proteomics__homeobox_tfs_2.png", width = 3.5, height=2)




### PCNA / Ki67 (MKI67) ----



plt |> dplyr::filter(grepl("^HMGB", protein_id)) |> dplyr::pull(protein_id)



lfc_cut <- log2(1.25)
plt <- plt |> dplyr::mutate(significant = adj.P.Val_cgc.cor < 0.01 & abs(logFC_cgc.cor) > lfc_cut)
ggplot(plt, aes(x= `logFC_cgc.cor`, y= -log(`adj.P.Val_cgc.cor`), col=significant, label=protein_id)) +
  
  geom_hline(yintercept=0,          col="black", lwd=theme_cellpress_lwd, lty=1) +
  geom_hline(yintercept=-log(0.01), col="black", lwd=theme_cellpress_lwd, lty=2) +
  #geom_hline(yintercept=-log(0.05), col="gray", lwd=theme_cellpress_lwd) +
  
  geom_vline(xintercept=0,        col="black", lwd=theme_cellpress_lwd, lty=2) +
  geom_vline(xintercept=-lfc_cut, col="black", lwd=theme_cellpress_lwd, lty=2) +
  geom_vline(xintercept= lfc_cut, col="black", lwd=theme_cellpress_lwd, lty=2) +
  
  geom_point(pch=19, cex=0.5) +
  ggrepel::geom_text_repel(data=subset(plt, protein_id %in% c(
    "PCNA", "PNCA", "MKI67", "CDC6", "PLK1", "MCM3" ,  "MCM4" ,  "MCM5" ,  "MCM7"  , "MCM2" ,  "MCM6",
    
    "TMPO","DEK","RFC4",
    "HMGB3", "HMGB1", "HMGB2"
    
    #"H6PD", "H2AC25;H2AC6;H2AC8", "H2BC11;H2BC21", "H1-0", "H1-4", "H1-5", "H1-2", "H2BC5", "H4C16", "H3-3B;H3C12;H3C13", "H1-1", "H2AC19;H2AC20", "H3-7", "H2AZ2", "H2AC21", "H1-10", "H2BC1", "H1-6", "H2BC19P;H2BC20P"
  )), size=theme_cellpress_size, col="black") +
  scale_color_manual(values=c(`TRUE`='red',`FALSE`='darkgray')) +
  theme_cellpress



### recursive Cor N/A adj ----


lfc_cut <- 0.5

plt <- metadata.proteomics.glass_od |>
  dplyr::mutate(significant = DPA__GLASS_OD__CGC__adj.P.Val < 0.01 &
                  abs(`DPA__GLASS_OD__CGC__logFC`) > lfc_cut) 
table(plt$significant)


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, starts_with("proteomics_")) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 137) 
    return(.)
  })()




data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id) |> 
  tibble::rownames_to_column("protein_id") |> 
  dplyr::filter(protein_id %in% (plt |> dplyr::filter(significant) |> dplyr::pull(Genes))) |> 
  tibble::column_to_rownames("protein_id")




dim(data)



cordata <- metadata |> 
  dplyr::inner_join(
    glass_od.metadata.array_samples |> 
      filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
      dplyr::filter(resection_id %in% metadata$resection_id) |> 
      assertr::verify(!duplicated(resection_id)),
    by=c('resection_id'='resection_id')
  ) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 118) 
    return(.)
  })() |> 
  dplyr::select(resection_id, proteomics_id, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::left_join(
    as.data.frame(t(data)) |> tibble::rownames_to_column('proteomics_id'), by=c('proteomics_id'='proteomics_id')
  ) |> 
  tibble::column_to_rownames('resection_id') |> 
  dplyr::mutate(proteomics_id = NULL)


cordata <- cor(cordata, use="pairwise.complete.obs", method="pearson") |> 
  as.data.frame() |> 
  dplyr::select(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  tibble::rownames_to_column('protein_id') |> 
  dplyr::filter(protein_id != "array_A_IDH_HG__A_IDH_LG_lr__lasso_fit") |> 
  dplyr::rename(`cor CGC` = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)




labels <- data |>
  tibble::rownames_to_column("gene_id") |> 
  dplyr::select("gene_id") |> 
  dplyr::mutate(FN1 = gene_id == "FINC") |> 
  dplyr::mutate(cycling = gene_id %in% c(
    "PCNA", "PNCA", "MKI67", "CDC6", "PLK1", "MCM3" ,  "MCM4" ,  "MCM5" ,  "MCM7"  , "MCM2" ,  "MCM6",
    "TMPO","DEK","RFC4",
    "HMGB3", "HMGB1", "HMGB2"
  )) |> 
  dplyr::mutate(GO_0002526_acute_inflammatory_response = gene_id %in% c(
    "BTK","CD6","PTGER3","ITIH4","CREB3L3","TFRC","FCGR2B","PTGS2","GSTP1","B4GALT1","EIF2AK1","OSM","RHBDD3","PIK3CG","EPHB6","TFR2","LIPA","VNN1","IL4","IL1A","ACVR1","FN1","PARK7","ASH1L","PLA2G2D","F3","TNFSF4","CNR1","TNFSF11","IL1B","C3","FFAR2","CCR7","IL22","APOL2","LBP","EPO","ASS1","F12","SELENOS","PPARG","CRP","APCS","ALOX5AP","SAA2","IL6ST","EDNRB","IL6","IL1RN","PRCP","MYLK3","TNFRSF11A","REG3G","AHSG","OSMR","PTGES","SAA4","FCGR1A","ADAM8","IL6R","NLRP3","ADORA1","DNASE1L3","NPY5R","KLKB1","IL31RA","MBL2","SERPINF2","SCN11A","MRGPRX1","REG3A","CEBPB","SAA1","IL20RB","NLRP6","A2M","NUPR1","ANO6","CD163","CTNNBIP1","FCER1A","F2","FUT7","CXCR2","EXT1","F8","SIGIRR","FFAR3","PLSCR1","ZP3","SERPINA3","TRPV1","SERPINA1","SPN","ELANE","C2CD4A","FCGR3A","HLA-E","C2CD4B","IGHG1","DNASE1","ORM2","ORM1","TNF","UGT1A1","INS","SAA2-SAA4","HP"
  )) |>
  dplyr::mutate(GO_0062023_collagen_containing_extracellular_matrix = gene_id %in% c(
    "DCN","SEMA3B","MARCO","SERPINB1","MYOC","TIMP2","VCAN","TNC","USH2A","COL9A2","LTBP1","ELN","LAMC3","COL23A1","LAMA3","ENTPD2","CBLN4","ITIH4","ITIH1","F7","LAMC2","COL11A1","WNT8A","GPC1","CDON","NTN1","COL17A1","FGFR2","PKM","FGF10","SPP2","ADAM11","NTN4","LMAN1","DLG1","GPC4","FBLN1","COL5A3","HSP90AA1","COL4A4","IMPG2","COL19A1","COL16A1","FCN1","ACHE","ADAMTS2","MMP2","NID2","LTBP4","ICAM1","P3H2","LAMB4","LAMB1","APOH","ANGPT2","CCDC80","CMA1","COL9A3","TGFB2","AMELY","HNRNPM","LGALS1","TIMP3","PDGFB","CHADL","CTSG","COCH","MMP9","BMP7","CTSZ","COL20A1","ANGPT4","LAMA1","MXRA5","SRPX","F9","TIMP1","SRPX2","CBLN1","ZP2","CTSH","SFRP1","IL7","FGL1","CLC","APLP1","TGFB1","COMP","WNT2","PTPRZ1","PCOLCE","SERPINE1","PLOD3","RARRES2","AEBP1","MEGF9","OGN",
    "ASPN","ECM2","AMBP","CXCL12","ACTA2","KAZALD1","LGALS3BP","COL1A1","VTN","SOD3","CTSC","HPX","APOA4","APOC3","MDK","VWF","WNT5B","MGP","COL12A1","COL9A1","SMOC2","IMPG1","LAMA4","ERBIN","LOX","SPARC","THBS4","KNG1","HRG","WNT5A","COL7A1","LOXL3","EFEMP1","FN1","TNR","ANGPTL1","PRG4","NID1","ZP4","MFAP2","P3H1","F3","SERPINC1","CTSD","MMP8","APOA1","CCN2","LTBP2","TGFB3","TNN","TGFBI",
    "CLU","A1BG","FMOD","PLG","ANXA11","COL10A1","SLPI","MATN4","NDP","F13A1","SERPINB6","COL21A1","AMELX","CFP","PZP","FGL2","LOXL1","APOE","NCAN","PXDN","GDF15","COL5A1","LAMA5","FIBCD1","ANGPTL6","F12","LGALS3","PODNL1","RTBDN","MATN3","EMILIN2","SERPINF1","ITGB4","MATN2","BCAN","HAPLN2","APCS","ANGPTL3","POSTN","LOXL2","ADAMDEC1","WNT2B","COL4A2","ADAMTS8","ANXA1","CTSL","ADAM19","AGT","LAMC1","SERPINE2","ANGPTL2","CCN3","IGFBPL1","TINAG","SULF1","THBS1",
    "EMILIN1","LOXL4","ANXA7","CILP","SEMA7A","MMRN1","FRAS1","FBN2","COL2A1","INHBE","LUM","FBLN5","PCSK6","LMAN1L","HAPLN3","MFGE8","TGFB1I1","CDH13","COL6A1","COL6A2","ADAMTS10","FCN3","HSPG2","CCN1","TINAGL1","DPT","HDGF","HMCN1","ECM1","ADAMTSL4","ANXA9","S100A8","S100A7","FLG","LEFTY2","COL8A1","AHSG","SFRP2","HAPLN1","CASK","GPC3","HMCN2","SERPING1","SERPINH1","NCAM1","ZP1","FREM2","ITIH2","DST","SPARCL1","ABI3BP","ANGPT1","ADAMTS1","ADAMTS5","ADAMTS3","PRG3","ACAN","COLEC12","ADAMTS4","LAD1","C1QC","CSTB","FCN2","AZGP1","COL26A1","MATN1","SDC3","CTSS","S100A9","COL6A3","IGFBP7","FBLN2","ADAMTS9","CRELD1","PF4","CPA3","TGM4","CLEC3B","ANXA5","EDIL3","EGFLAM","RELL2","SHH","COL1A2","CTSB","SBSPON","CTHRC1","FREM1","VWA2","HTRA1","ADAMTS15","FBN1","SERPINB8","MFAP4","HSP90B1","SERPINB12","NAV2","GREM1","SERPINF2","KRT1","ANGPTL4",
    "SOST","LTBP3","TNXB","COL3A1","NPNT","LGALS9","COL4A3","THBS3","COL22A1","SDC2","MUC17","SERPINB9","CDH2","COL24A1","FGG","FGA","FGB","LGALS4","COL8A2","ANGPTL7","LAMB2","LRRC15","TPSAB1","EFEMP2","COL6A5","EGFL7","ADAMTS20","MMRN2","C1QA","DAG1","CSPG4","CTSF","PODN","FGFBP3","ZG16","NPPA","A2M","CLEC14A","CD151","CALR","GPC5","VWA1","SSC5D","F2","ADIPOQ","OTOL1","BGN","MXRA7","SNORC","ANXA2","COL18A1","FREM3","GPC6","EMILIN3","EFNA5","THBS2","PRG2","C17orf58","ANGPTL5","BCAM","COL4A1","HAPLN4","AMTN","THSD4","COL14A1","EYS","COL4A5","AGRN","PLSCR1","ZP3","IFNA2","SERPINA5","COL25A1","VWC2","PRELP","MMP23B","SERPINA3","S100A4","PRTN3","LAMA2","COL27A1","LAMB3","ANXA4","ANXA6","SERPINA1","TPSB2","COL13A1","SPN","ELANE","COL4A6","MFAP5","PSAP","S100A10","HRNR","S100A6","SMOC1","EGFL6","L1CAM","TGM2","COL11A2","COL5A2","COL15A1","PRSS1","VIT","DEFA1","COL6A6","COLQ","GPC2","ANG","COL28A1","ORM2","ORM1","DEFA1B","MARCOL","GH1","SPON1","ANXA8","RBP3","GDF10","MMP28"
  )) |> 
  dplyr::mutate(`GO_0072562_blood_microparticle` = gene_id %in% c(
    "CFH","SLC4A1","PON1","ITGA2B","CP","ITIH4","ITIH1","GRIPAP1","TFRC","CD5L","ACTB","AFM","PSMC5","TF","APOL1","TGFB1","AMBP","ENG","PFN1","LGALS3BP","VTN","HSPA8","HPX","APOA4","C9","KNG1","HRG","BCHE","FN1","CFHR3","SLC2A1","SERPINC1","APOA1","CLU","A1BG","PLG","C4BPA","DNPEP","F13A1","C3","HSPA2","PZP","INTS11","APOE","JCHAIN","APCS","AGT","CIB2","SDCBP","TMPRSS13","FCN3","OAZ3","ACTA1","EIF2A","AHSG","GC","MSN","STOM","GSN","SERPING1","ITIH2","C8A","APOA2","C1QC","ACTC1","FCN2","ACTG2","ALB","ANXA5","YWHAZ","ACSM1","SERPINF2","KRT1","ANGPTL4","FGG","FGA","FGB","HSPA6","A2M","C8G","ZBTB38","CPN2","F2","C1S","ACTG1","PROS1","KDM4D","POTEE","HBA2","ZNF177","SERPINA3","HBG2","POTEF","HSPA1B","HSPA1A","HSPA1L","PRSS1","HBA1","IGKV4-1","IGLV1-47","IGLV3-25","IGLV3-21","IGLC2","IGLC3","IGHA2","IGHG4","IGHG2","IGHA1","IGHG1","IGHM","IGHV3-7","IGHV3-13","IGHV3-23","CLIC1","HBE1","HBD","C4B","ORM2","ORM1","IGKV3-20","IGKV1D-33","IGKV1-17","IGKV3-11","IGKV1-33","IGKV1-39","IGKV2D-28","IGKV2-30","IGKV1-5","CFB","CFHR1","IGKV3-15","C4A","HBB","IGKV2D-40","HP","HPR","IGKV1D-12"
  )) |>
  dplyr::mutate(GO_0002250_adaptive_immune_response = c(
    gene_id %in% c(
      "TMEM98", "CD79B", "MASP2", "CD4", "BTK", "HFE", "FYN", "ERCC1", "KDM5D", "CD6", "UFL1", "WAS", "SLC11A1", "CD74", "C8B", "SLAMF7", "BTN3A1", "TNFRSF1B", "C6", "PARP3", "TNFRSF17", "FOXP3", "RC3H2", "TRAF3IP2", "PRDM1", "PRKCQ", "CD84", "TP53BP1", "PRKCZ", "TFE3", "BCL3", "RORA", "RAB27A", "TFRC", "FCGR2B", "PVR", "TBX21", "TXK", "MLH1", "UNG", "PAG1", "RAP1GAP", "NFKB2", "IL4R", "LAMP3", "RAPGEF3", "CEACAM1", "STX7", "RIF1", "ARG2", "MEF2C", "PTPRC", "LAT2", "P2RX7", "LAG3", "ICAM1", "CD209", "MCOLN1", "RAPGEF4", "UNC13D", "MSH2", "TREM2", "SIRT1", "JAK2", "IL12RB1", "KLHL22", "LGALS1", "CSF2RB", "IL2RB", "EP300", "CD40", "SLA2", "JAG1", "SMAD7", "RNF125", "TLR8", "BMX", "CD40LG", "TNFSF13B", "PYCARD", "IL21R", "CSK", "CTSH", "RIPK2", "NBN", "FGL1", "RELB", "FCER2", "CLEC4M", "IL4I1", "LILRA1", "IL27RA", "CLC", "EBI3", "PRKD2", "TGFB1", "CD79A", "TYK2", "LILRB5", "JAK3", "PIK3CG", "EPHB6", "AHR", "C5", "EXOSC3", "GATA3", "PPP3CB", "C1QBP", "SUPT6H", "IL2", "NSD2", "CTSC", "CRTAM", "HSPA8", "UNC93B1", "KMT5B", "HPX", "CD81", "TCIRG1", "CD69",
      "IL23A", "IFNG", "PTPN6", "CLEC4A", "AICDA", "BTN3A3", "ULBP1", "IL17A", "IL17F", "RNF8", "BACH2", "CCR6", "TFEB", "DUSP22", "C7", "ITK", "IL12B", "BTNL8", "IL4", "C9", "BCL6", "CD86", "ZAP70", "LOXL3", "IL1R1", "IL1RL1", "IL18R1", "MSH6", "RNF19B", "MAD2L2", "PLA2G4A", "SLAMF1", "CD160", "CR2", "CD46", "MPL", "TNFSF4", "TNFAIP3", "ARG1", "CD274", "IFNA6", "IFNA8", "TNFSF18", "CLU", "PTK2B", "ADCY7", "KCNJ8", "CD80", "CCR2", "SASH3", "LAX1", "CD244", "LY9", "SHLD2", "PMS2", "PKN1", "NCKAP1L", "IL13RA2", "C4BPA", "C4BPB", "IL9R", "TREM1", "IRF1", "IL1B", "CD70", "C3", "STAT5A", "ZBTB1", "TRAF2", "IGLL1", "EIF2AK4", "RIPK3", "FOXJ1", "CRACR2A", "NECTIN2", "TRPM4", "ULBP2", "ULBP3", "LILRB2", "PRR7", "RFTN1", "NDFIP1", "PRKAA1", "RAP1GAP2", "JCHAIN", "CLEC10A", "CTNNBL1", "EPHB2", "KMT5C", "SWAP70", "VTCN1", "RSAD2", "IL6ST", "KLRD1", "KLRC1", "ANXA1", "CTSL", "HAVCR2", "AKIRIN2", "MAP3K7", "TEC", "RC3H1", "IL6", "DBNL", "MYO1G", "HLX", "IL10", "TLR4", "IL33", "SIT1", "IFNA21", "IRF4", "IL18BP", "STAT4", "LEF1", "C1RL", "CD27", "SLC15A4", "SKAP1", "SECTM1", "TNFRSF11A", "BRD4", "SIGLEC10", "BCL10", "XCL1", "FCGR2A", "RORC", "SUSD4", "DUSP10", "HSPD1", "NFKBIZ", "IL9", "FBXO38", "TNFRSF21", "EBAG9", "IFNA5", "IFNA16", "IFNK", "NOTCH1", "SERPING1", "FCGR1A", "CD226", "IL18", "ADAM17", "DCLRE1C", "CAMK4", "MR1", "CYRIB",
      "CD8A", "MCOLN2", "SAMSN1", "RAET1L", "BATF", "CXCL13", "C8A", "PAXIP1", "TNFRSF14", "CD1D", "CD1A", "CD1C", "CD1B", "CD1E", "HMHB1", "FCER1G", "C1QC", "C1R", "TNFRSF13C", "ICOSLG", "AIRE", "CD3G", "ZBTB7B", "IL6R", "PCYT1A", "ALOX15", "TNFSF13", "JAK1", "IL23R", "NLRP3", "SLAMF6", "FCGR3B", "FCMR", "FCAMR", "SANBR", "CTSS", "FZD5", "EOMES", "FCRL4", "TRAT1", "CTLA4", "ZC3H12A", "PRKCD", "RNF168", "ERAP1", "ERAP2", "RAET1E", "SYK", "SVEP1", "MARCHF8", "MBL2", "TSC1", "HPRT1", "JAM3", "C2", "RAG1", "PRKCB", "CLEC4D", "B2M", "STAT6", "PDIA3", "PHB1", "NOD2", "CD3D", "NFKBID", "LAIR1", "OTUB1", "FADD", "CX3CR1", "TAP1", "STAT3", "IL7R", "IL12A", "INPP5D", "GPR183", "APLF", "ALCAM", "YWHAG", "SOCS5", "JUNB", "FGA", "FGB", "PIK3CD", "IFNB1", "SHLD1", "CD8B", "MALT1", "CLEC7A", "KLHL6", "THEMIS", "CCL19", "MYD88", "SLC22A13", "C1QA", "CD7", "EXO1", "LIG4", "IL20RB", "HRAS", "ADGRE1", "TRAF6", "CLCF1", "ZNF683", "ATAD5", "C8G", "IFNW1", "CD19", "IL17RA", "CD28", "HLA-DQB1", "FCER1A", "FUT7", "PRF1", "NLRP10", "C1S", "CLEC4G", "HMCES", "PSG9", "ASCL2", "SH2D1A", "SOCS3", "IFNE", "IRF7", "BTLA", "BTN3A2", "IFNA10", "LILRB4", "CARD9", "C17orf99", "NCR3LG1", "ZP3", "IFNA2", "PDCD1", "IGHV1OR15-9", "HMGB1", "HLA-DRB1", "SEMA4A", "CD55", "HLA-DQA1", "ADA", "ARID5A", "IL27", "SPN", "GZMM", "MPEG1", "PDCD1LG2", "CR1L", "IFNA1", "ENTPD7", "CLEC4C", "HLA-DRB5", "SH2D1B", "MTOR", "CD247", "OPA1", "CD3E", "CR1", "RAET1G", "FCGR3A", "LIME1", "HLA-DOA", "BRD2", "HLA-DMA", "TAP2", "HLA-DRA", "AGER", "MICB", "HLA-C", "LILRB3", "GNL1", "HLA-E", "HLA-G", "HLA-F", "TRIM27", "CFI", "KLRC2", "CLEC6A", "HLA-A", "IGKJ1", "IGKV4-1", "IGKV5-2", "IGKV6-21", "IGKV2D-26", "IGKV3D-20", "IGKV6D-41", "IGKV3D-11", "IGKV1D-42", "IGLV4-69", "IGLV8-61", "IGLV4-60", "IGLV6-57", "IGLV11-55", "IGLV10-54", "IGLV5-52", "IGLV1-51", "IGLV1-50", "IGLV5-48", "IGLV1-47", "IGLV7-46", "IGLV5-45", "IGLV1-44", "IGLV7-43", "IGLV1-40", "IGLV5-37", "IGLV1-36", "IGLV2-33", "IGLV3-32", "IGLV3-27", "IGLV3-25", "IGLV2-23", "IGLV3-22", "IGLV3-21", "IGLV3-19", "IGLV2-18", "IGLV3-16", "IGLV2-14", "IGLV3-12", "IGLV2-11", "IGLV3-10", "IGLV3-9", "IGLV4-3", "IGLV3-1", "IGLJ1", "IGLC2", "IGLC3", "TRGC1", "TRGJ1", "TRGV11",
      "TRGV10", "TRGV9", "TRGV8", "TRGV5", "TRGV4", "TRGV3", "TRGV1", "TRBV6-1", "TRBV7-1", "TRBV4-1", "TRBV6-4", "TRBV7-3", "TRBV5-3", "TRBV9", "TRBV10-1", "TRBV11-1", "TRBV6-5", "TRBV6-6", "TRBV5-5", "TRBV7-6", "TRBV5-6", "TRBV5-7", "TRBV5-1", "TRBV4-2", "TRBV19", "TRBV20-1", "TRBV23-1", "TRBV24-1", "TRBV27", "TRBV28", "TRBJ2-1", "TRBJ2-3", "TRBJ2-7", "TRAV2", "TRAV3", "TRAV4", "TRAV5", "TRAV6", "TRAV7", "TRAV8-1", "TRAV9-1", "TRAV10", "TRAV12-1", "TRAV8-2", "TRAV8-3", "TRAV13-1", "TRAV12-2", "TRAV8-4", "TRAV13-2", "TRAV14DV4", "TRAV9-2", "TRAV12-3", "TRAV8-6", "TRAV16", "TRAV17", "TRAV18", "TRAV19", "TRAV20", "TRAV21", "TRAV22", "TRAV23DV6", "TRDV1", "TRAV24", "TRAV25", "TRAV26-1", "TRAV27", "TRAV29DV5", "TRAV26-2", "TRAV34", "TRAV36DV7", "TRAV38-1", "TRAV38-2DV8", "TRAV39", "TRAV40", "TRAV41", "TRDV2", "TRDJ1", "TRAJ42", "TRAJ31", "TRAJ3", "IGHA2", "IGHE", "IGHG4", "IGHG2", "IGHA1", "IGHG1", "IGHG3", "IGHM", "IGHJ1", "IGHV6-1", "IGHV1-2", "IGHV1-3", "IGHV2-5", "IGHV3-7", "IGHV3-11", "IGHV3-13", "IGHV3-15", "IGHV3-16", "IGHV1-18", "IGHV3-20", "IGHV3-21", "IGHV3-23", "IGHV1-24", "IGHV2-26", "IGHV4-28", "IGHV3-33", "IGHV4-34", "IGHV3-35", "IGHV3-38", "IGHV4-39", "IGHV1-45", "IGHV1-46", "IGHV3-48", "IGHV3-49", "IGHV5-51", "IGHV3-53", "IGHV1-58", "IGHV4-61", "IGHV3-66", "IGHV1-69", "IGHV2-70D", "IGHV3-73", "IGHV7-81", "DENND1B", "SIPA1", "LAT", "TREX1", "KLRK1", "EMP2", "IFNA7", "SCART1", "IGLV9-49", "EXOSC6", "IGHV3-64", "HLA-DPB1", "TRDD1", "IGKV3D-15", "IGHV4-59", "C4B", "IGHV3-74", "IGKV6D-21", "IGHV3-72", "TRBV2", "LTA", "TRGC2", "IFNA14", "IGKV3D-7", "TRBV10-2", "TRBV5-4", "HLA-DPA1", "IGHV4-31",
      "IGHV3-43", "HLA-DQB2", "TNF", "TRBV29-1", "TRGV2", "IGHV3OR16-10", "IFNA13", "HLA-B", "IFNA17", "IGHD1-1", "IFNA4", "TRBV30", "HLA-DQA2", "TRBV3-1", "IGKV2D-30", "TNFSF12", "IGKV1D-8", "IGKV1-6", "IGKV1-37", "IGKV3-20", "IGKV1D-33", "IGKV1-17", "TNFRSF13B", "IGKV1-8", "IGKV1-16", "HLA-DOB", "IGKV1D-16", "IGKV2-24", "IGKV3-11", "IGKV2D-24", "IGKV1-9", "IGKV1-33", "IGKV1-39", "IGKV2D-28", "HLA-DMB", "IGKV1D-43", "IGKV1D-17", "IGKV3-7", "IGKV2-30", "IGKV2D-29", "IGKV1-12", "IGKV1-5", "CFB", "IGKV2-28", "IGKV3-15", "LILRA6", "IGKV1-27", "C4A", "TARM1", "IGKV1D-37", "IGKV2D-40", "IGKV1D-39", "TRBV6-7", "SHLD3", "TRBV7-7", "TRBV7-4", "TRBV6-8", "LYN", "CD8B2", "IGHV8-51-1", "IGLL5", "TRAV1-1", "TRAV1-2", "TRDV3", "TRAV30", "IGHV4OR15-8", "IGHV2OR16-5", "IGHV3OR15-7", "IGHV3OR16-17", "OTUD7B", "IGHV3OR16-12", "IGHV3OR16-9", "IGHV1OR15-1", "IGHV3-30", "IGHV3OR16-8", "IGHV3OR16-13", "IGKV2-40", "IGHV2-70", "TRBV12-3", "TRBV12-5", "TRBV16", "TRBV14", "TRBV10-3", "ORAI1", "TRBV13", "TRBV18", "IGKV1D-13", "TRBV11-3", "IGHV4-4", "TRBV12-4", "IGHV1OR21-1", "TRBV17", "TRBV7-9", "IGLV2-8", "IGKV1D-12", "IGHV1-69D", "TRBJ1-4", "IGHV1-69-2", "IGHV7-4-1", "TRBJ1-3", "TRBJ1-5", "TRBJ1-1", "TRBJ1-2", "TRBD1", "TRBV25-1", "IGHV3-64D", "IGHV5-10-1", "TRBJ1-6", "TRBV7-2", "TRBV6-2", "IL9R" 
    )
  )) |> 
  dplyr::mutate(FN1 = NULL, cycling = NULL, `GO_0072562_blood_microparticle` = NULL, `GO_0002526_acute_inflammatory_response` = NULL) |> 
  tibble::column_to_rownames("gene_id")


labels_continuous <- data |>
  tibble::rownames_to_column("gene_id") |> 
  dplyr::select("gene_id") |> 
  dplyr::mutate(fraction_NAs = rowSums(is.na(data))) |> 
  dplyr::left_join(cordata, by=c('gene_id'='protein_id')) |> 
  tibble::column_to_rownames("gene_id")




pplt <- data |>
  tibble::rownames_to_column('__hugo_symbol__') |>
  dplyr::filter(!duplicated(.data$`__hugo_symbol__`)) |>
  tibble::column_to_rownames('__hugo_symbol__')


pplt <- pplt |>
  base::as.matrix() |>
  base::t() |> 
  stats::cor(use="pairwise.complete.obs") # copes with N/A values



h <- stats::hclust(stats::as.dist(1 - stats::cor(pplt)), method = "ward.D2" ) # recursive cor-based cluastering !!!


o <- h$labels[h$order] |>
  base::rev()

#write.csv(data.frame(o), file="output/tables/proteomics_rcur.txt")




ph <- ggdendro::ggdendrogram(h, rotate = TRUE, theme_dendro = FALSE) +
  ggdendro::theme_dendro()


pplt <- pplt |>
  base::as.data.frame() |>
  dplyr::select(dplyr::all_of(o)) |>
  base::t() |>
  base::as.data.frame() |>
  dplyr::select(dplyr::all_of(o)) |>
  base::t() |>
  base::as.matrix()


o.join <- base::data.frame(name = o, i = 1:length(o))


plt.expanded2 <- reshape2::melt(pplt) |>
  dplyr::rename(y = .data$`Var1`) |>
  dplyr::rename(x = .data$`Var2`) |>
  dplyr::mutate(x = as.factor(.data$`x`)) |>
  dplyr::mutate(y = as.factor(.data$`y`)) |>
  dplyr::left_join(o.join |> dplyr::rename(x.order = .data$`i`), by = c("x" = "name")) |>
  dplyr::left_join(o.join |> dplyr::mutate(i = dplyr::n() - .data$i + 1) |> dplyr::rename(y.order = .data$i), by = c("y" = "name"))

base::rm(o.join)


font_scale = 4.5
legend_scale = 1.2


p1 <- ggplot2::ggplot(plt.expanded2, ggplot2::aes(
  x = .data$x.order,
  y = .data$y.order,
  radius = ((abs(.data$value) * 0.7) + 0.3) / 2 - 0.05,
  # [0.3 , 0.8] + 0.2 smoothened from lwd/border
  fill = .data$value,
  col = .data$value,
  label = .data$x
)
) +
  ggplot2::geom_tile(col = "gray", fill = "white", lwd = 0.15) +
  ggplot2::scale_fill_gradientn(colours = col2(200), na.value = "grey50", limits = c(-1, 1), guide = "none") + # guide = "colourbar",
  ggplot2::scale_color_gradientn(colours = col2(200), na.value = "grey50", limits = c(-1, 1), guide = "none") +
  recursiveCorPlot::geom_circle(radius.fixed = T) + # from THIS repo
  ggplot2::scale_x_discrete(labels = NULL, breaks = NULL) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.y = ggplot2::element_text(size = font_scale, angle = 0, hjust = 1, vjust = 0.5), # used to be [3,6] reduce font size here, should become argument
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5, color = "gray80"),
    text = ggplot2::element_text(size = 13),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank()
  ) +
  ggplot2::labs(y = NULL, x = NULL, main = NULL) +
  ggplot2::coord_fixed() +
  ggplot2::scale_y_continuous(name = NULL, breaks = base::length(o):1, labels = o)





pplt <- base::data.frame(gid = o, i = 1:length(o)) |>
  dplyr::left_join(labels |> tibble::rownames_to_column("gid"), by = c("gid" = "gid")) |>
  reshape2::melt(id.vars = c("gid", "i")) |>
  dplyr::mutate(variable = factor(.data$variable, levels = base::rev(base::colnames(labels))))


p2 <- ggplot2::ggplot(pplt, ggplot2::aes(
  x = .data$i,
  y = .data$variable,
  fill = .data$value,
  label = .data$gid
)) +
  ggplot2::geom_tile(col = "white", lwd = 0.15) +
  # scale_x_discrete(position = "bottom")  +
  ggplot2::scale_x_discrete(labels = NULL, breaks = NULL) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5),
    axis.text.y = ggplot2::element_text(size = font_scale, angle = 0, hjust = 1, vjust = 0.5),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank()
  ) +
  ggplot2::guides(fill = "none") +
  ggplot2::coord_fixed(ratio = legend_scale) + # used to be 2.75
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray98"))




pplt <- base::data.frame(gid = o, i = 1:length(o)) |>
  dplyr::left_join(labels_continuous |> tibble::rownames_to_column("gid"), by = c("gid" = "gid")) |>
  dplyr::mutate(fraction_NAs = 1.0 * fraction_NAs) |> 
  reshape2::melt(id.vars = c("gid", "i")) |>
  dplyr::mutate(variable = factor(.data$variable, levels = base::rev(base::colnames(labels_continuous))))

pplt <- rbind(
  pplt,
  pplt |> dplyr::mutate(value=0)
)



p3 <- ggplot2::ggplot(pplt, ggplot2::aes(
  x = reorder(gid, .data$i),
  y = .data$value,
  #fill = .data$value,
  label = .data$gid
)) +
  ggplot2::geom_line(col = "black", lwd = 0.15) +
  facet_grid(rows = vars(.data$variable), scales = "free") + 
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    axis.ticks.x=element_blank(), 
    axis.ticks.y=element_blank()
  ) +
  ggplot2::guides(fill = "none") +
  ggplot2::labs(x = NULL, y = NULL) 




patchwork::wrap_plots(
  #A = p3,
  B = p2,
  C = p1
  #D = (ph + patchwork::plot_spacer()),
  #   design = 'A
  # B
  # C'
) + 
  patchwork::plot_annotation(caption = 'caption') +
  patchwork::plot_layout(ncol=1)


ggsave("output/figures/vis_differential_proteomics_corrplot.pdf", width=3.4 , height=3.4*2, dpi=600)



#### GSEA ----


metadata.proteomics.glass_od |> 
  dplyr::filter(!is.na(DPA__GLASS_OD__CGC__adj.P.Val)) |>
  dplyr::filter(DPA__GLASS_OD__CGC__adj.P.Val < 0.01) |>
  dplyr::filter(abs(DPA__GLASS_OD__CGC__logFC) > 0.5) |>
  dplyr::filter(DPA__GLASS_OD__CGC__logFC > 0) |> 
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::pull(Genes)


metadata.proteomics.glass_od |> 
  dplyr::filter(!is.na(DPA__GLASS_OD__CGC__adj.P.Val)) |>
  dplyr::filter(DPA__GLASS_OD__CGC__adj.P.Val < 0.01) |>
  dplyr::filter(abs(DPA__GLASS_OD__CGC__logFC) > 0.5) |>
  dplyr::filter(DPA__GLASS_OD__CGC__logFC < 0) |> 
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::pull(Genes)



write.table(
  plt |>
    dplyr::filter(!is.na(adj.P.Val_cgc)) |>
    dplyr::pull(protein_id) |> 
    stringr::str_split(";", n = Inf, simplify = FALSE) |> 
    unlist()
  ,
  "output/tables/GSEA_background_proteomics_cgc_cor.txt",  quote = F, row.names = F, col.names=F)



write.table(
  pplt |> 
    dplyr::filter(variable == 'cor CGC' & value != 0) |> 
    dplyr::filter(value > 0) |> 
    dplyr::pull(gid) |> 
    stringr::str_split(";", n = Inf, simplify = FALSE) |> 
    unlist()
  ,
  "output/tables/GSEA_proteins_up_proteomics_cgc_cor.txt",  quote = F, row.names = F, col.names=F)




write.table(
  pplt |> 
    dplyr::filter(variable == 'cor CGC' & value != 0) |> 
    dplyr::filter(value < 0) |> 
    dplyr::pull(gid) |> 
    stringr::str_split(";", n = Inf, simplify = FALSE) |> 
    unlist()
  ,
  "output/tables/GSEA_proteins_down_proteomics_cgc_cor.txt",  quote = F, row.names = F, col.names=F)



## nega: ----


plt |>
  dplyr::filter(significant & logFC_cgc < 0) |> 
  dplyr::arrange(P.Value_cgc) |> 
  dplyr::pull('protein_id')



# GLASS-OD: grade x prim-rec naive ----


plt <- metadata.proteomics.glass_od


ggplot(plt, aes(x=`DPA__GLASS_OD__prim-rec__naive__t`, y=`DPA__GLASS_OD__grade__naive__t`)) +

  geom_hline(yintercept=0, col="red", lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="red", lwd=theme_cellpress_lwd) +
  
  geom_point(size=theme_cellpress_size/3) +
  
  geom_hline(yintercept=0, col="red", alpha=0.5, lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="red", alpha=0.5, lwd=theme_cellpress_lwd) +
  
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)), col="1", cor.coef.name ="R", size=theme_cellpress_size, family=theme_cellpress_font_family) +

  theme_cellpress +
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +
  theme(aspect.ratio=1)



# GLASS-OD: grade x prim-rec patcor ----



plt <- metadata.proteomics.glass_od |> 
  dplyr::filter(!is.na(`DPA__GLASS_OD__prim-rec__pat_corrected__adj.P.Val`)) |> 
  dplyr::filter(!is.na(`DPA__GLASS_OD__grade__pat_corrected__adj.P.Val`)) |> 
  dplyr::mutate(significant = (`DPA__GLASS_OD__prim-rec__pat_corrected__adj.P.Val` < 0.01) | (DPA__GLASS_OD__grade__pat_corrected__adj.P.Val < 0.01))


ggplot(plt, aes(y=`DPA__GLASS_OD__prim-rec__pat_corrected__t`,
                x=`DPA__GLASS_OD__grade__pat_corrected__t`,
                col=significant,
                label=gsub("_HUMAN$", "", Protein.Names))) +
  
  geom_hline(yintercept=0, col="red", lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="red", lwd=theme_cellpress_lwd) +
  
  geom_point(size=theme_cellpress_size/3) +
  
  geom_hline(yintercept=0, col="red", alpha=0.5, lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="red", alpha=0.5, lwd=theme_cellpress_lwd) +
  
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)), col="1", cor.coef.name ="R", size=theme_cellpress_size, family=theme_cellpress_font_family) +
  
  ggrepel::geom_text_repel(data=subset(plt, significant), col="darkblue") +
  
  theme_cellpress +
  
  scale_color_manual(values=c(`TRUE`='red',`FALSE`='black')) +
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +
  theme(aspect.ratio=1)









ggplot(plt, aes(y=`DPA__GLASS_OD__prim-rec__pat_corrected__t`,
                x=`DPA__GLASS_OD__CGC__t`,
                col=significant,
                label=gsub("_HUMAN$", "", Protein.Names))) +
  
  geom_hline(yintercept=0, col="red", lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="red", lwd=theme_cellpress_lwd) +
  
  geom_point(size=theme_cellpress_size/3) +
  
  geom_hline(yintercept=0, col="red", alpha=0.5, lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="red", alpha=0.5, lwd=theme_cellpress_lwd)


## ggcorrplot style ----


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(proteomics_id %in% c("GLODprot_12_1048", "GLODprot_42_1169", "GLODprot_56_1230", "GLODprot_56_1232") == F) |>   # no pass qc
  
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, patient_id, starts_with("proteomics_"), starts_with("resection_"), contains("_excluded")) |> 
  dplyr::mutate(condition = factor(ifelse(resection_number == 1,"primary","recurrence"), levels=c("primary","recurrence"))) |> 
  assertr::verify(!duplicated(resection_id)) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 137) 
    return(.)
  })()


data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id)
data$nas <- rowSums(is.na(data))
data <- data |> 
  dplyr::filter(nas < 60) |> 
  dplyr::mutate(nas=NULL) |> 
  tibble::rownames_to_column('gene_ids')



data2 <- data |> 
  dplyr::left_join(
    plt |>  dplyr::select(Genes, significant, DPA__GLASS_OD__grade__pat_corrected__t),
    by=c('gene_ids'='Genes')) |>
  dplyr::filter(significant) |> 
  dplyr::mutate(significant = NULL) |> 
  dplyr::arrange(DPA__GLASS_OD__grade__pat_corrected__t) |> 
  dplyr::mutate(DPA__GLASS_OD__grade__pat_corrected__t = NULL) |> 
  tidyr::separate_rows(gene_ids, sep = ";") |> 
  tibble::column_to_rownames('gene_ids')


c = cor(na.omit(t(data2)))
g = ggcorrplot(c)
g

ggsave("output/figures/vis_differential_proteomics__ggcorr_grade_pr.pdf", width = 4, height = 2.8)



matching_meth <-  data.mvalues.probes |>
  dplyr::select(probe_id, UCSC_RefGene_Name,
                UCSC_RefGene_Group,
                DMP__g2_g3__pp_nc_PC1__t) |> 
  tidyr::separate_rows(
    UCSC_RefGene_Name,
    UCSC_RefGene_Group,
    sep = ";"
  ) |> 
  dplyr::filter(UCSC_RefGene_Name != "") |> 
  dplyr::distinct(probe_id, UCSC_RefGene_Name, UCSC_RefGene_Group, DMP__g2_g3__pp_nc_PC1__t, .keep_all = TRUE)    #UCSC_RefGene_Group %in% c("TSS200","TSS1500","1stExon","5'UTR"))



plt.panel = g$data |> dplyr::select(x, x.order) |>
  unique()  |> 
  tibble::as_tibble() |> 
  dplyr::rename(gene_id = x)



matching_meth2 <- matching_meth |> 
  dplyr::filter(UCSC_RefGene_Name %in% plt.panel$gene_id) |> 
  dplyr::mutate(probe_id = NULL) |> 
  dplyr::mutate(UCSC_RefGene_Group = ifelse(UCSC_RefGene_Group %in% c("TSS1500", "TSS200"),"TSS1500/200", UCSC_RefGene_Group)) |> 
  dplyr::group_by(UCSC_RefGene_Name, UCSC_RefGene_Group) |> 
  dplyr::summarise(
    across(where(is.numeric), mean, .names = "{.col}"),
    .groups = "drop"
  )  |> 
  dplyr::filter(!is.na(DMP__g2_g3__pp_nc_PC1__t))





matching_meth3 <-   plt.panel |>
  dplyr::left_join(
    rbind(matching_meth2, matching_meth2 |> dplyr::mutate(DMP__g2_g3__pp_nc_PC1__t = 0)),
    by=c('gene_id' = 'UCSC_RefGene_Name')
  )


complete_grid <- expand.grid(gene_id = unique(matching_meth3$gene_id),
                             UCSC_RefGene_Group = unique(matching_meth3$UCSC_RefGene_Group),
                             stringsAsFactors = FALSE) |> 
  dplyr::left_join(matching_meth3, by = c("gene_id", "UCSC_RefGene_Group")) |> 
  dplyr::select(-x.order) |> 
  dplyr::left_join(matching_meth3 |> dplyr::select(gene_id, x.order) |> dplyr::distinct(), by=c('gene_id'='gene_id')) |> 
  dplyr::filter(!is.na( UCSC_RefGene_Group))  |> 
  dplyr::filter(UCSC_RefGene_Group %in% c("3'UTR", "ExonBnd") == F) # not really relevant, only few CpGs and lot of space




ggplot(complete_grid,
       aes(x = DMP__g2_g3__pp_nc_PC1__t,
           y = reorder(gene_id, -x.order))
       ) + 
  facet_wrap(~UCSC_RefGene_Group, ncol=8) +
  geom_vline(aes(xintercept=0), lty=1, lwd=theme_cellpress_lwd, col="gray80") +
  geom_line(lwd=theme_cellpress_lwd) +
  geom_point(data=subset(complete_grid, 
                         is.na(DMP__g2_g3__pp_nc_PC1__t)),
             aes(x=0),
             col=mixcol("red","white",0.5),
             cex=theme_cellpress_size/10,
             pch=4) + 
  theme_cellpress +
  labs(y=NULL, x="Mean t-score methylation Grade 2 and Grade 3 oligodendroglioma")


ggsave("output/figures/vis_differential_proteomics__grade_pr_methylation.pdf", width=2.4, height=3.1)





plt.panel2 <-   plt.panel |> 
  dplyr::left_join(
    plt |> 
      #dplyr::mutate(significant_grade = DPA__GLASS_OD__grade__pat_corrected__adj.P.Val < 0.01) |> 
      #dplyr::mutate(significant_pr = `DPA__GLASS_OD__prim-rec__pat_corrected__adj.P.Val` < 0.01) |> 
      dplyr::select(Genes,
                    `DPA__GLASS_OD__grade__pat_corrected__adj.P.Val`,
                    `DPA__GLASS_OD__prim-rec__pat_corrected__adj.P.Val`,
                         #DPA__GLASS_OD__CGC__t,
                         `DPA__GLASS_OD__grade__pat_corrected__t`,
                         `DPA__GLASS_OD__prim-rec__pat_corrected__t`
    ),
    by=c('gene_id'='Genes')
  ) |> 
  tidyr::pivot_longer(
    cols = dplyr::starts_with("DPA__GLASS_OD"),
    names_to = c("variable", ".value"),
    names_pattern = "DPA__GLASS_OD__(.*)__pat_corrected__(.*)"
  )


#  tidyr::pivot_longer(cols=dplyr::starts_with("DPA__"))




plt.panel2 <- rbind(plt.panel2,
                    plt.panel2 |> dplyr::mutate(t = 0))
  #dplyr::mutate(name = gsub("DPA__GLASS_OD__", "", name)) |> 
  #dplyr::mutate(name = gsub("__pat_corrected", "", name)) |> 
  #dplyr::mutate(name = gsub("__t$", "", name))



ggplot(plt.panel2,
       aes(x = t,
           y = reorder(gene_id, -x.order),
           col = adj.P.Val < 0.01)
) + 
  facet_wrap(~variable, ncol=2) +
  geom_vline(aes(xintercept=0), lty=1, lwd=theme_cellpress_lwd, col="gray80") +
  geom_line(lwd=theme_cellpress_lwd) +
  geom_point(data=subset(plt.panel2, 
                         is.na(t)),
             aes(x=0),
             col=mixcol("red","white",0.5),
             cex=theme_cellpress_size/10,
             pch=4) + 
  theme_cellpress +
  scale_color_manual(values = c('TRUE'='darkgreen', 'FALSE'=mixcol('red','black',0.25))) +
  labs(y=NULL, x="Difference protein expression (t-score)", col="significant")



ggsave("output/figures/vis_differential_proteomics__grade_pr_cgc__proteomics.pdf", width=1.8, height=3.35)


## concept 2 ----



metadata.proteomics.glass_od |>
  dplyr::filter(DPA__GLASS_OD__CGC__P.Value < 0.01 ) |>
  dplyr::filter(abs(DPA__GLASS_OD__CGC__logFC) > 0.5) |>
  tidyr::separate_rows(Genes, sep = ";") |> 
  dplyr::pull(Genes)


stats.a <- data.frame(
  Genes  = c(
    metadata.proteomics.glass_od |> 
      dplyr::filter(!is.na(`DPA__GLASS_OD__prim-rec__pat_corrected__adj.P.Val`)) |> 
      dplyr::filter(`DPA__GLASS_OD__prim-rec__pat_corrected__adj.P.Val` < 0.01) |> 
      dplyr::pull(Genes),
    metadata.proteomics.glass_od |> 
      dplyr::filter(!is.na(`DPA__GLASS_OD__grade__pat_corrected__adj.P.Val`)) |> 
      dplyr::filter(DPA__GLASS_OD__grade__pat_corrected__adj.P.Val < 0.01) |> 
      dplyr::pull(Genes)
  )) |> 
  dplyr::distinct() |> 
  dplyr::pull(Genes)

stats.b = metadata.proteomics.glass_od |>
                             dplyr::filter(DPA__GLASS_OD__CGC__P.Value < 0.01 ) |>
                             dplyr::filter(abs(DPA__GLASS_OD__CGC__logFC) > 0.5) |>
                             tidyr::separate_rows(Genes, sep = ";") |> 
                             dplyr::pull(Genes)






plt <- data.frame(
  Genes  = c(
    metadata.proteomics.glass_od |> 
      dplyr::filter(!is.na(`DPA__GLASS_OD__prim-rec__pat_corrected__adj.P.Val`)) |> 
      dplyr::filter(`DPA__GLASS_OD__prim-rec__pat_corrected__adj.P.Val` < 0.01) |> 
      dplyr::pull(Genes),
    metadata.proteomics.glass_od |> 
      dplyr::filter(!is.na(`DPA__GLASS_OD__grade__pat_corrected__adj.P.Val`)) |> 
      dplyr::filter(DPA__GLASS_OD__grade__pat_corrected__adj.P.Val < 0.01) |> 
      dplyr::pull(Genes)
  )) |> 
  dplyr::distinct() |> 
  dplyr::left_join(
    metadata.proteomics.glass_od |>
      dplyr::select(Genes, DPA__GLASS_OD__CGC__t) |>
      dplyr::mutate(panel = ifelse(DPA__GLASS_OD__CGC__t > 0 , "Proteomics up", "Proteomics down")) |> 
      dplyr::mutate(panel = factor(panel, levels=c("Proteomics up", "Proteomics down"))) |> 
      dplyr::rename(rank = DPA__GLASS_OD__CGC__t),
    by=c('Genes'='Genes')
  ) |> 
  dplyr::left_join(
    metadata.proteomics.glass_od |> 
      dplyr::filter(           !is.na(`DPA__GLASS_OD__prim-rec__pat_corrected__adj.P.Val`)) |> 
      dplyr::mutate(pr_significant = (`DPA__GLASS_OD__prim-rec__pat_corrected__adj.P.Val` < 0.01)) |> 
      dplyr::mutate(pr_t =            `DPA__GLASS_OD__prim-rec__pat_corrected__t`) |> 
      dplyr::select(Genes, pr_significant, pr_t),
    by=c('Genes'='Genes')
  ) |> 
  dplyr::left_join(
    metadata.proteomics.glass_od |> 
      dplyr::filter(              !is.na(`DPA__GLASS_OD__grade__pat_corrected__adj.P.Val`)) |> 
      dplyr::mutate(grade_significant = (`DPA__GLASS_OD__grade__pat_corrected__adj.P.Val` < 0.01)) |> 
      dplyr::mutate(grade_t =            `DPA__GLASS_OD__grade__pat_corrected__t`) |> 
      dplyr::select(Genes, grade_significant, grade_t),
    by=c('Genes'='Genes')
  ) |> 
  tidyr::pivot_longer(
    cols = dplyr::ends_with("_significant") | dplyr::ends_with("_t"),
    names_to = c("comparison", ".value"),
    names_pattern = "(.*)_(significant|t)"
  ) |> 
  (function(.) {
    return(dplyr::bind_rows(. |> dplyr::mutate(type = "data"),
                            . |> dplyr::mutate(type = "reference point") 
                              |> dplyr::mutate(`t` = 0)))
  })()



p1 = ggplot(plt, aes(x=t, y=reorder(Genes, rank), col=significant)) +
  facet_grid(cols = vars(comparison), rows=vars(panel), scales = "free_y", space = "free_y") +
  geom_vline(aes(xintercept=0), lty=1, lwd=theme_cellpress_lwd, col="gray80") +
  geom_line(lwd=theme_cellpress_lwd) +
  scale_color_manual(values=c('TRUE'='black','FALSE'='gray70')) +
  theme_cellpress +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(y=NULL) +
  theme(panel.spacing = unit(0.125, "cm"))
p1






plt.meth_grade <- data.mvalues.probes |>
  dplyr::select(probe_id, UCSC_RefGene_Name,
                UCSC_RefGene_Group,
                DMP__g2_g3__pp_nc_PC1__t) |> 
  tidyr::separate_rows(
    UCSC_RefGene_Name,
    UCSC_RefGene_Group,
    sep = ";"
  ) |> 
  dplyr::filter(UCSC_RefGene_Name != "") |> 
  dplyr::distinct(probe_id, UCSC_RefGene_Name, UCSC_RefGene_Group, DMP__g2_g3__pp_nc_PC1__t, .keep_all = TRUE) |> #UCSC_RefGene_Group %in% c("TSS200","TSS1500","1stExon","5'UTR"))
  dplyr::mutate(UCSC_RefGene_Group = ifelse(UCSC_RefGene_Group %in% c("TSS1500", "TSS200"),"TSS1500/200", UCSC_RefGene_Group)) |> 
  dplyr::filter(!is.na(DMP__g2_g3__pp_nc_PC1__t)) |> 
  dplyr::mutate(probe_id = NULL) |> 
  dplyr::group_by(UCSC_RefGene_Name, UCSC_RefGene_Group) |> 
  dplyr::summarise(
    across(where(is.numeric), mean, .names = "{.col}"),
    .groups = "drop"
  )






complete_grid <- expand.grid(UCSC_RefGene_Name = plt |> dplyr::arrange(rank) |> dplyr::select(Genes) |> dplyr::distinct() |> dplyr::pull(Genes),
                             UCSC_RefGene_Group = unique(plt.meth_grade$UCSC_RefGene_Group),
                             stringsAsFactors = FALSE) |> 
  dplyr::left_join(plt.meth_grade, by = c("UCSC_RefGene_Name", "UCSC_RefGene_Group")) |> 
  dplyr::rename(Genes = UCSC_RefGene_Name) |> 
  dplyr::left_join(plt |> dplyr::select(Genes, rank, panel) |> dplyr::distinct(), by=c('Genes'='Genes')) |> 
  dplyr::filter(!is.na( UCSC_RefGene_Group))  |> 
  dplyr::filter(UCSC_RefGene_Group %in% c("3'UTR", "ExonBnd") == F) |>  # not really relevant, only few CpGs and lot of space
  (function(.) {
    return(dplyr::bind_rows(. |> dplyr::mutate(type = "data"),
                            . |> dplyr::mutate(type = "reference point") |> 
                              dplyr::mutate(`DMP__g2_g3__pp_nc_PC1__t` = 0)))
  })() |> 
  dplyr::mutate( UCSC_RefGene_Group = factor( UCSC_RefGene_Group, levels=c( "TSS1500/200", "1stExon", "5'UTR", "Body")))



p2 = ggplot(complete_grid, aes(x=DMP__g2_g3__pp_nc_PC1__t, y = reorder(Genes, rank))) +
  facet_grid(cols = vars(UCSC_RefGene_Group), rows=vars(panel), scales = "free_y", space = "free_y") +
  geom_vline(aes(xintercept=0), lty=1, lwd=theme_cellpress_lwd, col="gray80") +
  geom_line(lwd=theme_cellpress_lwd) +
  geom_point(data=subset(complete_grid, is.na(DMP__g2_g3__pp_nc_PC1__t)),
             aes(x=0),
             col=mixcol("red","white",0.5),
             cex=theme_cellpress_size/10,
             pch=4) + 
  theme_cellpress + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(x=NULL,y="Mean t-score between Grade 2 - Grade 3") +
  theme(panel.spacing = unit(0.125, "cm"))
p2



library(patchwork)
(p1 + p2) + plot_layout(
  widths = c(1.3, 2.5)
) #/ (plot_spacer() + p3)
ggsave("output/figures/vis_differential_proteomics__master1.pdf", width=3.45, height=1.8)



p3 = ggplot(complete_grid |> dplyr::filter(!is.na(DMP__g2_g3__pp_nc_PC1__t) & type == "data"), aes(x=panel, y=DMP__g2_g3__pp_nc_PC1__t)) +
  facet_grid(cols = vars(UCSC_RefGene_Group)) +
  ggbeeswarm::geom_quasirandom(size=theme_cellpress_size/8, width = 0.25) +
  ggpubr::stat_compare_means(aes(label = paste0("p = ", ..p.format..)),
                             method = "t.test",
                             label.x.npc=0.1,
                             label.y.npc=0.035,
                             show_guide  = FALSE,
                             size=theme_cellpress_size,
                             family=theme_cellpress_font_family
  ) +
  labs(x=NULL) +
  theme_cellpress + theme(
    axis.text.x = element_text(
      angle = 90, 
      hjust = 1.0,
      vjust = 0.5
    )
  ) +
  theme(panel.spacing = unit(0.25, "cm"))
p3

ggsave("output/figures/vis_differential_proteomics__master2.pdf", width=2, height=1.96)



plt.0 <- complete_grid |> 
  dplyr::filter(UCSC_RefGene_Group == "TSS1500/200") |>
  dplyr::filter(type == "reference point") |>
  dplyr::select(Genes, rank, `panel`) |> 
  dplyr::mutate(GO_0062023_collagen_containing_extracellular_matrix = Genes %in% c(
    "DCN","SEMA3B","MARCO","SERPINB1","MYOC","TIMP2","VCAN","TNC","USH2A","COL9A2","LTBP1","ELN","LAMC3","COL23A1","LAMA3","ENTPD2","CBLN4","ITIH4","ITIH1","F7","LAMC2","COL11A1","WNT8A","GPC1","CDON","NTN1","COL17A1","FGFR2","PKM","FGF10","SPP2","ADAM11","NTN4","LMAN1","DLG1","GPC4","FBLN1","COL5A3","HSP90AA1","COL4A4","IMPG2","COL19A1","COL16A1","FCN1","ACHE","ADAMTS2","MMP2","NID2","LTBP4","ICAM1","P3H2","LAMB4","LAMB1","APOH","ANGPT2","CCDC80","CMA1","COL9A3","TGFB2","AMELY","HNRNPM","LGALS1","TIMP3","PDGFB","CHADL","CTSG","COCH","MMP9","BMP7","CTSZ","COL20A1","ANGPT4","LAMA1","MXRA5","SRPX","F9","TIMP1","SRPX2","CBLN1","ZP2","CTSH","SFRP1","IL7","FGL1","CLC","APLP1","TGFB1","COMP","WNT2","PTPRZ1","PCOLCE","SERPINE1","PLOD3","RARRES2","AEBP1","MEGF9","OGN",
    "ASPN","ECM2","AMBP","CXCL12","ACTA2","KAZALD1","LGALS3BP","COL1A1","VTN","SOD3","CTSC","HPX","APOA4","APOC3","MDK","VWF","WNT5B","MGP","COL12A1","COL9A1","SMOC2","IMPG1","LAMA4","ERBIN","LOX","SPARC","THBS4","KNG1","HRG","WNT5A","COL7A1","LOXL3","EFEMP1","FN1","TNR","ANGPTL1","PRG4","NID1","ZP4","MFAP2","P3H1","F3","SERPINC1","CTSD","MMP8","APOA1","CCN2","LTBP2","TGFB3","TNN","TGFBI",
    "CLU","A1BG","FMOD","PLG","ANXA11","COL10A1","SLPI","MATN4","NDP","F13A1","SERPINB6","COL21A1","AMELX","CFP","PZP","FGL2","LOXL1","APOE","NCAN","PXDN","GDF15","COL5A1","LAMA5","FIBCD1","ANGPTL6","F12","LGALS3","PODNL1","RTBDN","MATN3","EMILIN2","SERPINF1","ITGB4","MATN2","BCAN","HAPLN2","APCS","ANGPTL3","POSTN","LOXL2","ADAMDEC1","WNT2B","COL4A2","ADAMTS8","ANXA1","CTSL","ADAM19","AGT","LAMC1","SERPINE2","ANGPTL2","CCN3","IGFBPL1","TINAG","SULF1","THBS1",
    "EMILIN1","LOXL4","ANXA7","CILP","SEMA7A","MMRN1","FRAS1","FBN2","COL2A1","INHBE","LUM","FBLN5","PCSK6","LMAN1L","HAPLN3","MFGE8","TGFB1I1","CDH13","COL6A1","COL6A2","ADAMTS10","FCN3","HSPG2","CCN1","TINAGL1","DPT","HDGF","HMCN1","ECM1","ADAMTSL4","ANXA9","S100A8","S100A7","FLG","LEFTY2","COL8A1","AHSG","SFRP2","HAPLN1","CASK","GPC3","HMCN2","SERPING1","SERPINH1","NCAM1","ZP1","FREM2","ITIH2","DST","SPARCL1","ABI3BP","ANGPT1","ADAMTS1","ADAMTS5","ADAMTS3","PRG3","ACAN","COLEC12","ADAMTS4","LAD1","C1QC","CSTB","FCN2","AZGP1","COL26A1","MATN1","SDC3","CTSS","S100A9","COL6A3","IGFBP7","FBLN2","ADAMTS9","CRELD1","PF4","CPA3","TGM4","CLEC3B","ANXA5","EDIL3","EGFLAM","RELL2","SHH","COL1A2","CTSB","SBSPON","CTHRC1","FREM1","VWA2","HTRA1","ADAMTS15","FBN1","SERPINB8","MFAP4","HSP90B1","SERPINB12","NAV2","GREM1","SERPINF2","KRT1","ANGPTL4",
    "SOST","LTBP3","TNXB","COL3A1","NPNT","LGALS9","COL4A3","THBS3","COL22A1","SDC2","MUC17","SERPINB9","CDH2","COL24A1","FGG","FGA","FGB","LGALS4","COL8A2","ANGPTL7","LAMB2","LRRC15","TPSAB1","EFEMP2","COL6A5","EGFL7","ADAMTS20","MMRN2","C1QA","DAG1","CSPG4","CTSF","PODN","FGFBP3","ZG16","NPPA","A2M","CLEC14A","CD151","CALR","GPC5","VWA1","SSC5D","F2","ADIPOQ","OTOL1","BGN","MXRA7","SNORC","ANXA2","COL18A1","FREM3","GPC6","EMILIN3","EFNA5","THBS2","PRG2","C17orf58","ANGPTL5","BCAM","COL4A1","HAPLN4","AMTN","THSD4","COL14A1","EYS","COL4A5","AGRN","PLSCR1","ZP3","IFNA2","SERPINA5","COL25A1","VWC2","PRELP","MMP23B","SERPINA3","S100A4","PRTN3","LAMA2","COL27A1","LAMB3","ANXA4","ANXA6","SERPINA1","TPSB2","COL13A1","SPN","ELANE","COL4A6","MFAP5","PSAP","S100A10","HRNR","S100A6","SMOC1","EGFL6","L1CAM","TGM2","COL11A2","COL5A2","COL15A1","PRSS1","VIT","DEFA1","COL6A6","COLQ","GPC2","ANG","COL28A1","ORM2","ORM1","DEFA1B","MARCOL","GH1","SPON1","ANXA8","RBP3","GDF10","MMP28"
  )) |> 
  dplyr::mutate(GO_0002250_adaptive_immune_response = c(
    Genes %in% c(
      "TMEM98", "CD79B", "MASP2", "CD4", "BTK", "HFE", "FYN", "ERCC1", "KDM5D", "CD6", "UFL1", "WAS", "SLC11A1", "CD74", "C8B", "SLAMF7", "BTN3A1", "TNFRSF1B", "C6", "PARP3", "TNFRSF17", "FOXP3", "RC3H2", "TRAF3IP2", "PRDM1", "PRKCQ", "CD84", "TP53BP1", "PRKCZ", "TFE3", "BCL3", "RORA", "RAB27A", "TFRC", "FCGR2B", "PVR", "TBX21", "TXK", "MLH1", "UNG", "PAG1", "RAP1GAP", "NFKB2", "IL4R", "LAMP3", "RAPGEF3", "CEACAM1", "STX7", "RIF1", "ARG2", "MEF2C", "PTPRC", "LAT2", "P2RX7", "LAG3", "ICAM1", "CD209", "MCOLN1", "RAPGEF4", "UNC13D", "MSH2", "TREM2", "SIRT1", "JAK2", "IL12RB1", "KLHL22", "LGALS1", "CSF2RB", "IL2RB", "EP300", "CD40", "SLA2", "JAG1", "SMAD7", "RNF125", "TLR8", "BMX", "CD40LG", "TNFSF13B", "PYCARD", "IL21R", "CSK", "CTSH", "RIPK2", "NBN", "FGL1", "RELB", "FCER2", "CLEC4M", "IL4I1", "LILRA1", "IL27RA", "CLC", "EBI3", "PRKD2", "TGFB1", "CD79A", "TYK2", "LILRB5", "JAK3", "PIK3CG", "EPHB6", "AHR", "C5", "EXOSC3", "GATA3", "PPP3CB", "C1QBP", "SUPT6H", "IL2", "NSD2", "CTSC", "CRTAM", "HSPA8", "UNC93B1", "KMT5B", "HPX", "CD81", "TCIRG1", "CD69",
      "IL23A", "IFNG", "PTPN6", "CLEC4A", "AICDA", "BTN3A3", "ULBP1", "IL17A", "IL17F", "RNF8", "BACH2", "CCR6", "TFEB", "DUSP22", "C7", "ITK", "IL12B", "BTNL8", "IL4", "C9", "BCL6", "CD86", "ZAP70", "LOXL3", "IL1R1", "IL1RL1", "IL18R1", "MSH6", "RNF19B", "MAD2L2", "PLA2G4A", "SLAMF1", "CD160", "CR2", "CD46", "MPL", "TNFSF4", "TNFAIP3", "ARG1", "CD274", "IFNA6", "IFNA8", "TNFSF18", "CLU", "PTK2B", "ADCY7", "KCNJ8", "CD80", "CCR2", "SASH3", "LAX1", "CD244", "LY9", "SHLD2", "PMS2", "PKN1", "NCKAP1L", "IL13RA2", "C4BPA", "C4BPB", "IL9R", "TREM1", "IRF1", "IL1B", "CD70", "C3", "STAT5A", "ZBTB1", "TRAF2", "IGLL1", "EIF2AK4", "RIPK3", "FOXJ1", "CRACR2A", "NECTIN2", "TRPM4", "ULBP2", "ULBP3", "LILRB2", "PRR7", "RFTN1", "NDFIP1", "PRKAA1", "RAP1GAP2", "JCHAIN", "CLEC10A", "CTNNBL1", "EPHB2", "KMT5C", "SWAP70", "VTCN1", "RSAD2", "IL6ST", "KLRD1", "KLRC1", "ANXA1", "CTSL", "HAVCR2", "AKIRIN2", "MAP3K7", "TEC", "RC3H1", "IL6", "DBNL", "MYO1G", "HLX", "IL10", "TLR4", "IL33", "SIT1", "IFNA21", "IRF4", "IL18BP", "STAT4", "LEF1", "C1RL", "CD27", "SLC15A4", "SKAP1", "SECTM1", "TNFRSF11A", "BRD4", "SIGLEC10", "BCL10", "XCL1", "FCGR2A", "RORC", "SUSD4", "DUSP10", "HSPD1", "NFKBIZ", "IL9", "FBXO38", "TNFRSF21", "EBAG9", "IFNA5", "IFNA16", "IFNK", "NOTCH1", "SERPING1", "FCGR1A", "CD226", "IL18", "ADAM17", "DCLRE1C", "CAMK4", "MR1", "CYRIB",
      "CD8A", "MCOLN2", "SAMSN1", "RAET1L", "BATF", "CXCL13", "C8A", "PAXIP1", "TNFRSF14", "CD1D", "CD1A", "CD1C", "CD1B", "CD1E", "HMHB1", "FCER1G", "C1QC", "C1R", "TNFRSF13C", "ICOSLG", "AIRE", "CD3G", "ZBTB7B", "IL6R", "PCYT1A", "ALOX15", "TNFSF13", "JAK1", "IL23R", "NLRP3", "SLAMF6", "FCGR3B", "FCMR", "FCAMR", "SANBR", "CTSS", "FZD5", "EOMES", "FCRL4", "TRAT1", "CTLA4", "ZC3H12A", "PRKCD", "RNF168", "ERAP1", "ERAP2", "RAET1E", "SYK", "SVEP1", "MARCHF8", "MBL2", "TSC1", "HPRT1", "JAM3", "C2", "RAG1", "PRKCB", "CLEC4D", "B2M", "STAT6", "PDIA3", "PHB1", "NOD2", "CD3D", "NFKBID", "LAIR1", "OTUB1", "FADD", "CX3CR1", "TAP1", "STAT3", "IL7R", "IL12A", "INPP5D", "GPR183", "APLF", "ALCAM", "YWHAG", "SOCS5", "JUNB", "FGA", "FGB", "PIK3CD", "IFNB1", "SHLD1", "CD8B", "MALT1", "CLEC7A", "KLHL6", "THEMIS", "CCL19", "MYD88", "SLC22A13", "C1QA", "CD7", "EXO1", "LIG4", "IL20RB", "HRAS", "ADGRE1", "TRAF6", "CLCF1", "ZNF683", "ATAD5", "C8G", "IFNW1", "CD19", "IL17RA", "CD28", "HLA-DQB1", "FCER1A", "FUT7", "PRF1", "NLRP10", "C1S", "CLEC4G", "HMCES", "PSG9", "ASCL2", "SH2D1A", "SOCS3", "IFNE", "IRF7", "BTLA", "BTN3A2", "IFNA10", "LILRB4", "CARD9", "C17orf99", "NCR3LG1", "ZP3", "IFNA2", "PDCD1", "IGHV1OR15-9", "HMGB1", "HLA-DRB1", "SEMA4A", "CD55", "HLA-DQA1", "ADA", "ARID5A", "IL27", "SPN", "GZMM", "MPEG1", "PDCD1LG2", "CR1L", "IFNA1", "ENTPD7", "CLEC4C", "HLA-DRB5", "SH2D1B", "MTOR", "CD247", "OPA1", "CD3E", "CR1", "RAET1G", "FCGR3A", "LIME1", "HLA-DOA", "BRD2", "HLA-DMA", "TAP2", "HLA-DRA", "AGER", "MICB", "HLA-C", "LILRB3", "GNL1", "HLA-E", "HLA-G", "HLA-F", "TRIM27", "CFI", "KLRC2", "CLEC6A", "HLA-A", "IGKJ1", "IGKV4-1", "IGKV5-2", "IGKV6-21", "IGKV2D-26", "IGKV3D-20", "IGKV6D-41", "IGKV3D-11", "IGKV1D-42", "IGLV4-69", "IGLV8-61", "IGLV4-60", "IGLV6-57", "IGLV11-55", "IGLV10-54", "IGLV5-52", "IGLV1-51", "IGLV1-50", "IGLV5-48", "IGLV1-47", "IGLV7-46", "IGLV5-45", "IGLV1-44", "IGLV7-43", "IGLV1-40", "IGLV5-37", "IGLV1-36", "IGLV2-33", "IGLV3-32", "IGLV3-27", "IGLV3-25", "IGLV2-23", "IGLV3-22", "IGLV3-21", "IGLV3-19", "IGLV2-18", "IGLV3-16", "IGLV2-14", "IGLV3-12", "IGLV2-11", "IGLV3-10", "IGLV3-9", "IGLV4-3", "IGLV3-1", "IGLJ1", "IGLC2", "IGLC3", "TRGC1", "TRGJ1", "TRGV11",
      "TRGV10", "TRGV9", "TRGV8", "TRGV5", "TRGV4", "TRGV3", "TRGV1", "TRBV6-1", "TRBV7-1", "TRBV4-1", "TRBV6-4", "TRBV7-3", "TRBV5-3", "TRBV9", "TRBV10-1", "TRBV11-1", "TRBV6-5", "TRBV6-6", "TRBV5-5", "TRBV7-6", "TRBV5-6", "TRBV5-7", "TRBV5-1", "TRBV4-2", "TRBV19", "TRBV20-1", "TRBV23-1", "TRBV24-1", "TRBV27", "TRBV28", "TRBJ2-1", "TRBJ2-3", "TRBJ2-7", "TRAV2", "TRAV3", "TRAV4", "TRAV5", "TRAV6", "TRAV7", "TRAV8-1", "TRAV9-1", "TRAV10", "TRAV12-1", "TRAV8-2", "TRAV8-3", "TRAV13-1", "TRAV12-2", "TRAV8-4", "TRAV13-2", "TRAV14DV4", "TRAV9-2", "TRAV12-3", "TRAV8-6", "TRAV16", "TRAV17", "TRAV18", "TRAV19", "TRAV20", "TRAV21", "TRAV22", "TRAV23DV6", "TRDV1", "TRAV24", "TRAV25", "TRAV26-1", "TRAV27", "TRAV29DV5", "TRAV26-2", "TRAV34", "TRAV36DV7", "TRAV38-1", "TRAV38-2DV8", "TRAV39", "TRAV40", "TRAV41", "TRDV2", "TRDJ1", "TRAJ42", "TRAJ31", "TRAJ3", "IGHA2", "IGHE", "IGHG4", "IGHG2", "IGHA1", "IGHG1", "IGHG3", "IGHM", "IGHJ1", "IGHV6-1", "IGHV1-2", "IGHV1-3", "IGHV2-5", "IGHV3-7", "IGHV3-11", "IGHV3-13", "IGHV3-15", "IGHV3-16", "IGHV1-18", "IGHV3-20", "IGHV3-21", "IGHV3-23", "IGHV1-24", "IGHV2-26", "IGHV4-28", "IGHV3-33", "IGHV4-34", "IGHV3-35", "IGHV3-38", "IGHV4-39", "IGHV1-45", "IGHV1-46", "IGHV3-48", "IGHV3-49", "IGHV5-51", "IGHV3-53", "IGHV1-58", "IGHV4-61", "IGHV3-66", "IGHV1-69", "IGHV2-70D", "IGHV3-73", "IGHV7-81", "DENND1B", "SIPA1", "LAT", "TREX1", "KLRK1", "EMP2", "IFNA7", "SCART1", "IGLV9-49", "EXOSC6", "IGHV3-64", "HLA-DPB1", "TRDD1", "IGKV3D-15", "IGHV4-59", "C4B", "IGHV3-74", "IGKV6D-21", "IGHV3-72", "TRBV2", "LTA", "TRGC2", "IFNA14", "IGKV3D-7", "TRBV10-2", "TRBV5-4", "HLA-DPA1", "IGHV4-31",
      "IGHV3-43", "HLA-DQB2", "TNF", "TRBV29-1", "TRGV2", "IGHV3OR16-10", "IFNA13", "HLA-B", "IFNA17", "IGHD1-1", "IFNA4", "TRBV30", "HLA-DQA2", "TRBV3-1", "IGKV2D-30", "TNFSF12", "IGKV1D-8", "IGKV1-6", "IGKV1-37", "IGKV3-20", "IGKV1D-33", "IGKV1-17", "TNFRSF13B", "IGKV1-8", "IGKV1-16", "HLA-DOB", "IGKV1D-16", "IGKV2-24", "IGKV3-11", "IGKV2D-24", "IGKV1-9", "IGKV1-33", "IGKV1-39", "IGKV2D-28", "HLA-DMB", "IGKV1D-43", "IGKV1D-17", "IGKV3-7", "IGKV2-30", "IGKV2D-29", "IGKV1-12", "IGKV1-5", "CFB", "IGKV2-28", "IGKV3-15", "LILRA6", "IGKV1-27", "C4A", "TARM1", "IGKV1D-37", "IGKV2D-40", "IGKV1D-39", "TRBV6-7", "SHLD3", "TRBV7-7", "TRBV7-4", "TRBV6-8", "LYN", "CD8B2", "IGHV8-51-1", "IGLL5", "TRAV1-1", "TRAV1-2", "TRDV3", "TRAV30", "IGHV4OR15-8", "IGHV2OR16-5", "IGHV3OR15-7", "IGHV3OR16-17", "OTUD7B", "IGHV3OR16-12", "IGHV3OR16-9", "IGHV1OR15-1", "IGHV3-30", "IGHV3OR16-8", "IGHV3OR16-13", "IGKV2-40", "IGHV2-70", "TRBV12-3", "TRBV12-5", "TRBV16", "TRBV14", "TRBV10-3", "ORAI1", "TRBV13", "TRBV18", "IGKV1D-13", "TRBV11-3", "IGHV4-4", "TRBV12-4", "IGHV1OR21-1", "TRBV17", "TRBV7-9", "IGLV2-8", "IGKV1D-12", "IGHV1-69D", "TRBJ1-4", "IGHV1-69-2", "IGHV7-4-1", "TRBJ1-3", "TRBJ1-5", "TRBJ1-1", "TRBJ1-2", "TRBD1", "TRBV25-1", "IGHV3-64D", "IGHV5-10-1", "TRBJ1-6", "TRBV7-2", "TRBV6-2", "IL9R" 
    )
  )) |> 
  tidyr::pivot_longer(cols = c(GO_0062023_collagen_containing_extracellular_matrix, GO_0002250_adaptive_immune_response))



ggplot(plt.0, aes(x=name, y=reorder(Genes, rank), col=value)) +
  facet_grid(rows=vars(panel), scales = "free_y", space = "free_y") +
  geom_point(data=subset(plt.0, value==F), pch='-', size=0.0) +
  geom_point(data=subset(plt.0, value==T), pch='-', size=3) +
  scale_color_manual(values=c(`TRUE`='red', `FALSE`="white")) +
  theme_cellpress +
  labs(y=NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


ggsave("output/figures/vis_differential_proteomics__master3.pdf", width=1, height=3.86)




# GLASS-OD: grade x prim-rec patcor col by genes from fig 2k ----


ids <- metadata.proteomics.glass_od |> 
  dplyr::select(Protein.Names, Genes) |> 
  dplyr::mutate(protein_id = paste0(Protein.Names, ";",Genes)) |> 
  dplyr::mutate(protein_id = sapply(strsplit(protein_id, ";"), function(x) {
    cleaned <- unique(stringr::str_replace(x, "_HUMAN$", ""))
    paste(cleaned, collapse = ";")
  })) |> 
  tidyr::separate_rows(protein_id, sep = ";") |>
  
  #dplyr::filter(protein_id %in% c("FN1", "FINC"))  |> # 1
  #dplyr::filter(grepl("^OR[0-9]", protein_id)) |>   # none
  #dplyr::filter(grepl("^(DAXX|ATRX|CIC|FUBP1|IDH1|IDH2|H3F3A|TCF12|TERT)", protein_id)) |>   # none
  #dplyr::filter(grepl("^(DAXX|ATRX|CIC|FUBP1|IDH1|IDH2|H3F3A|TCF12|TERT)", protein_id)) |>   # none
  dplyr::filter(protein_id %in% c(
    #genes_polycomb_suz12
    #genes_polycomb_eed,
    #genes_polycomb_h3k27
    genes_polycomb_prc2
  )) |> 

  dplyr::pull(Genes) |> 
  unique()

ids


plt <- metadata.proteomics.glass_od |> 
  dplyr::filter(!is.na(`DPA__GLASS_OD__prim-rec__pat_corrected__t`)) |> # excldued due to too many N/A's
  dplyr::mutate(col = ifelse(Genes %in% ids,"yes","no"))
  #dplyr::mutate(col = ifelse(grepl("mitochondrial", First.Protein.Description), "yes","no"))
plt |> 
  dplyr::filter(col == "yes") |> 
  dplyr::select(Genes, Protein.Names)
#table(plt$col)


ggplot(plt, aes(x=`DPA__GLASS_OD__prim-rec__pat_corrected__t`, y=`DPA__GLASS_OD__grade__pat_corrected__t`, label=`Genes`)) +
  
  geom_hline(yintercept=0, col="red", lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="red", lwd=theme_cellpress_lwd) +
  
  geom_point(data=subset(plt, col != "yes"), size=theme_cellpress_size/3, col="black") +
  geom_point(data=subset(plt, col == "yes"), size=theme_cellpress_size/3, col="red") +
  
  geom_hline(yintercept=0, col="red", alpha=0.5, lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="red", alpha=0.5, lwd=theme_cellpress_lwd) +
  
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)), col="1", cor.coef.name ="R", size=theme_cellpress_size, family=theme_cellpress_font_family) +
  
  
  ggrepel::geom_text_repel(data= subset(plt, col == "yes"), col="red", size=theme_cellpress_size, 
                           segment.size=theme_cellpress_lwd, family = theme_cellpress_font_family,
                           force=1,
                           max.overlaps=15,
                           #nudge_y = -1,
                           nudge_x = -2.75
  ) +
  
  theme_cellpress +
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +
  theme(aspect.ratio=1)


metadata.proteomics.glass_od  |>
  dplyr::filter(`DPA__GLASS_OD__prim-rec__pat_corrected__t` < -2.5)  |>
  dplyr::filter(`DPA__GLASS_OD__grade__pat_corrected__t` < -3)



# GLASS-OD x GLASS-NL ----
## CGC x CGC ----


plt.glass_od <- metadata.proteomics.glass_od |> 
  dplyr::filter(!is.na(DPA__GLASS_OD__CGC__t)) |> # excldued due to too many N/A's
  dplyr::mutate(protein_id = paste0(Protein.Names, ";",Genes)) |> 
  dplyr::mutate(protein_id = sapply(strsplit(protein_id, ";"), function(x) {
    cleaned <- unique(stringr::str_replace(x, "_HUMAN$", ""))
    paste(cleaned, collapse = ";")
  })) |> 

  tidyr::separate_rows(protein_id, sep = ";") |> 
  dplyr::mutate(protein_id =gsub("_HUMAN$","", protein_id)) |> 

  dplyr::mutate(protein_id = ifelse(protein_id == "H1-0", "H1F0", protein_id)) |> 
  
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-0", "H1F0", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-1", "HIST1H1A", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-10", "H1FX", protein_id)) |> 
  
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-2", "HIST1H1C", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-4", "HIST1H1E", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-5", "HIST1H1B", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-6", "HIST1H1T", protein_id)) |> 
  
  dplyr::mutate(protein_id = ifelse(protein_id == "H2AC19", "HIST2H2AA", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2AC20", "HIST2H2AC", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2AC21", "HIST2H2AB", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2AC25", "HIST3H2A", protein_id)) |> 
  
  
  dplyr::mutate(protein_id = ifelse(protein_id == "H2AC6", "HIST1H2AC", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2AC8", "HIST1H2AE", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2AZ2", "H2AFV", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2BC1", "HIST1H2BA", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2BC11", "HIST1H2BJ", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2BC21", "HIST2H2BE", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2BC19P", "HIST2H2BD", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2BC20P", "HIST2H2BC", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2BC5", "HIST1H2BD", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H3-3B", "H3F3B", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H3C12", "HIST1H3J", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H3C13", "HIST2H3D", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H3-7", "HIST2H3PS2", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H4C16", "HIST4H4", protein_id)) |> 
  
  dplyr::select(protein_id, starts_with("DPA__GLASS_OD__CGC")) 


plt.glass_nl <- metadata.proteomics.glass_nl |> 
  dplyr::filter(!is.na(DPA__GLASS_NL__CGC__t)) |> # excldued due to too many N/A's
  dplyr::select(gene_id, starts_with("DPA__GLASS_NL__CGC"))




plt <- dplyr::inner_join(
  plt.glass_od,
  plt.glass_nl, by=c('protein_id'='gene_id'))


print(dim(plt))


plt <- plt |> 
  dplyr::mutate(col = dplyr::case_when(
    grepl("^HOX[ABCD]", protein_id) ~ "HOX",
    T ~ "other"
  )) |> 
  dplyr::mutate(cycling = protein_id %in% c(
    "PCNA", "PNCA", "MKI67", "CDC6", "PLK1", "MCM3" , "MIB","MIB1",'MIB-1' , "MCM4" ,  "MCM5" ,  "MCM7"  , "MCM2" ,  "MCM6", "MCM3","MCM7",
    "TMPO","DEK","RFC4",
    "HMGB3", "HMGB1", "HMGB2",
    "BUB3"
  ))  |> 
  dplyr::mutate(is_fibronectin_fn1 = protein_id %in% c("FN1", "FINC", "FIBRONETIN", "FIBRONETIN1", "FN", "CIG", "GFND2", "LETS", "MSF")) |> 
  dplyr::mutate(cycling = protein_id %in% c(
    "PCNA", "PNCA", "MKI67", "CDC6", "PLK1", "MCM3" , "MIB","MIB1",'MIB-1' , "MCM4" ,  "MCM5" ,  "MCM7"  , "MCM2" ,  "MCM6", "MCM3","MCM7",
    "TMPO","DEK","RFC4",
    "HMGB3", "HMGB1", "HMGB2",
    "BUB3"
  )) |> 
  dplyr::mutate(GO_0002250_adaptive_immune_response = c(
    protein_id %in% c(
      "TMEM98", "CD79B", "MASP2", "CD4", "BTK", "HFE", "FYN", "ERCC1", "KDM5D", "CD6", "UFL1", "WAS", "SLC11A1", "CD74", "C8B", "SLAMF7", "BTN3A1", "TNFRSF1B", "C6", "PARP3", "TNFRSF17", "FOXP3", "RC3H2", "TRAF3IP2", "PRDM1", "PRKCQ", "CD84", "TP53BP1", "PRKCZ", "TFE3", "BCL3", "RORA", "RAB27A", "TFRC", "FCGR2B", "PVR", "TBX21", "TXK", "MLH1", "UNG", "PAG1", "RAP1GAP", "NFKB2", "IL4R", "LAMP3", "RAPGEF3", "CEACAM1", "STX7", "RIF1", "ARG2", "MEF2C", "PTPRC", "LAT2", "P2RX7", "LAG3", "ICAM1", "CD209", "MCOLN1", "RAPGEF4", "UNC13D", "MSH2", "TREM2", "SIRT1", "JAK2", "IL12RB1", "KLHL22", "LGALS1", "CSF2RB", "IL2RB", "EP300", "CD40", "SLA2", "JAG1", "SMAD7", "RNF125", "TLR8", "BMX", "CD40LG", "TNFSF13B", "PYCARD", "IL21R", "CSK", "CTSH", "RIPK2", "NBN", "FGL1", "RELB", "FCER2", "CLEC4M", "IL4I1", "LILRA1", "IL27RA", "CLC", "EBI3", "PRKD2", "TGFB1", "CD79A", "TYK2", "LILRB5", "JAK3", "PIK3CG", "EPHB6", "AHR", "C5", "EXOSC3", "GATA3", "PPP3CB", "C1QBP", "SUPT6H", "IL2", "NSD2", "CTSC", "CRTAM", "HSPA8", "UNC93B1", "KMT5B", "HPX", "CD81", "TCIRG1", "CD69",
      "IL23A", "IFNG", "PTPN6", "CLEC4A", "AICDA", "BTN3A3", "ULBP1", "IL17A", "IL17F", "RNF8", "BACH2", "CCR6", "TFEB", "DUSP22", "C7", "ITK", "IL12B", "BTNL8", "IL4", "C9", "BCL6", "CD86", "ZAP70", "LOXL3", "IL1R1", "IL1RL1", "IL18R1", "MSH6", "RNF19B", "MAD2L2", "PLA2G4A", "SLAMF1", "CD160", "CR2", "CD46", "MPL", "TNFSF4", "TNFAIP3", "ARG1", "CD274", "IFNA6", "IFNA8", "TNFSF18", "CLU", "PTK2B", "ADCY7", "KCNJ8", "CD80", "CCR2", "SASH3", "LAX1", "CD244", "LY9", "SHLD2", "PMS2", "PKN1", "NCKAP1L", "IL13RA2", "C4BPA", "C4BPB", "IL9R", "TREM1", "IRF1", "IL1B", "CD70", "C3", "STAT5A", "ZBTB1", "TRAF2", "IGLL1", "EIF2AK4", "RIPK3", "FOXJ1", "CRACR2A", "NECTIN2", "TRPM4", "ULBP2", "ULBP3", "LILRB2", "PRR7", "RFTN1", "NDFIP1", "PRKAA1", "RAP1GAP2", "JCHAIN", "CLEC10A", "CTNNBL1", "EPHB2", "KMT5C", "SWAP70", "VTCN1", "RSAD2", "IL6ST", "KLRD1", "KLRC1", "ANXA1", "CTSL", "HAVCR2", "AKIRIN2", "MAP3K7", "TEC", "RC3H1", "IL6", "DBNL", "MYO1G", "HLX", "IL10", "TLR4", "IL33", "SIT1", "IFNA21", "IRF4", "IL18BP", "STAT4", "LEF1", "C1RL", "CD27", "SLC15A4", "SKAP1", "SECTM1", "TNFRSF11A", "BRD4", "SIGLEC10", "BCL10", "XCL1", "FCGR2A", "RORC", "SUSD4", "DUSP10", "HSPD1", "NFKBIZ", "IL9", "FBXO38", "TNFRSF21", "EBAG9", "IFNA5", "IFNA16", "IFNK", "NOTCH1", "SERPING1", "FCGR1A", "CD226", "IL18", "ADAM17", "DCLRE1C", "CAMK4", "MR1", "CYRIB",
      "CD8A", "MCOLN2", "SAMSN1", "RAET1L", "BATF", "CXCL13", "C8A", "PAXIP1", "TNFRSF14", "CD1D", "CD1A", "CD1C", "CD1B", "CD1E", "HMHB1", "FCER1G", "C1QC", "C1R", "TNFRSF13C", "ICOSLG", "AIRE", "CD3G", "ZBTB7B", "IL6R", "PCYT1A", "ALOX15", "TNFSF13", "JAK1", "IL23R", "NLRP3", "SLAMF6", "FCGR3B", "FCMR", "FCAMR", "SANBR", "CTSS", "FZD5", "EOMES", "FCRL4", "TRAT1", "CTLA4", "ZC3H12A", "PRKCD", "RNF168", "ERAP1", "ERAP2", "RAET1E", "SYK", "SVEP1", "MARCHF8", "MBL2", "TSC1", "HPRT1", "JAM3", "C2", "RAG1", "PRKCB", "CLEC4D", "B2M", "STAT6", "PDIA3", "PHB1", "NOD2", "CD3D", "NFKBID", "LAIR1", "OTUB1", "FADD", "CX3CR1", "TAP1", "STAT3", "IL7R", "IL12A", "INPP5D", "GPR183", "APLF", "ALCAM", "YWHAG", "SOCS5", "JUNB", "FGA", "FGB", "PIK3CD", "IFNB1", "SHLD1", "CD8B", "MALT1", "CLEC7A", "KLHL6", "THEMIS", "CCL19", "MYD88", "SLC22A13", "C1QA", "CD7", "EXO1", "LIG4", "IL20RB", "HRAS", "ADGRE1", "TRAF6", "CLCF1", "ZNF683", "ATAD5", "C8G", "IFNW1", "CD19", "IL17RA", "CD28", "HLA-DQB1", "FCER1A", "FUT7", "PRF1", "NLRP10", "C1S", "CLEC4G", "HMCES", "PSG9", "ASCL2", "SH2D1A", "SOCS3", "IFNE", "IRF7", "BTLA", "BTN3A2", "IFNA10", "LILRB4", "CARD9", "C17orf99", "NCR3LG1", "ZP3", "IFNA2", "PDCD1", "IGHV1OR15-9", "HMGB1", "HLA-DRB1", "SEMA4A", "CD55", "HLA-DQA1", "ADA", "ARID5A", "IL27", "SPN", "GZMM", "MPEG1", "PDCD1LG2", "CR1L", "IFNA1", "ENTPD7", "CLEC4C", "HLA-DRB5", "SH2D1B", "MTOR", "CD247", "OPA1", "CD3E", "CR1", "RAET1G", "FCGR3A", "LIME1", "HLA-DOA", "BRD2", "HLA-DMA", "TAP2", "HLA-DRA", "AGER", "MICB", "HLA-C", "LILRB3", "GNL1", "HLA-E", "HLA-G", "HLA-F", "TRIM27", "CFI", "KLRC2", "CLEC6A", "HLA-A", "IGKJ1", "IGKV4-1", "IGKV5-2", "IGKV6-21", "IGKV2D-26", "IGKV3D-20", "IGKV6D-41", "IGKV3D-11", "IGKV1D-42", "IGLV4-69", "IGLV8-61", "IGLV4-60", "IGLV6-57", "IGLV11-55", "IGLV10-54", "IGLV5-52", "IGLV1-51", "IGLV1-50", "IGLV5-48", "IGLV1-47", "IGLV7-46", "IGLV5-45", "IGLV1-44", "IGLV7-43", "IGLV1-40", "IGLV5-37", "IGLV1-36", "IGLV2-33", "IGLV3-32", "IGLV3-27", "IGLV3-25", "IGLV2-23", "IGLV3-22", "IGLV3-21", "IGLV3-19", "IGLV2-18", "IGLV3-16", "IGLV2-14", "IGLV3-12", "IGLV2-11", "IGLV3-10", "IGLV3-9", "IGLV4-3", "IGLV3-1", "IGLJ1", "IGLC2", "IGLC3", "TRGC1", "TRGJ1", "TRGV11",
      "TRGV10", "TRGV9", "TRGV8", "TRGV5", "TRGV4", "TRGV3", "TRGV1", "TRBV6-1", "TRBV7-1", "TRBV4-1", "TRBV6-4", "TRBV7-3", "TRBV5-3", "TRBV9", "TRBV10-1", "TRBV11-1", "TRBV6-5", "TRBV6-6", "TRBV5-5", "TRBV7-6", "TRBV5-6", "TRBV5-7", "TRBV5-1", "TRBV4-2", "TRBV19", "TRBV20-1", "TRBV23-1", "TRBV24-1", "TRBV27", "TRBV28", "TRBJ2-1", "TRBJ2-3", "TRBJ2-7", "TRAV2", "TRAV3", "TRAV4", "TRAV5", "TRAV6", "TRAV7", "TRAV8-1", "TRAV9-1", "TRAV10", "TRAV12-1", "TRAV8-2", "TRAV8-3", "TRAV13-1", "TRAV12-2", "TRAV8-4", "TRAV13-2", "TRAV14DV4", "TRAV9-2", "TRAV12-3", "TRAV8-6", "TRAV16", "TRAV17", "TRAV18", "TRAV19", "TRAV20", "TRAV21", "TRAV22", "TRAV23DV6", "TRDV1", "TRAV24", "TRAV25", "TRAV26-1", "TRAV27", "TRAV29DV5", "TRAV26-2", "TRAV34", "TRAV36DV7", "TRAV38-1", "TRAV38-2DV8", "TRAV39", "TRAV40", "TRAV41", "TRDV2", "TRDJ1", "TRAJ42", "TRAJ31", "TRAJ3", "IGHA2", "IGHE", "IGHG4", "IGHG2", "IGHA1", "IGHG1", "IGHG3", "IGHM", "IGHJ1", "IGHV6-1", "IGHV1-2", "IGHV1-3", "IGHV2-5", "IGHV3-7", "IGHV3-11", "IGHV3-13", "IGHV3-15", "IGHV3-16", "IGHV1-18", "IGHV3-20", "IGHV3-21", "IGHV3-23", "IGHV1-24", "IGHV2-26", "IGHV4-28", "IGHV3-33", "IGHV4-34", "IGHV3-35", "IGHV3-38", "IGHV4-39", "IGHV1-45", "IGHV1-46", "IGHV3-48", "IGHV3-49", "IGHV5-51", "IGHV3-53", "IGHV1-58", "IGHV4-61", "IGHV3-66", "IGHV1-69", "IGHV2-70D", "IGHV3-73", "IGHV7-81", "DENND1B", "SIPA1", "LAT", "TREX1", "KLRK1", "EMP2", "IFNA7", "SCART1", "IGLV9-49", "EXOSC6", "IGHV3-64", "HLA-DPB1", "TRDD1", "IGKV3D-15", "IGHV4-59", "C4B", "IGHV3-74", "IGKV6D-21", "IGHV3-72", "TRBV2", "LTA", "TRGC2", "IFNA14", "IGKV3D-7", "TRBV10-2", "TRBV5-4", "HLA-DPA1", "IGHV4-31",
      "IGHV3-43", "HLA-DQB2", "TNF", "TRBV29-1", "TRGV2", "IGHV3OR16-10", "IFNA13", "HLA-B", "IFNA17", "IGHD1-1", "IFNA4", "TRBV30", "HLA-DQA2", "TRBV3-1", "IGKV2D-30", "TNFSF12", "IGKV1D-8", "IGKV1-6", "IGKV1-37", "IGKV3-20", "IGKV1D-33", "IGKV1-17", "TNFRSF13B", "IGKV1-8", "IGKV1-16", "HLA-DOB", "IGKV1D-16", "IGKV2-24", "IGKV3-11", "IGKV2D-24", "IGKV1-9", "IGKV1-33", "IGKV1-39", "IGKV2D-28", "HLA-DMB", "IGKV1D-43", "IGKV1D-17", "IGKV3-7", "IGKV2-30", "IGKV2D-29", "IGKV1-12", "IGKV1-5", "CFB", "IGKV2-28", "IGKV3-15", "LILRA6", "IGKV1-27", "C4A", "TARM1", "IGKV1D-37", "IGKV2D-40", "IGKV1D-39", "TRBV6-7", "SHLD3", "TRBV7-7", "TRBV7-4", "TRBV6-8", "LYN", "CD8B2", "IGHV8-51-1", "IGLL5", "TRAV1-1", "TRAV1-2", "TRDV3", "TRAV30", "IGHV4OR15-8", "IGHV2OR16-5", "IGHV3OR15-7", "IGHV3OR16-17", "OTUD7B", "IGHV3OR16-12", "IGHV3OR16-9", "IGHV1OR15-1", "IGHV3-30", "IGHV3OR16-8", "IGHV3OR16-13", "IGKV2-40", "IGHV2-70", "TRBV12-3", "TRBV12-5", "TRBV16", "TRBV14", "TRBV10-3", "ORAI1", "TRBV13", "TRBV18", "IGKV1D-13", "TRBV11-3", "IGHV4-4", "TRBV12-4", "IGHV1OR21-1", "TRBV17", "TRBV7-9", "IGLV2-8", "IGKV1D-12", "IGHV1-69D", "TRBJ1-4", "IGHV1-69-2", "IGHV7-4-1", "TRBJ1-3", "TRBJ1-5", "TRBJ1-1", "TRBJ1-2", "TRBD1", "TRBV25-1", "IGHV3-64D", "IGHV5-10-1", "TRBJ1-6", "TRBV7-2", "TRBV6-2", "IL9R" 
    )
  )) |> 
  dplyr::mutate(GO_0002526_acute_inflammatory_response = protein_id %in% c(
    "BTK","CD6","PTGER3","ITIH4","CREB3L3","TFRC","FCGR2B","PTGS2","GSTP1","B4GALT1","EIF2AK1","OSM","RHBDD3","PIK3CG","EPHB6","TFR2","LIPA","VNN1","IL4","IL1A","ACVR1","FN1","PARK7","ASH1L","PLA2G2D","F3","TNFSF4","CNR1","TNFSF11","IL1B","C3","FFAR2","CCR7","IL22","APOL2","LBP","EPO","ASS1","F12","SELENOS","PPARG","CRP","APCS","ALOX5AP","SAA2","IL6ST","EDNRB","IL6","IL1RN","PRCP","MYLK3","TNFRSF11A","REG3G","AHSG","OSMR","PTGES","SAA4","FCGR1A","ADAM8","IL6R","NLRP3","ADORA1","DNASE1L3","NPY5R","KLKB1","IL31RA","MBL2","SERPINF2","SCN11A","MRGPRX1","REG3A","CEBPB","SAA1","IL20RB","NLRP6","A2M","NUPR1","ANO6","CD163","CTNNBIP1","FCER1A","F2","FUT7","CXCR2","EXT1","F8","SIGIRR","FFAR3","PLSCR1","ZP3","SERPINA3","TRPV1","SERPINA1","SPN","ELANE","C2CD4A","FCGR3A","HLA-E","C2CD4B","IGHG1","DNASE1","ORM2","ORM1","TNF","UGT1A1","INS","SAA2-SAA4","HP"
  )) |>
  dplyr::mutate(GO_0062023_collagen_containing_extracellular_matrix = protein_id %in% c(
    "DCN","SEMA3B","MARCO","SERPINB1","MYOC","TIMP2","VCAN","TNC","USH2A","COL9A2","LTBP1","ELN","LAMC3","COL23A1","LAMA3","ENTPD2","CBLN4","ITIH4","ITIH1","F7","LAMC2","COL11A1","WNT8A","GPC1","CDON","NTN1","COL17A1","FGFR2","PKM","FGF10","SPP2","ADAM11","NTN4","LMAN1","DLG1","GPC4","FBLN1","COL5A3","HSP90AA1","COL4A4","IMPG2","COL19A1","COL16A1","FCN1","ACHE","ADAMTS2","MMP2","NID2","LTBP4","ICAM1","P3H2","LAMB4","LAMB1","APOH","ANGPT2","CCDC80","CMA1","COL9A3","TGFB2","AMELY","HNRNPM","LGALS1","TIMP3","PDGFB","CHADL","CTSG","COCH","MMP9","BMP7","CTSZ","COL20A1","ANGPT4","LAMA1","MXRA5","SRPX","F9","TIMP1","SRPX2","CBLN1","ZP2","CTSH","SFRP1","IL7","FGL1","CLC","APLP1","TGFB1","COMP","WNT2","PTPRZ1","PCOLCE","SERPINE1","PLOD3","RARRES2","AEBP1","MEGF9","OGN",
    "ASPN","ECM2","AMBP","CXCL12","ACTA2","KAZALD1","LGALS3BP","COL1A1","VTN","SOD3","CTSC","HPX","APOA4","APOC3","MDK","VWF","WNT5B","MGP","COL12A1","COL9A1","SMOC2","IMPG1","LAMA4","ERBIN","LOX","SPARC","THBS4","KNG1","HRG","WNT5A","COL7A1","LOXL3","EFEMP1","FN1","TNR","ANGPTL1","PRG4","NID1","ZP4","MFAP2","P3H1","F3","SERPINC1","CTSD","MMP8","APOA1","CCN2","LTBP2","TGFB3","TNN","TGFBI",
    "CLU","A1BG","FMOD","PLG","ANXA11","COL10A1","SLPI","MATN4","NDP","F13A1","SERPINB6","COL21A1","AMELX","CFP","PZP","FGL2","LOXL1","APOE","NCAN","PXDN","GDF15","COL5A1","LAMA5","FIBCD1","ANGPTL6","F12","LGALS3","PODNL1","RTBDN","MATN3","EMILIN2","SERPINF1","ITGB4","MATN2","BCAN","HAPLN2","APCS","ANGPTL3","POSTN","LOXL2","ADAMDEC1","WNT2B","COL4A2","ADAMTS8","ANXA1","CTSL","ADAM19","AGT","LAMC1","SERPINE2","ANGPTL2","CCN3","IGFBPL1","TINAG","SULF1","THBS1",
    "EMILIN1","LOXL4","ANXA7","CILP","SEMA7A","MMRN1","FRAS1","FBN2","COL2A1","INHBE","LUM","FBLN5","PCSK6","LMAN1L","HAPLN3","MFGE8","TGFB1I1","CDH13","COL6A1","COL6A2","ADAMTS10","FCN3","HSPG2","CCN1","TINAGL1","DPT","HDGF","HMCN1","ECM1","ADAMTSL4","ANXA9","S100A8","S100A7","FLG","LEFTY2","COL8A1","AHSG","SFRP2","HAPLN1","CASK","GPC3","HMCN2","SERPING1","SERPINH1","NCAM1","ZP1","FREM2","ITIH2","DST","SPARCL1","ABI3BP","ANGPT1","ADAMTS1","ADAMTS5","ADAMTS3","PRG3","ACAN","COLEC12","ADAMTS4","LAD1","C1QC","CSTB","FCN2","AZGP1","COL26A1","MATN1","SDC3","CTSS","S100A9","COL6A3","IGFBP7","FBLN2","ADAMTS9","CRELD1","PF4","CPA3","TGM4","CLEC3B","ANXA5","EDIL3","EGFLAM","RELL2","SHH","COL1A2","CTSB","SBSPON","CTHRC1","FREM1","VWA2","HTRA1","ADAMTS15","FBN1","SERPINB8","MFAP4","HSP90B1","SERPINB12","NAV2","GREM1","SERPINF2","KRT1","ANGPTL4",
    "SOST","LTBP3","TNXB","COL3A1","NPNT","LGALS9","COL4A3","THBS3","COL22A1","SDC2","MUC17","SERPINB9","CDH2","COL24A1","FGG","FGA","FGB","LGALS4","COL8A2","ANGPTL7","LAMB2","LRRC15","TPSAB1","EFEMP2","COL6A5","EGFL7","ADAMTS20","MMRN2","C1QA","DAG1","CSPG4","CTSF","PODN","FGFBP3","ZG16","NPPA","A2M","CLEC14A","CD151","CALR","GPC5","VWA1","SSC5D","F2","ADIPOQ","OTOL1","BGN","MXRA7","SNORC","ANXA2","COL18A1","FREM3","GPC6","EMILIN3","EFNA5","THBS2","PRG2","C17orf58","ANGPTL5","BCAM","COL4A1","HAPLN4","AMTN","THSD4","COL14A1","EYS","COL4A5","AGRN","PLSCR1","ZP3","IFNA2","SERPINA5","COL25A1","VWC2","PRELP","MMP23B","SERPINA3","S100A4","PRTN3","LAMA2","COL27A1","LAMB3","ANXA4","ANXA6","SERPINA1","TPSB2","COL13A1","SPN","ELANE","COL4A6","MFAP5","PSAP","S100A10","HRNR","S100A6","SMOC1","EGFL6","L1CAM","TGM2","COL11A2","COL5A2","COL15A1","PRSS1","VIT","DEFA1","COL6A6","COLQ","GPC2","ANG","COL28A1","ORM2","ORM1","DEFA1B","MARCOL","GH1","SPON1","ANXA8","RBP3","GDF10","MMP28"
  )) |> 
  dplyr::mutate(`GO_0072562_blood_microparticle` = protein_id %in% c(
    "CFH","SLC4A1","PON1","ITGA2B","CP","ITIH4","ITIH1","GRIPAP1","TFRC","CD5L","ACTB","AFM","PSMC5","TF","APOL1","TGFB1","AMBP","ENG","PFN1","LGALS3BP","VTN","HSPA8","HPX","APOA4","C9","KNG1","HRG","BCHE","FN1","CFHR3","SLC2A1","SERPINC1","APOA1","CLU","A1BG","PLG","C4BPA","DNPEP","F13A1","C3","HSPA2","PZP","INTS11","APOE","JCHAIN","APCS","AGT","CIB2","SDCBP","TMPRSS13","FCN3","OAZ3","ACTA1","EIF2A","AHSG","GC","MSN","STOM","GSN","SERPING1","ITIH2","C8A","APOA2","C1QC","ACTC1","FCN2","ACTG2","ALB","ANXA5","YWHAZ","ACSM1","SERPINF2","KRT1","ANGPTL4","FGG","FGA","FGB","HSPA6","A2M","C8G","ZBTB38","CPN2","F2","C1S","ACTG1","PROS1","KDM4D","POTEE","HBA2","ZNF177","SERPINA3","HBG2","POTEF","HSPA1B","HSPA1A","HSPA1L","PRSS1","HBA1","IGKV4-1","IGLV1-47","IGLV3-25","IGLV3-21","IGLC2","IGLC3","IGHA2","IGHG4","IGHG2","IGHA1","IGHG1","IGHM","IGHV3-7","IGHV3-13","IGHV3-23","CLIC1","HBE1","HBD","C4B","ORM2","ORM1","IGKV3-20","IGKV1D-33","IGKV1-17","IGKV3-11","IGKV1-33","IGKV1-39","IGKV2D-28","IGKV2-30","IGKV1-5","CFB","CFHR1","IGKV3-15","C4A","HBB","IGKV2D-40","HP","HPR","IGKV1D-12"
  )) |> 
  dplyr::mutate(homeobox_tf = protein_id %in% c(genes_homeobox_nonhox,
                                                genes_polycomb_eed_homeobox,
                                                genes_polycomb_h3k27_homeobox,
                                                genes_polycomb_prc2_homeobox,
                                                genes_polycomb_suz12_homeobox
  )) |> 
  dplyr::mutate(od = protein_id %in% c("MAG", "MBP", "SOX10", "NKX2.2", "CNP")) |> 
  dplyr::mutate(cycling2 = protein_id %in% c("MCM7", "MCM5", "MCM3", "MCM2", "MCM4", 
                                             "CCNB2", "PCNA","CCNB1", "PLK1", "BUB1","MKI67","CCND3","MYBL2","E2F1", "CCNE2", "CCND2", "CCNE1", "CCND1", "RPM1", "RFC2","PRIM2","RRM2",
                                             "MCM6"
  ))




# 
# plt |> 
#   dplyr::filter(grepl("^CNP", protein_id)) |> 
#   dplyr::pull(protein_id)




plt <- plt  |> 
  dplyr::mutate(col = dplyr::case_when(
    cycling ~ "Proliferation",
    GO_0062023_collagen_containing_extracellular_matrix ~ "GO:0062023 ECM",
    GO_0002250_adaptive_immune_response ~ "GO:002250: Adpt. immune resp.",
    T ~ "Other"
  ))
plt$col |> table()



ggplot(plt, aes(x=DPA__GLASS_OD__CGC__t, y=DPA__GLASS_NL__CGC__t, col=col, label=protein_id)) + 
  geom_hline(yintercept=0, col="red", lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="red", lwd=theme_cellpress_lwd) +
  
  geom_point(data=subset(plt, col == "Other"), size=theme_cellpress_size/3 / 5) +
  geom_point(data=subset(plt, col != "Other"), size=theme_cellpress_size/3 / 5) +
  
  geom_hline(yintercept=0, col="red", alpha=0.1, lwd=theme_cellpress_lwd) +
  geom_vline(xintercept=0, col="red", alpha=0.1, lwd=theme_cellpress_lwd) +
  
  ggrepel::geom_text_repel(data= subset(plt, col == "Proliferation"), col="black", size=theme_cellpress_size, 
                           segment.size=theme_cellpress_lwd, family = theme_cellpress_font_family,
                           force=1,
                           max.overlaps=15,
                           #nudge_y = -1,
                           nudge_x = -2.75
  ) +
  scale_color_manual(values = c(`Other`= 'darkgray',
                                `Proliferation`='blue',
                                `GO:0062023 ECM` = 'red',
                                `GO:002250: Adpt. immune resp.` = 'darkgreen')) +
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)), col="1", cor.coef.name ="R", size=theme_cellpress_size, family=theme_cellpress_font_family) +
  labs(col = NULL, x="t-statistic proteomics GLASS-OD", y="t-statistic proteomics GLASS-NL") +
  
  theme_cellpress +
  
  theme(plot.background = element_rect(fill="white", colour=NA)) +
  theme(aspect.ratio=1)



ggsave(paste0("output/figures/vis_differential_proteomics__GLASS-OD_CGC__x__GLASS-NL_CGC.png"),
       width=2.1 - 0.4, height=3.9 -0.4, dpi=1200)




# GLASS-OD: prot x dna meth ----


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



#plt.dna.meth <- gene_enrichment_0.GencodeCompV12_NAME
plt.dna.meth <- out.per.gene.GencodeCompV12_NAME
plt.proteomics <- metadata.proteomics.glass_od



plt.dna.meth |> dplyr::filter(grepl("^(HIST)", gene))
plt.proteomics |> 
  dplyr::mutate(Genes = ifelse(Genes == "H1-0", "H1F0", Genes)) |> 
  
  dplyr::mutate(Genes = ifelse(Genes == "H1-0", "H1F0", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H1-1", "HIST1H1A", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H1-10", "H1FX", Genes)) |> 
  
  dplyr::mutate(Genes = ifelse(Genes == "H1-2", "HIST1H1C", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H1-4", "HIST1H1E", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H1-5", "HIST1H1B", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H1-6", "HIST1H1T", Genes)) |> 
  
  dplyr::mutate(Genes = ifelse(Genes == "H2AC19", "HIST2H2AA", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H2AC20", "HIST2H2AC", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H2AC21", "HIST2H2AB", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H2AC25", "HIST3H2A", Genes)) |> 
  
  
  dplyr::mutate(Genes = ifelse(Genes == "H2AC6", "HIST1H2AC", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H2AC8", "HIST1H2AE", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H2AZ2", "H2AFV", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H2BC1", "HIST1H2BA", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H2BC11", "HIST1H2BJ", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H2BC21", "HIST2H2BE", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H2BC19P", "HIST2H2BD", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H2BC20P", "HIST2H2BC", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H2BC5", "HIST1H2BD", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H3-3B", "H3F3B", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H3C12", "HIST1H3J", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H3C13", "HIST2H3D", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H3-7", "HIST2H3PS2", Genes)) |> 
  dplyr::mutate(Genes = ifelse(Genes == "H4C16", "HIST4H4", Genes)) |> 
  
  dplyr::filter(grepl("^(HIST)", Protein.Ids) | grepl("^(HIST)", Genes))



join_ids <- plt.proteomics |> 
  dplyr::select(Protein.Names, Genes) |> 
  dplyr::mutate(protein_id = paste0(Protein.Names, ";",Genes)) |> 
  dplyr::mutate(protein_id = sapply(strsplit(protein_id, ";"), function(x) {
    cleaned <- unique(stringr::str_replace(x, "_HUMAN$", ""))
    paste(cleaned, collapse = ";")
  })) |> 
  tidyr::separate_rows(protein_id, sep = ";") |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-0", "H1F0", protein_id)) |> 
  
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-0", "H1F0", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-1", "HIST1H1A", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-10", "H1FX", protein_id)) |> 
  
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-2", "HIST1H1C", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-4", "HIST1H1E", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-5", "HIST1H1B", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H1-6", "HIST1H1T", protein_id)) |> 
  
  dplyr::mutate(protein_id = ifelse(protein_id == "H2AC19", "HIST2H2AA", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2AC20", "HIST2H2AC", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2AC21", "HIST2H2AB", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2AC25", "HIST3H2A", protein_id)) |> 
  
  
  dplyr::mutate(protein_id = ifelse(protein_id == "H2AC6", "HIST1H2AC", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2AC8", "HIST1H2AE", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2AZ2", "H2AFV", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2BC1", "HIST1H2BA", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2BC11", "HIST1H2BJ", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2BC21", "HIST2H2BE", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2BC19P", "HIST2H2BD", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2BC20P", "HIST2H2BC", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H2BC5", "HIST1H2BD", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H3-3B", "H3F3B", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H3C12", "HIST1H3J", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H3C13", "HIST2H3D", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H3-7", "HIST2H3PS2", protein_id)) |> 
  dplyr::mutate(protein_id = ifelse(protein_id == "H4C16", "HIST4H4", protein_id)) 

join_ids <- join_ids |> 
  dplyr::left_join(
    plt.dna.meth |>
      dplyr::mutate(gene_meth = gene) |> 
      dplyr::select(gene, gene_meth)
    , by=c('protein_id'='gene'), suffix=c('','')
  )

join_ids <- join_ids |> dplyr::filter(!is.na(gene_meth))
join_ids <- join_ids |> 
  dplyr::filter(protein_id != "MYL12B") |> #unclear between A and B
  dplyr::filter(protein_id != "MAP1LC3B2") |> #unclear dupli
  dplyr::filter(protein_id != "PLSCR1") |> #unclear dupli
  dplyr::filter(protein_id != "MGAM") |> #unclear dupli
  dplyr::filter(protein_id != "SLC25A20") |> #unclear dupli
  dplyr::filter(protein_id != "EIF1B") |>  #unclear dupli
  dplyr::filter(protein_id != "POLR1A") |>  #unclear dupli
  dplyr::filter(protein_id != "CALM2") |>  #unclear dupli
  dplyr::filter(protein_id != "CALM3") |>  #unclear dupli
  dplyr::filter(protein_id != "IFITM2") |>  #unclear dupli
  dplyr::filter(protein_id != "IFITM3")|>  #unclear dupli
  dplyr::filter(protein_id != "TF") |> 
  dplyr::filter(protein_id %in% c("NTN1", "PF4V1","HMGN3","HSPA1B","BGN") == F) |> 
  dplyr::filter(protein_id %in% c("RAE1", "RBMS3","SSTR2", "FASN","RANBP2") == F) |> 
  dplyr::filter(protein_id %in% c("METAP2", "UBE2D4", "ACTG1","ARF3", "SMS", "MAGOHB") == F) |> 
  dplyr::filter(protein_id %in% c("RPL39P5", "ADAMTSL3", "PTRH2", "AMDHD2", "RABL2B", "AAA1", "DPAGT1", "CSN3") == F)|> 
  dplyr::filter(protein_id %in% c("POLR1B", "HDC", "RBM38", "LPCAT1", "CLDND1", "CTDSPL2", "PLSCR3", "PTH") == F) |> 
  dplyr::filter(protein_id %in% c("RALBP1", "ALG10B", "PREP", "PITRM1", "GPRIN1", "CADPS2", "HEATR5A", "SPRYD4", "SYNM") == F)|> 
  dplyr::filter(protein_id %in% c("HIST1H2AC","HIST1H2BJ","HIST1H3J", "HIST1H2AE") == F)


join_ids |> 
  dplyr::filter(Genes %in% (join_ids |> dplyr::filter(duplicated(Genes)) |>  dplyr::pull(Genes))) |> as.data.frame()

stopifnot(!duplicated(join_ids$Genes))



joined <- join_ids |> dplyr::select(Genes, gene_meth) |> 
  dplyr::left_join(plt.proteomics, by=c('Genes'='Genes')) |> 
  dplyr::left_join(plt.dna.meth, by=c('gene_meth'='gene'))


join_ids |> 
  dplyr::filter(grepl("^(HIST|H1|H2|H3)", protein_id))



ggplot(joined, aes(DPA__GLASS_OD__grade__pat_corrected__t, y=DMP__g2_g3__pp_nc_PC1__t__mean)) +
  geom_point(size = theme_cellpress_size/3) +
  geom_point(data=subset(joined, grepl("^HIST", gene_meth)), size = theme_cellpress_size/3, col='red') + 
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)), col="1", cor.coef.name ="R", size=theme_cellpress_size, family=theme_cellpress_font_family, col="red") + 
  #ylim(-13,9) +
  labs(x="Change in protein (grade)",y="Change in methylation (grade)") +
  theme_cellpress


ggplot(joined, aes(DPA__GLASS_OD__grade__pat_corrected__t, y=DMP__g2_g3__pp_nc_PC1__t__median)) +
  geom_point(size = theme_cellpress_size/3) +
  geom_point(data=subset(joined, grepl("^HIST", gene_meth)), size = theme_cellpress_size/3, col='red') + 
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)), col="1", cor.coef.name ="R", size=theme_cellpress_size, family=theme_cellpress_font_family, col="red") + 
  #ylim(-13,9) +
  labs(x="Change in protein (grade)",y="Change in methylation (grade)") +
  theme_cellpress

# intersections ----
## ALX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^ALX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()


## BARX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^BARX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## DBX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^DBX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## DLX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^DLX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## EMX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^EMX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## EN ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^EN[12]", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## EVX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^EVX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## GBX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^GBX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## HMX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^HMX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## HOXA ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^HOXA", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## HOT/HOXC ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^HO(T|XC)", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## HOXB ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^HOXB", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## HOXD ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^HOXD", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## IRX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^IRX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## ISL ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^ISL", Genes)) |>
  dplyr::pull(Genes) |> 
  table()


## LHX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^LHX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()


## LMX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^LMX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()


## MKX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^MKX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()



## MNX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^MNX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()


## MSX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^MSX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## NKX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^NKX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## NKX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^NKX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## SIX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^SIX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()


## TBX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^TBX", Genes)) |>
  dplyr::pull(Genes) |> 
  table()



## TBX ----
metadata.proteomics.glass_od |>
  tidyr::separate_rows(Genes, sep=";") |>
  dplyr::filter(grepl("^OR[0-9]+", Genes)) |>
  dplyr::pull(Genes) |> 
  table()

## OR ----
metadata.proteomics.glass_od |>
   tidyr::separate_rows(Genes, sep=";") |>
   dplyr::filter(grepl("^OR[0-9]+", Genes)) |>
   dplyr::pull(Genes) |> 
   table()

length(plt.meth_grade |> dplyr::filter(grepl("^OR[0-9]", UCSC_RefGene_Name)) |> dplyr::pull(UCSC_RefGene_Name) |> table())
