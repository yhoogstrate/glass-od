#!/usr/bin/env R

# libs ----

#install.packages('Seurat')
library(Seurat)


# load s[c/n]ATAC data ----


oligo_data_all <- readRDS("data/GSM6261343_oligodendroglioma_atac_seurat_object.rds")
# oligo_data_all |> Features() # bins (w/ names)
oligo_data_all$subtype |> table()
oligo_data_all$grade |> table() # who grades
oligo_data_all |> Idents() |> table() # cell types of matching RNA data
oligo_data_all$labelling |> table() # cell types of matching RNA data

oligo_celltypes <- data.frame("type"=   c("Astro-like", "Astrocytes","Endothelial","Gradient","Microglia","Mixed","Neurons","Oligodendrocytes","OPC-like","Pericytes","RA"),
                              "tumor" = c(T,            F,            F,            T,         F,          T,     F,         F,                 T,         F,          T))



astro_data_all <- readRDS("data/GSM6261342_astrocytoma_atac_seurat_object.rds")
# astro_data_all |> Features() # bins (w/ names)
astro_data_all |> Idents() |> table() # cell types of matching RNA data
astro_data_all$grade |> table() # who grades
astro_data_all$subtype |> table()
astro_data_all$labelling |> table() # cell types of matching RNA data

astro_celltypes <- data.frame("type"=   c("Astro-like","Astrocytes","Cycling-like","Endothelial","Excluded","Gradient","Microglia","Mixed","Neurons","Oligodendrocytes","OPC-like","Pericytes","RA","T-Cells"),
                              "tumor" = c(T,            F,           T,             F,            F,         T,         F,          T,      F,        F,                 T,         F,          T,   F))



length(oligo_data_all |> Features())
length(astro_data_all |> Features())

length(intersect(oligo_data_all |> Features(), astro_data_all |> Features()))




# sandbox code ----

# sample identifiers?
# https://pubmed.ncbi.nlm.nih.gov/37883975/
#gsub("_[ACTG]+-1$","",colnames(counts)) |> table()

# IDH_NCH536	Grade 2 = OD
# IDH_NCH2111	Grade 3 = OD *
# IDH_NCH781	Grade 3 = OD
# IDH_NCH6341	Grade 2 = OD *
# IDH_ACB_AD_809	Grade 2 = OD *
# IDH_ACB_AD_540	Grade 3 = OD
# IDH_NCH6702	Grade 3 = OD
# IDH_ACB_AD_883	Grade 3 = OD *


# OE0145-IDH_ACB_AD_809 OE0145-IDH_ACB_AD_883    OE0145-IDH_NCH2111    OE0145-IDH_NCH6341 
#    3509                   292                  2514                  1502 




# select tumor cells only ----


oligo_data_tumorcells <- oligo_data_all[,oligo_data_all$labelling %in% c(oligo_celltypes |>  dplyr::filter(tumor == T) |> dplyr::pull(type))]
dim(oligo_data_tumorcells)
table(as.character(oligo_data_tumorcells$labelling))
rm(oligo_data_all, oligo_celltypes)
gc()



astro_data_tumorcells <- astro_data_all[,astro_data_all$labelling %in% c(astro_celltypes |>  dplyr::filter(tumor == T) |> dplyr::pull(type))]
dim(astro_data_tumorcells)
table(as.character(astro_data_tumorcells$labelling))
rm(astro_data_all, astro_celltypes)
gc()



# calc per bin sum over all cells from oligo's only ----


oligo_data_pseudobulked <- data.frame(
  ATAC_oligo_bin = rownames(oligo_data_tumorcells),
  ATAC_oligo_chr = gsub("(chr.+)-[0-9]+-[0-9]+$","\\1",rownames(oligo_data_tumorcells)),
  ATAC_oligo_start = as.numeric(gsub("chr.+-([0-9]+)-.+$","\\1",rownames(oligo_data_tumorcells))),
  ATAC_oligo_end = as.numeric(gsub("chr.+-[0-9]+-([0-9]+)$","\\1",rownames(oligo_data_tumorcells))),
  ATAC_oligo_counts_per_bin = rowSums(oligo_data_tumorcells[["ATAC"]]$counts),
  ATAC_oligo_data_per_bin = sparseMatrixStats::rowMeans2(oligo_data_tumorcells[["ATAC"]]$data)
) |>
  dplyr::mutate(ATAC_oligo_size = ATAC_oligo_end - ATAC_oligo_start) |> 
  dplyr::mutate(ATAC_oligo_counts_per_bin_per_base = ATAC_oligo_counts_per_bin / ATAC_oligo_size) |> 
  dplyr::mutate(ATAC_oligo_data_per_bin_per_base = ATAC_oligo_data_per_bin / ATAC_oligo_size) |> 
  tibble::remove_rownames() |> 
  
  dplyr::filter(ATAC_oligo_size > 35) # odd outliers

rm(oligo_data_tumorcells)
gc()





astro_data_pseudobulked <- data.frame(
  ATAC_astro_bin = rownames(astro_data_tumorcells),
  ATAC_astro_chr = gsub("(chr.+)-[0-9]+-[0-9]+$","\\1",rownames(astro_data_tumorcells)),
  ATAC_astro_start = as.numeric(gsub("chr.+-([0-9]+)-.+$","\\1",rownames(astro_data_tumorcells))),
  ATAC_astro_end = as.numeric(gsub("chr.+-[0-9]+-([0-9]+)$","\\1",rownames(astro_data_tumorcells))),
  ATAC_astro_counts_per_bin = rowSums(astro_data_tumorcells[["ATAC"]]$counts),
  ATAC_astro_data_per_bin = sparseMatrixStats::rowMeans2(astro_data_tumorcells[["ATAC"]]$data)
) |>
  dplyr::mutate(ATAC_astro_size = ATAC_astro_end - ATAC_astro_start) |> 
  dplyr::mutate(ATAC_astro_counts_per_bin_per_base = ATAC_astro_counts_per_bin / ATAC_astro_size) |> 
  dplyr::mutate(ATAC_astro_data_per_bin_per_base = ATAC_astro_data_per_bin / ATAC_astro_size) |> 
  tibble::remove_rownames() |> 
  
  dplyr::filter(ATAC_astro_size > 35) # odd outliers

rm(astro_data_tumorcells)
gc()





# more sandbox ----

# ggplot(bin_data, aes(x=i, y=counts_per_bin_per_base)) +
#   geom_point(pch=19, cex=0.01)
# 
# 
# ggplot(bin_data, aes(x=counts_per_bin_per_base, y=size)) +
#   geom_point(pch=19, cex=0.01) +
#   ylim(0, 500)
# 
# 
# bin_data |> 
#   #head() |> 
#   dplyr::filter(ATAC_size > 40) |> # these small bins have odd high base counts
#   ggplot(aes(x=ATAC_counts_per_bin_per_base, y=ATAC_size)) +
#   geom_point(pch=19, cex=0.01) +
#   ylim(0, 5000)
# 


# bin_data |> 
#   #head() |> 
#   dplyr::filter(ATAC_size > 40) |> # these small bins have odd high base counts
#   ggplot(aes(x=ATAC_data_per_bin_per_base, y=ATAC_size)) +
#   geom_point(pch=19, cex=0.01) +
#   ylim(0, 5000)


# # https://stackoverflow.com/questions/75629990/lookup-table-in-r-by-matching-ranges
# exp <- bin_data |>  
#   dplyr::mutate(i = NULL) 
#   #GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)


# export ----



saveRDS(oligo_data_pseudobulked, file="cache/analysis_ATAC_oligo.Rds")
saveRDS(astro_data_pseudobulked, file="cache/analysis_ATAC_astro.Rds")




