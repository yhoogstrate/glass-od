#!/usr/bin/env R
#' combined reading and normalization of GLASS-NL and GLASS-OD data

# load data ----


library(ggplot2)
library(minfi)

library(IlluminaHumanMethylation450kmanifest)
data(IlluminaHumanMethylation450kmanifest)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19) # BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#data(IlluminaHumanMethylation450kanno.ilmn12.hg19)


if(!exists('youri_gg_theme')) {
  source('scripts/youri_gg_theme.R')
}


if(!exists('tcga_laml.metadata.array_samples')) {
  source('scripts/load_TCGA-LAML_metadata.R')
}


source('scripts/load_functions.R')



metadata <- tcga_laml.metadata.array_samples |> 
  dplyr::filter(array_type == "450k")





# integrate & make m-values ----

targets <- metadata |> 
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","",array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL, array_channel_red = NULL) |> 
  dplyr::mutate(Sample_Name = array_sentrix_id) |>
  dplyr::mutate(Array = gsub("^.+_","", array_sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1",array_sentrix_id))


RGSet <- minfi::read.metharray.exp(targets = targets, force = T) #red/green channel together


detP <- minfi::detectionP(RGSet, type = "m+u")
proc <- minfi::preprocessNoob(RGSet, offset = 0, dyeCorr = TRUE, verbose = TRUE, dyeMethod="reference") 
rm(RGSet)
gc()


### m-values for all samples ----


mvalue <- minfi::ratioConvert(proc, what = "M")

stopifnot(rownames(mvalue) == rownames(detP))
stopifnot(colnames(mvalue) == colnames(detP))


mvalue <-  mvalue |> 
  assays() |> 
  purrr::pluck('listData') |> 
  purrr::pluck("M")


stopifnot(dim(mvalue) == dim(detP))


#mvalue.mask <- mvalue |> 
#  magrittr::multiply_by(ifelse(detP > 0.01 , NA, 1))

mvalue.mask <- ifelse(detP > 0.01 , NA, 1) |> 
  data.table::as.data.table(keep.rownames = "probe_id") |> 
  dplyr::filter(probe_id %in% (
    metadata.cg_probes.450k |> 
      dplyr::filter(exclude == F) |> 
      dplyr::pull(probe_id)
  )) |> 
  tibble::column_to_rownames('probe_id')

dim(mvalue.mask)


mvalue <- mvalue |> 
  data.table::as.data.table(keep.rownames = "probe_id") |> 
  dplyr::filter(probe_id %in% (
    metadata.cg_probes.450k |> 
      dplyr::filter(exclude == F) |> 
      dplyr::pull(probe_id)
  )) |> 
  tibble::column_to_rownames('probe_id')

dim(mvalue)




stopifnot(sum(is.na(mvalue)) == 0)
stopifnot(sum(is.na(mvalue.mask)) > 0)


stopifnot(targets$array_sentrix_id == colnames(mvalue))

# cleanup 

rm(detP, proc)
gc()


saveRDS(mvalue, "cache/mvalues.TCGA-LAML.Rds")
saveRDS(mvalue.mask, "cache/mvalues.TCGA-LAML.detP_mask.Rds")

rm(mvalue, mvalue.mask)
gc()




# test analysis ----

sel <- tcga_laml.metadata.array_samples |>  dplyr::filter(IDH & array_type == "450k") |> dplyr::pull(array_sentrix_id)
#sel %in% colnames(mval)

mval <- readRDS(file="cache/mvalues.TCGA-LAML.Rds") |> 
  dplyr::select(all_of(sel))

pc <- prcomp(t(mval), scale.=F)

# pc$rotation$PC

plt <- tcga_laml.metadata.array_samples |> 
  dplyr::left_join(
    pc$x|> 
      as.data.frame() |> 
      tibble::rownames_to_column('array_sentrix_id') ,
    by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')
  )


ggplot(plt, aes(x=PC1, y=PC2, col=IDH)) +
  geom_point()



pc_glassod <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_prcomp.Rds")



pc_aml <- pc$rotation[,1:15] |> 
  as.data.frame() |> 
  dplyr::rename_with( ~ paste0("aml_", .x)) |> 
  tibble::rownames_to_column('probe_id')

pc_oli <- pc_glassod$rotation  |> 
  as.data.frame() |> 
  dplyr::select(PC2) |> 
  dplyr::rename_with( ~ paste0("oli_", .x)) |> 
  tibble::rownames_to_column('probe_id')


plt_pc <- pc_aml |> 
  dplyr::inner_join(pc_oli, by=c('probe_id'='probe_id'))

dim(plt_pc)

ggplot(plt_pc, aes(x=aml_PC2, y=oli_PC1)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=aml_PC2, y=oli_PC2)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=aml_PC2, y=oli_PC3)) + geom_point(pch=19,cex=0.05, alpha=0.1)# sick
ggplot(plt_pc, aes(x=aml_PC2, y=oli_PC4)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=aml_PC2, y=oli_PC5)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=aml_PC2, y=oli_PC6)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=aml_PC2, y=oli_PC7)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=aml_PC2, y=oli_PC8)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=aml_PC2, y=oli_PC9)) + geom_point(pch=19,cex=0.05, alpha=0.1)


ggplot(plt_pc, aes(x=oli_PC2, y=aml_PC1)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=oli_PC2, y=aml_PC2)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=oli_PC2, y=aml_PC3)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=oli_PC2, y=aml_PC4)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=oli_PC2, y=aml_PC5)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=oli_PC2, y=aml_PC6)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=oli_PC2, y=aml_PC7)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=oli_PC2, y=aml_PC8)) + geom_point(pch=19,cex=0.05, alpha=0.1)
ggplot(plt_pc, aes(x=oli_PC2, y=aml_PC9)) + geom_point(pch=19,cex=0.05, alpha=0.1)


c = cor(plt_pc |>  tibble::column_to_rownames("probe_id"))
corrplot::corrplot(c, order="hclust")






