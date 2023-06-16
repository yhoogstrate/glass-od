#!/usr/bin/env R


# load data ----


library(ggplot2)
library(minfi)


if(!exists('youri_gg_theme')) {
  source('scripts/youri_gg_theme.R')
}


source('scripts/load_metadata.R')


# obtain patients ----

patients <- glass_od.metadata.patients |> 
  dplyr::filter(is.na(reason_excluded))

resections <- glass_od.metadata.patients |> 
  dplyr::filter(is.na(reason_excluded))

stopifnot(resections$patient_id %in% patients$patient_id)




# load all samples, qc and corr w/ qc stats ----


array_samples <- glass_od.metadata.idats |> 
  dplyr::filter(is.na(reason_excluded_patient)) |> 
  dplyr::filter(is.na(reason_excluded_resection)) |> 
  dplyr::filter(is.na(reason_excluded_resection_isolation)) |> 
  dplyr::filter(is.na(reason_excluded_array_sample))



rm(RGSet)
gc()
targets <- array_samples |> 
  dplyr::rename(Sample_Name = resection_isolation_id) |>
  dplyr::mutate(Array = gsub("^.+_","",sentrix_id)) |> 
  dplyr::rename(Slide = methylation_array_chip_id) |> 
  dplyr::mutate(Basename = gsub("_Grn.idat$","", channel_green)) |> 
  dplyr::select(Sample_Name, sentrix_id, Array,Slide, Basename)

#RGSet <- read.metharray.exp(targets = targets, force = T) #red/green channel together
#saveRDS(RGSet, "cache/RGSet.Rds")
#RGSet <- readRDS("cache/RGSet.Rds")


#detP <- detectionP(RGSet, type = "m+u")

proc <- preprocessNoob(RGSet, offset = 0, dyeCorr = TRUE, verbose = TRUE, dyeMethod="reference") 
rm(RGSet)
gc()

## m-values ----


mvalue <- ratioConvert(proc, what = "M")

stopifnot(rownames(mvalue) == rownames(detP))
stopifnot(colnames(mvalue) == colnames(detP))


mvalue <-  mvalue |> 
  assays() |> 
  purrr::pluck('listData') |> 
  purrr::pluck("M") |> 
  magrittr::multiply_by(ifelse(detP > 0.01 , NA, 1)) |> 
  data.table::as.data.table(keep.rownames = "probeID") |> 
  dplyr::filter(probeID %in% (
    read.delim("~/mnt/neuro-genomic-1-ro/catnon/Methylation - EPIC arrays/EPIC.hg38.manifest.tsv.gz") |> 
      dplyr::filter(MASK_general == F) |> 
      dplyr::pull(probeID)
  ))


#saveRDS(mvalue, "cache/mvalues.Rds")


# cleanup 

rm(detP, proc)
gc()


# pca ----

mvalue <- readRDS("cache/mvalues.Rds") |> 
  tibble::column_to_rownames("probeID")


pca <- mvalue |> 
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |>  # this synthax, oh my
  dplyr::filter(mad > 1.0) |> 
  dplyr::mutate(mad = NULL) %>% # remove the sd to obtain original vst matrix
  t() |> 
  prcomp() |> 
  purrr::pluck('x') %>%  # take coordinates
  as.data.frame(stringsAsFactor=F) %>% # transform back from matrix to data.frame 
  dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)



# add first 10 to qc
# qc plots detP x covars ----
# 
# plt <- glass_od.metadata.idats |> 
#   tibble::column_to_rownames('sentrix_id') |> 
#   dplyr::select(contains("detP") | starts_with("qc_"))

cplt <- glass_od.metadata.idats  |>
  dplyr::mutate(`qc_BISULFITE_CONVERSION_I_Beta_I-3_Beta_larger_0.12_0.18`=NULL) |> # praktisch allemaal NA
  dplyr::mutate(`qc_BISULFITE_CONVERSION_I_Beta_I-6_Beta_larger_0.2_0.3`=NULL) |> # praktisch allemaal NA
  tibble::column_to_rownames('sentrix_id') |> 
  dplyr::select(contains("detP") | starts_with("qc_")) |> 
  dplyr::filter(!is.na(`qc_BISULFITE_CONVERSION_I_Beta_I-1_Beta_larger_0.1_0.15`)) |> 
  dplyr::filter(!is.na(`qc_BISULFITE_CONVERSION_I_Beta_I-5_Beta_larger_0.2_0.3`)) |> 
  dplyr::filter(!is.na(`qc_BISULFITE_CONVERSION_II_Beta_II-1_Beta_larger_0.2_0.3`)) |> 
  dplyr::rename_with( ~ gsub("_larger.+$","", .x)) |> 
  dplyr::rename_with( ~ gsub("_smaller.+$","", .x)) |> 
  tibble::rownames_to_column('sentrix_id') |> 
  dplyr::left_join(
    pca |> tibble::rownames_to_column('sentrix_id'),
    by=c('sentrix_id'='sentrix_id')
  ) |> 
  tibble::column_to_rownames('sentrix_id') |> 
  dplyr::filter(!is.na(PC1)) |> 
  dplyr::mutate(PC1 = -1 * PC1) |> 
  dplyr::mutate(percentage.detP.signi = -1 * percentage.detP.signi)


corrplot::corrplot(abs(cor(cplt, method="spearman")),tl.cex=0.4, order="hclust")


ggplot(cplt, aes(x=percentage.detP.signi, y=PC1)) +
  geom_point()



plt <- cplt |> 
  tibble::rownames_to_column('sentrix_id') |> 
  dplyr::left_join(
    glass_od.metadata.idats |>  dplyr::select(sentrix_id, resection_id),
    by=c('sentrix_id'='sentrix_id')
  ) |> 
  dplyr::mutate(percentage.detP.signi = -1 * percentage.detP.signi) |> 
  dplyr::mutate(
    col = case_when(
      abs(percentage.detP.signi) > 5 ~ "detP too high (>5%)",
      resection_id %in% c("0006-R1", "0019-R1") ~ "obfuscating results when included in n=20 test set",
      T ~ "-"
    )
  )


ggplot(plt, aes(x=PC1, y=log(percentage.detP.signi + 1), label=resection_id, color = col)) +
  geom_point() +
  ggrepel::geom_text_repel() +
  youri_gg_theme


