#!/usr/bin/env R


# load data ----


library(ggplot2)
library(minfi)

library(IlluminaHumanMethylationEPICmanifest)



if(!exists('youri_gg_theme')) {
  source('scripts/youri_gg_theme.R')
}



if(!exists('glass_od.metadata.idats')) {
  source('scripts/load_metadata.R')
}




# replicates better than original? ----

#' messy code and not working anymore because samples already
#' excluded before generating the m-values 
#' it was just to see if the replicates were indeed better
#' only for 0003-R3 is was unclear from detP & PC1 stats
#' but CNV profile was clear

## load all samples, qc and corr w/ qc stats ----

mvalue <- readRDS("cache/mvalues.Rds") |> 
  tibble::column_to_rownames("probeID")




sum(array_samples$sentrix_id %in% colnames(mvalue)==F) # should be 0


replicated <- glass_od.metadata.idats |>
  dplyr::filter(duplicated(resection_id)) |>
  dplyr::pull(resection_id)


glass_od.metadata.idats |> 
  dplyr::filter(resection_id %in% replicated) |> 
  dplyr::select(resection_id, resection_isolation_id, 
              reason_excluded_patient,
              reason_excluded_resection,
              reason_excluded_resection_isolation,
              reason_excluded_array_sample
              ) |> 
  dplyr::arrange(resection_id, resection_isolation_id) |> 
  View()


array_samples <- glass_od.metadata.idats |> 
  dplyr::filter(is.na(reason_excluded_patient)) |> 
  dplyr::filter(is.na(reason_excluded_resection))
  #dplyr::filter(is.na(reason_excluded_resection_isolation)) 
  #dplyr::filter(is.na(reason_excluded_array_sample))

sum(duplicated(array_samples$resection_id))

targets <- array_samples |> 
  dplyr::rename(Sample_Name = resection_isolation_id) |>
  dplyr::mutate(Array = gsub("^.+_","",sentrix_id)) |> 
  dplyr::rename(Slide = methylation_array_chip_id) |> 
  dplyr::mutate(Basename = gsub("_Grn.idat$","", channel_green)) |> 
  dplyr::select(Sample_Name, sentrix_id, Array,Slide, Basename)


## pca ----


pca.raw <- mvalue |> 
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |>  # this synthax, oh my
  dplyr::filter(mad > 1.0) |> 
  dplyr::mutate(mad = NULL)  |>  # remove the sd to obtain original vst matrix
  t() |> 
  prcomp()

#saveRDS(pca.raw, "cache/analysis_supervised_qc_pca.Rds")
#pca.raw <- readRDS("cache/analysis_supervised_qc_pca.Rds")


pca2 <- pca.raw |> 
  purrr::pluck('x')  |>   # take coordinates
  as.data.frame(stringsAsFactor=F) |>  # transform back from matrix to data.frame 
  dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)




# add first 10 to qc
## qc plots detP x covars ----
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
    pca2 |> tibble::rownames_to_column('sentrix_id'),
    by=c('sentrix_id'='sentrix_id')
  ) |> 
  tibble::column_to_rownames('sentrix_id') |> 
  dplyr::filter(!is.na(PC1)) |> 
  #dplyr::mutate(PC1 = -1 * PC1) |> 
  dplyr::mutate(percentage.detP.signi = -1 * percentage.detP.signi)


corrplot::corrplot(abs(cor(cplt, method="spearman")),tl.cex=0.4, order="hclust")


ggplot(cplt, aes(x=percentage.detP.signi, y=PC1)) +
  geom_point()



plt <- cplt |> 
  tibble::rownames_to_column('sentrix_id') |> 
  dplyr::left_join(
    glass_od.metadata.idats |>  dplyr::select(sentrix_id, resection_id, resection_isolation_id),
    by=c('sentrix_id'='sentrix_id')
  ) |> 
  dplyr::mutate(percentage.detP.signi = -1 * percentage.detP.signi) |> 
  dplyr::mutate(
    col = dplyr::case_when(
      abs(percentage.detP.signi) > 5 ~ "detP too high (>5%)",
      resection_id %in% c("0006-R1", "0019-R1") ~ "obfuscating results when included in n=20 test set",
      T ~ "-"
    )
  ) |> 
  dplyr::mutate(resection = gsub("^.+\\-","", resection_id))

plt <- plt |> 
       dplyr::mutate(y = PC1 + (runif(nrow(plt)) * 10))

ggplot(plt 
         #dplyr::mutate(col = resection_id %in% replicated)
       , aes(x=PC1,
             y=percentage.detP.signi
             , 
                label=resection_isolation_id, color = col,
                group=resection_id)) +
  geom_point() +
  #geom_line() +
  ggrepel::geom_text_repel( col="black") +
  youri_gg_theme


plt |> 
  dplyr::filter(resection_id == "0003-R3") |> 
  dplyr::pull(sentrix_id)



# QC-check all of the samples ----


mvalue <- readRDS("cache/mvalues.Rds") |> 
  tibble::column_to_rownames("probeID") |> 
  (function(.) {
    assertthat::assert_that(ncol(.) == 199)
    return(.)
  })()


# pca on samples of all quality
pca.raw <- mvalue |> 
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |>  # this synthax, oh my
  dplyr::filter(mad > 1.0) |> 
  dplyr::mutate(mad = NULL)  |>  # remove the sd to obtain original vst matrix
  t() |> 
  prcomp()



plt <- pca.raw |> 
  purrr::pluck('x')  |>   # take coordinates
  as.data.frame(stringsAsFactor=F) |>  # transform back from matrix to data.frame 
  dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) |> 
  tibble::rownames_to_column('sentrix_id') |> 
  dplyr::left_join(
    glass_od.metadata.idats |>  dplyr::select(sentrix_id, resection_id, resection_isolation_id, patient_id),
    by=c('sentrix_id'='sentrix_id')
  ) |> 
  dplyr::left_join(
    
    glass_od.metadata.idats  |>
      dplyr::mutate(`qc_BISULFITE_CONVERSION_I_Beta_I-3_Beta_larger_0.12_0.18`=NULL) |> # praktisch allemaal NA
      dplyr::mutate(`qc_BISULFITE_CONVERSION_I_Beta_I-6_Beta_larger_0.2_0.3`=NULL) |> # praktisch allemaal NA
      tibble::column_to_rownames('sentrix_id') |> 
      dplyr::select(contains("detP") | starts_with("qc_")) |> 
      dplyr::filter(!is.na(`qc_BISULFITE_CONVERSION_I_Beta_I-1_Beta_larger_0.1_0.15`)) |> 
      dplyr::filter(!is.na(`qc_BISULFITE_CONVERSION_I_Beta_I-5_Beta_larger_0.2_0.3`)) |> 
      dplyr::filter(!is.na(`qc_BISULFITE_CONVERSION_II_Beta_II-1_Beta_larger_0.2_0.3`)) |> 
      dplyr::rename_with( ~ gsub("_larger.+$","", .x)) |> 
      dplyr::rename_with( ~ gsub("_smaller.+$","", .x)) |> 
      tibble::rownames_to_column('sentrix_id'),
    
    by=c('sentrix_id'='sentrix_id'), suffix=c('','')
  ) |> 
  dplyr::mutate(
    col = dplyr::case_when(
      abs(percentage.detP.signi) > 5 ~ "detP too high (>5%)",
      resection_id %in% c("0006-R1", "0019-R1") ~ "obfuscating results when included in n=20 test set",
      T ~ "-"
    )
  )




PC1.cutoff <- 425
ggplot(plt, aes(x = PC1, y = percentage.detP.signi, label = resection_isolation_id, color = col, group = patient_id)) +
  geom_vline(xintercept = PC1.cutoff, col = "blue", alpha = 0.5, lty = 2) +
  # geom_line(col="black") +
  geom_point() +
  ggrepel::geom_text_repel(data = plt |> dplyr::filter(PC1 <= PC1.cutoff), col = "black") +
  # ggrepel::geom_text_repel(data = plt |> dplyr::filter(PC1 > PC1.cutoff), col="black") +
  youri_gg_theme



qc.outliers <- plt |>
  dplyr::filter(PC1 > PC1.cutoff) |>
  assertr::verify(c("0019-R2", "0020-R1", "0042-R1", "0019-R1", "0062-R3") %in% resection_isolation_id) |>
  assertr::verify(c("0031-R3", "0060-R1", "0054-R3", "0099-R2") %in% resection_isolation_id == F)


out <- plt |>
  dplyr::select(sentrix_id, PC1) |> 
  dplyr::rename(qc.pca.comp1 = PC1) |> 
  dplyr::mutate(qc.pca.outlier = sentrix_id %in% qc.outliers$sentrix_id)


saveRDS(out, file="cache/unsupervised_qc_qc.outliers.Rds")




