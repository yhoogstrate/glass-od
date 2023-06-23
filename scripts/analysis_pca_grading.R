#!/usr/bin/env R


# load data ----


library(ggplot2)


if(!exists('youri_gg_theme')) {
  source('scripts/youri_gg_theme.R')
}



if(!exists('glass_od.metadata.idats')) {
  source('scripts/load_metadata.R')
}



source('scripts/load_mvalues.R')


# load all samples, qc and corr w/ qc stats ----


metadata <- glass_od.metadata.idats |> 
  dplyr::filter(is.na(reason_excluded_patient)) |> 
  dplyr::filter(is.na(reason_excluded_resection)) |> 
  dplyr::filter(is.na(reason_excluded_resection_isolation)) |> 
  dplyr::filter(is.na(reason_excluded_array_sample)) |>  # those are typically the replicates
  assertr::verify(!duplicated(resection_id)) |> 
  assertr::verify(sentrix_id != "204808700074_R04C01")  # 0003-R3-repA



# load m-values ----


data <- glass_od.data.mvalues |> 
    dplyr::select(all_of(metadata$sentrix_id))


#data.full = data
#data = data.full |> 
#  dplyr::slice_head(n=2000)


pca.raw <- data |>
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |>  # this synthax, oh my
  dplyr::filter(mad > 1.0) |>
  dplyr::mutate(mad = NULL)  |>  # remove the sd to obtain original vst matrix
  t() |>
  prcomp() |> 
  purrr::pluck('x')  |>   # take coordinates
  as.data.frame(stringsAsFactor=F) |>   # transform back from matrix to data.frame
  tibble::rownames_to_column('sentrix_id')


## find most variable features ---



plt <- metadata |>
  dplyr::left_join(pca.raw |>  dplyr::select(sentrix_id, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10),
                   by=c('sentrix_id'='sentrix_id'),
                   suffix=c('','')
                   )

plt <- rbind(
  # plt |> 
  #   dplyr::mutate(panel = "primary vs. non primary") |> 
  #   dplyr::mutate(col = ifelse(resection_number == 1, "primary", "recurrence"))
  # ,
  plt |> 
    dplyr::mutate(panel = "grading") |> 
    dplyr::mutate(col = paste0("g",resection_tumor_grade)) |> 
    dplyr::mutate(col = ifelse(col == "g4","g3", col))
)


ggplot(plt, aes(x=PC4, y=PC2, col=col)) + 
  facet_grid(cols = vars(panel), scales = "free", space="free")+ 
  geom_point() +
  theme_bw()






