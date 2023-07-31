#!/usr/bin/env R 


# all hq ----


data.mvalues.hq_samples <- readRDS("cache/mvalues.HQ_samples.Rds") 
data.mvalues.mask.hq_samples <- readRDS("cache/mvalues.HQ_samples.detP_mask.Rds")

stopifnot(colnames(data.mvalues.hq_samples) == colnames(data.mvalues.mask.hq_samples))
stopifnot(rownames(data.mvalues.hq_samples) == rownames(data.mvalues.mask.hq_samples))



data.mvalues.good_probes <- data.mvalues.mask.hq_samples |> 
  is.na() |> 
  rowSums() |> 
  as.data.frame() |> 
  dplyr::rename(n_na = 1) |> 
  dplyr::filter(n_na == 0) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (694299))
    return(.)
  })() |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::pull(probe_id)




