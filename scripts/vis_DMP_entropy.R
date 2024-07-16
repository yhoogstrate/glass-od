#!/usr/bin/env R


# https://www.yacinemahdid.com/shannon-entropy-from-theory-to-python/


if(!exists('data.mvalues.probes')) {
  source('scripts/load_mvalues_hq_samples.R')
}


## windowed - size: 100 probes ----


data <- data.mvalues.probes |> 
  dplyr::filter(!is.na(DMP__g2_g3__pp_nc_PC1__t)) |> 
  dplyr::filter(CHR_hg38 == "chr2") |> 
  dplyr::arrange(Start_hg38) |> 
  dplyr::mutate(data = DMP__g2_g3__pp_nc_PC1__t < 0) |> 
  dplyr::pull(data)


windows <- data.frame(data=data) |> 
  dplyr::mutate(i = (1:dplyr::n()) - 1) |> 
  dplyr::mutate(window = round(i / 100) ) |> 


plot(data, pch=19,cex=0.1)



