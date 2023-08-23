#!/usr/bin/env R


data.probes.alzheimer <- read.delim("data/Alzheimer-DMPs-10.1186_s13148-019-0672-7-table-1.txt") |> 
  dplyr::rename(probe_id = DMPs) |> 
  dplyr::rename(chr = Chromosome) |> 
  dplyr::mutate(chr = paste0("chr",chr)) |> 
  dplyr::rename(pos = Genomic.location) |> 
  assertr::verify(chr %in% c(paste0("chr",1:22))) |> 
  assertr::verify(is.numeric(pos) & pos > 0) |> 
  assertr::verify(is.numeric(Beta..difference) & abs(Beta..difference) >= 0 & abs(Beta..difference) <= 1)  |> 
  assertr::verify(is.numeric(FDR.p.value) & FDR.p.value >= 0 & FDR.p.value <= 0.05) |> 
  tibble::tibble()

