#!/usr/bin/env R


# BMC Clinical genetics 2019 ----

data.probes.alzheimer <- read.delim("data/Alzheimer-DMPs-10.1186_s13148-019-0672-7-table-1.txt") |> 
  dplyr::rename(probe_id = DMPs) |> 
  dplyr::rename(chr = Chromosome) |> 
  dplyr::mutate(chr = paste0("chr",chr)) |> 
  dplyr::rename(pos = Genomic.location) |> 
  assertr::verify(chr %in% c(paste0("chr",1:22))) |> 
  assertr::verify(is.numeric(pos) & pos > 0) |> 
  assertr::verify(is.numeric(Beta..difference) & abs(Beta..difference) >= 0 & abs(Beta..difference) <= 1)  |> 
  assertr::verify(is.numeric(FDR.p.value) & FDR.p.value >= 0 & FDR.p.value <= 0.05) |> 
  dplyr::mutate(label = paste0(GeneID1,";",GeneID2)) |> 
  tibble::tibble()


# Nat Com 2022 ----

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203332
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197305


