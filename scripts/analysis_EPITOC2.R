#!/usr/bin/env R

# load m-values ----


source('scripts/load_beta.values_hq_samples.R')



# load EPITOC2 probes ----


load('data/epiTOC2/dataETOC2.Rd')
source('data/epiTOC2/epiTOC2.R')



### run epitoc ----



epiTOC2.data.out <- epiTOC2(as.matrix(data.beta.values.hq_samples))
stopifnot(names(epiTOC2.data.out$tnsc) == names(epiTOC2.data.out$tnsc2))
stopifnot(names(epiTOC2.data.out$tnsc) == names(epiTOC2.data.out$pcgtAge))
stopifnot(names(epiTOC2.data.out$tnsc) == names(epiTOC2.data.out$hypoSC))


out <- data.frame(
  tnsc = epiTOC2.data.out$tnsc,
  tnsc2 = epiTOC2.data.out$tnsc2,
  pcgtAge = epiTOC2.data.out$pcgtAge,
  hypoSC = epiTOC2.data.out$hypoSC
) |> 
  dplyr::rename_with( ~ paste0("array_epiTOC2_", .x)) |> 
  tibble::rownames_to_column('array_sentrix_id')



saveRDS(out, file="cache/analysis_EPITOC2.Rds")



rm(epiTOC2.data.out, out)
gc()


