#!/usr/bin/env R

# load m-values ----


source('scripts/load_beta.values_hq_samples.R')



# load EPITOC2 probes ----


load('data/epiTOC2/dataETOC2.Rd')
source('data/epiTOC2/epiTOC2.R')



### run epitoc ----


dat <- data.beta.values.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.beta.values.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 693017)
    assertthat::assert_that(ncol(.) == 510)
    return(.)
  })()


epiTOC2.data.out <- epiTOC2(as.matrix(dat))
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


