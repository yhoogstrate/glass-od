#!/usr/bin/env R


# load data ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')
#source('scripts/load_palette.R')
#source('scripts/load_themes.R')



if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}


if(!exists('data.intensities.combined.hq_samples')) {
  source('scripts/load_intensities_hq_samples.R')
}


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}

#' too few primary - recurrences to reasonably conduct analysis

# analyses: GLASS-OD WHO Grade ----



metadata <- glass_od.metadata.array_samples |> 
  dplyr::filter(array_methylation_bins_1p19q_purity > 0.1) |> # discard those with near flat CNV profile
  filter_OD_validation_idats(CONST_N_OD_VALIDATION_INCLUDED_SAMPLES) |> 

  dplyr::filter(resection_number > 0) |>  # some unknown resections from oligosarcoma study, could be replicates?
    
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 84)
    return(.)
  })() |> 
  
  dplyr::filter(!is.na(resection_tumor_grade)) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 83)
    return(.)
  })() |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(ifelse(is.paired, paste0("p",patient_id), "a_remainder"))) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2", "Grade3")) |> 
  dplyr::mutate(gr.status = factor(gr.status, levels=c("Grade2","Grade3"))) |> 
  assertr::verify(array_sentrix_id %in% colnames(data.mvalues.hq_samples))



# metadata$array_sentrix_id %in% colnames(data.mvalues.hq_samples)



data <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
    assertthat::assert_that(ncol(.) == 83)
    return(.)
  })()




design <- model.matrix(~factor(patient) + factor(gr.status), data=metadata)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit, trend=T)
stats <- limma::topTable(fit,
                            n=nrow(data),
                            coef="factor(gr.status)Grade3",
                            sort.by = "none",
                            adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design, fit)


sum(stats$P.Value < 0.01)
sum(stats$adj.P.Val < 0.01)



saveRDS(stats, file="cache/analysis_differential__g2_g3__partial_paired_nc__validationset__stats.Rds")
rm(stats)




