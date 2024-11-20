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


if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}



# GLASS-NL: primary - recurrent  ----


metadata.glass_nl.prim_rec <- glass_nl.metadata.array_samples |> 
  dplyr::mutate(patiend_id = paste0("p", patient_id)) |>  # avoid numbers because regression models may go nuts
  filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(178) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id, "remainder")))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary", "recurrence"),levels=c("primary","recurrence"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 178) 
    return(.)
  })()



data.glass_nl.prim_rec <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.glass_nl.prim_rec$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



## naive ----



design.glass_nl.prim_rec_naive <- model.matrix(~factor(patient) + factor(pr.status), data=metadata.glass_nl.prim_rec)
fit.glass_nl.prim_rec_naive <- limma::lmFit(data.glass_nl.prim_rec, design.glass_nl.prim_rec_naive)
fit.glass_nl.prim_rec_naive <- limma::eBayes(fit.glass_nl.prim_rec_naive, trend=T)
stats.glass_nl.prim_rec_naive <- limma::topTable(fit.glass_nl.prim_rec_naive,
                                                 n=nrow(data.glass_nl.prim_rec),
                                                 coef="factor(pr.status)recurrence",
                                                 sort.by = "none",
                                                 adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.glass_nl.prim_rec_naive)


sum(stats.glass_nl.prim_rec_naive$P.Value < 0.01)
sum(stats.glass_nl.prim_rec_naive$adj.P.Val < 0.01)



saveRDS(stats.glass_nl.prim_rec_naive, file="cache/analysis_differential__GLASS-NL__primary_recurrence__partial_paired_nc__stats.Rds")




rm(fit.glass_nl.prim_rec_naive, stats.glass_nl.prim_rec_naive)




## quality (intr. PC1) ----




design.glass_nl.prim_rec_qual <- model.matrix(~PC1 + factor(patient) + factor(pr.status), data=metadata.glass_nl.prim_rec)
fit.glass_nl.prim_rec_qual <- limma::lmFit(data.glass_nl.prim_rec, design.glass_nl.prim_rec_qual)
fit.glass_nl.prim_rec_qual <- limma::eBayes(fit.glass_nl.prim_rec_qual, trend=T)
stats.glass_nl.prim_rec_qual <- limma::topTable(fit.glass_nl.prim_rec_qual,
                                                n=nrow(data.glass_nl.prim_rec),
                                                coef="factor(pr.status)recurrence",
                                                sort.by = "none",
                                                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.glass_nl.prim_rec_qual)



sum(stats.glass_nl.prim_rec_qual$P.Value < 0.01)
sum(stats.glass_nl.prim_rec_qual$adj.P.Val < 0.01)



saveRDS(stats.glass_nl.prim_rec_qual, file="cache/analysis_differential__GLASS-NL__primary_recurrence__PC1__partial_paired_nc__stats.Rds")




rm(fit.glass_nl.prim_rec_qual, stats.glass_nl.prim_rec_qual)



## clean-up ----


rm(metadata.glass_nl.prim_rec)
rm(data.glass_nl.prim_rec)



# GLASS-NL: (2) vs. (g3 & g4) ----


metadata.glass_nl.g2_vs_g3_4.pp.nc <- glass_nl.metadata.array_samples |> 
  dplyr::mutate(patiend_id = paste0("p", patient_id)) |>  # avoid numbers because regression models may go nuts
  filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES) |> 
  filter_first_G2_and_last_G3_G4(135) |> 
  dplyr::filter(!is.na(WHO_Classification2021)) |> # 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder")))) |> 
  dplyr::mutate(resection_tumor_grade = as.numeric(gsub("^.+(.)$", "\\1", WHO_Classification2021))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3,4)) |> 
  dplyr::mutate(gr.status = as.factor(ifelse(resection_tumor_grade == 2, "Grade2", "Grade3,4"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 135) 
    return(.)
  })()


data.glass_nl.g2_vs_g3_4.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.glass_nl.g2_vs_g3_4.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



## naive ----



design.glass_nl.g2_vs_g3_4.naive.pp.nc <- model.matrix(~factor(patient) + factor(gr.status), data=metadata.glass_nl.g2_vs_g3_4.pp.nc)
fit.glass_nl.g2_vs_g3_4.naive.pp.nc <- limma::lmFit(data.glass_nl.g2_vs_g3_4.pp.nc, design.glass_nl.g2_vs_g3_4.naive.pp.nc)
fit.glass_nl.g2_vs_g3_4.naive.pp.nc <- limma::eBayes(fit.glass_nl.g2_vs_g3_4.naive.pp.nc, trend=T)
stats.glass_nl.g2_vs_g3_4.naive.pp.nc <- limma::topTable(fit.glass_nl.g2_vs_g3_4.naive.pp.nc,
                                                            n=nrow(data.glass_nl.g2_vs_g3_4.pp.nc),
                                                            coef="factor(gr.status)Grade3,4",
                                                            sort.by = "none",
                                                            adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

sum(stats.glass_nl.g2_vs_g3_4.naive.pp.nc$P.Value < 0.01)
sum(stats.glass_nl.g2_vs_g3_4.naive.pp.nc$adj.P.Val < 0.01)


rm(design.glass_nl.g2_vs_g3_4.naive.pp.nc)


saveRDS(stats.glass_nl.g2_vs_g3_4.naive.pp.nc, file="cache/analysis_differential__GLASS-NL__g2_g3_4__partial_paired_nc__stats.Rds")


rm(fit.glass_nl.g2_vs_g3_4.naive.pp.nc, stats.glass_nl.g2_vs_g3_4.naive.pp.nc)





## quality (intr. PC1) ----


design.glass_nl.g2_vs_g3_4.quality.pp.nc <- model.matrix(~PC1 + factor(patient) + factor(gr.status), data=metadata.glass_nl.g2_vs_g3_4.pp.nc)
fit.glass_nl.g2_vs_g3_4.quality.pp.nc <- limma::lmFit(data.glass_nl.g2_vs_g3_4.pp.nc, design.glass_nl.g2_vs_g3_4.quality.pp.nc)
fit.glass_nl.g2_vs_g3_4.quality.pp.nc <- limma::eBayes(fit.glass_nl.g2_vs_g3_4.quality.pp.nc, trend=T)
stats.glass_nl.g2_vs_g3_4.quality.pp.nc <- limma::topTable(fit.glass_nl.g2_vs_g3_4.quality.pp.nc,
                                                              n=nrow(data.glass_nl.g2_vs_g3_4.pp.nc),
                                                              coef="factor(gr.status)Grade3,4",
                                                              sort.by = "none",
                                                              adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

sum(stats.glass_nl.g2_vs_g3_4.quality.pp.nc$P.Value < 0.01)
sum(stats.glass_nl.g2_vs_g3_4.quality.pp.nc$adj.P.Val < 0.01)


rm(design.glass_nl.g2_vs_g3_4.quality.pp.nc)


saveRDS(stats.glass_nl.g2_vs_g3_4.quality.pp.nc, file="cache/analysis_differential__GLASS-NL__g2_g3_4__PC1__partial_paired_nc__stats.Rds")


rm(fit.glass_nl.g2_vs_g3_4.quality.pp.nc, stats.glass_nl.g2_vs_g3_4.quality.pp.nc)


## clean-up ----


rm(metadata.glass_nl.g2_vs_g3_4.pp.nc)
rm(data.glass_nl.g2_vs_g3_4.pp.nc)



# GLASS-NL: (g2 & g3) vs. (g4) ----


metadata.glass_nl.g2_3_vs_g4.pp.nc <- glass_nl.metadata.array_samples |> 
  dplyr::mutate(patiend_id = paste0("p", patient_id)) |>  # avoid numbers because regression models may go nuts
  filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES) |> 
  filter_first_G2_G3_and_last_G4(133) |> 
  dplyr::filter(!is.na(WHO_Classification2021)) |> # 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder")))) |> 
  dplyr::mutate(resection_tumor_grade = as.numeric(gsub("^.+(.)$", "\\1", WHO_Classification2021))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3,4)) |> 
  dplyr::mutate(gr.status = factor(ifelse(resection_tumor_grade == 4, "Grade4", "Grade2,3"), levels = c("Grade2,3", "Grade4"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 133) 
    return(.)
  })()


data.glass_nl.g2_3_vs_g4.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.glass_nl.g2_3_vs_g4.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


## naive ----


design.glass_nl.g2_3_vs_g4.naive.pp.nc <- model.matrix(~factor(patient) + factor(gr.status), data=metadata.glass_nl.g2_3_vs_g4.pp.nc)
fit.glass_nl.g2_3_vs_g4.naive.pp.nc <- limma::lmFit(data.glass_nl.g2_3_vs_g4.pp.nc, design.glass_nl.g2_3_vs_g4.naive.pp.nc)
fit.glass_nl.g2_3_vs_g4.naive.pp.nc <- limma::eBayes(fit.glass_nl.g2_3_vs_g4.naive.pp.nc, trend=T)
stats.glass_nl.g2_3_vs_g4.naive.pp.nc <- limma::topTable(fit.glass_nl.g2_3_vs_g4.naive.pp.nc,
                                                         n=nrow(data.glass_nl.g2_3_vs_g4.pp.nc),
                                                         coef="factor(gr.status)Grade4",
                                                         sort.by = "none",
                                                         adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 


sum(stats.glass_nl.g2_3_vs_g4.naive.pp.nc$P.Value < 0.01)
sum(stats.glass_nl.g2_3_vs_g4.naive.pp.nc$adj.P.Val < 0.01)


rm(design.glass_nl.g2_3_vs_g4.naive.pp.nc)


saveRDS(stats.glass_nl.g2_3_vs_g4.naive.pp.nc, file="cache/analysis_differential__GLASS-NL__g2_3_g4__partial_paired_nc__stats.Rds")


rm(fit.glass_nl.g2_3_vs_g4.naive.pp.nc, stats.glass_nl.g2_3_vs_g4.naive.pp.nc)




## quality (intr. PC1) ----


design.glass_nl.g2_3_vs_g4.quality.pp.nc <- model.matrix(~PC1 + factor(patient) + factor(gr.status), data=metadata.glass_nl.g2_3_vs_g4.pp.nc)
fit.glass_nl.g2_3_vs_g4.quality.pp.nc <- limma::lmFit(data.glass_nl.g2_3_vs_g4.pp.nc, design.glass_nl.g2_3_vs_g4.quality.pp.nc)
fit.glass_nl.g2_3_vs_g4.quality.pp.nc <- limma::eBayes(fit.glass_nl.g2_3_vs_g4.quality.pp.nc, trend=T)
stats.glass_nl.g2_3_vs_g4.quality.pp.nc <- limma::topTable(fit.glass_nl.g2_3_vs_g4.quality.pp.nc,
                                                           n=nrow(data.glass_nl.g2_3_vs_g4.pp.nc),
                                                           coef="factor(gr.status)Grade4",
                                                           sort.by = "none",
                                                           adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

sum(stats.glass_nl.g2_3_vs_g4.quality.pp.nc$P.Value < 0.01)
sum(stats.glass_nl.g2_3_vs_g4.quality.pp.nc$adj.P.Val < 0.01)


rm(design.glass_nl.g2_3_vs_g4.quality.pp.nc)


saveRDS(stats.glass_nl.g2_3_vs_g4.quality.pp.nc, file="cache/analysis_differential__GLASS-NL__g2_3_g4__PC1__partial_paired_nc__stats.Rds")


rm(fit.glass_nl.g2_3_vs_g4.quality.pp.nc, stats.glass_nl.g2_3_vs_g4.quality.pp.nc)



## clean-up ----

rm(metadata.glass_nl.g2_3_vs_g4.pp.nc)
rm(data.glass_nl.g2_3_vs_g4.pp.nc)





