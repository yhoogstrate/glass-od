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



# analyses: GLASS-OD primary - last recurrence ----
## data: example ----


metadata.example <- data.frame(pid = (1:100) %% 50) |> 
  dplyr::arrange(pid) |> 
  dplyr::mutate(pid = factor(paste0("p", pid))) |> 
  dplyr::mutate(sid = factor(paste0("s",1:100))) |> 
  dplyr::mutate(condition = factor(c(rep("c1",75),rep("c2",25))))


replag <- function(c) {
  d <- c()
  
  x <- 1:length(c)
  y2 <- (x)*2
  y1 <- y2 - 1
  
  d[y1] <- c[x]
  d[y2] <- c[x]
  
  return(d)
}

reppos <- function(c, pl) {
  d = c()
  
  for(k in 1:length(c)) {
    p = levels(pl)[k]
    
    d[which(pl == p)] <- c[k]
  }
  
  return(d)
}


#reppos(runif(3), factor(c("p1","p4","p6","p6","p1","p4")))





metadata.example <- data.frame(
  condition = c("c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                
                "c2","c2",
                "c2","c2",
                "c2","c2",
                "c2","c2",
                "c2","c2")
) |> 
  dplyr::mutate(pat = factor(replag(paste0("p",1:35)))) |> 
  dplyr::arrange(condition, pat)



data.example <- data.frame(
  `nodiff_A_001` = c(runif(45), runif(25)) + reppos(runif(35), metadata.example$pat),
  `nodiff_A_002` = c(runif(45), runif(25)) + reppos(runif(35), metadata.example$pat),
  `nodiff_A_003` = c(runif(45), runif(25)) + reppos(runif(35), metadata.example$pat),
  `nodiff_A_004` = c(runif(45), runif(25)) + reppos(runif(35), metadata.example$pat),
  `nodiff_A_005` = c(runif(45), runif(25)) + reppos(runif(35), metadata.example$pat),
  
  `nodiff_B_001` = c(runif(45), runif(25)) + reppos(rep(0, 35), metadata.example$pat),
  `nodiff_B_002` = c(runif(45), runif(25)) + reppos(rep(0, 35), metadata.example$pat),
  `nodiff_B_003` = c(runif(45), runif(25)) + reppos(rep(0, 35), metadata.example$pat),
  `nodiff_B_004` = c(runif(45), runif(25)) + reppos(rep(0, 35), metadata.example$pat),
  `nodiff_B_005` = c(runif(45), runif(25)) + reppos(rep(0, 35), metadata.example$pat),
  
  `nodiff_high_batch_A_001` = c(runif(45), runif(25)) + reppos(runif(35,10,20), metadata.example$pat),
  `nodiff_high_batch_A_002` = c(runif(45), runif(25)) + reppos(runif(35,10,20), metadata.example$pat),
  `nodiff_high_batch_A_003` = c(runif(45), runif(25)) + reppos(runif(35,10,20), metadata.example$pat),
  `nodiff_high_batch_A_004` = c(runif(45), runif(25)) + reppos(runif(35,10,20), metadata.example$pat),
  `nodiff_high_batch_A_005` = c(runif(45), runif(25)) + reppos(runif(35,10,20), metadata.example$pat),
  
  `C1_high_batch_A_001` = c(runif(45,1,2), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_002` = c(runif(45,1,2), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_003` = c(runif(45,1.5,3), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_004` = c(runif(45,1.5,3), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_005` = c(runif(45,10,12), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_006` = c(runif(45,7,8), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_007` = c(runif(45,0,5), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_008` = c(runif(45,0,5), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_009` = c(runif(45,0,5), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_010` = c(runif(45,0,5), runif(25)) + reppos(runif(35), metadata.example$pat),
  
  
  `C1_low_batch_A_001` = c(runif(45,0,0.1), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_low_batch_A_002` = c(runif(45,0,0.1), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_low_batch_A_003` = c(runif(45,0,0.1), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_low_batch_A_004` = c(runif(45,0,0.1), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_low_batch_A_005` = c(runif(45,0,0.1), runif(25)) + reppos(runif(35), metadata.example$pat)
  
) |> 
  t() |> 
  as.data.frame()


### test: indeed p.values and not f-test ----


design.example <- model.matrix(~pat + condition, data=metadata.example)
fit.example <- limma::lmFit(as.matrix(data.example), design.example)
fit.example <- limma::eBayes(fit.example,trend=T)
stats.example <- limma::topTable(fit.example, n=nrow(fit.example)) # adjust="BH", sort="p")

pval = as.data.frame(fit.example$p.value) |> 
  tibble::rownames_to_column('testid') 

ggplot(pval, aes(x=testid, y=conditionc2)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))


coef <- as.data.frame(fit.example$coefficients) |> 
  tibble::rownames_to_column('testid') 

ggplot(coef, aes(x=testid, y=conditionc2)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))


plot(coef$conditionc2, -log(pval$conditionc2))



## data: partially paired [w/o FFPE/frozen batch correct] ----


metadata.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(180) |> 
  
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"remainder")))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 
  dplyr::select(array_sentrix_id, patient, pr.status) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 180) 
    return(.)
  })()


data.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



design.pp.nc <- model.matrix(~factor(patient) + factor(pr.status), data=metadata.pp.nc)
fit.pp.nc <- limma::lmFit(data.pp.nc, design.pp.nc)
fit.pp.nc <- limma::eBayes(fit.pp.nc, trend=T)
stats.pp.nc <- limma::topTable(fit.pp.nc,
                               n=nrow(data.pp.nc),
                               coef="factor(pr.status)recurrence",
                               sort.by = "none",
                               adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.pp.nc)


sum(stats.pp.nc$P.Value < 0.01)
sum(stats.pp.nc$adj.P.Val < 0.01)



#saveRDS(fit.pp.nc, file="cache/analysis_differential__primary_recurrence__partial_paired_nc__fit.Rds")
saveRDS(stats.pp.nc, file="cache/analysis_differential__primary_recurrence__partial_paired_nc__stats.Rds")


rm(fit.pp.nc, stats.pp.nc)


## data: partially paired + quality [w/o FFPE/frozen batch correct] ----


metadata.pp.prim_rec.quality.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(180) |> 
  
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"remainder")))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 

  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 180) 
    return(.)
  })()


data.pp.prim_rec.quality.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.pp.prim_rec.quality.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_PROBES_UNMASKED_AND_DETP)) 
    return(.)
  })()



design.pp.prim_rec.quality.nc <- model.matrix(~array_PC1 + factor(patient) + factor(pr.status), data=metadata.pp.prim_rec.quality.nc)
fit.pp.prim_rec.quality.nc <- limma::lmFit(data.pp.prim_rec.quality.nc, design.pp.prim_rec.quality.nc)
fit.pp.prim_rec.quality.nc <- limma::eBayes(fit.pp.prim_rec.quality.nc, trend=T)
stats.pp.prim_rec.quality.nc <- limma::topTable(fit.pp.prim_rec.quality.nc,
                                                n=nrow(data.pp.prim_rec.quality.nc),
                                                coef="factor(pr.status)recurrence",
                                                sort.by = "none",
                                                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.pp.prim_rec.quality.nc)


sum(stats.pp.prim_rec.quality.nc$P.Value < 0.01)
sum(stats.pp.prim_rec.quality.nc$adj.P.Val < 0.01)



saveRDS(stats.pp.prim_rec.quality.nc, file="cache/analysis_differential__primary_recurrence__PC1__partial_paired_nc__stats.Rds")




rm(fit.pp.prim_rec.quality.nc, stats.pp.prim_rec.quality.nc)



## data: partially paired INTENSITY [w/o FFPE/frozen batch correct] ----


metadata.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(180) |> 
  
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired, patient_id, "remainder")))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1, "primary", "recurrence"),levels=c("primary","recurrence"))) |> 
  dplyr::select(array_sentrix_id, patient, pr.status) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 180) 
    return(.)
  })()




data.pp.nc <- data.intensities.combined.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.intensities.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()





design.pp.nc <- model.matrix(~factor(patient) + factor(pr.status), data=metadata.pp.nc)
fit.pp.nc <- limma::lmFit(data.pp.nc, design.pp.nc)
fit.pp.nc <- limma::eBayes(fit.pp.nc, trend=T)
stats.pp.nc <- limma::topTable(fit.pp.nc,
                               n=nrow(data.pp.nc),
                               coef="factor(pr.status)recurrence",
                               sort.by = "none",
                               adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.pp.nc)


sum(stats.pp.nc$P.Value < 0.01)
sum(stats.pp.nc$adj.P.Val < 0.01)



saveRDS(stats.pp.nc, file="cache/analysis_differential_intensities__primary_recurrence__partial_paired_nc__stats.Rds")


rm(fit.pp.nc, stats.pp.nc, data.pp.nc)





## data: partially paired INTENSITY [w/o FFPE/frozen batch correct] ----


metadata.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(180) |> 
  
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"remainder")))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 
  dplyr::select(array_sentrix_id, patient, pr.status) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 180) 
    return(.)
  })()




data.pp.nc <- data.intensities.combined.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.intensities.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



data.pp.nc <- exp(data.pp.nc) - 1


design.pp.nc <- model.matrix(~factor(patient) + factor(pr.status), data=metadata.pp.nc)
fit.pp.nc <- limma::lmFit(data.pp.nc, design.pp.nc)
fit.pp.nc <- limma::eBayes(fit.pp.nc, trend=T)
stats.pp.nc <- limma::topTable(fit.pp.nc,
                               n=nrow(data.pp.nc),
                               coef="factor(pr.status)recurrence",
                               sort.by = "none",
                               adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.pp.nc)


sum(stats.pp.nc$P.Value < 0.01)
sum(stats.pp.nc$adj.P.Val < 0.01)



#saveRDS(fit.pp.nc, file="cache/analysis_differential_intensities__primary_recurrence__partial_paired_nc__fit.Rds")
saveRDS(stats.pp.nc, file="cache/analysis_differential_intensities__primary_recurrence__partial_paired_nc__non_log__stats.Rds")


rm(fit.pp.nc, stats.pp.nc, data.pp.nc)




## x data: full paired only ----



metadata.fp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(179) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::filter(dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",patient_id))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 
  dplyr::select(array_sentrix_id, patient, pr.status, batch) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (138)) # number of perfect pairs is smaller 
    return(.)
  })()


tmp <- rbind(
  metadata.fp |> 
    dplyr::filter(pr.status == "primary") |> 
    dplyr::mutate(patient = as.character(patient)) |> 
    dplyr::mutate(patient = sample(patient, replace=F)) |> 
    assertr::verify(!duplicated(patient))
  ,
  metadata.fp |> 
    dplyr::filter(pr.status == "recurrence") |> 
    dplyr::mutate(patient = as.character(patient)) |> 
    dplyr::mutate(patient = sample(patient, replace=F)) |> 
    assertr::verify(!duplicated(patient))
) |> 
  #dplyr::mutate(patient_shuffled = as.factor(patient)) |> 
  dplyr::select(array_sentrix_id)

metadata.fp <- metadata.fp |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))
rm(tmp)



data.fp <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.fp$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_PROBES_UNMASKED_AND_DETP)) 
    return(.)
  })()


design.fp <- model.matrix(~factor(patient) + factor(pr.status), data=metadata.fp)
fit.fp <- limma::lmFit(data.fp, design.fp)
fit.fp <- limma::eBayes(fit.fp, trend=T)
stats.fp <- limma::topTable(fit.fp,
                n=nrow(data.fp),
                coef="factor(pr.status)recurrence",
                sort.by = "none",
                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 


rm(design.fp)


sum(stats.fp$P.Value < 0.01)
sum(stats.fp$adj.P.Val < 0.01)



saveRDS(fit.fp, file="cache/analysis_differential__primary_recurrence__full_paired__fit.Rds")
saveRDS(stats.fp, file="cache/analysis_differential__primary_recurrence__full_paired__stats.Rds")



rm(fit.fp, stats.fp)








## x data: partially paired [w/ FFPE/frozen batch correct] ----


metadata.pp.cor <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(211) |> 
  filter_primaries_and_last_recurrences(179) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  #dplyr::mutate(batch = factor(ifelse(isolation_person_name == "USA / Duke", "Batch US [FF]", "Batch EU [FFPE]"))) |> 
  dplyr::mutate(patient_correct = 
                  dplyr::case_when(
                    is.paired ~ paste0("p",patient_id),
                    !is.paired & isolation_person_name == "USA / Duke" ~ "remainder - batch US [FF]",
                    !is.paired & isolation_person_name != "USA / Duke" ~ "remainder - batch EU [FFPE]",
                    T ~ "error"
                  )) |>
  assertr::verify(patient_correct != "error") |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 
  
  dplyr::select(array_sentrix_id, patient_correct, pr.status) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (179)) 
    return(.)
  })()



data.pp.cor <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.pp.cor$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_PROBES_UNMASKED_AND_DETP)) 
    return(.)
  })()



design.pp.cor <- model.matrix(~ factor(patient_correct) + factor(pr.status), data=metadata.pp.cor)
fit.pp.cor <- limma::lmFit(data.pp.cor, design.pp.cor)
fit.pp.cor <- limma::eBayes(fit.pp.cor, trend=T)
stats.pp.cor <- limma::topTable(fit.pp.cor,
                                n=nrow(data.pp.cor),
                                coef="factor(pr.status)recurrence",
                                sort.by = "none",
                                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 


rm(design.pp.cor)


sum(stats.pp.cor$P.Value < 0.01)
sum(stats.pp.cor$adj.P.Val < 0.01)



saveRDS(fit.pp.cor, file="cache/analysis_differential__primary_recurrence__partial_paired_cor__fit.Rds")
saveRDS(stats.pp.cor, file="cache/analysis_differential__primary_recurrence__partial_paired_cor__stats.Rds")


rm(fit.pp.cor, stats.pp.cor)



## x data: unpaired [w/o FFPE/frozen batch correct] ----


metadata.up.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(211) |> 
  filter_primaries_and_last_recurrences(179) |> 
  
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 
  
  dplyr::select(array_sentrix_id, pr.status) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (179)) 
    return(.)
  })()



data.up.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.up.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_PROBES_UNMASKED_AND_DETP)) 
    return(.)
  })()



design.up.nc <- model.matrix(~ factor(pr.status), data=metadata.up.nc)
fit.up.nc <- limma::lmFit(data.up.nc, design.up.nc)
fit.up.nc <- limma::eBayes(fit.up.nc, trend=T)
stats.up.nc <- limma::topTable(fit.up.nc,
                                n=nrow(data.up.nc),
                                coef="factor(pr.status)recurrence",
                                sort.by = "none",
                                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 


rm(design.up.nc)


sum(stats.up.nc$P.Value < 0.01)
sum(stats.up.nc$adj.P.Val < 0.01)



saveRDS(fit.up.nc, file="cache/analysis_differential__primary_recurrence__unpaired_nc__fit.Rds")
saveRDS(stats.up.nc, file="cache/analysis_differential__primary_recurrence__unpaired_nc__stats.Rds")


rm(fit.up.cor, stats.up.cor)




## x data: unpaired [w/ FFPE/frozen batch correct] ----

metadata.up.cor <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(211) |> 
  filter_primaries_and_last_recurrences(179) |> 
  
  dplyr::mutate(batch = factor(ifelse(isolation_person_name == "USA / Duke", "Batch US [FF]", "Batch EU [FFPE]"))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 
  
  dplyr::select(array_sentrix_id, batch, pr.status) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (179)) 
    return(.)
  })()



data.up.cor <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.up.cor$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_PROBES_UNMASKED_AND_DETP)) 
    return(.)
  })()



design.up.cor <- model.matrix(~ factor(batch) + factor(pr.status), data=metadata.up.cor)
fit.up.cor <- limma::lmFit(data.up.cor, design.up.cor)
fit.up.cor <- limma::eBayes(fit.up.cor, trend=T)
stats.up.cor <- limma::topTable(fit.up.cor,
                                n=nrow(data.up.cor),
                                coef="factor(pr.status)recurrence",
                                sort.by = "none",
                                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 


rm(design.up.cor)


sum(stats.up.cor$P.Value < 0.01)
sum(stats.up.cor$adj.P.Val < 0.01)



saveRDS(fit.up.cor, file="cache/analysis_differential__primary_recurrence__unpaired_cor__fit.Rds")
saveRDS(stats.up.cor, file="cache/analysis_differential__primary_recurrence__unpaired_cor__stats.Rds")


rm(fit.up.cor, stats.up.cor)







## x plots ----



plt <- readRDS(file="cache/analysis_differential__primary_recurrence__partial_paired_nc__stats.Rds") |> 
  dplyr::mutate(significant = adj.P.Val < 0.01 & abs(logFC) > 0.5)


ggplot(plt, aes(x=logFC , y=-log10(adj.P.Val), col=significant)) +
  geom_point(pch=19,cex=0.01) +
  labs(title = "DMPs primary - recurrence") +
  theme_cellpress

ggsave("/tmp/volcano_pp.png", width=8.5/2,height=8.5/2)



plt <- readRDS(file="cache/analysis_differential__g2_g3__partial_paired_nc__stats.Rds") |> 
  dplyr::mutate(significant = adj.P.Val < 0.01 & abs(logFC) > 0.5)


ggplot(plt, aes(x=logFC , y=-log10(adj.P.Val), col=significant)) +
  geom_point(pch=19,cex=0.01) +
  labs(title = "DMPs WHO Grade 2 - WHO Grade 3") +
  theme_cellpress

ggsave("/tmp/volcano_g3g2.png", width=8.5/2,height=8.5/2)


#saveRDS(stats.pp.nc, file="cache/analysis_differential__primary_recurrence__partial_paired_nc__stats.Rds")
#saveRDS(stats.pp.cor, file="cache/analysis_differential__primary_recurrence__partial_paired_cor__stats.Rds")
#saveRDS(stats.up.nc, file="cache/analysis_differential__primary_recurrence__unpaired_nc__stats.Rds")
#saveRDS(stats.up.cor, file="cache/analysis_differential__primary_recurrence__unpaired_cor__stats.Rds")





# analyses: GLASS-OD g2 - g3 ----
#' @todo ASK what to do if first resection is G3 while last is G2 !!!


## data: partially paired [w/o FFPE/frozen batch correct] ----


metadata.g2g3.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_first_G2_and_last_G3(156) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder")))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2", "Grade3")) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 156) 
    return(.)
  })()



data.g2g3.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.g2g3.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_PROBES_UNMASKED_AND_DETP)) 
    return(.)
  })()





design.g2g3.pp.nc <- model.matrix(~factor(patient) + factor(gr.status), data=metadata.g2g3.pp.nc)
fit.g2g3.pp.nc <- limma::lmFit(data.g2g3.pp.nc, design.g2g3.pp.nc)
fit.g2g3.pp.nc <- limma::eBayes(fit.g2g3.pp.nc, trend=T)
stats.g2g3.pp.nc <- limma::topTable(fit.g2g3.pp.nc,
                                    n=nrow(data.g2g3.pp.nc),
                                    coef="factor(gr.status)Grade3",
                                    sort.by = "none",
                                    adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.g2g3.pp.nc)


sum(stats.g2g3.pp.nc$P.Value < 0.01)
sum(stats.g2g3.pp.nc$adj.P.Val < 0.01)



#saveRDS(fit.g2g3.pp.nc, file="cache/analysis_differential__g2_g3__partial_paired_nc__fit.Rds")
saveRDS(stats.g2g3.pp.nc, file="cache/analysis_differential__g2_g3__partial_paired_nc__stats.Rds")


rm(fit.g2g3.pp.nc, stats.g2g3.pp.nc)



## data: partially paired + quality [w/o FFPE/frozen batch correct] ----


metadata.g2g3.quality.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_first_G2_and_last_G3(156) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder")))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2", "Grade3")) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 156) 
    return(.)
  })()



data.g2g3.quality.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.g2g3.quality.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_PROBES_UNMASKED_AND_DETP)) 
    return(.)
  })()




design.g2g3.quality.pp.nc <- model.matrix(~array_PC1 + factor(patient) + factor(gr.status), data=metadata.g2g3.quality.pp.nc)
fit.g2g3.quality.pp.nc <- limma::lmFit(data.g2g3.quality.pp.nc, design.g2g3.quality.pp.nc)
fit.g2g3.quality.pp.nc <- limma::eBayes(fit.g2g3.quality.pp.nc, trend=T)
stats.g2g3.quality.pp.nc <- limma::topTable(fit.g2g3.quality.pp.nc,
                                            n=nrow(data.g2g3.quality.pp.nc),
                                            coef="factor(gr.status)Grade3",
                                            sort.by = "none",
                                            adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

sum(stats.g2g3.quality.pp.nc$P.Value < 0.01)
sum(stats.g2g3.quality.pp.nc$adj.P.Val < 0.01)


rm(design.g2g3.quality.pp.nc)



saveRDS(stats.g2g3.quality.pp.nc, file="cache/analysis_differential__g2_g3__PC1__partial_paired_nc__stats.Rds")


rm(fit.g2g3.quality.pp.nc, stats.g2g3.quality.pp.nc)





## data: partially paired INTENSITIES [w/o FFPE/frozen batch correct] ----


metadata.g2g3.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_first_G2_and_last_G3(156) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder")))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2", "Grade3")) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 156) 
    return(.)
  })()



data.g2g3.pp.nc <- data.intensities.combined.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.intensities.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.g2g3.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



design.g2g3.pp.nc <- model.matrix(~factor(patient) + factor(gr.status), data=metadata.g2g3.pp.nc)
fit.g2g3.pp.nc <- limma::lmFit(data.g2g3.pp.nc, design.g2g3.pp.nc)
fit.g2g3.pp.nc <- limma::eBayes(fit.g2g3.pp.nc, trend=T)
stats.g2g3.pp.nc <- limma::topTable(fit.g2g3.pp.nc,
                                    n=nrow(data.g2g3.pp.nc),
                                    coef="factor(gr.status)Grade3",
                                    sort.by = "none",
                                    adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.g2g3.pp.nc)


sum(stats.g2g3.pp.nc$P.Value < 0.01)
sum(stats.g2g3.pp.nc$adj.P.Val < 0.01)



#saveRDS(fit.g2g3.pp.nc, file="cache/analysis_differential_intensities__g2_g3__partial_paired_nc__fit.Rds")
saveRDS(stats.g2g3.pp.nc, file="cache/analysis_differential_intensities__g2_g3__partial_paired_nc__stats.Rds")


rm(fit.g2g3.pp.nc, stats.g2g3.pp.nc, data.g2g3.pp.nc)


## data: partially paired INTENSITIES nonlog [w/o FFPE/frozen batch correct] ----


metadata.g2g3.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_first_G2_and_last_G3(156) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder")))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2", "Grade3")) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 156) 
    return(.)
  })()



data.g2g3.pp.nc <- data.intensities.combined.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.intensities.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.g2g3.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



data.g2g3.pp.nc <- exp(data.g2g3.pp.nc) - 1



design.g2g3.pp.nc <- model.matrix(~factor(patient) + factor(gr.status), data=metadata.g2g3.pp.nc)
fit.g2g3.pp.nc <- limma::lmFit(data.g2g3.pp.nc, design.g2g3.pp.nc)
fit.g2g3.pp.nc <- limma::eBayes(fit.g2g3.pp.nc, trend=T)
stats.g2g3.pp.nc <- limma::topTable(fit.g2g3.pp.nc,
                                    n=nrow(data.g2g3.pp.nc),
                                    coef="factor(gr.status)Grade3",
                                    sort.by = "none",
                                    adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.g2g3.pp.nc)


sum(stats.g2g3.pp.nc$P.Value < 0.01)
sum(stats.g2g3.pp.nc$adj.P.Val < 0.01)



#saveRDS(fit.g2g3.pp.nc, file="cache/analysis_differential_intensities__g2_g3__partial_paired_nc__fit.Rds")
saveRDS(stats.g2g3.pp.nc, file="cache/analysis_differential_intensities__g2_g3__partial_paired_nc__non_log__stats.Rds")


rm(fit.g2g3.pp.nc, stats.g2g3.pp.nc, data.g2g3.pp.nc)


## data: only non-FFPE samples ----


## data: partially paired [w/o FFPE/frozen batch correct] ----


metadata.g2g3.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_first_G2_and_last_G3(156) |> 
  dplyr::filter(isolation_material == "tissue") |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder")))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2", "Grade3")) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 37) 
    return(.)
  })() |> 
  dplyr::mutate(array_PC1 = scale(array_PC1, scale=F, center=T)) |> 
  dplyr::mutate(array_PC2 = scale(array_PC2, scale=F, center=T))



data.g2g3.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.g2g3.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_PROBES_UNMASKED_AND_DETP)) 
    return(.)
  })()





design.g2g3.pp.nc <- model.matrix(~factor(patient) + 
                                    array_PC1 + 
                                    #array_PC2
                                  array_A_IDH_HG__A_IDH_LG_lr__lasso_fit
                                    #+  
                                    #factor(gr.status)
                                  , data=metadata.g2g3.pp.nc)
fit.g2g3.pp.nc <- limma::lmFit(data.g2g3.pp.nc, design.g2g3.pp.nc)
fit.g2g3.pp.nc <- limma::eBayes(fit.g2g3.pp.nc, trend=T)
stats.g2g3.pp.nc <- limma::topTable(fit.g2g3.pp.nc,
                                    n=nrow(data.g2g3.pp.nc),
                                    #coef="factor(gr.status)Grade3",
                                    #coef="array_PC1",
                                    #coef="array_PC2",
                                    coef="array_A_IDH_HG__A_IDH_LG_lr__lasso_fit",
                                    sort.by = "none",
                                    adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.g2g3.pp.nc)


sum(stats.g2g3.pp.nc$P.Value < 0.01)
sum(stats.g2g3.pp.nc$adj.P.Val < 0.01)

ggplot(stats.g2g3.pp.nc, aes(x= logFC, y= -log(`P.Value`))) +
  
  geom_hline(yintercept=0,          col="black", lwd=theme_nature_lwd, lty=1) +
#  geom_hline(yintercept=-log(0.01), col="black", lwd=theme_nature_lwd, lty=2) +
  #geom_hline(yintercept=-log(0.05), col="gray", lwd=theme_nature_lwd) +
  
  geom_vline(xintercept=0,        col="black", lwd=theme_nature_lwd, lty=2) +

  geom_point(pch=19, cex=0.5) +
  theme_nature







#saveRDS(stats.g2g3.pp.nc, file="cache/analysis_differential__g2_g3__partial_paired_nc__stats.Rds")


rm(fit.g2g3.pp.nc, stats.g2g3.pp.nc)





## x data: full paired only ----


metadata.g2g3.fp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(211) |> 
  filter_first_G2_and_last_G3(140) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::filter(dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",patient_id))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2","Grade3")) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (72)) 
    return(.)
  })()



data.g2g3.fp <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.g2g3.fp$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_PROBES_UNMASKED_AND_DETP)) 
    return(.)
  })()



design.g2g3.fp <- model.matrix(~factor(patient) + factor(gr.status), data=metadata.g2g3.fp)
fit.g2g3.fp <- limma::lmFit(data.g2g3.fp, design.g2g3.fp)
fit.g2g3.fp <- limma::eBayes(fit.g2g3.fp, trend=T)
stats.g2g3.fp <- limma::topTable(fit.g2g3.fp,
                                 n=nrow(data.g2g3.fp),
                                 coef="factor(gr.status)Grade3",
                                 sort.by = "none",
                                 adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.g2g3.fp)


sum(stats.g2g3.fp$P.Value < 0.01)
sum(stats.g2g3.fp$adj.P.Val < 0.01)



saveRDS(fit.g2g3.fp, file="cache/analysis_differential__g2_g3__full_paired__fit.Rds")
saveRDS(stats.g2g3.fp, file="cache/analysis_differential__g2_g3__full_paired__stats.Rds")




## x data: partially paired [w/ FFPE/frozen batch correct] ----


metadata.g2g3.pp.cor <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(211) |> 
  filter_first_G2_and_last_G3(140) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient_correct = 
                  dplyr::case_when(
                    is.paired ~ paste0("p",patient_id),
                    !is.paired & isolation_person_name == "USA / Duke" ~ "remainder - batch US [FF]",
                    !is.paired & isolation_person_name != "USA / Duke" ~ "remainder - batch EU [FFPE]",
                    T ~ "error"
                  )) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2","Grade3")) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (140)) 
    return(.)
  })()



data.g2g3.pp.cor <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.g2g3.pp.cor$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_PROBES_UNMASKED_AND_DETP)) 
    return(.)
  })()





design.g2g3.pp.cor <- model.matrix(~factor(patient_correct) + factor(gr.status), data=metadata.g2g3.pp.cor)
fit.g2g3.pp.cor <- limma::lmFit(data.g2g3.pp.cor, design.g2g3.pp.cor)
fit.g2g3.pp.cor <- limma::eBayes(fit.g2g3.pp.cor, trend=T)
stats.g2g3.pp.cor <- limma::topTable(fit.g2g3.pp.cor,
                                    n=nrow(data.g2g3.pp.cor),
                                    coef="factor(gr.status)Grade3",
                                    sort.by = "none",
                                    adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.g2g3.pp.cor)


sum(stats.g2g3.pp.cor$P.Value < 0.01)
sum(stats.g2g3.pp.cor$adj.P.Val < 0.01)



saveRDS(fit.g2g3.pp.cor, file="cache/analysis_differential__g2_g3__partial_paired_cor__fit.Rds")
saveRDS(stats.g2g3.pp.cor, file="cache/analysis_differential__g2_g3__partial_paired_cor__stats.Rds")



## x data: unpaired [w/o FFPE/frozen batch correct] ----


metadata.g2g3.up.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(211) |> 
  filter_first_G2_and_last_G3(140) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  #dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder")))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2","Grade3")) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (140)) 
    return(.)
  })()



data.g2g3.up.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.g2g3.up.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_PROBES_UNMASKED_AND_DETP)) 
    return(.)
  })()




design.g2g3.up.nc <- model.matrix(~ factor(gr.status), data=metadata.g2g3.up.nc)
fit.g2g3.up.nc <- limma::lmFit(data.g2g3.up.nc, design.g2g3.up.nc)
fit.g2g3.up.nc <- limma::eBayes(fit.g2g3.up.nc, trend=T)
stats.g2g3.up.nc <- limma::topTable(fit.g2g3.up.nc,
                                    n=nrow(data.g2g3.up.nc),
                                    coef="factor(gr.status)Grade3",
                                    sort.by = "none",
                                    adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.g2g3.up.nc)


sum(stats.g2g3.up.nc$P.Value < 0.01)
sum(stats.g2g3.up.nc$adj.P.Val < 0.01)



saveRDS(fit.g2g3.pp.nc, file="cache/analysis_differential__g2_g3__unpaired_nc__fit.Rds")
saveRDS(stats.g2g3.pp.nc, file="cache/analysis_differential__g2_g3__unpaired_nc__stats.Rds")



## x data: unpaired [w/ FFPE/frozen batch correct] ----


metadata.g2g3.up.cor <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(211) |> 
  filter_first_G2_and_last_G3(140) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(batch = factor(ifelse(isolation_person_name == "USA / Duke", "Batch US [FF]", "Batch EU [FFPE]"))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2","Grade3")) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (140)) 
    return(.)
  })()



data.g2g3.up.cor <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.g2g3.up.cor$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_PROBES_UNMASKED_AND_DETP)) 
    return(.)
  })()





design.g2g3.up.cor <- model.matrix(~factor(batch) + factor(gr.status), data=metadata.g2g3.up.cor)
fit.g2g3.up.cor <- limma::lmFit(data.g2g3.up.cor, design.g2g3.up.cor)
fit.g2g3.up.cor <- limma::eBayes(fit.g2g3.up.cor, trend=T)
stats.g2g3.up.cor <- limma::topTable(fit.g2g3.up.cor,
                                     n=nrow(data.g2g3.up.cor),
                                     coef="factor(gr.status)Grade3",
                                     sort.by = "none",
                                     adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.g2g3.up.cor)


sum(stats.g2g3.up.cor$P.Value < 0.01)
sum(stats.g2g3.up.cor$adj.P.Val < 0.01)



saveRDS(fit.g2g3.up.cor, file="cache/analysis_differential__g2_g3__unpaired_cor__fit.Rds")
saveRDS(stats.g2g3.up.cor, file="cache/analysis_differential__g2_g3__unpaired_cor__stats.Rds")




## x plots ----


plt.a <- readRDS("cache/analysis_differential__g2_g3__partial_paired_nc__stats.Rds") |> 
  dplyr::rename_with( ~ paste0(.x, "__g2_g3__partial_paired_nc"), .cols=!matches("^probe_id$",perl = T))

plt.b <- readRDS("cache/analysis_differential__primary_recurrence__partial_paired_nc__stats.Rds") |> 
  dplyr::rename_with( ~ paste0(.x, "__primary_recurrence__partial_paired_nc"), .cols=!matches("^probe_id$",perl = T))

plt.c <- readRDS("cache/analysis_differential__ad_co__stats.Rds") |> 
  dplyr::rename_with( ~ paste0(.x, "__ad"), .cols=!matches("^probe_id$",perl = T))




plt <- dplyr::left_join(plt.a, plt.b, by=c('probe_id'='probe_id'), suffix=c('','')) |> 
  dplyr::left_join( metadata.cg_probes.epic , by=c('probe_id'='probe_id'), suffix=c('','') )




ggplot(plt, aes(x=t__g2_g3__partial_paired_nc,
                y=t__primary_recurrence__partial_paired_nc,
                col=glass_nl_prim_rec__deep_significant)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  geom_point(data=subset(plt, glass_nl_prim_rec__deep_significant == F),pch=19, cex=0.0015, alpha=0.035) +
  geom_point(data=subset(plt, glass_nl_prim_rec__deep_significant == T),pch=19, cex=0.0015, alpha=0.4) +
  #geom_bin2d(bins = 350) + 
  #geom_density_2d(h=0.8) +
  #stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", colour="white") +
  #scale_fill_continuous(type = "viridis") +
  #ggplot2::scale_fill_gradientn(colours = mixcol("gray90",col3(2)[1],0:100/100), na.value = "grey50")  +
  scale_color_manual(values=c('#c04040', 'darkblue')) +
  labs(x = "GLASS-OD: t-score (Grade 2 ~ Grade 3)", y="GLASS-OD: t-score (primary ~ recurrence)", col="GLASS-NL: deep significant") +
  theme_cellpress +
  xlim(-10,10) +
  ylim(-6,6)


ggsave("/home/r361003/volcano_pp_x_g2_g3.png", width=8.5/2,height=3.5)


### AD ----


plt <- plt.a |> 
  dplyr::left_join(plt.b, by=c('probe_id'='probe_id'), suffix=c('','')) |> 
  dplyr::left_join(plt.c, by=c('probe_id'='probe_id'), suffix=c('','')) |> 
  dplyr::left_join( metadata.cg_probes.epic , by=c('probe_id'='probe_id'), suffix=c('','') )  |> 
  dplyr::left_join(data.probes.alzheimer, by=c('probe_id'='probe_id'),suffix=c('','')) |> 
  #dplyr::mutate(col = ifelse(is.na(Beta..difference),"","Paper-reported")) |> 
  dplyr::mutate(genesUniq = gsub("PHOX|SHOX","", genesUniq)) |> 
  dplyr::mutate(genesUniq = gsub("PAXIP1|PAXBP","", genesUniq)) |>
  dplyr::mutate(genesUniq = gsub("TBXAS|TBXT|TBXA","", genesUniq)) |>
  dplyr::mutate(genesUniq = gsub("QSOX","", genesUniq)) |>
  dplyr::mutate(col = ifelse(grepl("HOX", genesUniq), "HOX probe", "non-HOX probe")) |> 
  #dplyr::mutate(col = ifelse(grepl("PAX3", genesUniq), "PAX probe", "non PAX probe"))
  #dplyr::mutate(col = ifelse(grepl("TBX", genesUniq), "TBX probe", "non TBX probe")) |> 
  #dplyr::mutate(col = ifelse(grepl("SOX", genesUniq), "SOX probe", "non SOX probe")) 


plt |> 
  dplyr::filter(col == "SOX probe") |> 
  dplyr::pull(genesUniq)


ggplot(plt, aes(x=t__g2_g3__partial_paired_nc,
                y=t__ad,
                col=col,
                label=genesUniq)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  geom_point(data= subset(plt, col != "HOX probe"), pch=19, cex=0.25, alpha=0.45) +
  geom_point(data= subset(plt, col == "HOX probe"), pch=19, cex=0.65, alpha=0.45) +
  #ggrepel::geom_text_repel(data= subset(plt, col == "X probe"), size = 2.5) +
  theme_cellpress +
  scale_color_manual(values=c("darkblue","red")) +
  labs(x = "T-score limma Alzheimer Disease ~ Control", y="T-score limma Oligodendroglioma Grade 3 ~ Grade 2")

cor(
  plt |> dplyr::filter(!is.na(t__g2_g3__partial_paired_nc) & !is.na(t__ad)) |>  dplyr::pull(t__g2_g3__partial_paired_nc),
  plt |> dplyr::filter(!is.na(t__g2_g3__partial_paired_nc) & !is.na(t__ad)) |>  dplyr::pull(t__ad)
  ,
  method = "spearman"
)

cor(
  plt |> dplyr::filter(!is.na(t__g2_g3__partial_paired_nc) & !is.na(t__ad)) |>  dplyr::pull(t__primary_recurrence__partial_paired_nc),
  plt |> dplyr::filter(!is.na(t__g2_g3__partial_paired_nc) & !is.na(t__ad)) |>  dplyr::pull(t__ad)
  ,
  method = "spearman"
)

cor(
  plt |> dplyr::filter(!is.na(t__g2_g3__partial_paired_nc) & !is.na(t__ad)) |>  dplyr::pull(t__primary_recurrence__partial_paired_nc),
  plt |> dplyr::filter(!is.na(t__g2_g3__partial_paired_nc) & !is.na(t__ad)) |>  dplyr::pull(t__g2_g3__partial_paired_nc)
  ,
  method = "spearman"
)


ggplot(plt, aes(x=t__primary_recurrence__partial_paired_nc,
                y=t__ad,
                col=col)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  geom_point(data= subset(plt, col == ""), pch=19, cex=0.25, alpha=0.15) +
  geom_point(data= subset(plt, col != ""), pch=19, cex=0.65, alpha=0.45) +
  labs(x = "t-score (Grade 2 ~ Grade 3)", y="t-score (AD ~ control)") +
  theme_cellpress




### Select these strange probes ----


sslope <- 1.7
lb <- 0.2
rb <- -2.2

lb <- 0.2
rb <- -2.2

top <- -1.7
bot <- -4.0

plt <- dplyr::left_join(plt.a, plt.b, by=c('probe_id'='probe_id'), suffix=c('','')) |> 
  dplyr::left_join( metadata.cg_probes.epic , by=c('probe_id'='probe_id'), suffix=c('','') )  |> 
  dplyr::mutate(col = dplyr::case_when(
    glass_nl_prim_rec__deep_significant == T ~ "GLASS-NL deep signi",
    
    
    (t__primary_recurrence__partial_paired_nc < top) &
    (t__primary_recurrence__partial_paired_nc > bot) & 
    (t__g2_g3__partial_paired_nc > ((t__primary_recurrence__partial_paired_nc - lb) / sslope)) & 
    (t__g2_g3__partial_paired_nc < ((t__primary_recurrence__partial_paired_nc - rb) / sslope))
      
      ~ "SEL",
    
    T ~ "other"
  )) |> 
  dplyr::mutate(probe_suffix = gsub("^.+(....)$","\\1", AlleleB_ProbeSeq))



ggplot(plt, aes(x=t__g2_g3__partial_paired_nc,
                y=t__primary_recurrence__partial_paired_nc,
                col=col)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=19, cex=0.001, alpha=0.15) + 
  
  geom_abline(intercept = lb, slope=sslope, col="darkgray", lty=2) +
  geom_abline(intercept = rb, slope=sslope, col="red", lty=2) +

  geom_hline(yintercept=top, col="purple", lty=2) +
  geom_hline(yintercept=bot, col="red", lty=2) +
  
  
  labs(x = "t-score (Grade 2 ~ Grade 3)", y="t-score (primary ~ recurrence)") +
  theme_bw() + #theme_cellpress +
  xlim(-10,10) +
  ylim(-6,6)



ggplot(plt, aes(y=col, x=as.factor(probe_suffix))) +
  geom_bar(stat="identity")
plot(table(plt |>  dplyr::select(col, probe_suffix)) / rowSums(table(plt |>  dplyr::select(col, probe_suffix))))



prbs <- plt |> 
  dplyr::filter(col == "SEL") |>
  dplyr::filter(!is.na(genesUniq)) |> 
  dplyr::pull(genesUniq) |> 
  stringr::str_split(";") |> 
  unlist() |> 
  table() |> 
  sort(decreasing=T)

prbs <- plt |> 
  dplyr::filter(col == "SEL") |>
  dplyr::filter(!is.na(genesUniq)) |> 
  dplyr::pull(probe_id) 


gplot <- plt |> 
  dplyr::filter(col == "SEL") |>
  dplyr::mutate(chr = factor(CpG_chrm, levels=gtools::mixedsort(unique(as.character(CpG_chrm))) ))

p1 = ggplot(gplot, aes(x=pos, y=logFC__primary_recurrence__partial_paired_nc, col=!grepl("COL",genesUniq))) +
  #facet_wrap(~chr, scales="free_x", ncol=23) +
  facet_grid(cols = vars(chr), scales = "free", space="free") + 
  #geom_vline(xintercept=47000000) +
  #geom_vline(xintercept=50000000) +
  geom_point(pch=19,cex=0.01)



gplot2 <- plt |> 
  dplyr::filter(col != "SEL") |>
  dplyr::mutate(chr = factor(CpG_chrm, levels=gtools::mixedsort(unique(as.character(CpG_chrm))) ))

p2 = ggplot(gplot2, aes(x=pos, y=logFC__primary_recurrence__partial_paired_nc, col=!grepl("COL",genesUniq))) +
  #facet_wrap(~chr, scales="free_x", ncol=23) +
  facet_grid(cols = vars(chr), scales = "free", space="free") + 
  #geom_vline(xintercept=47000000) +
  #geom_vline(xintercept=50000000) +
  geom_point(pch=19,cex=0.001,alpha=0.5)


p1 / p2



### mk PCA of these probes specifically ----



m <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) 

d <- data.mvalues.hq_samples |> 
  dplyr::select(m$array_sentrix_id) |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes & probe_id %in% prbs) |> 
  tibble::column_to_rownames('probe_id')

dim(m)
dim(d)


pc <- prcomp(t(d))
head(pc$x)

factoextra::fviz_eig(pc)

prim_rec_secondary_effect <- pc |>
  purrr::pluck('x') |> 
  as.data.frame() |>
  dplyr::select(PC1) |> 
  dplyr::mutate(inverse = sum(pc$rotation[,1] < 0) > sum(pc$rotation[,1] > 0)) |>
  dplyr::mutate(PC1 = ifelse(inverse , -PC1, PC1), inverse = NULL ) |>
  dplyr::rename(prim_rec_secondary_effect_PC = PC1) |> 
  tibble::rownames_to_column('array_sentrix_id')


saveRDS(prim_rec_secondary_effect, file="cache/analysis_differential__prim_rec_secondary_effect_PC.Rds")



### plt & figure out stuff ----


plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(211) |> 
  #filter_primaries_and_last_recurrences(179) |> 
  dplyr::left_join(prim_rec_secondary_effect, by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')) |> 
  dplyr::mutate(is.prim = resection_number == 1) |> 
  dplyr::mutate(time_since_ok = Sys.Date() - resection_date) |> 
  dplyr::mutate(batch = ifelse(isolation_person_name == "USA / Duke", "USA", "EU")) |> 
  dplyr::mutate(col = paste0("Real: ", isolation_material, "    MNP: ", array_mnp_QC_v12.8_predicted_sample_type)) 




ggplot(plt, aes(x=-time_between_resection_and_array, y=prim_rec_secondary_effect_PC, col=col, label=isolation_id)) +
  geom_point() +
  ggrepel::geom_text_repel(data=subset(plt, col=="EU-ffpe"), col="black", size=2.5) 




ggplot(plt, aes(x=patient_id, y=prim_rec_secondary_effect_PC, col=is.prim)) +
  theme_cellpress +
  geom_point()



plt2 <- plt |> 
  dplyr::filter(!is.na(prim_rec_secondary_effect_PC))
plot(plt2$array_PC1, plt2$prim_rec_secondary_effect_PC, col=as.numeric(as.factor(plt2$batch)))
plot(plt2$array_PC2, plt2$prim_rec_secondary_effect_PC, col=as.numeric(as.factor(plt2$batch)))
plot(plt2$array_PC3, plt2$prim_rec_secondary_effect_PC, col=as.numeric(as.factor(plt2$batch)))
plot(plt2$array_PC4, plt2$prim_rec_secondary_effect_PC, col=as.numeric(as.factor(plt2$batch)))
plot(plt2$array_PC5, plt2$prim_rec_secondary_effect_PC, col=as.numeric(as.factor(plt2$batch)))
plot(plt2$array_PC6, plt2$prim_rec_secondary_effect_PC, col=as.numeric(as.factor(plt2$batch)))





# analyses: Treatment(s) ----
## data ----



metadata.trtmnt.quality.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(!is.na(resection_treatment_status_chemo)) |> 
  dplyr::filter(!is.na(resection_treatment_status_radio)) |> 
  filter_first_G2_and_last_G3(155) |> 
  
  dplyr::mutate(in_ffpe = isolation_material == "ffpe") |> 
  dplyr::filter(!is.na(in_ffpe)) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder")))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2", "Grade3")) |>
  dplyr::mutate(TMZ = grepl("TMZ", resection_treatment_status_summary)) |> 
  dplyr::mutate(PCV = grepl("PCV", resection_treatment_status_summary)) |> 
  
  dplyr::mutate(chemo = ifelse(resection_treatment_status_chemo, "Yes", "No")) |> 
  dplyr::mutate(TMZ = ifelse(TMZ, "Yes", "No")) |> 
  dplyr::mutate(PCV = ifelse(PCV, "Yes", "No")) |> 
  dplyr::mutate(radio = ifelse(resection_treatment_status_radio, "Yes", "No")) |> 
  dplyr::mutate(treated = ifelse(resection_treatment_status_chemo | resection_treatment_status_radio, "Yes", "No")) |> 
  
  dplyr::mutate(chemo = factor(chemo, levels=c("No", "Yes"))) |> 
  dplyr::mutate(radio = factor(radio, levels=c("No", "Yes"))) |> 
  dplyr::mutate(TMZ = factor(TMZ, levels=c("No", "Yes"))) |> 
  dplyr::mutate(PCV = factor(PCV, levels=c("No", "Yes"))) |> 
  dplyr::mutate(treated = factor(treated, levels=c("No", "Yes"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 149) 
    return(.)
  })() 


data.trtmnt.quality.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.trtmnt.quality.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()




## test: PC1 + pat + trt1 + trt2 ----



design.trtmnt.quality.pp.nc <- model.matrix(~
                                              array_PC1 +
                                              factor(patient) +
                                              factor(gr.status) *
                                              factor(TMZ)
                                              #factor(treated)
                                              #factor(radio)
                                            , data=metadata.trtmnt.quality.pp.nc)
fit.trtmnt.quality.pp.nc <- limma::lmFit(data.trtmnt.quality.pp.nc, design.trtmnt.quality.pp.nc)
fit.trtmnt.quality.pp.nc <- limma::eBayes(fit.trtmnt.quality.pp.nc, trend=T)



colnames(fit.trtmnt.quality.pp.nc$coefficients)


tmp <- limma::topTable(fit.trtmnt.quality.pp.nc,
                       n=nrow(data.trtmnt.quality.pp.nc),
                       coef="factor(gr.status)Grade3",
                       sort.by = "none",
                       adjust.method="fdr") |> 
  dplyr::select(t, adj.P.Val) |> 
  tibble::rownames_to_column('probe_id') 

sum(tmp$adj.P.Val < 0.01)




tmp <- limma::topTable(fit.trtmnt.quality.pp.nc,
                       n=nrow(data.trtmnt.quality.pp.nc),
                       coef="factor(TMZ)Yes",
                       sort.by = "none",
                       adjust.method="fdr") |> 
  dplyr::select(t, adj.P.Val) |> 
  tibble::rownames_to_column('probe_id') 

sum(tmp$adj.P.Val < 0.01)



tmp <- limma::topTable(fit.trtmnt.quality.pp.nc,
                       n=nrow(data.trtmnt.quality.pp.nc),
                       coef="factor(gr.status)Grade3:factor(TMZ)Yes",
                       sort.by = "none",
                       adjust.method="fdr") |> 
  dplyr::select(t, adj.P.Val) |> 
  tibble::rownames_to_column('probe_id') 


sum(tmp$adj.P.Val < 0.01)




## test: PC1 + pat + grade + trt1 + trt2 ----


## grade in non-treated ----



metadata.trtmnt.quality.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(!is.na(resection_treatment_status_chemo)) |> 
  dplyr::filter(!is.na(resection_treatment_status_radio)) |> 
  filter_first_G2_and_last_G3(155) |> 
  
  #dplyr::filter(!is.na(resection_treatment_status_summary != "no") & resection_treatment_status_summary != "no") |> 
  #dplyr::filter(resection_tumor_grade == 2) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder")))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2", "Grade3")) |>
  dplyr::mutate(TMZ = grepl("TMZ", resection_treatment_status_summary)) |> 
  dplyr::mutate(PCV = grepl("PCV", resection_treatment_status_summary)) |> 
  
  dplyr::mutate(chemo = ifelse(resection_treatment_status_chemo, "Yes", "No")) |> 
  dplyr::mutate(TMZ = ifelse(TMZ, "Yes", "No")) |> 
  dplyr::mutate(PCV = ifelse(PCV, "Yes", "No")) |> 
  dplyr::mutate(radio = ifelse(resection_treatment_status_radio, "Yes", "No")) |> 
  dplyr::mutate(treated = ifelse(resection_treatment_status_chemo | resection_treatment_status_radio, "Yes", "No")) |> 
  
  dplyr::mutate(condition = ifelse(gr.status == "Grade3" & chemo == "Yes", "G3TRT","other")) |> 
  
  dplyr::mutate(chemo = factor(chemo, levels=c("No", "Yes"))) |> 
  dplyr::mutate(radio = factor(radio, levels=c("No", "Yes"))) |> 
  dplyr::mutate(TMZ = factor(TMZ, levels=c("No", "Yes"))) |> 
  dplyr::mutate(PCV = factor(PCV, levels=c("No", "Yes"))) |> 
  dplyr::mutate(treated = factor(treated, levels=c("No", "Yes"))) |> 
  dplyr::mutate(condition = factor(condition, levels=c("other", "G3TRT"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 155) 
    return(.)
  })() 


data.trtmnt.quality.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.trtmnt.quality.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


table(metadata.trtmnt.quality.pp.nc$condition)


design.trtmnt.quality.pp.nc <- model.matrix(~
                                              array_PC1 +
                                              #factor(in_ffpe) +
                                              factor(patient) +
                                              #factor(gr.status)
                                              factor(condition)
                                              #factor(radio)
                                            , data=metadata.trtmnt.quality.pp.nc)
fit.trtmnt.quality.pp.nc <- limma::lmFit(data.trtmnt.quality.pp.nc, design.trtmnt.quality.pp.nc)
fit.trtmnt.quality.pp.nc <- limma::eBayes(fit.trtmnt.quality.pp.nc, trend=T)


colnames(fit.trtmnt.quality.pp.nc)


tmp <- limma::topTable(fit.trtmnt.quality.pp.nc,
     n=nrow(data.trtmnt.quality.pp.nc),
     coef="factor(condition)G3TRT",
     sort.by = "none",
     adjust.method="fdr") |> 
  dplyr::select(t, adj.P.Val) |> 
  dplyr::rename_with( ~ paste0("DMP_trt__grade_", .x)) |>
  tibble::rownames_to_column('probe_id') 


sum(tmp$DMP_trt__grade_adj.P.Val < 0.01)




## chemo ----


design.trtmnt.quality.pp.nc <- model.matrix(~
                                              #array_PC1 +
                                              factor(patient) +
                                              factor(chemo)
                                            , data=metadata.trtmnt.quality.pp.nc)
fit.trtmnt.quality.pp.nc <- limma::lmFit(data.trtmnt.quality.pp.nc, design.trtmnt.quality.pp.nc)
fit.trtmnt.quality.pp.nc <- limma::eBayes(fit.trtmnt.quality.pp.nc, trend=T)



tmp <- limma::topTable(fit.trtmnt.quality.pp.nc,
                       n=nrow(data.trtmnt.quality.pp.nc),
                       coef="factor(chemo)Yes",
                       sort.by = "none",
                       adjust.method="fdr") |> 
  dplyr::select(t, adj.P.Val) |> 
  dplyr::rename_with( ~ paste0("DMP_trt__chemo_", .x)) |>
  tibble::rownames_to_column('probe_id') 


sum(tmp$DMP_trt__chemo_adj.P.Val < 0.01)


# plot


plt <- data.mvalues.probes |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'))

ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP_trt__chemo_t)) + 
  geom_point(pch=19,cex=0.01, alpha=0.1)


ggplot(plt, aes(x=DMP__g2_g3__pp_nc_PC1__t, y=DMP_trt__chemo_t)) + 
  geom_point(pch=19,cex=0.01, alpha=0.1)


ggplot(plt, aes(x=DMP__PCs__pp_nc__PC1_t, y=DMP_trt__chemo_t)) +  # DMP__g2_g3__pp_nc_PC1__t
  geom_point(pch=19,cex=0.01, alpha=0.1)



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t)) + 
  geom_point(pch=19,cex=0.01, alpha=0.1)


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP_trt__chemo_t)) + 
  geom_point(pch=19,cex=0.01, alpha=0.1) +
  scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-8, 8), oob = scales::squish)




## chemo + grade ----


design.trtmnt.quality.pp.nc <- model.matrix(~
                                              array_PC1 +
                                              #factor(in_ffpe) +
                                              factor(patient) +
                                              factor(gr.status) + factor(chemo)
                                            #factor(radio)
                                            , data=metadata.trtmnt.quality.pp.nc)
fit.trtmnt.quality.pp.nc <- limma::lmFit(data.trtmnt.quality.pp.nc, design.trtmnt.quality.pp.nc)
fit.trtmnt.quality.pp.nc <- limma::eBayes(fit.trtmnt.quality.pp.nc, trend=T)


tmp <- limma::topTable(fit.trtmnt.quality.pp.nc,
                       n=nrow(data.trtmnt.quality.pp.nc),
                       coef="factor(chemo)Yes",
                       sort.by = "none",
                       adjust.method="fdr") |> 
  dplyr::select(t, adj.P.Val) |> 
  tibble::rownames_to_column('probe_id') 

sum(tmp$adj.P.Val < 0.01)


tmp <- limma::topTable(fit.trtmnt.quality.pp.nc,
                       n=nrow(data.trtmnt.quality.pp.nc),
                       coef="factor(gr.status)Grade3",
                       sort.by = "none",
                       adjust.method="fdr") |> 
  dplyr::select(t, adj.P.Val) |> 
  tibble::rownames_to_column('probe_id') 

sum(tmp$adj.P.Val < 0.01)




## TMZ ----


design.trtmnt.quality.pp.nc <- model.matrix(~
                                              array_PC1 +
                                              #factor(in_ffpe) +
                                              factor(patient) +
                                              factor(gr.status) +
                                              factor(TMZ)
                                            , data=metadata.trtmnt.quality.pp.nc)
fit.trtmnt.quality.pp.nc <- limma::lmFit(data.trtmnt.quality.pp.nc, design.trtmnt.quality.pp.nc)
fit.trtmnt.quality.pp.nc <- limma::eBayes(fit.trtmnt.quality.pp.nc, trend=T)


tmp <- limma::topTable(fit.trtmnt.quality.pp.nc,
                       n=nrow(data.trtmnt.quality.pp.nc),
                       coef="factor(TMZ)Yes",
                       sort.by = "none",
                       adjust.method="fdr") |> 
  dplyr::select(t, adj.P.Val) |> 
  tibble::rownames_to_column('probe_id') 

sum(tmp$adj.P.Val < 0.01)


tmp <- limma::topTable(fit.trtmnt.quality.pp.nc,
                       n=nrow(data.trtmnt.quality.pp.nc),
                       coef="factor(gr.status)Grade3",
                       sort.by = "none",
                       adjust.method="fdr") |> 
  dplyr::select(t, adj.P.Val) |> 
  tibble::rownames_to_column('probe_id') 

sum(tmp$adj.P.Val < 0.01)



## PCV ----


design.trtmnt.quality.pp.nc <- model.matrix(~
                                              array_PC1 +
                                              #factor(in_ffpe) +
                                              factor(patient) +
                                              factor(gr.status) +
                                              factor(PCV) +
                                              factor(TMZ)
                                            , data=metadata.trtmnt.quality.pp.nc)
fit.trtmnt.quality.pp.nc <- limma::lmFit(data.trtmnt.quality.pp.nc, design.trtmnt.quality.pp.nc)
fit.trtmnt.quality.pp.nc <- limma::eBayes(fit.trtmnt.quality.pp.nc, trend=T)


tmp <- limma::topTable(fit.trtmnt.quality.pp.nc,
                       n=nrow(data.trtmnt.quality.pp.nc),
                       coef="factor(PCV)Yes",
                       sort.by = "none",
                       adjust.method="fdr") |> 
  dplyr::select(t, adj.P.Val) |> 
  tibble::rownames_to_column('probe_id') 

sum(tmp$adj.P.Val < 0.01)

tmp <- limma::topTable(fit.trtmnt.quality.pp.nc,
                       n=nrow(data.trtmnt.quality.pp.nc),
                       coef="factor(TMZ)Yes",
                       sort.by = "none",
                       adjust.method="fdr") |> 
  dplyr::select(t, adj.P.Val) |> 
  tibble::rownames_to_column('probe_id') 

sum(tmp$adj.P.Val < 0.01)


tmp <- limma::topTable(fit.trtmnt.quality.pp.nc,
                       n=nrow(data.trtmnt.quality.pp.nc),
                       coef="factor(gr.status)Grade3",
                       sort.by = "none",
                       adjust.method="fdr") |> 
  dplyr::select(t, adj.P.Val) |> 
  tibble::rownames_to_column('probe_id') 

sum(tmp$adj.P.Val < 0.01)




## radio ----


design.trtmnt.quality.pp.nc <- model.matrix(~
                                              array_PC1 +
                                              #factor(in_ffpe) +
                                              factor(patient) +
                                              #factor(gr.status) +
                                              #factor(chemo)
                                              factor(radio)
                                            , data=metadata.trtmnt.quality.pp.nc)
fit.trtmnt.quality.pp.nc <- limma::lmFit(data.trtmnt.quality.pp.nc, design.trtmnt.quality.pp.nc)
fit.trtmnt.quality.pp.nc <- limma::eBayes(fit.trtmnt.quality.pp.nc, trend=T)



tmp <- limma::topTable(fit.trtmnt.quality.pp.nc,
                                                     n=nrow(data.trtmnt.quality.pp.nc),
                                                     coef="factor(radio)Yes",
                                                     sort.by = "none",
                                                     adjust.method="fdr") |> 
  dplyr::select(t, adj.P.Val) |> 
  dplyr::rename_with( ~ paste0("DMP_trt__radio_", .x)) |>
  tibble::rownames_to_column('probe_id') 


sum(tmp$DMP_trt__radio_adj.P.Val < 0.01)


## tables ----


metadata.trtmnt.quality.pp.nc |> 
  dplyr::select(chemo, gr.status) |> 
  table()


metadata.trtmnt.quality.pp.nc |> 
  dplyr::select(radio, gr.status) |> 
  table()


metadata.trtmnt.quality.pp.nc |> 
  dplyr::select(chemo, radio) |> 
  table()


# analyses: GLASS-OD PC1-8 ----


metadata.PCs.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  assertr::verify(!is.na(array_PC1)) |> 
  assertr::verify(is.numeric(array_PC1)) |> 
  dplyr::mutate(array_PC1 = scale(array_PC1)) |> 
  dplyr::mutate(array_PC2 = scale(array_PC2)) |> 
  dplyr::mutate(array_PC3 = scale(array_PC3)) |> 
  dplyr::mutate(array_PC4 = scale(array_PC4)) |> 
  dplyr::mutate(array_PC5 = scale(array_PC5)) |> 
  dplyr::mutate(array_PC6 = scale(array_PC6)) |> 
  dplyr::mutate(array_PC7 = scale(array_PC7)) |> 
  dplyr::mutate(array_PC8 = scale(array_PC8)) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |>  # also pairs with n = 3 & n= 4
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_",ifelse(is.paired,patient_id,"remainder")))) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_GLASS_OD_INCLUDED_SAMPLES)) 
    return(.)
  })()




data.PCs.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.PCs.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()




design.PCs.pp.nc <- model.matrix(~factor(patient) + array_PC1 + array_PC2 +
                                   array_PC3 + array_PC4 + array_PC5 + 
                                   array_PC6 + array_PC7 + array_PC8 , data=metadata.PCs.pp.nc)
fit.PCs.pp.nc <- limma::lmFit(data.PCs.pp.nc, design.PCs.pp.nc)
fit.PCs.pp.nc <- limma::eBayes(fit.PCs.pp.nc, trend=T)

stats.PCs.PC1.pp.nc <- limma::topTable(fit.PCs.pp.nc,
                                       n=nrow(data.PCs.pp.nc),
                                       coef="array_PC1",
                                       sort.by = "none",
                                       adjust.method="fdr") |> 
  dplyr::rename_with(~paste0("PC1_", .x)) |> 
  tibble::rownames_to_column('probe_id') 

stats.PCs.PC2.pp.nc <- limma::topTable(fit.PCs.pp.nc,
                                       n=nrow(data.PCs.pp.nc),
                                       coef="array_PC2",
                                       sort.by = "none",
                                       adjust.method="fdr") |> 
  dplyr::rename_with(~paste0("PC2_", .x)) |> 
  tibble::rownames_to_column('probe_id') 

stats.PCs.PC3.pp.nc <- limma::topTable(fit.PCs.pp.nc,
                                       n=nrow(data.PCs.pp.nc),
                                       coef="array_PC3",
                                       sort.by = "none",
                                       adjust.method="fdr") |> 
  dplyr::rename_with(~paste0("PC3_", .x)) |> 
  tibble::rownames_to_column('probe_id') 


stats.PCs.PC4.pp.nc <- limma::topTable(fit.PCs.pp.nc,
                                       n=nrow(data.PCs.pp.nc),
                                       coef="array_PC4",
                                       sort.by = "none",
                                       adjust.method="fdr") |> 
  dplyr::rename_with(~paste0("PC4_", .x)) |> 
  tibble::rownames_to_column('probe_id') 


stats.PCs.PC5.pp.nc <- limma::topTable(fit.PCs.pp.nc,
                                       n=nrow(data.PCs.pp.nc),
                                       coef="array_PC5",
                                       sort.by = "none",
                                       adjust.method="fdr") |> 
  dplyr::rename_with(~paste0("PC5_", .x)) |> 
  tibble::rownames_to_column('probe_id') 


stats.PCs.PC6.pp.nc <- limma::topTable(fit.PCs.pp.nc,
                                       n=nrow(data.PCs.pp.nc),
                                       coef="array_PC6",
                                       sort.by = "none",
                                       adjust.method="fdr") |> 
  dplyr::rename_with(~paste0("PC6_", .x)) |> 
  tibble::rownames_to_column('probe_id') 


stats.PCs.PC7.pp.nc <- limma::topTable(fit.PCs.pp.nc,
                                       n=nrow(data.PCs.pp.nc),
                                       coef="array_PC7",
                                       sort.by = "none",
                                       adjust.method="fdr") |> 
  dplyr::rename_with(~paste0("PC7_", .x)) |> 
  tibble::rownames_to_column('probe_id') 


stats.PCs.PC8.pp.nc <- limma::topTable(fit.PCs.pp.nc,
                                       n=nrow(data.PCs.pp.nc),
                                       coef="array_PC8",
                                       sort.by = "none",
                                       adjust.method="fdr") |> 
  dplyr::rename_with(~paste0("PC8_", .x)) |> 
  tibble::rownames_to_column('probe_id') 


stats.PCs.pp.nc <- stats.PCs.PC1.pp.nc |> 
  dplyr::left_join(stats.PCs.PC2.pp.nc, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(stats.PCs.PC3.pp.nc, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(stats.PCs.PC4.pp.nc, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(stats.PCs.PC5.pp.nc, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(stats.PCs.PC6.pp.nc, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(stats.PCs.PC7.pp.nc, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(stats.PCs.PC8.pp.nc, by=c('probe_id'='probe_id'))



rm(design.PCs.pp.nc)
rm(stats.PCs.PC1.pp.nc, stats.PCs.PC2.pp.nc, stats.PCs.PC3.pp.nc, stats.PCs.PC4.pp.nc,
   stats.PCs.PC5.pp.nc, stats.PCs.PC6.pp.nc, stats.PCs.PC7.pp.nc, stats.PCs.PC8.pp.nc)


colnames(stats.PCs.pp.nc)



saveRDS(stats.PCs.pp.nc, file="cache/analysis_differential__PCs__partial_paired_nc__stats.Rds")



# analysis: QC metrics ----

# array_percentage.detP.signi # log


metadata.qc.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |>  # also pairs with n = 3 & n= 4
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_",ifelse(is.paired,patient_id,"remainder")))) |> 
  
  dplyr::mutate( `array_qc_BISULFITE_CONVERSION_I_Beta_I-3_Beta_larger_0.12_0.18` = NULL) |> 
  dplyr::mutate( `array_qc_BISULFITE_CONVERSION_I_Beta_I-6_Beta_larger_0.2_0.3` = NULL) |> 
  dplyr::mutate(pct_detP_signi = scale(log(array_percentage.detP.signi))) |> 
  dplyr::rename( array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000 = `array_qc_STAINING_Biotin_(High)_Grn_smaller_6000_2000` ) |> 
  dplyr::rename( array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000 = `array_qc_STAINING_Biotin_(Bkg)_Grn_larger_500_1000` ) |> 
  dplyr::rename( array_qc_STAINING_DNP_High_Red_smaller_9000_3000 = `array_qc_STAINING_DNP_(High)_Red_smaller_9000_3000` ) |> 
  dplyr::rename( array_qc_STAINING_DNP_Bkg_Red_larger_750_1500 = `array_qc_STAINING_DNP_(Bkg)_Red_larger_750_1500` ) |> 
  dplyr::rename( array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000 = `array_qc_EXTENSION_Extension_(C)_Grn_smaller_20000_10000` ) |> 
  dplyr::rename( array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000 = `array_qc_EXTENSION_Extension_(G)_Grn_smaller_20000_10000` ) |> 
  dplyr::rename( array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000 = `array_qc_EXTENSION_Extension_(A)_Red_smaller_30000_15000` ) |> 
  dplyr::rename( array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000 = `array_qc_EXTENSION_Extension_(T)_Red_smaller_30000_15000` ) |> 
  dplyr::rename( array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000 = `array_qc_HYBRIDIZATION_Hyb_(High)_Grn_smaller_16000_12000` ) |> 
  dplyr::rename( array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000 = `array_qc_HYBRIDIZATION_Hyb_(Medium)_Grn_smaller_8000_6000` ) |> 
  dplyr::rename( array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000 = `array_qc_HYBRIDIZATION_Hyb_(Low)_Grn_smaller_4000_3000` ) |> 
  dplyr::rename( array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985 = `array_qc_HYBRIDIZATION_Hyb_(Correlation)_Grn_smaller_NA_0.985` ) |> 
  dplyr::rename( array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15 = `array_qc_BISULFITE_CONVERSION_I_Beta_I-1_Beta_larger_0.1_0.15` ) |> 
  dplyr::rename( array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12 = `array_qc_BISULFITE_CONVERSION_I_Beta_I-2_Beta_larger_0.08_0.12` ) |> 
  dplyr::rename( array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15 = `array_qc_BISULFITE_CONVERSION_I_Beta_I-4_Beta_larger_0.1_0.15` ) |> 
  dplyr::rename( array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3 = `array_qc_BISULFITE_CONVERSION_I_Beta_I-5_Beta_larger_0.2_0.3` ) |> 
  dplyr::rename( array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3 = `array_qc_BISULFITE_CONVERSION_II_Beta_II-1_Beta_larger_0.2_0.3` ) |> 
  dplyr::rename( array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12 = `array_qc_BISULFITE_CONVERSION_II_Beta_II-2_Beta_larger_0.08_0.12` ) |> 
  dplyr::rename( array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24 = `array_qc_BISULFITE_CONVERSION_II_Beta_II-3_Beta_larger_0.16_0.24` ) |> 
  dplyr::rename( array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12 = `array_qc_BISULFITE_CONVERSION_II_Beta_II-4_Beta_larger_0.08_0.12` ) |> 
  dplyr::rename( array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2 = `array_qc_SPECIFICITY_I_GT_Mismatch_1_(PM)_Grn_smaller_NA_1.2` ) |> 
  dplyr::rename( array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1 = `array_qc_SPECIFICITY_I_GT_Mismatch_1_(MM)_Grn_larger_NA_0.1` ) |> 
  dplyr::rename( array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1 = `array_qc_SPECIFICITY_I_GT_Mismatch_3_(MM)_Grn_larger_NA_0.1` ) |> 
  dplyr::rename( array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2 = `array_qc_SPECIFICITY_I_GT_Mismatch_4_(PM)_Red_smaller_NA_2` ) |> 
  dplyr::rename( array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5 = `array_qc_SPECIFICITY_I_GT_Mismatch_5_(PM)_Red_smaller_NA_0.5` ) |> 
  dplyr::rename( array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5 = `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5` ) |> 
  dplyr::rename( array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3 = `array_qc_SPECIFICITY_I_GT_Mismatch_5_(MM)_Red_larger_NA_0.3` ) |> 
  dplyr::rename( array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5 = `array_qc_SPECIFICITY_I_GT_Mismatch_6_(MM)_Red_larger_NA_0.5` ) |> 
  dplyr::rename( array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6 = `array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1.6` ) |> 
  dplyr::rename( array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2 = `array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0.2` ) |> 
  dplyr::rename( array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15 = `array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0.15` ) |> 
  dplyr::rename( array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15 = `array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0.15` ) |> 
  dplyr::rename( array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4 = `array_qc_NON-POLYMORPHIC_NP_(A)_Grn_larger_NA_0.4` ) |> 
  dplyr::rename( array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4 = `array_qc_NON-POLYMORPHIC_NP_(T)_Grn_larger_NA_0.4` ) |> 
  dplyr::rename( array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5 = `array_qc_NON-POLYMORPHIC_NP_(C)_Grn_smaller_NA_1.5` ) |> 
  dplyr::rename( array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1 = `array_qc_NON-POLYMORPHIC_NP_(G)_Grn_smaller_NA_1` ) |> 
  dplyr::rename( array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5 = `array_qc_NON-POLYMORPHIC_NP_(A)_Red_smaller_NA_1.5` ) |> 
  dplyr::rename( array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5 = `array_qc_NON-POLYMORPHIC_NP_(T)_Red_smaller_NA_1.5` ) |> 
  dplyr::rename( array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4 = `array_qc_NON-POLYMORPHIC_NP_(C)_Red_larger_NA_0.4` ) |> 
  dplyr::rename( array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6 = `array_qc_NON-POLYMORPHIC_NP_(G)_Red_larger_NA_0.6` )


data.qc.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.qc.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()





## 1. pct_detP_signi ----

design.pct_detP_signi.pp.nc <- model.matrix(~factor(patient) + pct_detP_signi, data=metadata.qc.pp.nc |> dplyr::mutate(pct_detP_signi=scale(pct_detP_signi)))
fit.pct_detP_signi.pp.nc <- limma::lmFit(data.qc.pp.nc, design.pct_detP_signi.pp.nc)
fit.pct_detP_signi.pp.nc <- limma::eBayes(fit.pct_detP_signi.pp.nc, trend=T)
stats.pct_detP_signi.pp.nc <- limma::topTable(fit.pct_detP_signi.pp.nc,
                                              n=nrow(data.qc.pp.nc),
                                              coef="pct_detP_signi",
                                              sort.by = "none",
                                              adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t, adj.P.Val)


saveRDS(stats.pct_detP_signi.pp.nc, file="cache/analysis_differential__pct_detP_signi__partial_paired_nc__stats.Rds")


rm(
  design.pct_detP_signi.pp.nc,
  fit.pct_detP_signi.pp.nc,
  stats.pct_detP_signi.pp.nc
)




## 2. array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000 ----

design.array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000.pp.nc <- model.matrix(~factor(patient) + array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000, data=metadata.qc.pp.nc)
fit.array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000.pp.nc)
fit.array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000.pp.nc <- limma::eBayes(fit.array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000.pp.nc, trend=T)
stats.array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000.pp.nc <- limma::topTable(fit.array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000.pp.nc,
                                                                                   n=nrow(data.qc.pp.nc),
                                                                                   coef="array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000",
                                                                                   sort.by = "none",
                                                                                   adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000.pp.nc, file="cache/analysis_differential__array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000.pp.nc,
  fit.array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000.pp.nc,
  stats.array_qc_STAINING_Biotin_High_Grn_smaller_6000_2000.pp.nc
)




## 3. array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000 ----

design.array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000.pp.nc <- model.matrix(~factor(patient) + array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000, data=metadata.qc.pp.nc)
fit.array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000.pp.nc)
fit.array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000.pp.nc <- limma::eBayes(fit.array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000.pp.nc, trend=T)
stats.array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000.pp.nc <- limma::topTable(fit.array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000.pp.nc,
                                                                                n=nrow(data.qc.pp.nc),
                                                                                coef="array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000",
                                                                                sort.by = "none",
                                                                                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000.pp.nc, file="cache/analysis_differential__array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000.pp.nc,
  fit.array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000.pp.nc,
  stats.array_qc_STAINING_Biotin_Bkg_Grn_larger_500_1000.pp.nc
)




## 4. array_qc_STAINING_DNP_High_Red_smaller_9000_3000 ----

design.array_qc_STAINING_DNP_High_Red_smaller_9000_3000.pp.nc <- model.matrix(~factor(patient) + array_qc_STAINING_DNP_High_Red_smaller_9000_3000, data=metadata.qc.pp.nc)
fit.array_qc_STAINING_DNP_High_Red_smaller_9000_3000.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_STAINING_DNP_High_Red_smaller_9000_3000.pp.nc)
fit.array_qc_STAINING_DNP_High_Red_smaller_9000_3000.pp.nc <- limma::eBayes(fit.array_qc_STAINING_DNP_High_Red_smaller_9000_3000.pp.nc, trend=T)
stats.array_qc_STAINING_DNP_High_Red_smaller_9000_3000.pp.nc <- limma::topTable(fit.array_qc_STAINING_DNP_High_Red_smaller_9000_3000.pp.nc,
                                                                                n=nrow(data.qc.pp.nc),
                                                                                coef="array_qc_STAINING_DNP_High_Red_smaller_9000_3000",
                                                                                sort.by = "none",
                                                                                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_STAINING_DNP_High_Red_smaller_9000_3000.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_STAINING_DNP_High_Red_smaller_9000_3000.pp.nc, file="cache/analysis_differential__array_qc_STAINING_DNP_High_Red_smaller_9000_3000__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_STAINING_DNP_High_Red_smaller_9000_3000.pp.nc,
  fit.array_qc_STAINING_DNP_High_Red_smaller_9000_3000.pp.nc,
  stats.array_qc_STAINING_DNP_High_Red_smaller_9000_3000.pp.nc
)




## 5. array_qc_STAINING_DNP_Bkg_Red_larger_750_1500 ----

design.array_qc_STAINING_DNP_Bkg_Red_larger_750_1500.pp.nc <- model.matrix(~factor(patient) + array_qc_STAINING_DNP_Bkg_Red_larger_750_1500, data=metadata.qc.pp.nc)
fit.array_qc_STAINING_DNP_Bkg_Red_larger_750_1500.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_STAINING_DNP_Bkg_Red_larger_750_1500.pp.nc)
fit.array_qc_STAINING_DNP_Bkg_Red_larger_750_1500.pp.nc <- limma::eBayes(fit.array_qc_STAINING_DNP_Bkg_Red_larger_750_1500.pp.nc, trend=T)
stats.array_qc_STAINING_DNP_Bkg_Red_larger_750_1500.pp.nc <- limma::topTable(fit.array_qc_STAINING_DNP_Bkg_Red_larger_750_1500.pp.nc,
                                                                             n=nrow(data.qc.pp.nc),
                                                                             coef="array_qc_STAINING_DNP_Bkg_Red_larger_750_1500",
                                                                             sort.by = "none",
                                                                             adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_STAINING_DNP_Bkg_Red_larger_750_1500.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_STAINING_DNP_Bkg_Red_larger_750_1500.pp.nc, file="cache/analysis_differential__array_qc_STAINING_DNP_Bkg_Red_larger_750_1500__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_STAINING_DNP_Bkg_Red_larger_750_1500.pp.nc,
  fit.array_qc_STAINING_DNP_Bkg_Red_larger_750_1500.pp.nc,
  stats.array_qc_STAINING_DNP_Bkg_Red_larger_750_1500.pp.nc
)




## 6. array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000 ----

design.array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000.pp.nc <- model.matrix(~factor(patient) + array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000, data=metadata.qc.pp.nc)
fit.array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000.pp.nc)
fit.array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000.pp.nc <- limma::eBayes(fit.array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000.pp.nc, trend=T)
stats.array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000.pp.nc <- limma::topTable(fit.array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000.pp.nc,
                                                                                      n=nrow(data.qc.pp.nc),
                                                                                      coef="array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000",
                                                                                      sort.by = "none",
                                                                                      adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000.pp.nc, file="cache/analysis_differential__array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000.pp.nc,
  fit.array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000.pp.nc,
  stats.array_qc_EXTENSION_Extension_C_Grn_smaller_20000_10000.pp.nc
)




## 7. array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000 ----

design.array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000.pp.nc <- model.matrix(~factor(patient) + array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000, data=metadata.qc.pp.nc)
fit.array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000.pp.nc)
fit.array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000.pp.nc <- limma::eBayes(fit.array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000.pp.nc, trend=T)
stats.array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000.pp.nc <- limma::topTable(fit.array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000.pp.nc,
                                                                                      n=nrow(data.qc.pp.nc),
                                                                                      coef="array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000",
                                                                                      sort.by = "none",
                                                                                      adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000.pp.nc, file="cache/analysis_differential__array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000.pp.nc,
  fit.array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000.pp.nc,
  stats.array_qc_EXTENSION_Extension_G_Grn_smaller_20000_10000.pp.nc
)




## 8. array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000 ----

design.array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000.pp.nc <- model.matrix(~factor(patient) + array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000, data=metadata.qc.pp.nc)
fit.array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000.pp.nc)
fit.array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000.pp.nc <- limma::eBayes(fit.array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000.pp.nc, trend=T)
stats.array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000.pp.nc <- limma::topTable(fit.array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000.pp.nc,
                                                                                      n=nrow(data.qc.pp.nc),
                                                                                      coef="array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000",
                                                                                      sort.by = "none",
                                                                                      adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000.pp.nc, file="cache/analysis_differential__array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000.pp.nc,
  fit.array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000.pp.nc,
  stats.array_qc_EXTENSION_Extension_A_Red_smaller_30000_15000.pp.nc
)




## 9. array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000 ----

design.array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000.pp.nc <- model.matrix(~factor(patient) + array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000, data=metadata.qc.pp.nc)
fit.array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000.pp.nc)
fit.array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000.pp.nc <- limma::eBayes(fit.array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000.pp.nc, trend=T)
stats.array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000.pp.nc <- limma::topTable(fit.array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000.pp.nc,
                                                                                      n=nrow(data.qc.pp.nc),
                                                                                      coef="array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000",
                                                                                      sort.by = "none",
                                                                                      adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000.pp.nc, file="cache/analysis_differential__array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000.pp.nc,
  fit.array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000.pp.nc,
  stats.array_qc_EXTENSION_Extension_T_Red_smaller_30000_15000.pp.nc
)




## 10. array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000 ----

design.array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000.pp.nc <- model.matrix(~factor(patient) + array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000, data=metadata.qc.pp.nc)
fit.array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000.pp.nc)
fit.array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000.pp.nc <- limma::eBayes(fit.array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000.pp.nc, trend=T)
stats.array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000.pp.nc <- limma::topTable(fit.array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000.pp.nc,
                                                                                       n=nrow(data.qc.pp.nc),
                                                                                       coef="array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000",
                                                                                       sort.by = "none",
                                                                                       adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000.pp.nc, file="cache/analysis_differential__array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000.pp.nc,
  fit.array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000.pp.nc,
  stats.array_qc_HYBRIDIZATION_Hyb_High_Grn_smaller_16000_12000.pp.nc
)




## 11. array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000 ----

design.array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000.pp.nc <- model.matrix(~factor(patient) + array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000, data=metadata.qc.pp.nc)
fit.array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000.pp.nc)
fit.array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000.pp.nc <- limma::eBayes(fit.array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000.pp.nc, trend=T)
stats.array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000.pp.nc <- limma::topTable(fit.array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000.pp.nc,
                                                                                       n=nrow(data.qc.pp.nc),
                                                                                       coef="array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000",
                                                                                       sort.by = "none",
                                                                                       adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000.pp.nc, file="cache/analysis_differential__array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000.pp.nc,
  fit.array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000.pp.nc,
  stats.array_qc_HYBRIDIZATION_Hyb_Medium_Grn_smaller_8000_6000.pp.nc
)




## 12. array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000 ----

design.array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000.pp.nc <- model.matrix(~factor(patient) + array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000, data=metadata.qc.pp.nc)
fit.array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000.pp.nc)
fit.array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000.pp.nc <- limma::eBayes(fit.array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000.pp.nc, trend=T)
stats.array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000.pp.nc <- limma::topTable(fit.array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000.pp.nc,
                                                                                    n=nrow(data.qc.pp.nc),
                                                                                    coef="array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000",
                                                                                    sort.by = "none",
                                                                                    adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000.pp.nc, file="cache/analysis_differential__array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000.pp.nc,
  fit.array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000.pp.nc,
  stats.array_qc_HYBRIDIZATION_Hyb_Low_Grn_smaller_4000_3000.pp.nc
)




## 13. array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985 ----

design.array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985.pp.nc <- model.matrix(~factor(patient) + array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985, data=metadata.qc.pp.nc)
fit.array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985.pp.nc)
fit.array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985.pp.nc <- limma::eBayes(fit.array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985.pp.nc, trend=T)
stats.array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985.pp.nc <- limma::topTable(fit.array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985.pp.nc,
                                                                                           n=nrow(data.qc.pp.nc),
                                                                                           coef="array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985",
                                                                                           sort.by = "none",
                                                                                           adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985.pp.nc, file="cache/analysis_differential__array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985.pp.nc,
  fit.array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985.pp.nc,
  stats.array_qc_HYBRIDIZATION_Hyb_Correlation_Grn_smaller_NA_0_985.pp.nc
)




## 14. array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15 ----

design.array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15.pp.nc <- model.matrix(~factor(patient) + array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15, data=metadata.qc.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15.pp.nc <- limma::eBayes(fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15.pp.nc, trend=T)
stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15.pp.nc <- limma::topTable(fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15.pp.nc,
                                                                                             n=nrow(data.qc.pp.nc),
                                                                                             coef="array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15",
                                                                                             sort.by = "none",
                                                                                             adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15.pp.nc, file="cache/analysis_differential__array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15.pp.nc,
  fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15.pp.nc,
  stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_1_Beta_larger_0_1_0_15.pp.nc
)




## 15. array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12 ----

design.array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12.pp.nc <- model.matrix(~factor(patient) + array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12, data=metadata.qc.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12.pp.nc <- limma::eBayes(fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12.pp.nc, trend=T)
stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12.pp.nc <- limma::topTable(fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12.pp.nc,
                                                                                              n=nrow(data.qc.pp.nc),
                                                                                              coef="array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12",
                                                                                              sort.by = "none",
                                                                                              adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12.pp.nc, file="cache/analysis_differential__array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12.pp.nc,
  fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12.pp.nc,
  stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_2_Beta_larger_0_08_0_12.pp.nc
)




## 16. array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15 ----

design.array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15.pp.nc <- model.matrix(~factor(patient) + array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15, data=metadata.qc.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15.pp.nc <- limma::eBayes(fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15.pp.nc, trend=T)
stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15.pp.nc <- limma::topTable(fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15.pp.nc,
                                                                                             n=nrow(data.qc.pp.nc),
                                                                                             coef="array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15",
                                                                                             sort.by = "none",
                                                                                             adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15.pp.nc, file="cache/analysis_differential__array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15.pp.nc,
  fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15.pp.nc,
  stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_4_Beta_larger_0_1_0_15.pp.nc
)




## 17. array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3 ----

design.array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3.pp.nc <- model.matrix(~factor(patient) + array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3, data=metadata.qc.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3.pp.nc <- limma::eBayes(fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3.pp.nc, trend=T)
stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3.pp.nc <- limma::topTable(fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3.pp.nc,
                                                                                            n=nrow(data.qc.pp.nc),
                                                                                            coef="array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3",
                                                                                            sort.by = "none",
                                                                                            adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3.pp.nc, file="cache/analysis_differential__array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3.pp.nc,
  fit.array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3.pp.nc,
  stats.array_qc_BISULFITE_CONVERSION_I_Beta_I_5_Beta_larger_0_2_0_3.pp.nc
)




## 18. array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3 ----

design.array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3.pp.nc <- model.matrix(~factor(patient) + array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3, data=metadata.qc.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3.pp.nc <- limma::eBayes(fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3.pp.nc, trend=T)
stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3.pp.nc <- limma::topTable(fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3.pp.nc,
                                                                                              n=nrow(data.qc.pp.nc),
                                                                                              coef="array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3",
                                                                                              sort.by = "none",
                                                                                              adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3.pp.nc, file="cache/analysis_differential__array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3.pp.nc,
  fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3.pp.nc,
  stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_1_Beta_larger_0_2_0_3.pp.nc
)




## 19. array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12 ----

design.array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12.pp.nc <- model.matrix(~factor(patient) + array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12, data=metadata.qc.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12.pp.nc <- limma::eBayes(fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12.pp.nc, trend=T)
stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12.pp.nc <- limma::topTable(fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12.pp.nc,
                                                                                                n=nrow(data.qc.pp.nc),
                                                                                                coef="array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12",
                                                                                                sort.by = "none",
                                                                                                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12.pp.nc, file="cache/analysis_differential__array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12.pp.nc,
  fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12.pp.nc,
  stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_2_Beta_larger_0_08_0_12.pp.nc
)




## 20. array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24 ----

design.array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24.pp.nc <- model.matrix(~factor(patient) + array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24, data=metadata.qc.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24.pp.nc <- limma::eBayes(fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24.pp.nc, trend=T)
stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24.pp.nc <- limma::topTable(fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24.pp.nc,
                                                                                                n=nrow(data.qc.pp.nc),
                                                                                                coef="array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24",
                                                                                                sort.by = "none",
                                                                                                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24.pp.nc, file="cache/analysis_differential__array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24.pp.nc,
  fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24.pp.nc,
  stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_3_Beta_larger_0_16_0_24.pp.nc
)




## 21. array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12 ----

design.array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12.pp.nc <- model.matrix(~factor(patient) + array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12, data=metadata.qc.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12.pp.nc)
fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12.pp.nc <- limma::eBayes(fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12.pp.nc, trend=T)
stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12.pp.nc <- limma::topTable(fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12.pp.nc,
                                                                                                n=nrow(data.qc.pp.nc),
                                                                                                coef="array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12",
                                                                                                sort.by = "none",
                                                                                                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12.pp.nc, file="cache/analysis_differential__array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12.pp.nc,
  fit.array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12.pp.nc,
  stats.array_qc_BISULFITE_CONVERSION_II_Beta_II_4_Beta_larger_0_08_0_12.pp.nc
)




## 22. array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2 ----

design.array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2.pp.nc <- model.matrix(~factor(patient) + array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2, data=metadata.qc.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2.pp.nc <- limma::eBayes(fit.array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2.pp.nc, trend=T)
stats.array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2.pp.nc <- limma::topTable(fit.array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2.pp.nc,
                                                                                          n=nrow(data.qc.pp.nc),
                                                                                          coef="array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2",
                                                                                          sort.by = "none",
                                                                                          adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2.pp.nc, file="cache/analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2.pp.nc,
  fit.array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2.pp.nc,
  stats.array_qc_SPECIFICITY_I_GT_Mismatch_1_PM_Grn_smaller_NA_1_2.pp.nc
)




## 23. array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1 ----

design.array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1.pp.nc <- model.matrix(~factor(patient) + array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1, data=metadata.qc.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1.pp.nc <- limma::eBayes(fit.array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1.pp.nc, trend=T)
stats.array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1.pp.nc <- limma::topTable(fit.array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1.pp.nc,
                                                                                         n=nrow(data.qc.pp.nc),
                                                                                         coef="array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1",
                                                                                         sort.by = "none",
                                                                                         adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1.pp.nc, file="cache/analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1.pp.nc,
  fit.array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1.pp.nc,
  stats.array_qc_SPECIFICITY_I_GT_Mismatch_1_MM_Grn_larger_NA_0_1.pp.nc
)




## 24. array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1 ----

design.array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1.pp.nc <- model.matrix(~factor(patient) + array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1, data=metadata.qc.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1.pp.nc <- limma::eBayes(fit.array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1.pp.nc, trend=T)
stats.array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1.pp.nc <- limma::topTable(fit.array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1.pp.nc,
                                                                                         n=nrow(data.qc.pp.nc),
                                                                                         coef="array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1",
                                                                                         sort.by = "none",
                                                                                         adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1.pp.nc, file="cache/analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1.pp.nc,
  fit.array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1.pp.nc,
  stats.array_qc_SPECIFICITY_I_GT_Mismatch_3_MM_Grn_larger_NA_0_1.pp.nc
)




## 25. array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2 ----

design.array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2.pp.nc <- model.matrix(~factor(patient) + array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2, data=metadata.qc.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2.pp.nc <- limma::eBayes(fit.array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2.pp.nc, trend=T)
stats.array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2.pp.nc <- limma::topTable(fit.array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2.pp.nc,
                                                                                        n=nrow(data.qc.pp.nc),
                                                                                        coef="array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2",
                                                                                        sort.by = "none",
                                                                                        adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2.pp.nc, file="cache/analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2.pp.nc,
  fit.array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2.pp.nc,
  stats.array_qc_SPECIFICITY_I_GT_Mismatch_4_PM_Red_smaller_NA_2.pp.nc
)




## 26. array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5 ----

design.array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5.pp.nc <- model.matrix(~factor(patient) + array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5, data=metadata.qc.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5.pp.nc <- limma::eBayes(fit.array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5.pp.nc, trend=T)
stats.array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5.pp.nc <- limma::topTable(fit.array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5.pp.nc,
                                                                                          n=nrow(data.qc.pp.nc),
                                                                                          coef="array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5",
                                                                                          sort.by = "none",
                                                                                          adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5.pp.nc, file="cache/analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5.pp.nc,
  fit.array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5.pp.nc,
  stats.array_qc_SPECIFICITY_I_GT_Mismatch_5_PM_Red_smaller_NA_0_5.pp.nc
)




## 27. array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5 ----

design.array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5.pp.nc <- model.matrix(~factor(patient) + array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5, data=metadata.qc.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5.pp.nc <- limma::eBayes(fit.array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5.pp.nc, trend=T)
stats.array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5.pp.nc <- limma::topTable(fit.array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5.pp.nc,
                                                                                          n=nrow(data.qc.pp.nc),
                                                                                          coef="array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5",
                                                                                          sort.by = "none",
                                                                                          adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t, adj.P.Val)

print(paste0("Restoration: ", sum(stats.array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5.pp.nc, file="cache/analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5.pp.nc,
  fit.array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5.pp.nc,
  stats.array_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5.pp.nc
)




## 28. array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3 ----

design.array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3.pp.nc <- model.matrix(~factor(patient) + array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3, data=metadata.qc.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3.pp.nc <- limma::eBayes(fit.array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3.pp.nc, trend=T)
stats.array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3.pp.nc <- limma::topTable(fit.array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3.pp.nc,
                                                                                         n=nrow(data.qc.pp.nc),
                                                                                         coef="array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3",
                                                                                         sort.by = "none",
                                                                                         adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3.pp.nc, file="cache/analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3.pp.nc,
  fit.array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3.pp.nc,
  stats.array_qc_SPECIFICITY_I_GT_Mismatch_5_MM_Red_larger_NA_0_3.pp.nc
)




## 29. array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5 ----

design.array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5.pp.nc <- model.matrix(~factor(patient) + array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5, data=metadata.qc.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5.pp.nc)
fit.array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5.pp.nc <- limma::eBayes(fit.array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5.pp.nc, trend=T)
stats.array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5.pp.nc <- limma::topTable(fit.array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5.pp.nc,
                                                                                         n=nrow(data.qc.pp.nc),
                                                                                         coef="array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5",
                                                                                         sort.by = "none",
                                                                                         adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5.pp.nc, file="cache/analysis_differential__array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5.pp.nc,
  fit.array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5.pp.nc,
  stats.array_qc_SPECIFICITY_I_GT_Mismatch_6_MM_Red_larger_NA_0_5.pp.nc
)




## 30. array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6 ----

design.array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6.pp.nc <- model.matrix(~factor(patient) + array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6, data=metadata.qc.pp.nc)
fit.array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6.pp.nc)
fit.array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6.pp.nc <- limma::eBayes(fit.array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6.pp.nc, trend=T)
stats.array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6.pp.nc <- limma::topTable(fit.array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6.pp.nc,
                                                                                        n=nrow(data.qc.pp.nc),
                                                                                        coef="array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6",
                                                                                        sort.by = "none",
                                                                                        adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6.pp.nc, file="cache/analysis_differential__array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6.pp.nc,
  fit.array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6.pp.nc,
  stats.array_qc_SPECIFICITY_II_Specificity_1_Red_smaller_NA_1_6.pp.nc
)




## 31. array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2 ----

design.array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2.pp.nc <- model.matrix(~factor(patient) + array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2, data=metadata.qc.pp.nc)
fit.array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2.pp.nc)
fit.array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2.pp.nc <- limma::eBayes(fit.array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2.pp.nc, trend=T)
stats.array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2.pp.nc <- limma::topTable(fit.array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2.pp.nc,
                                                                                       n=nrow(data.qc.pp.nc),
                                                                                       coef="array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2",
                                                                                       sort.by = "none",
                                                                                       adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2.pp.nc, file="cache/analysis_differential__array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2.pp.nc,
  fit.array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2.pp.nc,
  stats.array_qc_SPECIFICITY_II_Specificity_1_Grn_larger_NA_0_2.pp.nc
)




## 32. array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15 ----

design.array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15.pp.nc <- model.matrix(~factor(patient) + array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15, data=metadata.qc.pp.nc)
fit.array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15.pp.nc)
fit.array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15.pp.nc <- limma::eBayes(fit.array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15.pp.nc, trend=T)
stats.array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15.pp.nc <- limma::topTable(fit.array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15.pp.nc,
                                                                                        n=nrow(data.qc.pp.nc),
                                                                                        coef="array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15",
                                                                                        sort.by = "none",
                                                                                        adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15.pp.nc, file="cache/analysis_differential__array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15.pp.nc,
  fit.array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15.pp.nc,
  stats.array_qc_SPECIFICITY_II_Specificity_2_Grn_larger_NA_0_15.pp.nc
)




## 33. array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15 ----

design.array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15.pp.nc <- model.matrix(~factor(patient) + array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15, data=metadata.qc.pp.nc)
fit.array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15.pp.nc)
fit.array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15.pp.nc <- limma::eBayes(fit.array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15.pp.nc, trend=T)
stats.array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15.pp.nc <- limma::topTable(fit.array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15.pp.nc,
                                                                                        n=nrow(data.qc.pp.nc),
                                                                                        coef="array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15",
                                                                                        sort.by = "none",
                                                                                        adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15.pp.nc, file="cache/analysis_differential__array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15.pp.nc,
  fit.array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15.pp.nc,
  stats.array_qc_SPECIFICITY_II_Specificity_3_Grn_larger_NA_0_15.pp.nc
)




## 34. array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4 ----

design.array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4.pp.nc <- model.matrix(~factor(patient) + array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4, data=metadata.qc.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4.pp.nc <- limma::eBayes(fit.array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4.pp.nc, trend=T)
stats.array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4.pp.nc <- limma::topTable(fit.array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4.pp.nc,
                                                                               n=nrow(data.qc.pp.nc),
                                                                               coef="array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4",
                                                                               sort.by = "none",
                                                                               adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4.pp.nc, file="cache/analysis_differential__array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4.pp.nc,
  fit.array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4.pp.nc,
  stats.array_qc_NON_POLYMORPHIC_NP_A_Grn_larger_NA_0_4.pp.nc
)




## 35. array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4 ----

design.array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4.pp.nc <- model.matrix(~factor(patient) + array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4, data=metadata.qc.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4.pp.nc <- limma::eBayes(fit.array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4.pp.nc, trend=T)
stats.array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4.pp.nc <- limma::topTable(fit.array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4.pp.nc,
                                                                               n=nrow(data.qc.pp.nc),
                                                                               coef="array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4",
                                                                               sort.by = "none",
                                                                               adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4.pp.nc, file="cache/analysis_differential__array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4.pp.nc,
  fit.array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4.pp.nc,
  stats.array_qc_NON_POLYMORPHIC_NP_T_Grn_larger_NA_0_4.pp.nc
)




## 36. array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5 ----

design.array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5.pp.nc <- model.matrix(~factor(patient) + array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5, data=metadata.qc.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5.pp.nc <- limma::eBayes(fit.array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5.pp.nc, trend=T)
stats.array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5.pp.nc <- limma::topTable(fit.array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5.pp.nc,
                                                                                n=nrow(data.qc.pp.nc),
                                                                                coef="array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5",
                                                                                sort.by = "none",
                                                                                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5.pp.nc, file="cache/analysis_differential__array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5.pp.nc,
  fit.array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5.pp.nc,
  stats.array_qc_NON_POLYMORPHIC_NP_C_Grn_smaller_NA_1_5.pp.nc
)




## 37. array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1 ----

design.array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1.pp.nc <- model.matrix(~factor(patient) + array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1, data=metadata.qc.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1.pp.nc <- limma::eBayes(fit.array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1.pp.nc, trend=T)
stats.array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1.pp.nc <- limma::topTable(fit.array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1.pp.nc,
                                                                              n=nrow(data.qc.pp.nc),
                                                                              coef="array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1",
                                                                              sort.by = "none",
                                                                              adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1.pp.nc, file="cache/analysis_differential__array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1.pp.nc,
  fit.array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1.pp.nc,
  stats.array_qc_NON_POLYMORPHIC_NP_G_Grn_smaller_NA_1.pp.nc
)




## 38. array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5 ----

design.array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5.pp.nc <- model.matrix(~factor(patient) + array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5, data=metadata.qc.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5.pp.nc <- limma::eBayes(fit.array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5.pp.nc, trend=T)
stats.array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5.pp.nc <- limma::topTable(fit.array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5.pp.nc,
                                                                                n=nrow(data.qc.pp.nc),
                                                                                coef="array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5",
                                                                                sort.by = "none",
                                                                                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5.pp.nc, file="cache/analysis_differential__array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5.pp.nc,
  fit.array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5.pp.nc,
  stats.array_qc_NON_POLYMORPHIC_NP_A_Red_smaller_NA_1_5.pp.nc
)




## 39. array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5 ----

design.array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5.pp.nc <- model.matrix(~factor(patient) + array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5, data=metadata.qc.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5.pp.nc <- limma::eBayes(fit.array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5.pp.nc, trend=T)
stats.array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5.pp.nc <- limma::topTable(fit.array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5.pp.nc,
                                                                                n=nrow(data.qc.pp.nc),
                                                                                coef="array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5",
                                                                                sort.by = "none",
                                                                                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5.pp.nc, file="cache/analysis_differential__array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5.pp.nc,
  fit.array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5.pp.nc,
  stats.array_qc_NON_POLYMORPHIC_NP_T_Red_smaller_NA_1_5.pp.nc
)




## 40. array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4 ----

design.array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4.pp.nc <- model.matrix(~factor(patient) + array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4, data=metadata.qc.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4.pp.nc <- limma::eBayes(fit.array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4.pp.nc, trend=T)
stats.array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4.pp.nc <- limma::topTable(fit.array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4.pp.nc,
                                                                               n=nrow(data.qc.pp.nc),
                                                                               coef="array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4",
                                                                               sort.by = "none",
                                                                               adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4.pp.nc, file="cache/analysis_differential__array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4.pp.nc,
  fit.array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4.pp.nc,
  stats.array_qc_NON_POLYMORPHIC_NP_C_Red_larger_NA_0_4.pp.nc
)




## 41. array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6 ----

design.array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6.pp.nc <- model.matrix(~factor(patient) + array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6, data=metadata.qc.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6.pp.nc <- limma::lmFit(data.qc.pp.nc, design.array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6.pp.nc)
fit.array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6.pp.nc <- limma::eBayes(fit.array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6.pp.nc, trend=T)
stats.array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6.pp.nc <- limma::topTable(fit.array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6.pp.nc,
                                                                               n=nrow(data.qc.pp.nc),
                                                                               coef="array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6",
                                                                               sort.by = "none",
                                                                               adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6.pp.nc, file="cache/analysis_differential__array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6__partial_paired_nc__stats.Rds")

rm(
  design.array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6.pp.nc,
  fit.array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6.pp.nc,
  stats.array_qc_NON_POLYMORPHIC_NP_G_Red_larger_NA_0_6.pp.nc
)






# analyses: GLASS-OD ewastools ----


metadata.ewas.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |>  # also pairs with n = 3 & n= 4
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_",ifelse(is.paired,patient_id,"remainder")))) |> 
  dplyr::select(array_sentrix_id, patient, contains("ewastools")) |> 
  assertr::verify(!is.na( array_ewastools_qc_Restoration)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Staining.Green)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Staining.Red)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Extension.Green)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Extension.Red)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Hybridization.High.Medium)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Hybridization.Medium.Low)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Target.Removal.1)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Target.Removal.2)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Bisulfite.Conversion.I.Green)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Bisulfite.Conversion.I.Red)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Bisulfite.Conversion.II)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Specificity.I.Green)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Specificity.I.Red)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Specificity.II)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Non.polymorphic.Green)) |> 
  assertr::verify(!is.na( array_ewastools_qc_Non.polymorphic.Red))


data.ewas.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.ewas.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



## 1. Restoration ----

design.ewas_Restoration.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Restoration, data=metadata.ewas.pp.nc)
fit.ewas_Restoration.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Restoration.pp.nc)
fit.ewas_Restoration.pp.nc <- limma::eBayes(fit.ewas_Restoration.pp.nc, trend=T)
stats.ewas_Restoration.pp.nc <- limma::topTable(fit.ewas_Restoration.pp.nc,
                                                n=nrow(data.ewas.pp.nc),
                                                coef="array_ewastools_qc_Restoration",
                                                sort.by = "none",
                                                adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)

print(paste0("Restoration: ", sum(stats.ewas_Restoration.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Restoration.pp.nc, file="cache/analysis_differential__ewastools_Restoration__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Restoration.pp.nc,
  fit.ewas_Restoration.pp.nc,
  stats.ewas_Restoration.pp.nc
)



## 2. Staining.Green ----

design.ewas_Staining.Green.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Staining.Green, data=metadata.ewas.pp.nc)
fit.ewas_Staining.Green.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Staining.Green.pp.nc)
fit.ewas_Staining.Green.pp.nc <- limma::eBayes(fit.ewas_Staining.Green.pp.nc, trend=T)
stats.ewas_Staining.Green.pp.nc <- limma::topTable(fit.ewas_Staining.Green.pp.nc,
                                                   n=nrow(data.ewas.pp.nc),
                                                   coef="array_ewastools_qc_Staining.Green",
                                                   sort.by = "none",
                                                   adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Staining.Green: ", sum(stats.ewas_Staining.Green.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Staining.Green.pp.nc, file="cache/analysis_differential__ewastools_Staining.Green__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Staining.Green.pp.nc,
  fit.ewas_Staining.Green.pp.nc,
  stats.ewas_Staining.Green.pp.nc
)


## 3. Staining.Red ----

design.ewas_Staining.Red.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Staining.Red, data=metadata.ewas.pp.nc)
fit.ewas_Staining.Red.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Staining.Red.pp.nc)
fit.ewas_Staining.Red.pp.nc <- limma::eBayes(fit.ewas_Staining.Red.pp.nc, trend=T)
stats.ewas_Staining.Red.pp.nc <- limma::topTable(fit.ewas_Staining.Red.pp.nc,
                                                 n=nrow(data.ewas.pp.nc),
                                                 coef="array_ewastools_qc_Staining.Red",
                                                 sort.by = "none",
                                                 adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Staining.Red: ", sum(stats.ewas_Staining.Red.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Staining.Red.pp.nc, file="cache/analysis_differential__ewastools_Staining.Red__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Staining.Red.pp.nc,
  fit.ewas_Staining.Red.pp.nc,
  stats.ewas_Staining.Red.pp.nc
)


## 4. Extension.Green ----

design.ewas_Extension.Green.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Extension.Green, data=metadata.ewas.pp.nc)
fit.ewas_Extension.Green.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Extension.Green.pp.nc)
fit.ewas_Extension.Green.pp.nc <- limma::eBayes(fit.ewas_Extension.Green.pp.nc, trend=T)
stats.ewas_Extension.Green.pp.nc <- limma::topTable(fit.ewas_Extension.Green.pp.nc,
                                                    n=nrow(data.ewas.pp.nc),
                                                    coef="array_ewastools_qc_Extension.Green",
                                                    sort.by = "none",
                                                    adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Extension.Green: ", sum(stats.ewas_Extension.Green.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Extension.Green.pp.nc, file="cache/analysis_differential__ewastools_Extension.Green__partial_paired_nc__stats.Rds")


rm(
  design.ewas_Extension.Green.pp.nc,
  fit.ewas_Extension.Green.pp.nc,
  stats.ewas_Extension.Green.pp.nc
)


## 5. Extension.Red ----

design.ewas_Extension.Red.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Extension.Red, data=metadata.ewas.pp.nc)
fit.ewas_Extension.Red.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Extension.Red.pp.nc)
fit.ewas_Extension.Red.pp.nc <- limma::eBayes(fit.ewas_Extension.Red.pp.nc, trend=T)
stats.ewas_Extension.Red.pp.nc <- limma::topTable(fit.ewas_Extension.Red.pp.nc,
                                                  n=nrow(data.ewas.pp.nc),
                                                  coef="array_ewastools_qc_Extension.Red",
                                                  sort.by = "none",
                                                  adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Extension.Red: ", sum(stats.ewas_Extension.Red.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Extension.Red.pp.nc, file="cache/analysis_differential__ewastools_Extension.Red__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Extension.Red.pp.nc,
  fit.ewas_Extension.Red.pp.nc,
  stats.ewas_Extension.Red.pp.nc
)


## 6. Hybridization.High.Medium ----

design.ewas_Hybridization.High.Medium.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Hybridization.High.Medium, data=metadata.ewas.pp.nc)
fit.ewas_Hybridization.High.Medium.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Hybridization.High.Medium.pp.nc)
fit.ewas_Hybridization.High.Medium.pp.nc <- limma::eBayes(fit.ewas_Hybridization.High.Medium.pp.nc, trend=T)
stats.ewas_Hybridization.High.Medium.pp.nc <- limma::topTable(fit.ewas_Hybridization.High.Medium.pp.nc,
                                                              n=nrow(data.ewas.pp.nc),
                                                              coef="array_ewastools_qc_Hybridization.High.Medium",
                                                              sort.by = "none",
                                                              adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Hybridization.High.Medium: ", sum(stats.ewas_Hybridization.High.Medium.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Hybridization.High.Medium.pp.nc, file="cache/analysis_differential__ewastools_Hybridization.High.Medium__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Hybridization.High.Medium.pp.nc,
  fit.ewas_Hybridization.High.Medium.pp.nc,
  stats.ewas_Hybridization.High.Medium.pp.nc
)


## 7. Hybridization.Medium.Low ----

design.ewas_Hybridization.Medium.Low.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Hybridization.Medium.Low, data=metadata.ewas.pp.nc)
fit.ewas_Hybridization.Medium.Low.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Hybridization.Medium.Low.pp.nc)
fit.ewas_Hybridization.Medium.Low.pp.nc <- limma::eBayes(fit.ewas_Hybridization.Medium.Low.pp.nc, trend=T)
stats.ewas_Hybridization.Medium.Low.pp.nc <- limma::topTable(fit.ewas_Hybridization.Medium.Low.pp.nc,
                                                             n=nrow(data.ewas.pp.nc),
                                                             coef="array_ewastools_qc_Hybridization.Medium.Low",
                                                             sort.by = "none",
                                                             adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Hybridization.Medium.Low: ", sum(stats.ewas_Hybridization.Medium.Low.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Hybridization.Medium.Low.pp.nc, file="cache/analysis_differential__ewastools_Hybridization.Medium.Low__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Hybridization.Medium.Low.pp.nc,
  fit.ewas_Hybridization.Medium.Low.pp.nc,
  stats.ewas_Hybridization.Medium.Low.pp.nc
)



## 8. Target.Removal.1 ----

design.ewas_Target.Removal.1.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Target.Removal.1, data=metadata.ewas.pp.nc)
fit.ewas_Target.Removal.1.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Target.Removal.1.pp.nc)
fit.ewas_Target.Removal.1.pp.nc <- limma::eBayes(fit.ewas_Target.Removal.1.pp.nc, trend=T)
stats.ewas_Target.Removal.1.pp.nc <- limma::topTable(fit.ewas_Target.Removal.1.pp.nc,
                                                     n=nrow(data.ewas.pp.nc),
                                                     coef="array_ewastools_qc_Target.Removal.1",
                                                     sort.by = "none",
                                                     adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Target.Removal.1: ", sum(stats.ewas_Target.Removal.1.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Target.Removal.1.pp.nc, file="cache/analysis_differential__ewastools_Target.Removal.1__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Target.Removal.1.pp.nc,
  fit.ewas_Target.Removal.1.pp.nc,
  stats.ewas_Target.Removal.1.pp.nc
)


## 9. Target.Removal.2 ----

design.ewas_Target.Removal.2.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Target.Removal.2, data=metadata.ewas.pp.nc)
fit.ewas_Target.Removal.2.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Target.Removal.2.pp.nc)
fit.ewas_Target.Removal.2.pp.nc <- limma::eBayes(fit.ewas_Target.Removal.2.pp.nc, trend=T)
stats.ewas_Target.Removal.2.pp.nc <- limma::topTable(fit.ewas_Target.Removal.2.pp.nc,
                                                     n=nrow(data.ewas.pp.nc),
                                                     coef="array_ewastools_qc_Target.Removal.2",
                                                     sort.by = "none",
                                                     adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Target.Removal.2: ", sum(stats.ewas_Target.Removal.2.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Target.Removal.2.pp.nc, file="cache/analysis_differential__ewastools_Target.Removal.2__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Target.Removal.2.pp.nc,
  fit.ewas_Target.Removal.2.pp.nc,
  stats.ewas_Target.Removal.2.pp.nc
)



## 10. Bisulfite.Conversion.I.Green ----

design.ewas_Bisulfite.Conversion.I.Green.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Bisulfite.Conversion.I.Green, data=metadata.ewas.pp.nc)
fit.ewas_Bisulfite.Conversion.I.Green.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Bisulfite.Conversion.I.Green.pp.nc)
fit.ewas_Bisulfite.Conversion.I.Green.pp.nc <- limma::eBayes(fit.ewas_Bisulfite.Conversion.I.Green.pp.nc, trend=T)
stats.ewas_Bisulfite.Conversion.I.Green.pp.nc <- limma::topTable(fit.ewas_Bisulfite.Conversion.I.Green.pp.nc,
                                                                 n=nrow(data.ewas.pp.nc),
                                                                 coef="array_ewastools_qc_Bisulfite.Conversion.I.Green",
                                                                 sort.by = "none",
                                                                 adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Bisulfite.Conversion.I.Green: ", sum(stats.ewas_Bisulfite.Conversion.I.Green.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Bisulfite.Conversion.I.Green.pp.nc, file="cache/analysis_differential__ewastools_Bisulfite.Conversion.I.Green__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Bisulfite.Conversion.I.Green.pp.nc,
  fit.ewas_Bisulfite.Conversion.I.Green.pp.nc,
  stats.ewas_Bisulfite.Conversion.I.Green.pp.nc
)


## 11. Bisulfite.Conversion.I.Red ----

design.ewas_Bisulfite.Conversion.I.Red.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Bisulfite.Conversion.I.Red, data=metadata.ewas.pp.nc)
fit.ewas_Bisulfite.Conversion.I.Red.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Bisulfite.Conversion.I.Red.pp.nc)
fit.ewas_Bisulfite.Conversion.I.Red.pp.nc <- limma::eBayes(fit.ewas_Bisulfite.Conversion.I.Red.pp.nc, trend=T)
stats.ewas_Bisulfite.Conversion.I.Red.pp.nc <- limma::topTable(fit.ewas_Bisulfite.Conversion.I.Red.pp.nc,
                                                               n=nrow(data.ewas.pp.nc),
                                                               coef="array_ewastools_qc_Bisulfite.Conversion.I.Red",
                                                               sort.by = "none",
                                                               adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Bisulfite.Conversion.I.Red: ", sum(stats.ewas_Bisulfite.Conversion.I.Red.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Bisulfite.Conversion.I.Red.pp.nc, file="cache/analysis_differential__ewastools_Bisulfite.Conversion.I.Red__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Bisulfite.Conversion.I.Red.pp.nc,
  fit.ewas_Bisulfite.Conversion.I.Red.pp.nc,
  stats.ewas_Bisulfite.Conversion.I.Red.pp.nc
)


## 12. Bisulfite.Conversion.II ----

design.ewas_Bisulfite.Conversion.II.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Bisulfite.Conversion.II, data=metadata.ewas.pp.nc)
fit.ewas_Bisulfite.Conversion.II.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Bisulfite.Conversion.II.pp.nc)
fit.ewas_Bisulfite.Conversion.II.pp.nc <- limma::eBayes(fit.ewas_Bisulfite.Conversion.II.pp.nc, trend=T)
stats.ewas_Bisulfite.Conversion.II.pp.nc <- limma::topTable(fit.ewas_Bisulfite.Conversion.II.pp.nc,
                                                            n=nrow(data.ewas.pp.nc),
                                                            coef="array_ewastools_qc_Bisulfite.Conversion.II",
                                                            sort.by = "none",
                                                            adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Bisulfite.Conversion.II: ", sum(stats.ewas_Bisulfite.Conversion.II.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Bisulfite.Conversion.II.pp.nc, file="cache/analysis_differential__ewastools_Bisulfite.Conversion.II__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Bisulfite.Conversion.II.pp.nc,
  fit.ewas_Bisulfite.Conversion.II.pp.nc,
  stats.ewas_Bisulfite.Conversion.II.pp.nc
)


## 13. Specificity.I.Green ----

design.ewas_Specificity.I.Green.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Specificity.I.Green, data=metadata.ewas.pp.nc)
fit.ewas_Specificity.I.Green.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Specificity.I.Green.pp.nc)
fit.ewas_Specificity.I.Green.pp.nc <- limma::eBayes(fit.ewas_Specificity.I.Green.pp.nc, trend=T)
stats.ewas_Specificity.I.Green.pp.nc <- limma::topTable(fit.ewas_Specificity.I.Green.pp.nc,
                                                        n=nrow(data.ewas.pp.nc),
                                                        coef="array_ewastools_qc_Specificity.I.Green",
                                                        sort.by = "none",
                                                        adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Specificity.I.Green: ", sum(stats.ewas_Specificity.I.Green.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Specificity.I.Green.pp.nc, file="cache/analysis_differential__ewastools_Specificity.I.Green__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Specificity.I.Green.pp.nc,
  fit.ewas_Specificity.I.Green.pp.nc,
  stats.ewas_Specificity.I.Green.pp.nc
)


## 14. Specificity.I.Red ----

design.ewas_Specificity.I.Red.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Specificity.I.Red, data=metadata.ewas.pp.nc)
fit.ewas_Specificity.I.Red.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Specificity.I.Red.pp.nc)
fit.ewas_Specificity.I.Red.pp.nc <- limma::eBayes(fit.ewas_Specificity.I.Red.pp.nc, trend=T)
stats.ewas_Specificity.I.Red.pp.nc <- limma::topTable(fit.ewas_Specificity.I.Red.pp.nc,
                                                      n=nrow(data.ewas.pp.nc),
                                                      coef="array_ewastools_qc_Specificity.I.Red",
                                                      sort.by = "none",
                                                      adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Specificity.I.Red: ", sum(stats.ewas_Specificity.I.Red.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Specificity.I.Red.pp.nc, file="cache/analysis_differential__ewastools_Specificity.I.Red__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Specificity.I.Red.pp.nc,
  fit.ewas_Specificity.I.Red.pp.nc,
  stats.ewas_Specificity.I.Red.pp.nc
)


## 15. Specificity.II ----

design.ewas_Specificity.II.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Specificity.II, data=metadata.ewas.pp.nc)
fit.ewas_Specificity.II.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Specificity.II.pp.nc)
fit.ewas_Specificity.II.pp.nc <- limma::eBayes(fit.ewas_Specificity.II.pp.nc, trend=T)
stats.ewas_Specificity.II.pp.nc <- limma::topTable(fit.ewas_Specificity.II.pp.nc,
                                                   n=nrow(data.ewas.pp.nc),
                                                   coef="array_ewastools_qc_Specificity.II",
                                                   sort.by = "none",
                                                   adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Specificity.II: ", sum(stats.ewas_Specificity.II.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Specificity.II.pp.nc, file="cache/analysis_differential__ewastools_Specificity.II__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Specificity.II.pp.nc,
  fit.ewas_Specificity.II.pp.nc,
  stats.ewas_Specificity.II.pp.nc
)


## 16. Non.polymorphic.Green ----

design.ewas_Non.polymorphic.Green.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Non.polymorphic.Green, data=metadata.ewas.pp.nc)
fit.ewas_Non.polymorphic.Green.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Non.polymorphic.Green.pp.nc)
fit.ewas_Non.polymorphic.Green.pp.nc <- limma::eBayes(fit.ewas_Non.polymorphic.Green.pp.nc, trend=T)
stats.ewas_Non.polymorphic.Green.pp.nc <- limma::topTable(fit.ewas_Non.polymorphic.Green.pp.nc,
                                                          n=nrow(data.ewas.pp.nc),
                                                          coef="array_ewastools_qc_Non.polymorphic.Green",
                                                          sort.by = "none",
                                                          adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Non.polymorphic.Green: ", sum(stats.ewas_Non.polymorphic.Green.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Non.polymorphic.Green.pp.nc, file="cache/analysis_differential__ewastools_Non.polymorphic.Green__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Non.polymorphic.Green.pp.nc,
  fit.ewas_Non.polymorphic.Green.pp.nc,
  stats.ewas_Non.polymorphic.Green.pp.nc
)


## 17. Non.polymorphic.Red ----

design.ewas_Non.polymorphic.Red.pp.nc <- model.matrix(~factor(patient) + array_ewastools_qc_Non.polymorphic.Red, data=metadata.ewas.pp.nc)
fit.ewas_Non.polymorphic.Red.pp.nc <- limma::lmFit(data.ewas.pp.nc, design.ewas_Non.polymorphic.Red.pp.nc)
fit.ewas_Non.polymorphic.Red.pp.nc <- limma::eBayes(fit.ewas_Non.polymorphic.Red.pp.nc, trend=T)
stats.ewas_Non.polymorphic.Red.pp.nc <- limma::topTable(fit.ewas_Non.polymorphic.Red.pp.nc,
                                                        n=nrow(data.ewas.pp.nc),
                                                        coef="array_ewastools_qc_Non.polymorphic.Red",
                                                        sort.by = "none",
                                                        adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')  |> 
  dplyr::select(probe_id, t)

print(paste0("Non.polymorphic.Red: ", sum(stats.ewas_Non.polymorphic.Red.pp.nc$adj.P.Val < 0.01)))
saveRDS(stats.ewas_Non.polymorphic.Red.pp.nc, file="cache/analysis_differential__ewastools_Non.polymorphic.Red__partial_paired_nc__stats.Rds")

rm(
  design.ewas_Non.polymorphic.Red.pp.nc,
  fit.ewas_Non.polymorphic.Red.pp.nc,
  stats.ewas_Non.polymorphic.Red.pp.nc
)




# analyses: GLASS-OD AcCGAP CGC & GLASSNL sig ----
## data: CGC partially paired [w/o FFPE/frozen batch correct] ----


metadata.AccGAP.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  assertr::verify(is.numeric(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  dplyr::mutate(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit = scale(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |>  # also pairs with n = 3 & n= 4
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_",ifelse(is.paired,patient_id,"remainder")))) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_GLASS_OD_INCLUDED_SAMPLES)) 
    return(.)
  })()



data.AccGAP.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.AccGAP.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()




design.AccGAP.pp.nc <- model.matrix(~factor(patient) + array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, data=metadata.AccGAP.pp.nc)
fit.AccGAP.pp.nc <- limma::lmFit(data.AccGAP.pp.nc, design.AccGAP.pp.nc)
fit.AccGAP.pp.nc <- limma::eBayes(fit.AccGAP.pp.nc, trend=T)
stats.AccGAP.pp.nc <- limma::topTable(fit.AccGAP.pp.nc,
                                      n=nrow(data.AccGAP.pp.nc),
                                      coef="array_A_IDH_HG__A_IDH_LG_lr__lasso_fit",
                                      sort.by = "none",
                                      adjust.method="fdr") |> 
    tibble::rownames_to_column('probe_id')


rm(design.AccGAP.pp.nc)


sum(stats.AccGAP.pp.nc$P.Value < 0.01)
sum(stats.AccGAP.pp.nc$adj.P.Val < 0.01)


#saveRDS(fit.AccGAP.pp.nc, file="cache/analysis_differential__AcCGAP__partial_paired_nc__fit.Rds")
saveRDS(stats.AccGAP.pp.nc, file="cache/analysis_differential__AcCGAP__partial_paired_nc__stats.Rds")


rm(fit.AccGAP.pp.nc, stats.AccGAP.pp.nc)







## data: GLASS-NL sig ----




metadata.GLASS_NL_sig.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  assertr::verify(!is.na(array_GLASS_NL_g2_g3_sig)) |> 
  assertr::verify(is.numeric(array_GLASS_NL_g2_g3_sig)) |> 
  dplyr::mutate(array_GLASS_NL_g2_g3_sig = scale(array_GLASS_NL_g2_g3_sig)) |> 
  #dplyr::mutate( = as.numeric(array_GLASS_NL_g2_g3_sig)) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |>  # also pairs with n = 3 & n= 4
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_",ifelse(is.paired,patient_id,"remainder")))) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_GLASS_OD_INCLUDED_SAMPLES)) 
    return(.)
  })()



data.GLASS_NL_sig.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.GLASS_NL_sig.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()




design.GLASS_NL_sig.pp.nc <- model.matrix(~factor(patient) + array_GLASS_NL_g2_g3_sig, data=metadata.GLASS_NL_sig.pp.nc)
fit.GLASS_NL_sig.pp.nc <- limma::lmFit(data.GLASS_NL_sig.pp.nc, design.GLASS_NL_sig.pp.nc)
fit.GLASS_NL_sig.pp.nc <- limma::eBayes(fit.GLASS_NL_sig.pp.nc, trend=T)
stats.GLASS_NL_sig.pp.nc <- limma::topTable(fit.GLASS_NL_sig.pp.nc,
                                            n=nrow(data.GLASS_NL_sig.pp.nc),
                                            coef="array_GLASS_NL_g2_g3_sig",
                                            sort.by = "none",
                                            adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')


rm(design.GLASS_NL_sig.pp.nc)


sum(stats.GLASS_NL_sig.pp.nc$P.Value < 0.01)
sum(stats.GLASS_NL_sig.pp.nc$adj.P.Val < 0.01)


saveRDS(stats.GLASS_NL_sig.pp.nc, file="cache/analysis_differential__GLASS_NL_sig__partial_paired_nc__stats.Rds")


rm(fit.GLASS_NL_sig.pp.nc, stats.GLASS_NL_sig.pp.nc)





## x plots: test pat + condition  ----


design.lgc <- model.matrix(~factor(patient) + as.numeric(A_IDH_HG__A_IDH_LG_lr__lasso_fit), data=metadata.pp)
fit.lgc <- limma::lmFit(data.pp, design.lgc)
fit.lgc <- limma::eBayes(fit.lgc, trend=T)
stats.lgc <- limma::topTable(fit.lgc,
                            n=nrow(data.pp),
                            coef="as.numeric(A_IDH_HG__A_IDH_LG_lr__lasso_fit)",
                            sort.by = "none",
                            adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')

rm(design.lgc, fit.lgc)


sum(stats.lgc$P.Value < 0.01)
sum(stats.lgc$adj.P.Val < 0.01)

plot(sort(stats.lgc$P.Value),type="l")



# analyses: FFPE & FF ----
## data: partially paired [w/o FFPE/frozen batch correct] ----


metadata.ffpe_or_ff.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(!is.na(isolation_material)) |>
  dplyr::mutate(isolation_material = factor(isolation_material, levels=c("ffpe", "tissue"))) |> 
  dplyr::mutate(ffpe_or_ff_time = ifelse(isolation_material == "ffpe", -time_between_resection_and_array, 0)) |> 
  dplyr::filter(!is.na(ffpe_or_ff_time)) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_",ifelse(is.paired, patient_id, "remainder")))) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 205)
    return(.)
  })()


data.ffpe_or_ff.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.ffpe_or_ff.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



design.ffpe_or_ff.pp.nc <- model.matrix(~factor(patient) + isolation_material, data=metadata.ffpe_or_ff.pp.nc)
fit.ffpe_or_ff.pp.nc <- limma::lmFit(data.ffpe_or_ff.pp.nc, design.ffpe_or_ff.pp.nc)
fit.ffpe_or_ff.pp.nc <- limma::eBayes(fit.ffpe_or_ff.pp.nc, trend=T)

stats.ffpe_or_ff.pp.nc <- limma::topTable(fit.ffpe_or_ff.pp.nc,
                                          n=nrow(data.ffpe_or_ff.pp.nc),
                                          coef="isolation_materialtissue",
                                          sort.by = "none",
                                          adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')



rm(design.ffpe_or_ff.pp.nc)


sum(stats.ffpe_or_ff.pp.nc$P.Value < 0.01)
sum(stats.ffpe_or_ff.pp.nc$adj.P.Val < 0.01)



#saveRDS(fit.ffpe_or_ff.pp.nc, file="cache/analysis_differential__ffpe_or_ff__partial_paired_nc__fit.Rds")
saveRDS(stats.ffpe_or_ff.pp.nc, file="cache/analysis_differential__ffpe_or_ff__partial_paired_nc__stats.Rds")


rm(fit.ffpe_or_ff.pp.nc, stats.ffpe_or_ff.pp.nc)





# analyses: FFPE-decay predictor ----
## data: partially paired [w/o FFPE/frozen batch correct] ----
### 1. time in FFPE ----


metadata.ffpe_decay.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(!is.na(isolation_material)) |> 
  #dplyr::filter(!is.na(time_between_resection_and_array)) |> # may remove NA's which are non-ffpe
  dplyr::mutate(ffpe_decay_time = ifelse(isolation_material == "ffpe", -time_between_resection_and_array, 0)) |> 
  dplyr::filter(!is.na(ffpe_decay_time)) |> 
  dplyr::mutate(ffpe_decay_time = scale(ffpe_decay_time, center=F)) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_",ifelse(is.paired, patient_id, "remainder")))) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 205) 
    return(.)
  })()


data.ffpe_decay.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.ffpe_decay.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



design.ffpe_decay.pp.nc <- model.matrix(~factor(patient) + ffpe_decay_time, data=metadata.ffpe_decay.pp.nc)
fit.ffpe_decay.pp.nc <- limma::lmFit(data.ffpe_decay.pp.nc, design.ffpe_decay.pp.nc)
fit.ffpe_decay.pp.nc <- limma::eBayes(fit.ffpe_decay.pp.nc, trend=T)

stats.ffpe_decay.pp.nc <- limma::topTable(fit.ffpe_decay.pp.nc,
                                          n=nrow(data.ffpe_decay.pp.nc),
                                          coef="ffpe_decay_time",
                                          sort.by = "none",
                                          adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')



rm(design.ffpe_decay.pp.nc)


sum(stats.ffpe_decay.pp.nc$P.Value < 0.01)
sum(stats.ffpe_decay.pp.nc$adj.P.Val < 0.01)



#saveRDS(fit.ffpe_decay.pp.nc, file="cache/analysis_differential__ffpe-decay-time__partial_paired_nc__fit.Rds")
saveRDS(stats.ffpe_decay.pp.nc, file="cache/analysis_differential__ffpe-decay-time__partial_paired_nc__stats.Rds")


rm(fit.ffpe_decay.pp.nc, stats.ffpe_decay.pp.nc, metadata.ffpe_decay.pp.nc)




### 2. time in freezer ----


metadata.freezer_decay.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(!is.na(isolation_material)) |> 
  dplyr::mutate(freezer_decay_time = ifelse(isolation_material == "tissue", -time_between_resection_and_array, 0)) |> 
  dplyr::filter(!is.na(freezer_decay_time)) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_",ifelse(is.paired, patient_id, "remainder")))) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 203) 
    return(.)
  })()



data.freezer_decay.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.freezer_decay.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



design.freezer_decay.pp.nc <- model.matrix(~factor(patient) + freezer_decay_time, data=metadata.freezer_decay.pp.nc)
fit.freezer_decay.pp.nc <- limma::lmFit(data.freezer_decay.pp.nc, design.freezer_decay.pp.nc)
fit.freezer_decay.pp.nc <- limma::eBayes(fit.freezer_decay.pp.nc, trend=T)

stats.freezer_decay.pp.nc <- limma::topTable(fit.freezer_decay.pp.nc,
                                             n=nrow(data.freezer_decay.pp.nc),
                                             coef="freezer_decay_time",
                                             sort.by = "none",
                                             adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t)



rm(design.freezer_decay.pp.nc)


sum(stats.freezer_decay.pp.nc$P.Value < 0.01)
sum(stats.freezer_decay.pp.nc$adj.P.Val < 0.01)



#saveRDS(fit.freezer_decay.pp.nc, file="cache/analysis_differential__freezer-decay-time__partial_paired_nc__fit.Rds")
saveRDS(stats.freezer_decay.pp.nc, file="cache/analysis_differential__freezer-decay-time__partial_paired_nc__stats.Rds")


rm(fit.freezer_decay.pp.nc,
   stats.freezer_decay.pp.nc,
   data.freezer_decay.pp.nc,
   metadata.freezer_decay.pp.nc)







### 3. time in freezer & ffpe multivariate ----


metadata.ffpe_freezer_multivar.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(!is.na(isolation_material)) |> 
  
  dplyr::mutate(ffpe_decay_time = ifelse(isolation_material == "ffpe", -time_between_resection_and_array, 0)) |> 
  #dplyr::mutate(ffpe_decay_time = -log((-ffpe_decay_time) + 1)) |> 
  
  dplyr::mutate(freezer_decay_time = ifelse(isolation_material == "tissue", -time_between_resection_and_array, 0)) |> 
  #dplyr::mutate(freezer_decay_time = -log((-freezer_decay_time) + 1)) |> 
  
  dplyr::filter(!is.na(ffpe_decay_time) & !is.na(freezer_decay_time)) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_",ifelse(is.paired, patient_id, "remainder")))) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 203) 
    return(.)
  })()



data.ffpe_freezer_multivar.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.ffpe_freezer_multivar.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()




design.ffpe_freezer_multivar.pp.nc <- model.matrix(~factor(patient) + 
                                                     ffpe_decay_time +
                                                     freezer_decay_time
                                                     , data=metadata.ffpe_freezer_multivar.pp.nc)
fit.ffpe_freezer_multivar.pp.nc <- limma::lmFit(data.ffpe_freezer_multivar.pp.nc, design.ffpe_freezer_multivar.pp.nc)
fit.ffpe_freezer_multivar.pp.nc <- limma::eBayes(fit.ffpe_freezer_multivar.pp.nc, trend=T)

stats.ffpe_freezer_multivar__ffpe.pp.nc <- limma::topTable(fit.ffpe_freezer_multivar.pp.nc,
                                                           n=nrow(data.ffpe_freezer_multivar.pp.nc),
                                                           coef="ffpe_decay_time",
                                                           sort.by = "none",
                                                           adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t) |> 
  dplyr::rename(t__ffpe = t)


sum(abs(stats.ffpe_freezer_multivar__ffpe.pp.nc$t) > 1.5)

stats.ffpe_freezer_multivar__freezer.pp.nc <- limma::topTable(fit.ffpe_freezer_multivar.pp.nc,
                                                           n=nrow(data.ffpe_freezer_multivar.pp.nc),
                                                           coef="freezer_decay_time",
                                                           sort.by = "none",
                                                           adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::select(probe_id, t) |> 
  dplyr::rename(t__freezer = t)


stats.ffpe_freezer_multivar.pp.nc <- stats.ffpe_freezer_multivar__ffpe.pp.nc |> 
  dplyr::left_join(stats.ffpe_freezer_multivar__freezer.pp.nc, by=c('probe_id'='probe_id'))



rm(design.ffpe_freezer_multivar.pp.nc)




saveRDS(stats.ffpe_freezer_multivar.pp.nc, file="cache/analysis_differential__ffpe_freezer_decay-time_multivar__partial_paired_nc__stats.Rds")



rm(fit.ffpe_freezer_multivar.pp.nc,
   stats.ffpe_freezer_multivar__ffpe.pp.nc,
   stats.ffpe_freezer_multivar__freezer.pp.nc,
   stats.ffpe_freezer_multivar.pp.nc,
   data.ffpe_freezer_multivar.pp.nc,
   metadata.ffpe_freezer_multivar.pp.nc
   )




### lm per gene outcome x CG or C ----


if(!exists('data.mvalues.probes')) {
  source('scripts/load_mvalues_hq_samples.R')
}

#tmp <- readRDS(file="cache/analysis_differential__ffpe-decay-time__partial_paired_nc__stats.Rds")


tmp <- data.mvalues.probes |> 
  dplyr::filter(!is.na(DMP__FFPE_decay_time__pp_nc__t)) |>
  dplyr::filter(!is.na(n_CG)) |> 
  dplyr::select(probe_id, DMP__FFPE_decay_time__pp_nc__t,  n_CG, n_CG_to_CR, n_CG_to_CA, n_independent_C) |> 
  dplyr::mutate(total_c = n_CG + n_independent_C)  |> 
  dplyr::mutate(indep_c_frac = n_independent_C / (50 - (2 * n_CG)))



stat <- stats::glm(tmp$DMP__FFPE_decay_time__pp_nc__t ~ tmp$n_CG + tmp$n_independent_C)
summary(stat)


stat <- stats::glm(tmp$DMP__FFPE_decay_time__pp_nc__t ~ tmp$n_CG + tmp$indep_c_frac)
summary(stat)



stat <- lm(tmp$DMP__FFPE_decay_time__pp_nc__t ~ tmp$n_independent_C )
summary(stat)


stat <- lm(tmp$DMP__FFPE_decay_time__pp_nc__t ~ tmp$indep_c_frac )
summary(stat)



stat <- lm(tmp$DMP__FFPE_decay_time__pp_nc__t ~  tmp$total_c)
summary(stat)


stat <- lm(tmp$DMP__FFPE_decay_time__pp_nc__t ~  tmp$n_CG )
summary(stat)

stat <- lm(tmp$DMP__FFPE_decay_time__pp_nc__t ~  tmp$n_CG_to_CR )
summary(stat)

stat <- lm(tmp$DMP__FFPE_decay_time__pp_nc__t ~  tmp$n_CG_to_CA )
summary(stat)



ggplot(tmp, aes(y=DMP__FFPE_decay_time__pp_nc__t, x=factor(n_CG))) + 
  geom_boxplot(width=0.75, outlier.shape=NA, outlier.color=NA, col="darkgray")


ggplot(tmp, aes(y=DMP__FFPE_decay_time__pp_nc__t, x=factor(n_independent_C))) + 
  geom_boxplot(width=0.75, outlier.shape=NA, outlier.color=NA, col="darkgray")


ggplot(tmp, aes(y=DMP__FFPE_decay_time__pp_nc__t, x=indep_c_frac)) + 
  geom_point(cex=0.1)
  #geom_boxplot(width=0.75, outlier.shape=NA, outlier.color=NA, col="darkgray")





## data: COMBINED INTENSITIES ----


metadata.ffpe_decay.intensities.total <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  
  dplyr::filter(!is.na(isolation_material)) |> 
  #dplyr::filter(!is.na(time_between_resection_and_array)) |> 
  dplyr::mutate(in_ffpe = isolation_material == "ffpe") |> 
  dplyr::mutate(ffpe_decay_time = ifelse(in_ffpe, time_between_resection_and_array, 0)) |> # don't change sign here (!!)
  #dplyr::filter(in_ffpe) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_",ifelse(is.paired,patient_id,"remainder")))) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 205) 
    return(.)
  })() 
#dplyr::mutate(ffpe_decay_time = scale(ffpe_decay_time, center=F, scale=T))


data.ffpe_decay.intensities.total <- data.intensities.combined.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.intensities.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.ffpe_decay.intensities.total$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()
#t() |> 
#scale(center = TRUE) |> 
#t()



design.ffpe_decay.intensities.total <- model.matrix(~factor(patient) + 
                                                      ffpe_decay_time, data=metadata.ffpe_decay.intensities.total)
fit.ffpe_decay.intensities.total <- limma::lmFit(data.ffpe_decay.intensities.total, design.ffpe_decay.intensities.total)
fit.ffpe_decay.intensities.total <- limma::eBayes(fit.ffpe_decay.intensities.total, trend=T)
stats.ffpe_decay.intensities.total <- limma::topTable(fit.ffpe_decay.intensities.total,
                                                      n=nrow(data.ffpe_decay.intensities.total),
                                                      coef="ffpe_decay_time",
                                                      sort.by = "none",
                                                      adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')



#plot(sort(metadata.ffpe_decay.intensities.total$ffpe_decay_time ), pch=16, type="l")
plot(sort(stats.ffpe_decay.intensities.total$t), pch=16, type="l")
abline(h=0, col="red")
#abline(h= -3.258864, col="red")



rm(design.ffpe_decay.intensities.total)


sum(stats.ffpe_decay.intensities.total$P.Value < 0.01)
sum(stats.ffpe_decay.intensities.total$adj.P.Val < 0.01)



saveRDS(stats.ffpe_decay.intensities.total, file="cache/analysis_differential_intensities__ffpe-decay-time__partial_paired_nc__stats.Rds")



rm(fit.ffpe_decay.intensities.total, stats.ffpe_decay.intensities.total, metadata.ffpe_decay.intensities.total)



# top3

# stats.ffpe_decay.pp.nc |> 
#   dplyr::arrange(P.Value) |> 
#   head(n=650000) |> 
#   tail(n=3)
# 
# 
# plt <- data.ffpe_decay.pp.nc |> 
#   as.data.frame() |> 
#   dplyr::select(metadata.ffpe_decay.pp.nc$array_sentrix_id) |> 
#   tibble::rownames_to_column("probe_id") |> 
#   dplyr::filter(probe_id  %in% c("cg16405174")) |> 
#   tibble::column_to_rownames("probe_id") |> 
#   t() |> 
#   as.data.frame() |> 
#   tibble::rownames_to_column("array_sentrix_id") |> 
#   dplyr::left_join(metadata.ffpe_decay.pp.nc, by=c('array_sentrix_id'='array_sentrix_id'))
# 
# 
# ggplot(plt, aes(x=ffpe_decay_time, y=cg16405174, label=resection_id, col=in_ffpe)) +  
#   ggrepel::geom_text_repel(size=2) +
#   geom_point(cex=0.5) +
#   ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho")
# 
# 
# 
# ggplot(data = cats, aes(x = life_years, y = happy_years, color = cat_type)) +
#   # facet_wrap(~ cat_type) +
#   stat_poly_line(formula = y ~ x) +
#   stat_poly_eq(formula = y ~ x, aes(label = after_stat(f.value))) +
#   #stat_poly_eq(formula = y ~ x, aes(label = after_stat(coef.ls))) +
#   geom_point()


# https://www.rdocumentation.org/packages/ggpmisc/versions/0.6.0/topics/stat_poly_eq

## data: METHYLATED INTENSITIES ----



metadata.ffpe_decay.intensities.methylated <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(!is.na(isolation_material)) |> 
  dplyr::mutate(ffpe_decay_time = ifelse(isolation_material == "ffpe", time_between_resection_and_array, 0)) |> # don't change sign here (!!)
  dplyr::filter(!is.na(ffpe_decay_time)) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_", ifelse(is.paired,patient_id, "remainder")))) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 205) 
    return(.)
  })()


data.ffpe_decay.intensities.methylated <- data.intensities.methylated.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.intensities.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.ffpe_decay.intensities.methylated$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



design.ffpe_decay.intensities.methylated <- model.matrix(~factor(patient) + ffpe_decay_time, data=metadata.ffpe_decay.intensities.methylated)
fit.ffpe_decay.intensities.methylated <- limma::lmFit(data.ffpe_decay.intensities.methylated, design.ffpe_decay.intensities.methylated)
fit.ffpe_decay.intensities.methylated <- limma::eBayes(fit.ffpe_decay.intensities.methylated, trend=T)
stats.ffpe_decay.intensities.methylated <- limma::topTable(fit.ffpe_decay.intensities.methylated,
                                                           n=nrow(data.ffpe_decay.intensities.methylated),
                                                           coef="ffpe_decay_time",
                                                           sort.by = "none",
                                                           adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')



rm(design.ffpe_decay.intensities.methylated)


plot(sort(stats.ffpe_decay.intensities.methylated$t), pch=16, type="l")
abline(h=0, col="red")


sum(stats.ffpe_decay.intensities.methylated$P.Value < 0.01)
sum(stats.ffpe_decay.intensities.methylated$adj.P.Val < 0.01)



saveRDS(stats.ffpe_decay.intensities.methylated, file="cache/analysis_differential_methylated_intensities__ffpe-decay-time__partial_paired_nc__stats.Rds")


rm(fit.ffpe_decay.intensities.methylated, stats.ffpe_decay.intensities.methylated, metadata.ffpe_decay.intensities.methylated)






## data: UNMETHYLATED INTENSITIES ----


metadata.ffpe_decay.intensities.unmethylated <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(!is.na(isolation_material)) |> 
  dplyr::mutate(ffpe_decay_time = ifelse(isolation_material == "ffpe", time_between_resection_and_array, 0)) |> # don't change sign here (!!)
  dplyr::filter(!is.na(ffpe_decay_time)) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_", ifelse(is.paired,patient_id, "remainder")))) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 205) 
    return(.)
  })()


data.ffpe_decay.intensities.unmethylated <- data.intensities.unmethylated.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.intensities.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.ffpe_decay.intensities.unmethylated$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



design.ffpe_decay.intensities.unmethylated <- model.matrix(~factor(patient) + ffpe_decay_time, data=metadata.ffpe_decay.intensities.unmethylated)
fit.ffpe_decay.intensities.unmethylated <- limma::lmFit(data.ffpe_decay.intensities.unmethylated, design.ffpe_decay.intensities.unmethylated)
fit.ffpe_decay.intensities.unmethylated <- limma::eBayes(fit.ffpe_decay.intensities.unmethylated, trend=T)
stats.ffpe_decay.intensities.unmethylated <- limma::topTable(fit.ffpe_decay.intensities.unmethylated,
                                                             n=nrow(data.ffpe_decay.intensities.unmethylated),
                                                             coef="ffpe_decay_time",
                                                             sort.by = "none",
                                                             adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')



rm(design.ffpe_decay.intensities.unmethylated)


sum(stats.ffpe_decay.intensities.unmethylated$P.Value < 0.01)
sum(stats.ffpe_decay.intensities.unmethylated$adj.P.Val < 0.01)


saveRDS(stats.ffpe_decay.intensities.unmethylated, file="cache/analysis_differential_unmethylated_intensities__ffpe-decay-time__partial_paired_nc__stats.Rds")


rm(fit.ffpe_decay.intensities.unmethylated, stats.ffpe_decay.intensities.unmethylated, metadata.ffpe_decay.intensities.unmethylated)





## data: unpaired [w/o FFPE/frozen batch correct] ----


metadata.ffpe_decay.up.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(!is.na(isolation_material)) |> 
  dplyr::mutate(ffpe_decay_time = ifelse(isolation_material == "ffpe", -time_between_resection_and_array, 0)) |> 
  dplyr::filter(!is.na(ffpe_decay_time)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 205) 
    return(.)
  })()


data.ffpe_decay.up.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.ffpe_decay.up.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


design.ffpe_decay.up.nc <- model.matrix(~ffpe_decay_time, data=metadata.ffpe_decay.up.nc)
fit.ffpe_decay.up.nc <- limma::lmFit(data.ffpe_decay.up.nc, design.ffpe_decay.up.nc)
fit.ffpe_decay.up.nc <- limma::eBayes(fit.ffpe_decay.up.nc, trend=T)

stats.ffpe_decay.up.nc <- limma::topTable(fit.ffpe_decay.up.nc,
                                          n=nrow(data.ffpe_decay.up.nc),
                                          coef="ffpe_decay_time",
                                          sort.by = "none",
                                          adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id')



rm(design.ffpe_decay.up.nc)


sum(stats.ffpe_decay.up.nc$P.Value < 0.01)
sum(stats.ffpe_decay.up.nc$adj.P.Val < 0.01)


saveRDS(stats.ffpe_decay.up.nc, file="cache/analysis_differential__ffpe-decay-time__unpaired_nc__stats.Rds")


rm(fit.ffpe_decay.up.nc, stats.ffpe_decay.up.nc, metadata.ffpe_decay.up.nc)





# analyses: epiGenetic clocks ----
## data: unpaired [w/o FFPE/frozen batch correct] ----


clocks <- glass_od.metadata.array_samples |> 
  dplyr::select(contains("epiTOC") | contains("dnaMethyAge_") | contains("_RepliTali")) |> 
  colnames()


for(clock in clocks) {
  print(clock)
  
  
  metadata.current_clock.up.nc <- glass_od.metadata.array_samples |> 
    filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
    data.table::copy() |> # odd hack needed, because setnames also affects the former "glass_od.metadata.array_samples" object...
    data.table::setnames(old = c(clock), new = c("array_current_clock")) |> 
    dplyr::filter(!is.na(array_current_clock)) |> 
    
    dplyr::mutate(array_current_clock = scale(array_current_clock)) |> 
    
    dplyr::group_by(patient_id) |> 
    dplyr::mutate(is.paired = dplyr::n() >= 2) |> 
    dplyr::ungroup() |> 
    dplyr::mutate(patient = as.factor(paste0("p_", ifelse(is.paired,patient_id, "remainder")))) |> 
    
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_INCLUDED_SAMPLES) 
      return(.)
    })()
  
  
  data.current_clock.up.nc <- data.mvalues.hq_samples |> 
    tibble::rownames_to_column('probe_id') |> 
    dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
    tibble::column_to_rownames('probe_id') |> 
    
    dplyr::select(metadata.current_clock.up.nc$array_sentrix_id) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
      return(.)
    })()
  
  
  
  
  design.current_clock.up.nc <- model.matrix(~factor(patient) + array_current_clock, data=metadata.current_clock.up.nc)
  fit.current_clock.up.nc <- limma::lmFit(data.current_clock.up.nc, design.current_clock.up.nc)
  fit.current_clock.up.nc <- limma::eBayes(fit.current_clock.up.nc, trend=T)
  
  stats.current_clock.up.nc <- limma::topTable(fit.current_clock.up.nc,
                                               n=nrow(data.current_clock.up.nc),
                                               coef="array_current_clock",
                                               sort.by = "none",
                                               adjust.method="fdr") |> 
    tibble::rownames_to_column('probe_id')
  
  
  
  rm(design.current_clock.up.nc)
  
  
  #sum(stats.current_clock.up.nc$P.Value < 0.01)
  #sum(stats.current_clock.up.nc$adj.P.Val < 0.01)
  
  
  
  saveRDS(stats.current_clock.up.nc, file=paste0("cache/analysis_differential__",gsub("array_","", clock),"__unpaired_nc__stats.Rds"))
  
  
  rm(fit.current_clock.up.nc,
     stats.current_clock.up.nc,
     clock,
     metadata.current_clock.up.nc)
}






# analyses: GLASS-OD A_IDH_HG - OLIGOSARC ----


metadata.hg_olsc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(211) |> 
  dplyr::filter(array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG", "OLIGOSARC_IDH")) |> 
  dplyr::filter(isolation_id != "0104-R2")   # actual astro


p1 = ggplot(metadata.hg_olsc, aes(x=array_mnp_predictBrain_v12.8_cal_class, y=array_methylation_bins_1p19q_purity, label=isolation_id)) +
  ggbeeswarm::geom_quasirandom(size=theme_cellpress_size/2) +
  theme_cellpress


p2 = ggplot(metadata.hg_olsc, aes(x=array_mnp_predictBrain_v12.8_cal_class, y=array_median.overall.methylation, label=isolation_id)) +
  ggbeeswarm::geom_quasirandom(size=theme_cellpress_size/2) +
  theme_cellpress


p3 = ggplot(metadata.hg_olsc, aes(x=array_mnp_predictBrain_v12.8_cal_class, y=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, label=isolation_id)) +
  ggbeeswarm::geom_quasirandom(size=theme_cellpress_size/2) +
  theme_cellpress


p4 = ggplot(metadata.hg_olsc, aes(x=array_mnp_predictBrain_v12.8_cal_class, y=array_qc.pca.comp1, label=isolation_id)) +
  ggbeeswarm::geom_quasirandom(size=theme_cellpress_size/2) +
  theme_cellpress


p1 + p2 + p3 + p4





data.hg_olsc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.hg_olsc$array_sentrix_id) 


design.hg_olsc <- model.matrix(~  factor(array_mnp_predictBrain_v12.8_cal_class), data=metadata.hg_olsc)

fit.hg_olsc <- limma::lmFit(data.hg_olsc, design.hg_olsc)
fit.hg_olsc <- limma::eBayes(fit.hg_olsc, trend=T)
stats.hg_olsc <- limma::topTable(fit.hg_olsc,
                            n=nrow(data.hg_olsc),
                            coef="factor(array_mnp_predictBrain_v12.8_cal_class)OLIGOSARC_IDH",
                            sort.by = "none",
                            adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.hg_olsc, fit.hg_olsc)


paste0("DMPs: ",sum(stats.hg_olsc$adj.P.Val < 0.01), " / ", nrow(stats.hg_olsc), " (padj < 0.01)")





#sum(stats.hg_olsc$P.Value < 0.01)
plot(sort(stats.hg_olsc$P.Value),type="l")



hg_olsc.borderline.probes <- stats.hg_olsc |> 
  dplyr::filter(P.Value < 0.01) |> 
  dplyr::left_join(
    data.mvalues.probes, by=c('probe_id'='probe_id'),suffix=c('','')
  ) |> 
  dplyr::mutate(x = ((CpG_beg + CpG_end) / 2 ) / 100000000) |> 
  dplyr::mutate(chr = factor(CpG_chrm, levels=gtools::mixedsort(unique(as.character(CpG_chrm))) ))


plt <- rbind(
  hg_olsc.borderline.probes |> 
    dplyr::mutate(col = logFC),
  hg_olsc.borderline.probes |> 
    dplyr::mutate(col = logFC) |> 
    dplyr::mutate(logFC = 0)
) 




ggplot(plt, aes(x = x, y= logFC, group=probe_id, col=col)) +
  facet_grid(cols = vars(chr), scales = "free", space="free") +
  geom_line(lwd=0.22, alpha=0.3) +
  theme_cellpress +
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(-3, 3), breaks=c(-3, 3), oob = scales::squish)


ggplot(plt |> dplyr::filter(chr == "chr2"), aes(x = x, y= logFC, group=probe_id, col=col)) +
  facet_grid(cols = vars(chr), scales = "free", space="free") +
  geom_line(lwd=0.22, alpha=0.3) +
  theme_cellpress +
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(-3, 3), breaks=c(-3, 3), oob = scales::squish)


ggplot(plt |> dplyr::filter(chr == "chr3"), aes(x = x, y= logFC, group=probe_id, col=col)) +
  facet_grid(cols = vars(chr), scales = "free", space="free") +
  geom_line(lwd=0.22, alpha=0.3) +
  theme_cellpress +
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(-3, 3), breaks=c(-3, 3), oob = scales::squish)



ggplot(plt |> dplyr::filter(chr == "chr4"), aes(x = x, y= logFC, group=probe_id, col=col)) +
  facet_grid(cols = vars(chr), scales = "free", space="free") +
  geom_line(lwd=0.22, alpha=0.3) +
  theme_cellpress +
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(-3, 3), breaks=c(-3, 3), oob = scales::squish)



ggplot(plt |> dplyr::filter(chr == "chr6"), aes(x = x, y= logFC, group=probe_id, col=col)) +
  facet_grid(cols = vars(chr), scales = "free", space="free") +
  geom_line(lwd=0.22, alpha=0.3) +
  theme_cellpress +
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(-3, 3), breaks=c(-3, 3), oob = scales::squish)






hg_olsc.borderline.probes |>
  dplyr::filter(chr == "chr5") |> 
  dplyr::arrange(-P.Value) |>  
  dplyr::select(probe_id,  logFC ,   AveExpr   ,      t,     P.Value, adj.P.Val, CpG_beg ,genesUniq) |> 
  View()




# plot: g2-g3 x lgc ----



plt <- stats.gr |> 
  dplyr::left_join(stats.lgc,by=c('probe_id'='probe_id'),suffx=c('.gr','.lgc'))


plot(sort(stats.lgc$P.Value),type="l")


## powerplot ----
#' take intersect, then order


plt.pre <- stats.gr |> 
  dplyr::left_join(stats.lgc,by=c('probe_id'='probe_id'),suffix=c('.gr','.lgc')) |> 
  dplyr::mutate(col = ifelse(
    probe_id %in% (data.mvalues.probes |> dplyr::filter(catnon_embryionic_development) |>  dplyr::pull(probe_id)),
    "embrionic development", "other"
  ))

plt <- data.frame(
  p.gr = plt.pre |> dplyr::arrange(P.Value.gr) |> dplyr::pull(P.Value.gr),
  p.lgc = plt.pre |> dplyr::arrange(P.Value.lgc) |> dplyr::pull(P.Value.lgc)
) |> 
  dplyr::mutate(x = 1:dplyr::n()) |> 
  dplyr::mutate(delta = p.gr - p.lgc)


#plot(plt$x, plt$delta, type="l")
plot(plt$x, plt$p.gr, type="l")
lines(plt$x, plt$p.lgc)


plot(plt.pre$logFC.gr, plt.pre$logFC.lgc, pch=19,cex=0.1)
abline(h=0, col="red")
abline(v=0, col="red")


ggplot(plt.pre, aes(x=logFC.lgc, y=-log(P.Value.lgc), col=col)) + 
  geom_point(data = plt.pre |> dplyr::filter(col =="other"), pch=19, cex=0.1, alpha=0.1) + 
  geom_point(data = plt.pre |> dplyr::filter(col !="other"), pch=19, cex=0.5,alpha=0.8) +
  geom_vline(xintercept=0, col="red", lty=2, lwd=0.5) +
  geom_hline(yintercept=0, col="red", lty=2, lwd=0.5) +
  theme_bw()



# analyses: GLASS-OD + GLASS-NL LG/HG ----
#' multi-dataset mixed model to find changes specific to
#' OD
#' AC
#' general grading
## test pat + HG_AC + HG_OD + overall_HG ----


metadata.od <- glass_od.metadata.array_samples |> 
    filter_GLASS_OD_idats(163) |> 
    filter_first_G2_and_last_G3(105) |> 
    dplyr::select(sentrix_id, resection_id, patient_id, LG_HG_status, A_IDH_HG__A_IDH_LG_lr__lasso_fit) |>  # has to be WHO since OD has no methylation based split
    dplyr::rename(sample_id = resection_id) |> 
    dplyr::mutate(dataset = "GLASS-OD") |> 
    dplyr::mutate(LG_HG_status = factor(LG_HG_status, levels=c("LG","HG")))
metadata.ac <- glass_nl.metadata.array_samples |> 
      filter_GLASS_NL_idats(218) |> 
      filter_first_G2_and_last_G3(130) |> 
      dplyr::select(sentrix_id, Sample_Name, patient_id, LG_HG_status, A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV) |>  # has to be WHO since OD has no methylation based split
      dplyr::rename(sample_id = Sample_Name) |> 
      dplyr::rename(`A_IDH_HG__A_IDH_LG_lr__lasso_fit` = `A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV`) |> 
      dplyr::mutate(dataset = "GLASS-NL") |> 
      dplyr::mutate(LG_HG_status = factor(LG_HG_status, levels=c("LG","HG")))

data.od <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.od$sentrix_id)

data.ac <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.ac$sentrix_id)


design.od <- model.matrix(~factor(patient_id) + factor(LG_HG_status), data=metadata.od)
design.ac <- model.matrix(~factor(patient_id) + factor(LG_HG_status), data=metadata.ac)

fit.od <- limma::eBayes(limma::lmFit(data.od, design.od),trend=T)
fit.ac <- limma::eBayes(limma::lmFit(data.ac, design.ac),trend=T)

stats.od <- limma::topTable(fit.od,
                              n=nrow(data.od),
                              coef="factor(LG_HG_status)HG",
                              sort.by = "none", 
                              adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')
stats.ac <- limma::topTable(fit.ac,
                              n=nrow(data.ac),
                              coef="factor(LG_HG_status)HG",
                              sort.by = "none",
                              adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')


stopifnot(stats.od$probe_id == stats.ac$probe_id)

plt <- stats.od |> 
  dplyr::left_join(stats.ac, by=c('probe_id'='probe_id'),suffix=c('.od','.ac')) |> 
  dplyr::mutate(col1 = probe_id %in% (stats.b.hg_od |> dplyr::filter(adj.P.Val < 0.01) |> dplyr::pull(probe_id)))|> 
  dplyr::mutate(col2 = probe_id %in% (stats.b.hg_ac |> dplyr::filter(adj.P.Val < 0.01) |> dplyr::pull(probe_id)))


plt <- plt |> 
  dplyr::mutate(col3 = probe_id %in% (stats.b.hg_ac$probe_id[stats.b.hg_ac$adj.P.Val < 0.01]))

ggplot(plt, aes(x=`t.od`,y=`t.ac`, col=col1)) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  geom_point(pch=19, alpha=0.15,cex=0.02) +
  theme_bw()


ggplot(plt, aes(x=`t.od`,y=`t.ac`, col=col2)) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  geom_point(pch=19, alpha=0.15,cex=0.02) +
  theme_bw()


ggplot(plt, aes(x=`t.od`,y=`t.ac`, col=col3)) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  geom_point(data=subset(plt, col3==F),pch=19, alpha=0.15,cex=0.02) +
  geom_point(data=subset(plt, col3==T),pch=19, alpha=0.65,cex=0.05) +
  theme_bw()


ggplot(plt, aes(x=`t.od`,y=`t.ac`, col=col2)) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  geom_point(data=subset(plt, col2==F),pch=19, alpha=0.15,cex=0.02) +
  geom_point(data=subset(plt, col2==T),pch=19, alpha=0.65,cex=0.05) +
  theme_bw()



# AD dataset ----


metadata.ad <- ad_bmc_clin_epi.metadata.array_samples |> 
  dplyr::mutate(ad.status = factor(ifelse(Type == "AD control", "Control","AlzheimerDisease"), levels=c( "Control","AlzheimerDisease"))) |> 
  #dplyr::filter(is.na(reason_excluded)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 38)
    return(.)
  })()



data.ad <- data.mvalues.alzheimer.dirty |> 
  dplyr::select(metadata.ad$DNAm_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 470392) #470691
    return(.)
  })()



design.ad <- model.matrix(~factor(ad.status), data=metadata.ad)
fit.ad <- limma::lmFit(data.ad, design.ad)
fit.ad <- limma::eBayes(fit.ad, trend=T)
stats.ad <- limma::topTable(fit.ad,
                            n=nrow(data.ad),
                            coef="factor(ad.status)AlzheimerDisease",
                            sort.by = "none",
                            adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 


rm(design.ad)


sum(stats.ad$P.Value < 0.01)
sum(stats.ad$adj.P.Val < 0.01)



#saveRDS(fit.ad, file="cache/analysis_differential__ad_co__fit.Rds")
saveRDS(stats.ad, file="cache/analysis_differential__ad_co__stats.Rds")


rm(fit.ad, stats.ad)



## check plot ----

plt <- stats.ad |> 
  dplyr::left_join(
    data.probes.alzheimer,
    by=c('probe_id'='probe_id'), suffix=c('','')
  )

# plot(plt$Beta..difference, plt$logFC)
# plot(plt$FDR.p.value, plt$adj.P.Val)





## x alzheimer quick test ----

plt <- stats.od |> 
  dplyr::left_join(data.probes.alzheimer, by=c('probe_id'='probe_id'),suffix=c('',''))

# ggplot(plt, aes(x=logFC, y=-adj.P.Val, col=is.na(FDR.p.value))) +
#   geom_point(data = subset(plt, is.na(FDR.p.value)),pch=19,cex=0.01) +
#   geom_point(data = subset(plt, !is.na(FDR.p.value)),pch=19,cex=0.15)

ggplot(plt, aes(x=Beta..difference, y=logFC, col=is.na(FDR.p.value), label=label)) +
  geom_point(data = subset(plt, is.na(FDR.p.value)),pch=19,cex=0.01) +
  geom_point(data = subset(plt, !is.na(FDR.p.value)),pch=19,cex=0.35)
#ggrepel::geom_text_repel(col="black",size=3)



plt <- stats.ac |> 
  dplyr::left_join(data.probes.alzheimer, by=c('probe_id'='probe_id'),suffix=c('',''))


# ggplot(plt, aes(x=logFC, y=-adj.P.Val, col=is.na(FDR.p.value))) +
#   geom_point(data = subset(plt, is.na(FDR.p.value)),pch=19,cex=0.01) +
#   geom_point(data = subset(plt, !is.na(FDR.p.value)),pch=19,cex=0.15)


# yes - concordance
ggplot(plt, aes(x=Beta..difference, y=logFC, col=is.na(FDR.p.value), label=label)) +
  geom_point(data = subset(plt, is.na(FDR.p.value)),pch=19,cex=0.01) +
  geom_point(data = subset(plt, !is.na(FDR.p.value)),pch=19,cex=0.35)
#ggrepel::geom_text_repel(col="black",size=3)






### subtests HG_AC specific ----


metadata <- rbind( metadata.od , metadata.ac) |> 
  dplyr::mutate(HG_OD = factor(ifelse(dataset == 'GLASS-OD' & LG_HG_status == "HG","HG_OD","non_HG_OD"),levels=c('non_HG_OD','HG_OD'))) |> 
  dplyr::mutate(HG_AC = factor(ifelse(dataset == 'GLASS-NL' & LG_HG_status == "HG","HG_AC","non_HG_AC"),levels=c('non_HG_AC','HG_AC'))) |> 
  dplyr::mutate(HG_stat = factor(dplyr::case_when(
    HG_OD == "HG_OD" ~ "HG_OD",
    HG_AC == "HG_AC" ~ "HG_AC",
    T ~ "LG"
  ), levels=c("LG","HG_AC","HG_OD")))

data <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata$sentrix_id)


#### A HG_AC -1 - HG_OD 1 ----

# HG_stat of HG_OD + HG_AC maakt niet veel uit
#model.matrix(~HG_OD + HG_AC, data=metadata) |> head()
#model.matrix(~HG_stat, data=metadata) |> head()


design <- model.matrix(~factor(patient_id) + factor(HG_OD) + factor(HG_AC), data=metadata)
fit <- limma::eBayes(limma::lmFit(data, design),trend=T)

#limma::makeContrasts(HG_OD - HG_AC - LG, levels=colnames(design))
contrast <- limma::makeContrasts(HG_statHG_OD - HG_statHG_AC - HG_statLG, levels=colnames(design))
fit2 <- limma::eBayes( limma::contrasts.fit( fit, contrast ) )
# 
# HG_statHG_OD - HG_statHG_AC
f = rownames(fit2$p.value)[fit2$p.value < 0.01]


#### B HG_AC 1 & HG_OD 1 ----

# HG_stat of HG_OD + HG_AC maakt niet veel uit
#model.matrix(~HG_OD + HG_AC, data=metadata) |> head()
#model.matrix(~HG_stat, data=metadata) |> head()


design <- model.matrix(~factor(patient_id) + factor(dataset) + factor(LG_HG_status), data=metadata)
fit <- limma::eBayes(limma::lmFit(data, design),trend=T)
stats.hg <- limma::topTable(fit,
                                 n=nrow(data),
                                 coef="factor(LG_HG_status)HG",
                                 sort.by = "none",
                                 adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')



design <- model.matrix(~factor(patient_id) + factor(HG_OD) + factor(HG_AC), data=metadata)
fit <- limma::eBayes(limma::lmFit(data, design),trend=T)
stats.a.hg_ac <- limma::topTable(fit,
                                 n=nrow(data),
                                 coef="factor(HG_AC)HG_AC",
                                 sort.by = "none",
                                 adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')
stats.a.hg_od <- limma::topTable(fit,
                                 n=nrow(data),
                                 coef="factor(HG_OD)HG_OD",
                                 sort.by = "none",
                                 adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')



design <- model.matrix(~factor(patient_id) + factor(LG_HG_status) + factor(HG_AC), data=metadata)
fit <- limma::eBayes(limma::lmFit(data, design),trend=T)
stats.b.hg_ac <- limma::topTable(fit,
                                 n=nrow(data),
                                 coef="factor(HG_AC)HG_AC",
                                 sort.by = "none",
                                 adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')
#plot(plt$t.ac, stats.b.hg_ac$t, pch=19,cex=0.01)


design <- model.matrix(~factor(patient_id) + factor(LG_HG_status) + factor(HG_OD), data=metadata)
fit <- limma::eBayes(limma::lmFit(data, design),trend=T)
stats.b.hg_od <- limma::topTable(fit,
                            n=nrow(data),
                            coef="factor(HG_OD)HG_OD",
                            sort.by = "none",
                            adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')
#plot(plt$t.ac, stats.b.hg_ac$t, pch=19,cex=0.01)
#plot(plt$t.od, stats.b.hg_od$t, pch=19,cex=0.01)


design <- model.matrix(~factor(patient_id) + factor(HG_AC), data=metadata)
fit <- limma::eBayes(limma::lmFit(data, design),trend=T)
stats.c.hg_ac <- limma::topTable(fit,
                                 n=nrow(data),
                                 coef="factor(HG_AC)HG_AC",
                                 sort.by = "none",
                                 adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')


design <- model.matrix(~factor(patient_id) + factor(HG_OD), data=metadata)
fit <- limma::eBayes(limma::lmFit(data, design),trend=T)
stats.c.hg_od <- limma::topTable(fit,
                                 n=nrow(data),
                                 coef="factor(HG_OD)HG_OD",
                                 sort.by = "none",
                                 adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')



plt <- stats.hg |>  dplyr::rename_with(~ paste0(.x, ".hg")) |>
  dplyr::left_join(
    stats.a.hg_ac |> dplyr::rename_with(~ paste0(.x, ".a.hg_ac")),
    by = c("probe_id.hg" = "probe_id.a.hg_ac")
  ) |>
  dplyr::left_join(
    stats.a.hg_od |> dplyr::rename_with(~ paste0(.x, ".a.hg_od")),
    by = c("probe_id.hg" = "probe_id.a.hg_od")
  ) |>
  dplyr::left_join(
    stats.b.hg_ac |> dplyr::rename_with(~ paste0(.x, ".b.hg_ac")),
    by = c("probe_id.hg" = "probe_id.b.hg_ac")
  ) |>
  dplyr::left_join(
    stats.b.hg_od |> dplyr::rename_with(~ paste0(.x, ".b.hg_od")),
    by = c("probe_id.hg" = "probe_id.b.hg_od")
  ) |>
  dplyr::left_join(
    stats.c.hg_ac |> dplyr::rename_with(~ paste0(.x, ".c.hg_ac")),
    by = c("probe_id.hg" = "probe_id.c.hg_ac")
  ) |>
  dplyr::left_join(
    stats.c.hg_od |> dplyr::rename_with(~ paste0(.x, ".c.hg_od")),
    by = c("probe_id.hg" = "probe_id.c.hg_od")
  )


c = plt |>
  dplyr::select(probe_id.hg, starts_with("t.")) |> 
  tibble::column_to_rownames('probe_id.hg') |> 
  as.matrix() |> 
  cor()

corrplot::corrplot(c, order="hclust")


sum(stats.hg$adj.P.Val < 0.01)

sum(stats.a.hg_od$adj.P.Val < 0.01)
sum(stats.a.hg_ac$adj.P.Val < 0.01)
plot(stats.a.hg_od$t, stats.a.hg_ac$t, pch=19,cex=0.01)

sum(stats.b.hg_od$adj.P.Val < 0.01)
sum(stats.b.hg_ac$adj.P.Val < 0.01)
plot(stats.b.hg_od$t, stats.b.hg_ac$t, pch=19,cex=0.01)

sum(stats.c.hg_od$adj.P.Val < 0.01)
sum(stats.c.hg_ac$adj.P.Val < 0.01)
plot(stats.c.hg_od$t, stats.c.hg_ac$t, pch=19,cex=0.01)


plt <- plt |> 
  dplyr::mutate(col1 = ifelse(adj.P.Val.a.hg_ac < 0.01, "signi AC", "AC-n")) |> 
  dplyr::mutate(col2 = ifelse(adj.P.Val.a.hg_od < 0.01, "signi OD", "OD-n")) |> 
  #dplyr::mutate(col1 = ifelse(P.Value.c.hg_ac < 0.01, "signi AC-only", "AC-n")) |> 
  #dplyr::mutate(col2 = ifelse(P.Value.c.hg_od < 0.01, "signi OD-only", "OD-n")) |> 
  dplyr::mutate(col3 = ifelse(adj.P.Val.b.hg_ac < 0.01, "signi AC-only", "AC-n")) |> 
  dplyr::mutate(col4 = ifelse(adj.P.Val.b.hg_od < 0.01, "signi OD-only", "OD-n")) |> 
  dplyr::mutate(col5 = paste0(col1, " - ", col2, " - ", col3, " - ", col4))


ggplot(plt, aes(x=t.a.hg_ac, y=t.a.hg_od, col=col5)) +
  geom_point(pch=19,cex=0.01,alpha=0.5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  xlim(-12,12) +
  ylim(-12,12) +
  theme_bw()

ggplot(plt, aes(x=t.b.hg_ac, y=t.b.hg_od, col=col5)) +
  geom_point(pch=19,cex=0.01,alpha=0.5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  xlim(-12,12) +
  ylim(-12,12) +
  theme_bw()

ggplot(plt, aes(x=t.c.hg_ac, y=t.c.hg_od, col=col5)) +
  geom_point(pch=19,cex=0.01,alpha=0.5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  xlim(-12,12) +
  ylim(-12,12) +
  theme_bw()


# sandbox ----


d = data.frame(
  status=c(
    rep("g2",25),
    rep("g2+tmz",4),
    rep("g3+tmz",18),
    rep("g3",17)
  )
) |> 
  dplyr::mutate(grade = gsub("^(..).+$","\\1",status)) |> 
  dplyr::mutate(tmz = grepl("tmz", status)) |> 
  dplyr::mutate(dat1 = runif(dplyr::n()) + (runif(dplyr::n()) * (grade == "g3") * 0.8)) |> 
  dplyr::mutate(dat1 = runif(dplyr::n()) + (runif(dplyr::n()) * (grade == "g3") * 0.8)) |> 
  dplyr::mutate(dat1 = runif(dplyr::n()) + (runif(dplyr::n()) * (grade == "g3") * 0.8) + (runif(dplyr::n()) * (status == "g2+tmz") * 0.4))


summary(lm(dat1 ~ grade , data=d))$coefficients
summary(lm(dat1 ~ tmz , data=d))$coefficients
summary(lm(dat1 ~ grade + tmz, data=d))$coefficients




