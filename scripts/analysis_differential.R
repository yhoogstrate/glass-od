#!/usr/bin/env R


# load data ----


library(ggplot2)
#library(minfi)


source('scripts/load_constants.R')
source('scripts/load_functions.R')
source('scripts/load_themes.R')


if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
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
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_SAMPLES_INCLUDED) |> 
  filter_primaries_and_last_recurrences(178) |> 
  
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"remainder")))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 
  dplyr::select(array_sentrix_id, patient, pr.status) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (178)) 
    return(.)
  })()


data.pp.nc <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_PROBES_UNMASKED_AND_DETP)) 
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


## data: partially paired INTENSITY [w/o FFPE/frozen batch correct] ----



metadata.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_SAMPLES_INCLUDED) |> 
  filter_primaries_and_last_recurrences(178) |> 
  
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"remainder")))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 
  dplyr::select(array_sentrix_id, patient, pr.status) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (178)) 
    return(.)
  })()



data.pp.nc <- data.intensities.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.intensities.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693060)) 
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



#saveRDS(fit.pp.nc, file="cache/analysis_differential_intensities__primary_recurrence__partial_paired_nc__fit.Rds")
saveRDS(stats.pp.nc, file="cache/analysis_differential_intensities__primary_recurrence__partial_paired_nc__stats.Rds")


rm(fit.pp.nc, stats.pp.nc)



## x data: full paired only ----



metadata.fp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_SAMPLES_INCLUDED) |> 
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
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_SAMPLES_INCLUDED) |> 
  filter_first_G2_and_last_G3(150) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder")))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2", "Grade3")) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (150)) 
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




## data: partially paired INTENSITIES [w/o FFPE/frozen batch correct] ----


metadata.g2g3.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_SAMPLES_INCLUDED) |> 
  filter_first_G2_and_last_G3(150) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder")))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2", "Grade3")) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (150)) 
    return(.)
  })()



data.g2g3.pp.nc <- data.intensities.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.intensities.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.g2g3.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693060)) 
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
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_SAMPLES_INCLUDED) 

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



# analyses: GLASS-OD AcCGAP ----
## data: partially paired [w/o FFPE/frozen batch correct] ----


metadata.AccGAP.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_SAMPLES_INCLUDED) |> 
  assertr::verify(!is.na(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  assertr::verify(is.numeric(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  dplyr::mutate(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit = scale(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |>  # also pairs with n = 3 & n= 4
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_",ifelse(is.paired,patient_id,"remainder")))) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (CONST_N_GLASS_OD_SAMPLES_INCLUDED)) 
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
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_SAMPLES_INCLUDED) |> 
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
    assertthat::assert_that(nrow(.) == (202)) 
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


metadata.ffpe_decay.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_SAMPLES_INCLUDED) |> 
  dplyr::filter(!is.na(isolation_material)) |> 
  dplyr::mutate(ffpe_decay_time = ifelse(isolation_material == "ffpe", -time_between_resection_and_array, 0)) |> 
  dplyr::filter(!is.na(ffpe_decay_time)) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_",ifelse(is.paired,patient_id,"remainder")))) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 202) 
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


rm(fit.ffpe_decay.pp.nc, stats.ffpe_decay.pp.nc)



## data: partially paired INTENSITIES [w/o FFPE/frozen batch correct] ----


metadata.ffpe_decay.pp.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_SAMPLES_INCLUDED) |> 
  dplyr::filter(!is.na(isolation_material)) |> 
  dplyr::mutate(ffpe_decay_time = ifelse(isolation_material == "ffpe", -time_between_resection_and_array, 0)) |> 
  dplyr::filter(!is.na(ffpe_decay_time)) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() >= 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p_",ifelse(is.paired,patient_id,"remainder")))) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (202)) 
    return(.)
  })()


data.ffpe_decay.pp.nc <- data.intensities.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.intensities.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.ffpe_decay.pp.nc$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693060)) 
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



#saveRDS(fit.ffpe_decay.pp.nc, file="cache/analysis_differential_intensities__ffpe-decay-time__partial_paired_nc__fit.Rds")
saveRDS(stats.ffpe_decay.pp.nc, file="cache/analysis_differential_intensities__ffpe-decay-time__partial_paired_nc__stats.Rds")


rm(fit.ffpe_decay.pp.nc, stats.ffpe_decay.pp.nc)



## data: unpaired [w/o FFPE/frozen batch correct] ----


metadata.ffpe_decay.up.nc <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_SAMPLES_INCLUDED) |> 
  dplyr::filter(!is.na(isolation_material)) |> 
  dplyr::mutate(ffpe_decay_time = ifelse(isolation_material == "ffpe", -time_between_resection_and_array, 0)) |> 
  dplyr::filter(!is.na(ffpe_decay_time)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 202) 
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



#saveRDS(fit.ffpe_decay.up.nc, file="cache/analysis_differential__ffpe-decay-time__unpaired_nc__fit.Rds")
saveRDS(stats.ffpe_decay.up.nc, file="cache/analysis_differential__ffpe-decay-time__unpaired_nc__stats.Rds")


rm(fit.ffpe_decay.up.nc, stats.ffpe_decay.up.nc)



# analyses: epiGenetic clocks ----
## data: unpaired [w/o FFPE/frozen batch correct] ----


clocks <- glass_od.metadata.array_samples |> 
  dplyr::select(contains("epiTOC") | contains("dnaMethyAge_")) |> 
  colnames()


for(clock in clocks) {
  print(clock)
  
  
  metadata.current_clock.up.nc <- glass_od.metadata.array_samples |> 
    filter_GLASS_OD_idats(CONST_N_GLASS_OD_SAMPLES_INCLUDED) |> 
    data.table::copy() |> # odd hack needed, because setnames also affects the former "glass_od.metadata.array_samples" object...
    data.table::setnames(old = c(clock), new = c("array_current_clock")) |> 
    dplyr::filter(!is.na(array_current_clock)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == CONST_N_GLASS_OD_SAMPLES_INCLUDED) 
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
  
  
  
  
  design.current_clock.up.nc <- model.matrix(~array_current_clock, data=metadata.current_clock.up.nc)
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
  
  
  rm(fit.current_clock.up.nc)
  rm(stats.current_clock.up.nc)
  rm(clock)
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
    assertthat::assert_that(nrow(.) == (470392)) #470691
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



saveRDS(fit.ad, file="cache/analysis_differential__ad_co__fit.Rds")
saveRDS(stats.ad, file="cache/analysis_differential__ad_co__stats.Rds")


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


