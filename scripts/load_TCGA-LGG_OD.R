#!/usr/bin/env R


# load data ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')





# a. clinical.tsv ----


tcga_lgg_od.metadata.array_samples.a <- read.delim("data/tcga-lgg/clinical.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) |>
  dplyr::filter(project.project_id == "TCGA-LGG") |> 
  dplyr::filter(diagnoses.classification_of_tumor == "primary") |>
  dplyr::rename(a_diagnoses.primary_diagnosis = diagnoses.primary_diagnosis) |> 
  dplyr::mutate(a_clincal.tsv__days_to_death = as.numeric(ifelse(demographic.days_to_death == "'--", NA, demographic.days_to_death))) |>
  dplyr::select(cases.submitter_id, a_diagnoses.primary_diagnosis, a_clincal.tsv__days_to_death) |> 
  dplyr::arrange(cases.submitter_id, a_diagnoses.primary_diagnosis) |> 
  dplyr::distinct()


# b. TCGAbiolinks::GDCquery_clinic ----

## Clinical  -----


if(!file.exists('cache/TCGA-LGG.TCGAbiolinks.Rds')) {
  clin_lgg <-  TCGAbiolinks::GDCquery_clinic(project = "TCGA-LGG", type = "clinical", save.csv = T)
  saveRDS(clin_lgg, file='cache/TCGA-LGG.TCGAbiolinks.Rds')
} else {
  clin_lgg <- readRDS('cache/TCGA-LGG.TCGAbiolinks.Rds')
}



tcga_lgg_od.metadata.array_samples.b <- rbind(clin_lgg) |> 
  dplyr::mutate(sites_of_involvement = NULL) |> 
  dplyr::mutate(Deceased = ifelse(vital_status == "Dead", 1,0)) |> 
  dplyr::mutate(Deceased = as.factor(Deceased)) |> 
  dplyr::mutate(Survival = ifelse(vital_status == "Alive", days_to_last_follow_up, days_to_death)) |> 
  dplyr::mutate(b_tcga_biolinks_survival = as.numeric(Survival)) |> 
  dplyr::mutate(Deceased = ifelse( (!is.na(Survival)) & is.na(Deceased) , 0, Deceased)) |> 
  dplyr::mutate(age_at_diagnosis = age_at_diagnosis / 365) |> 
  dplyr::mutate(b_tcgabiolinks_deceased = ifelse(Deceased == 2, 1, 0)) |> 
  dplyr::rename(b_tcga_biolinks_primary_diagnosis = primary_diagnosis) |> 
  dplyr::select(submitter_id, 
                tumor_grade,
                b_tcga_biolinks_primary_diagnosis,
                b_tcgabiolinks_deceased, 
                b_tcga_biolinks_survival,
                prior_treatment, 
                age_at_diagnosis,
                gender
                )

rm(clin_lgg)


# c. santoesha catnon folder paper suppl ----


tcga_lgg_od.metadata.array_samples.c <- readxl::read_xlsx('data/tcga-lgg/TCGA methylation subtypes.xlsx',
                              sheet = 'S1A. TCGA discovery dataset', skip = 1) |> 
  dplyr::filter(`IDH/codel subtype` == 'IDHmut-codel')  |> 
  dplyr::select(Case, 
                `IDH/codel subtype`,
                Grade, 
                `Vital status (1=dead)`,
                `Survival (months)`) |> 
  dplyr::mutate(c_xlsx_survival_event = as.numeric(ifelse(`Vital status (1=dead)` == "NA", NA, `Vital status (1=dead)`)), `Vital status (1=dead)` = NULL) |> 
  dplyr::mutate(`Survival (months)` = as.numeric(ifelse(`Survival (months)` == "NA", NA , `Survival (months)`))) |> 
  dplyr::mutate(c_xlsx_survival = round(`Survival (months)` * CONST_DAYS_PER_MONTH), `Survival (months)` = NULL)



# d. nationwidechildrens.org_clinical_patient_lgg.txt ----


tcga_lgg_od.metadata.array_samples.d  <- read.delim('data/tcga-lgg/nationwidechildrens.org_clinical_patient_lgg.txt') |> 
  dplyr::filter(bcr_patient_uuid %in% c('bcr_patient_uuid', 'CDE_ID:') == F) |> 
  dplyr::select(bcr_patient_barcode, tumor_grade, histologic_diagnosis, death_days_to, last_contact_days_to) |> 
  dplyr::rename(d_histologic_diagnosis = histologic_diagnosis) |> 
  dplyr::mutate(d_survival_nationwidechildrens.org = as.numeric(dplyr::case_when(
    death_days_to == "[Not Applicable]" & last_contact_days_to == "[Not Available]" ~ NA,
    death_days_to == "[Not Applicable]" & last_contact_days_to != "[Not Available]" ~ last_contact_days_to,
    death_days_to != "[Not Applicable]" & last_contact_days_to == "[Not Available]" ~ death_days_to
  ))) |> 
  dplyr::mutate(d_survival_nationwidechildrens.org_event = as.numeric(dplyr::case_when(
    death_days_to == "[Not Applicable]" & last_contact_days_to == "[Not Available]" ~ NA,
    death_days_to == "[Not Applicable]" & last_contact_days_to != "[Not Available]" ~ 0,
    death_days_to != "[Not Applicable]" & last_contact_days_to == "[Not Available]" ~ 1
  ))) |> 
  dplyr::mutate(
    death_days_to = NULL,  last_contact_days_to = NULL
  )


# merge ----

# a is redundant t b but lacks other data
# d includes only entries already (or superseded) in c

tcga_lgg_od.metadata.array_samples <- tcga_lgg_od.metadata.array_samples.c |> 
  dplyr::left_join(tcga_lgg_od.metadata.array_samples.b, by=c('Case'='submitter_id')) |> 

  dplyr::select(
    Case,
    
    prior_treatment,
    age_at_diagnosis,
    gender,
    
    `IDH/codel subtype`,
    b_tcga_biolinks_primary_diagnosis,
    
    
    Grade, 
    tumor_grade,
    
    b_tcgabiolinks_deceased,
    c_xlsx_survival_event,

    b_tcga_biolinks_survival,
    c_xlsx_survival
  ) |> 
  dplyr::rename(patient_id = Case)


rm(tcga_lgg_od.metadata.array_samples.a)
rm(tcga_lgg_od.metadata.array_samples.b)
rm(tcga_lgg_od.metadata.array_samples.c)
rm(tcga_lgg_od.metadata.array_samples.d)
gc()


# fix grade ----


tcga_lgg_od.metadata.array_samples <- tcga_lgg_od.metadata.array_samples |> 
  dplyr::mutate(tumor_grade = ifelse(is.na(tumor_grade), gsub("G4","G3", Grade), tumor_grade) , Grade =  NULL) |> 
  dplyr::mutate(survival = dplyr::case_when(
    b_tcgabiolinks_deceased == 1 & is.na(b_tcga_biolinks_survival) ~ c_xlsx_survival, # biolinks entry does not make sense
    b_tcgabiolinks_deceased == 1 & !is.na(b_tcga_biolinks_survival) ~ b_tcga_biolinks_survival, # biolinks entry does not make sense
    !is.na(c_xlsx_survival) & c_xlsx_survival > 0 ~ c_xlsx_survival,
    T ~ NA
  )) |> 
  dplyr::mutate(survival_event = dplyr::case_when(
    b_tcgabiolinks_deceased == 1 & is.na(b_tcga_biolinks_survival) ~ c_xlsx_survival_event, # biolinks entry does not make sense
    b_tcgabiolinks_deceased == 1 & !is.na(b_tcga_biolinks_survival) ~ b_tcgabiolinks_deceased, # biolinks entry does not make sense
    !is.na(c_xlsx_survival) & c_xlsx_survival > 0 ~ c_xlsx_survival_event,
    T ~ NA
  )) |> 
  dplyr::mutate(b_tcgabiolinks_deceased = NULL, b_tcga_biolinks_survival = NULL) |> 
  dplyr::mutate(c_xlsx_survival_event = NULL, c_xlsx_survival = NULL)


# 
# tcga_lgg_od.metadata.array_samples |> 
#   View()
# 


# add CGC ----



tmp <- readRDS(file="cache/analysis_A_IDH_HG__A_IDH_LG_lr__lasso_fit__TCGA-LGG_450k.Rds") |> 
  dplyr::rename(sample_id_primary_resection = array_sentrix_id) |> 
  dplyr::mutate(patient_id = gsub("\\-[0-9A-Z]{3}$","", sample_id_primary_resection))


tcga_lgg_od.metadata.array_samples <- tcga_lgg_od.metadata.array_samples |> 
  assertr::verify(tmp$patient_id %in% patient_id) |> 
  dplyr::left_join(tmp, by=c('patient_id'='patient_id'), suffix = c('',''))



rm(tmp)



# load idats ----



if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}



probes.850k <-  data.mvalues.hq_samples |>  
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  dplyr::pull('probe_id')


tmp <- list.files(path =  "/data/cognition/data/DNA_methylation/idat/TCGA-LGG/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |> 
  data.frame(array_filename = _) |>
  dplyr::mutate(array_filename = paste0("/data/cognition/data/DNA_methylation/idat/TCGA-LGG/", array_filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (2 * 534))
    return(.)
  })() |> 
  assertr::verify(file.exists(array_filename)) |>
  dplyr::mutate(tcga_sample_id = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", array_filename)) |>
  
  dplyr::filter(grepl("01[A-Z]$", tcga_sample_id)) |> 
  
  dplyr::mutate(array_channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", array_filename)) |>
  tidyr::pivot_wider(id_cols = tcga_sample_id, names_from = array_channel, values_from = c(array_filename)) |>
  dplyr::mutate(tcga_patient_id = gsub("^(.+\\-.+)-.+$", "\\1", tcga_sample_id)) |>
  dplyr::rename(array_channel_green = Grn) |>
  dplyr::rename(array_channel_red = Red) |>
  dplyr::mutate(array_channel_green_filesize = file.info(array_channel_green)$size) |> 
  dplyr::mutate(array_channel_red_filesize   = file.info(array_channel_red)$size) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 516)
    return(.)
  })() |>
  assertr::verify(!is.na(array_channel_green)) |>
  assertr::verify(!is.na(array_channel_red)) |> 
  assertr::verify(!duplicated(tcga_sample_id)) |> 
  assertr::verify(!duplicated(tcga_patient_id)) |> 
  
  dplyr::filter(tcga_patient_id %in% tcga_lgg_od.metadata.array_samples$patient_id) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 169)
    return(.)
  })() |> 
  
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL) |> 
  dplyr::mutate(Sample_Name = tcga_sample_id)


if(!file.exists("cache/mvalues.tcga-lgg.Rds")) {
  RGSet <- minfi::read.metharray.exp(targets = tmp, force = T, verbose=F) #red/green channel together
  rm(tmp)
  proc <- minfi::preprocessNoob(RGSet, offset = 0, dyeCorr = TRUE, verbose = TRUE, dyeMethod="single")  # dyeMeth
  
  
  mvalue <- minfi::ratioConvert(proc, what = "M") |> 
    assays() |> 
    purrr::pluck('listData') |> 
    purrr::pluck("M") |> 
    data.table::as.data.table(keep.rownames = "probe_id") |> 
    dplyr::filter(probe_id %in% probes.850k)
  
  
  
  saveRDS(mvalue, file=paste0("cache/mvalues.tcga-lgg.Rds"))
  
  
  rm(proc, RGSet, mvalue)
  gc()
}


tcga_lgg.mvalues <- readRDS("cache/mvalues.tcga-lgg.Rds")

