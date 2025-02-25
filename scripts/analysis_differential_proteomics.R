#!/usr/bin/env

# libs ----


library(limma)
library(ggplot2)
library(patchwork)
library(EnhancedVolcano)

library(recursiveCorPlot)


# load data ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')


if(!exists('glass_od.metadata.proteomics')) {
  source('scripts/load_proteomics.R')
}





# limma N/A ----




pp <- function(measurements, patient_ids, conditions, coeff) {
  stopifnot(length(patient_ids) == length(conditions))
  #stopifnot(is.factor(conditions))
  stopifnot(length(patient_ids) == length(measurements))
  
  metadata <- data.frame(
    patient_id = patient_ids[!is.na(measurements)],
    condition  = conditions[  !is.na(measurements)]
  ) |> 
    dplyr::mutate(patient_id = as.character(patient_id)) |> 
    dplyr::group_by(patient_id) |> 
    dplyr::mutate(is.paired = dplyr::n() == 2) |> 
    dplyr::ungroup() |> 
    
    dplyr::mutate(patient = as.factor(ifelse(is.paired, paste0("p",patient_id), "a_remainder"))) |> 
    dplyr::mutate(patient_id = NULL)
  
  measurements <- tibble::tibble(t(tibble::tibble(measurements[!is.na(measurements)])))

  tryCatch({
    design <- model.matrix(~factor(patient) + condition, data=metadata)
    fit <- limma::lmFit(measurements, design)
    fit <- limma::eBayes(fit, trend=T)
    
    #print(colnames(fit$coefficients))
    
    stats <- limma::topTable(fit,
                             n=nrow(measurements),
                             coef=coeff,
                             sort.by = "none",
                             adjust.method="fdr") |> 
      tibble::rownames_to_column('probe_id') |> 
      dplyr::mutate(probe_id = NULL) |> 
      dplyr::mutate(n = nrow(metadata))
    
    return(stats)
  },
  error = function(w) { return(NA)}, # in case there are not sufficient samples per condition (too many NA's)
  warning = function(w) { return(NA)}
  )
  
}





## correlation w/ KI67 ----





ki67.proteomics <-  data.proteomics.glass_od  |>
  as.data.frame() |> 
  tibble::rownames_to_column('gene') |> 
  dplyr::filter(grepl("KI67", gene)) |> 
  tibble::column_to_rownames('gene') |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('proteomics_id') |> 
  dplyr::filter(proteomics_id != "GLODprot_03_1521") |> # doesn't matter which one
  dplyr::filter(proteomics_id != "GLODprot_09_1522") |> # no value
  dplyr::filter(proteomics_id != "GLODprot_11_1524")
dim(ki67.proteomics)


ki67.proteomics <-  ki67.proteomics |> 
  dplyr::left_join(
    glass_od.metadata.proteomics |>  dplyr::select(proteomics_id, resection_id, proteomics_notes),
    by=c('proteomics_id'='proteomics_id')
  )
dim(ki67.proteomics)
stopifnot(duplicated(ki67.proteomics$resection_id) == 0)








ki67.if <- glass_od.metadata.array_samples |>
  dplyr::select(resection_id,isolation_id, contains("67")) |>
  dplyr::filter(!is.na(staining_KI67_md5sum) & isolation_id %in% c("0047-RX", "0047-RY") == F)

sum(duplicated(ki67.if$resection_id))

dups <- ki67.if |>
  dplyr::filter(duplicated(resection_id)) |> 
  dplyr::pull(resection_id)

ki67.if <- ki67.if |> dplyr::filter((resection_id %in% dups == F) | (resection_id %in% dups == T & !is.na(array_PC67)))

sum(duplicated(ki67.if$resection_id))




ki67.joined = ki67.proteomics |> dplyr::inner_join(ki67.if, by=c('resection_id'='resection_id'))
sum(duplicated(ki67.joined$resection_id))


ggplot(ki67.joined, aes(x=MKI67, y=staining_KI67_pos_per_detected_cells)) +
  ggpubr::stat_cor(method = "spearman") +
  geom_point() +
  theme_nature



ggplot(ki67.joined, aes(x=MKI67, y=staining_KI67_lr_pos_neg_cells)) +
  ggpubr::stat_cor(method = "spearman") +
  geom_point() +
  theme_nature




ggplot(ki67.joined, aes(x=MKI67, y=staining_KI67_pos_neg_cell_density)) +
  ggpubr::stat_cor(method = "spearman") +
  geom_point() +
  theme_nature




ggplot(ki67.joined, aes(x=MKI67, y=staining_KI67_pos_per_area_um2)) +
  ggpubr::stat_cor(method = "spearman") +
  geom_point() +
  theme_nature




cor(exp(ki67.joined$MKI67),
    ki67.joined$staining_KI67_pos_per_area_um2,
    use="complete.obs", method="spearman")



# 0.4589544



## WHO Grade ----
### regular ----

View(glass_od.metadata.proteomics)

metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::filter(!is.na(resection_tumor_grade)) |>  # 0053-R2 pathology was not decisve
  dplyr::select(resection_id, patient_id, starts_with("proteomics_"), starts_with("resection_"), contains("_excluded")) |> 
  dplyr::mutate(condition = factor(paste0("g", resection_tumor_grade), levels=c("g2", "g3"))) |> 
  
  # dplyr::group_by(patient_id) |> 
  # dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  # dplyr::ungroup() |> 
  # dplyr::mutate(patient = as.factor(ifelse(is.paired,paste0("p",patient_id),"a_remainder"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 139) 
    return(.)
  })()




data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id)
data$nas <- rowSums(is.na(data))
data <- data |> 
  dplyr::filter(nas < 60) |> 
  dplyr::mutate(nas=NULL)


stopifnot(metadata$proteomics_id == colnames(data))



design <- model.matrix(~factor(condition), data=metadata)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit, trend=T)
stats.grade <- limma::topTable(fit,
                               n=nrow(data),
                               coef="factor(condition)g3",
                               sort.by = "none",
                               adjust.method="fdr") |> 
  tibble::rownames_to_column('protein_id') 

sum(stats.grade$adj.P.Val < 0.01)


rm(design, fit, data, metadata)



### pat corrected ----



metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::filter(!is.na(resection_tumor_grade)) |>  # 0053-R2 pathology was not decisve
  
  dplyr::select(proteomics_id, patient_id, resection_tumor_grade) |> 
  dplyr::filter(!is.na(resection_tumor_grade)) |> 
  dplyr::mutate(resection_tumor_grade = factor(ifelse(resection_tumor_grade == 2, "Grade2", "Grade3"), levels=c("Grade2", "Grade3"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 139) 
    return(.)
  })()


data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id) 



stats.grade.pat.corrected <- data |> 
  tibble::rownames_to_column("protein_id") |> 
  dplyr::select("protein_id") |> 
  dplyr::mutate(tmp = pbapply::pbapply(data, 1, pp, metadata$patient_id, metadata$resection_tumor_grade, "conditionGrade3")) |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(adj.P.Val = NULL) |> 
  dplyr::mutate(adj.P.Val = p.adjust(P.Value, "fdr")) |> 
  tibble::column_to_rownames("protein_id")


sum((stats.grade.pat.corrected |> dplyr::filter(!is.na(P.Value)) |> dplyr::pull(P.Value) ) < 0.01)
sum((stats.grade.pat.corrected |> dplyr::filter(!is.na(adj.P.Val)) |> dplyr::pull(adj.P.Val) ) < 0.01)


rm(metadata, data)



## prim-rec ----
### regular ----


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, patient_id, starts_with("proteomics_"), starts_with("resection_"), contains("_excluded")) |> 
  dplyr::mutate(condition = factor(ifelse(resection_number == 1,"primary","recurrence"), levels=c("primary","recurrence"))) |> 

  # dplyr::group_by(patient_id) |> 
  # dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  # dplyr::ungroup() |> 
  # dplyr::mutate(patient = as.factor(ifelse(is.paired,paste0("p",patient_id),"a_remainder"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 140) 
    return(.)
  })()




data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id)
data$nas <- rowSums(is.na(data))
data <- data |> 
  dplyr::filter(nas < 60) |> 
  dplyr::mutate(nas=NULL)

stopifnot(metadata$proteomics_id == colnames(data))




design <- model.matrix(~ factor(condition), data=metadata)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit, trend=T)
stats.time <- limma::topTable(fit,
                         n=nrow(data),
                         coef="factor(condition)recurrence",
                         sort.by = "none",
                         adjust.method="fdr") |> 
  tibble::rownames_to_column('protein_id') 

sum(stats.time$adj.P.Val < 0.01)


rm(design, fit, data, metadata)




### pat corrected ----



metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 

  dplyr::select(proteomics_id, patient_id, resection_number) |> 
  dplyr::mutate(primrec = factor(ifelse(resection_number == 1,"primary","recurrence"), levels=c("primary","recurrence"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 140) 
    return(.)
  })()



data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id) 



stats.time.pat.corrected <- data |> 
  tibble::rownames_to_column("protein_id") |> 
  dplyr::select("protein_id") |> 
  dplyr::mutate(tmp = pbapply::pbapply(data, 1, pp, metadata$patient_id, metadata$primrec, "conditionrecurrence")) |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(adj.P.Val = NULL) |> 
  dplyr::mutate(adj.P.Val = p.adjust(P.Value, "fdr")) |> 
  tibble::column_to_rownames("protein_id")


sum((stats.time.pat.corrected |> dplyr::filter(!is.na(P.Value)) |> dplyr::pull(P.Value) ) < 0.01)
sum((stats.time.pat.corrected |> dplyr::filter(!is.na(adj.P.Val)) |> dplyr::pull(adj.P.Val) ) < 0.01)



rm(metadata, data)




## CGC ----
### regular ----



metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, starts_with("proteomics_")) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 140) 
    return(.)
  })()

metadata <- metadata |> 
  dplyr::inner_join(
    glass_od.metadata.array_samples |> 
      filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
      dplyr::filter(resection_id %in% metadata$resection_id) |> 
      assertr::verify(!duplicated(resection_id)),
    by=c('resection_id'='resection_id')
  ) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 118) 
    return(.)
  })() |> 

  # dplyr::group_by(patient_id) |>
  # dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  # dplyr::ungroup() |> 
  # dplyr::mutate(patient = as.factor(ifelse(is.paired,paste0("p",patient_id),"a_remainder"))) |> 
   
  dplyr::mutate(CGC = scale(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit))




data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id)
data$nas <- rowSums(is.na(data))
data <- data |> 
  dplyr::filter(nas < 60) |> 
  dplyr::mutate(nas=NULL)

stopifnot(metadata$proteomics_id == colnames(data))





design <- model.matrix(~CGC, data=metadata)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit, trend=T)
stats.cgc <- limma::topTable(fit,
                              n=nrow(data),
                              coef="CGC",
                              sort.by = "none",
                              adjust.method="fdr") |> 
  tibble::rownames_to_column('protein_id') 



rm(design, fit, data, metadata)




# ### pat corrected ----
# # not useful, purpose is to find those proteins that explain
# # CGC best, even if high in 2 resections of same  patient
# 
# metadata <- glass_od.metadata.proteomics |> 
#   dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
#   dplyr::filter(is.na(resection_reason_excluded)) |> 
#   dplyr::filter(is.na(patient_reason_excluded)) |> 
#   dplyr::filter(patient_study_name == "GLASS-OD") |> 
#   dplyr::select(resection_id, starts_with("proteomics_")) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == 140) 
#     return(.)
#   })()
# 
# metadata <- metadata |>
#   dplyr::inner_join(
#     glass_od.metadata.array_samples |> 
#       filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
#       dplyr::filter(resection_id %in% metadata$resection_id) |> 
#       assertr::verify(!duplicated(resection_id)),
#     by=c('resection_id'='resection_id')
#   ) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(nrow(.) == 119) 
#     return(.)
#   })() |> 
#   
#   dplyr::mutate(CGC = scale(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit))
# 
# 
# 
# 
# data <- data.proteomics.glass_od |> 
#   dplyr::select(metadata$proteomics_id)
# 
# 
# 
# stats.cgc.pat.corrected <- data |> 
#   tibble::rownames_to_column("protein_id") |> 
#   dplyr::select("protein_id") |> 
#   dplyr::mutate(tmp = pbapply::pbapply(data, 1, pp, metadata$patient_id, metadata$CGC, "condition")) |> 
#   tidyr::unnest(tmp) |> 
#   dplyr::mutate(adj.P.Val = NULL) |> 
#   dplyr::mutate(adj.P.Val = p.adjust(P.Value, "fdr")) |> 
#   tibble::column_to_rownames("protein_id")
# 
# 
# sum((stats.cgc.pat.corrected |> dplyr::filter(!is.na(P.Value)) |> dplyr::pull(P.Value) ) < 0.01)
# sum((stats.cgc.pat.corrected |> dplyr::filter(!is.na(adj.P.Val)) |> dplyr::pull(adj.P.Val) ) < 0.01)
# 
# 
# 
# rm(metadata, data)
# 



## plt tmp ----
### cell cycling ----


clng <- c("AAAS","AATF","ABCB1","ABL1","ABRAXAS1","ABRAXAS2","ACTB","ACTL6A","ACTL6B","ACTR2","ACTR3","ACTR5",
          "ACTR8","ACVR1","ACVR1B","ADAM17","ADAMTS1","ADARB1","ADCYAP1","AFAP1L2","AGO4","AHCTF1","AHR","AICDA",
          "AIF1","AJUBA","AKAP8","AKAP8L","AKT1","AKT2","ALKBH4","ALMS1","ALOX15B","AMBRA1","ANAPC1","ANAPC10",
          "ANAPC11","ANAPC13","ANAPC15","ANAPC16","ANAPC2","ANAPC4","ANAPC5","ANAPC7","ANGEL2","ANK3","ANKFN1",
          "ANKK1","ANKLE1","ANKLE2","ANKRD17","ANKRD31","ANKRD53","ANLN","ANXA1","ANXA11","APBB1","APBB2","APC",
          "APEX2","APP","APPL1","APPL2","ARF1","ARF6","ARHGEF10","ARHGEF2","ARID1A","ARID1B","ARID2","ARL2","ARL3",
          "ARL8A","ARL8B","ARPP19","ASAH2","ASCL1","ASNS","ASPM","ASZ1","ATAD5","ATF2","ATF5","ATM","ATP2B4","ATR","ATRIP","ATRX","AUNIP","AURKA","AURKAIP1","AURKB","AURKC","AVEN","AVPI1","AXIN2","AZI2","BABAM1","BABAM2","BAG6","BAK1","BANF1","BANP","BAP1","BARD1","BAX","BAZ1B","BBS4","BCAT1","BCCIP","BCL2","BCL2L1","BCL2L11","BCL6","BCL7A","BCL7B","BCL7C","BCR","BECN1","BEX2","BEX4","BID","BIN1","BIN3","BIRC2","BIRC3","BIRC5","BIRC6","BIRC7","BLCAP","BLM","BMAL1","BMP2","BMP4","BMP7","BOD1","BOD1L1","BOD1L2","BOLL","BOP1","BORA","BRCA1","BRCA2","BRCC3","BRD4","BRD7","BRD8","BRDT","BRINP1","BRINP2","BRINP3","BRIP1","BRME1","BROX","BRSK1","BRSK2","BTBD18","BTC","BTG1","BTG2","BTG3","BTG4","BTN2A2","BTRC","BUB1","BUB1B","BUB3","C10orf90","C11orf80","C14orf39","C1orf146","C2CD3",
          "C6orf89","C9orf78","CABLES1","CABLES2","CACNB4","CACUL1","CALM1",
          "CALM2","CALM3","CALR","CAMK1","CAMK2A","CAMSAP3","CAPN3","CASP2","CASP3","CATSPERZ","CAV2",
          "CCAR1","CCAR2","CCDC102B","CCDC124","CCDC42","CCDC57","CCDC61","CCDC69","CCDC8","CCL2","CCNA1","CCNA2","CCNB1","CCNB1IP1","CCNB2","CCNB3","CCNC","CCND1","CCND2","CCND3","CCNDBP1","CCNE1","CCNE2","CCNF","CCNG1","CCNG2","CCNH","CCNI","CCNI2","CCNJ","CCNJL","CCNK","CCNL1","CCNL2","CCNO","CCNP","CCNQ","CCNT1","CCNT2","CCNY","CCNYL1","CCP110","CCPG1","CCSAP","CD28","CD2AP","CDC123",
          "CDC14A","CDC14B","CDC14C","CDC16","CDC20","CDC23","CDC25A","CDC25B","CDC25C","CDC26","CDC27","CDC34","CDC37","CDC42","CDC45","CDC5L","CDC6","CDC7","CDC73","CDCA2","CDCA3","CDCA5","CDCA8","CDK1","CDK10","CDK11A","CDK11B","CDK12","CDK13","CDK14","CDK15","CDK16","CDK17","CDK18","CDK19","CDK2","CDK20","CDK2AP1","CDK2AP2","CDK3","CDK4","CDK5","CDK5R1","CDK5R2","CDK5RAP1","CDK5RAP2","CDK5RAP3","CDK6","CDK7","CDK8","CDK9","CDKL1","CDKL2","CDKL3","CDKL4","CDKL5","CDKN1A","CDKN1B","CDKN1C","CDKN2A","CDKN2B","CDKN2C","CDKN2D","CDKN3","CDT1","CEBPA","CECR2","CENATAC","CENPA","CENPC","CENPE","CENPF","CENPH","CENPI","CENPJ","CENPK","CENPL","CENPM","CENPN","CENPO","CENPP","CENPQ","CENPS",
          "CENPT","CENPU","CENPV","CENPW","CENPX","CEP120","CEP126","CEP131","CEP135","CEP152","CEP164","CEP192","CEP250","CEP295","CEP295NL","CEP44","CEP55","CEP63","CEP68","CEP72","CEP76","CEP85","CEP97","CETN1","CETN2","CETN3","CFL1","CGREF1","CGRRF1","CHAF1A","CHAF1B","CHAMP1","CHD3","CHEK1","CHEK2","CHFR","CHMP1A","CHMP1B","CHMP2A","CHMP2B","CHMP3","CHMP4A","CHMP4B","CHMP4C","CHMP5","CHMP6","CHMP7","CHORDC1","CHTF18","CHTF8","CIAO1","CIAO2A","CIAO2B","CIB1","CINP","CIT","CITED2","CKAP2","CKAP5","CKS1B","CKS2","CLASP1","CLASP2","CLIP1","CLOCK","CLSPN","CLTA","CLTC","CLTCL1","CNPPD1","CNTD1","CNTLN","CNTRL","CNTROB","COPS5","CPSF3","CRADD","CREB3","CREBL2","CRLF3","CRNN","CROCC","CRY1","CSAG1","CSNK1A1","CSNK1D","CSNK2A1","CSNK2A2","CSNK2A3","CSPP1","CTBP1","CTC1","CTCF",
          "CTCFL","CTDNEP1","CTDP1","CTDSP1","CTDSP2",
          "CTDSPL","CTNNB1","CUL1","CUL2","CUL3","CUL4A","CUL4B","CUL5","CUL7","CUL9","CUZD1","CXCR5","CYLD","CYP1A1","CYP26B1","CYP27B1","DAB2IP","DACH1","DACT1","DAPK3",
          "DAZL","DBF4","DBF4B","DCAF13","DCDC1","DCLRE1A","DCTN1","DCTN2","DCTN3","DCTN6","DCUN1D3","DDB1","DDIAS","DDIT3","DDR2","DDRGK1","DDX11","DDX39B","DDX3X","DDX4","DEUP1","DGKZ","DIAPH3","DIRAS3","DIS3L2","DLG1","DLGAP5","DMAP1","DMC1","DMRT1","DMRTC2","DMTF1","DNA2","DNM2","DNMT3L","DONSON","DOT1L","DPF1","DPF2","DPF3","DR1","DRD2","DRD3","DRG1","DSCC1","DSN1","DTL","DTX3L","DUSP1","DUSP13B","DUSP3","DUX4","DYNC1H1","DYNC1LI1","DYNLT1","DYNLT3","DYRK3","E2F1","E2F2","E2F3","E2F4","E2F6","E2F7","E2F8","E4F1","ECD","ECRG4","ECT2","EDN1","EDN3","EDNRA","EFHC1","EGF","EGFR","EHMT2","EID1","EIF2AK4","EIF4E","EIF4EBP1","EIF4G1","EIF4G2","EME1","EME2","EML1","EML3","EML4","ENKD1","ENSA","ENTR1","EP300",
          "EP400","EPB41","EPB41L2","EPC1","EPC2","EPGN","EPM2A","EPS8","ERCC2","ERCC3","ERCC4","ERCC6","ERCC6L","EREG","ERH","ESCO1","ESCO2","ESPL1","ESX1","ETAA1","ETS1","EVI2B","EVI5","EXD1","EXO1","EXOC1","EXOC2","EXOC3","EXOC4","EXOC5","EXOC6","EXOC6B","EXOC7","EXOC8","EZH2","EZR","FAM107A","FAM110A","FAM32A","FAM83D","FAM9A","FAM9B","
          FAM9C","FANCA","FANCD2","FANCI","FANCM","FAP","FBXL15","FBXL7","FBXO31","FBXO4","FBXO43","FBXO5","FBXO6","FBXO7","FBXW11","FBXW5","FBXW7","FEM1B","FEN1","FES","FGF10","FGF2","FGF8","FGFR1","FGFR2","FHL1","FIGN","FIGNL1","FIRRM","FKBP6","FLCN","FLNA","FMN2","FOSL1","FOXA1","FOXC1","FOXE3","FOXG1","FOXJ2","FOXJ3","FOXM1","FOXN3","FOXO4","FSD1","FZD3","FZD9","FZR1","GADD45A","GADD45B","GADD45G","GADD45GIP1","GAK","GAS1","GAS2","GATA3","GATA6","GBF1","GEM","GEN1","GFI1B","GIGYF2","GINS1","GINS3","GIPC1","GIT1","GJA1","GJC2","GLI1","GML","GMNC","GMNN","GNAI1","GNAI2","GNAI3","GOLGA2","GOLGA8B","GOLGA8F","GOLGA8S","GPER1","GPNMB","GPR132","GPR15LG","GPR3","GPSM1","GPSM2","GRK5","GSPT1","GSPT2","GTF2B","GTF2H1","GTPBP4","GTSE1","H1-8","H2AX","HACE1","HASPIN","HAUS1","HAUS2","HAUS3",
          "HAUS4","HAUS5","HAUS6","HAUS7","HAUS8","HBP1","HCFC1","HDAC3",
          "HDAC8","HECA","HECW2","HELLS","HEPACAM","HEPACAM2","HERC5","HES1","HEXIM1","HEXIM2","HFM1","HGF","HHEX","HINFP","HIPK2","HJURP",
          "HLA-G","HMCN1","HMG20B","HMGA2","HNRNPU","HORMAD1","HORMAD2","HOXA13","HOXC9","HOXD10","HPGD","HRAS","HSF1","HSF2BP","HSP90AB1","HSPA1A","HSPA1B","HSPA2","HTRA2","HTT","HUS1","HUS1B","HYAL1","ID2","ID4","IER3","IFFO1","IFNW1","IGF1","IGF1R","IGF2","IHO1","IK","IKZF1","IL10","IL1A","IL1B","INCA1",
          "INCENP","ING1","ING2","ING3","ING4","ING5","INHA","INHBA","INIP","INO80","INO80B","INO80C","INO80D","INO80E","INPPL1","INS","INSC","INSM1","INSM2","INSR","INTS13","INTS3","INTS7","IPO5","IPO7","IQGAP1","IQGAP2","IQGAP3","IRF1","IST1","ITGB1","ITGB3BP","JADE1","JADE2","JADE3","JTB","JUN","JUNB","JUND","KANK2","KASH5","KAT14","KAT2A","KAT2B","KAT5","KAT7","KATNA1","KATNB1","KCNA5","KCNH5","KCTD11","KCTD19","KDM8","KHDRBS1","KIAA0753","KIAA1614","KIF11","KIF13A","KIF14","KIF15","KIF18A","KIF18B","KIF20A","KIF20B","KIF22","KIF23","KIF25","KIF2A","KIF2B","KIF2C","KIF3A",
          "KIF3B","KIF4A","KIF4B","KIFC1","KIFC2","KIZ","KLF4","KLHDC3","KLHDC8B","KLHL13","KLHL18","KLHL21","KLHL22","KLHL42","KLHL9","KLK10","KLLN","KMT2E","KMT5A","KNL1","KNSTRN","KNTC1","KPNB1","KRT18","L3MBTL1","LATS1","LATS2","LCMT1","LEF1","LEP","LFNG","LGMN","LIF","LIG1","LIG3","LIG4","LIMK2","LIN54","LIN9","LIPA","LLGL1","LLGL2","LMLN","LMNA","LPIN1","LRP5","LRP6","LRRCC1","LSM10","LSM11","LSM14A","LZTS1","LZTS2","M1AP","MACROH2A1","MAD1L1","MAD2L1","MAD2L1BP","MAD2L2","MADD","MAEA","MAEL","MAGEA4","MAJIN","MAP10","MAP1S","MAP2K6","MAP3K11","MAP3K20","MAP3K7","MAP3K8","MAP4","MAP9","MAPK1","MAPK12","MAPK13","MAPK14","MAPK15","MAPK3","MAPK4","MAPK6","MAPK7","
          MAPRE1","MAPRE2","MAPRE3","MARF1","MARK4","MARVELD1","MASTL","MAU2","MBIP","MBLAC1","MBTD1","MBTPS1","MBTPS2","MCIDAS","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MCM8","MCMBP","MCMDC2","MCPH1","MCRS1","MCTS1","MDC1","MDM1","MDM2","MDM4","MEAF6","MECOM","MECP2","MED1","MEI1","MEI4","MEIKIN","MEIOB","MEIOC","MEIOSIN","MEIS2","MELK","MEN1","MEPCE","METTL13","METTL3","MICAL3","MIIP","MIS12","MIS18A","MIS18BP1","MISP","MITD1","MKI67","MLF1","MLH1","MLH3","MMS19","MN1","MNAT1","MND1","MNS1","MNT","MOK","MORF4L1","MORF4L2","MOS","MOV10L1","MPLKIP","MRE11","MRFAP1L2","MRGBP","MRGPRX2","MRNIP","MRPL41","MSH2","MSH4","MSH5","MSH6","MSX1","MSX2","MTA3","MTBP","MTCL1","MUC1","MUS81","MX2","MYB","MYBBP1A","MYBL1","MYBL2","MYC","MYH10","MYH9","MYO16","MYO19","MYOCD","MYOG","MZT1","NAA10","NAA50","NAA60","NABP1","NABP2","NAE1","NANOGP8","NANOS2","NANOS3","NASP","NAT10","NBN","NCAPD2","NCAPD3","NCAPG","NCAPG2","NCAPH","NCAPH2","NCOR1",
          "NDC1","NDC80","NDE1","NDEL1","NDP","NEDD1","NEDD9","NEK1","NEK10","NEK11","NEK2","NEK3","NEK4","NEK6","NEK7","NEK9","NES","NEUROG1","NF2","NFIA","NFIB","NFRKB","NHERF1","NIN","NIPBL","NKX3-1","NLE1","NLRP2B","NLRP5","NME6","NOLC1","None","None","None","None","None","None","None","None","None","None","NOP53","NOX5","NPAT","NPM1","NPM2",
          "NPPC","NPR2","NR2E1","NR2F2","NR3C1","NR4A1","NRDE2","NSFL1C","NSL1","NSMCE2","NSUN2","NTMT1",
          "NUBP1","NUDC","NUDT15","NUDT16","NUDT6","NUF2","NUGGC","NUMA1","NUP214","NUP37","NUP43","NUP62","NUPR1","NUPR2","NUSAP1","OBSL1","ODF2","OIP5","OOEP","OPN1LW","OPN1MW","OR1A2","OR2A4","ORC1","ORC4","OSGIN2","OVOL1","OVOL2","P3H4","PABIR1","PAF1","PAFAH1B1","PAGR1","PAK4","PARD3","PARD3B","PARD6A","PARD6B","PARD6G","PARP3","PARP9","PAX6","PAXIP1","PBK","PBRM1","PBX1","PCID2","PCLAF","PCM1","PCNA",
          "PCNP","PCNT","PDCD2L","PDCD6IP","PDE3A","PDE4DIP","PDGFB","PDGFRB","PDS5A","PDS5B","PDXP","PELO","PER2","PES1","PHACTR4","PHB2","PHF10","PHF13","PHF8","PHGDH","PHIP","PHOX2B","PIAS1","PIBF1","PIDD1","PIK3C3","PIK3R4","PIM1","PIM2","PIM3","PIMREG","PIN1","PINX1","PIWIL1","PIWIL2","PIWIL3","PIWIL4","PKD1","PKD2","PKHD1","PKIA","PKMYT1","PKN2","PKP4","PLAGL1","PLCB1","PLCG2","PLD6","PLEC","PLK1","PLK2","PLK3","PLK4","PLK5","PLRG1","PLSCR1","PMF1-BGLAP","PMF1","PML","PNPT1","POC1B","POC5","POGZ","POLA1","POLDIP2","POLE","POU4F1","PPM1A","PPM1D",
          "PPM1G","PPME1","PPP1CA","PPP1CB","PPP1CC","PPP1R10","PPP1R12A","PPP1R13B","PPP1R15A","PPP1R1C","PPP1R35","PPP1R9B","PPP2CA","PPP2CB","PPP2R1A","PPP2R2D","PPP2R3B","PPP2R5B","PPP2R5C","PPP3CA","PPP5C","PPP6C","PRAP1","PRC1","PRCC","PRDM11","PRDM5","PRDM7","PRDM9","PRICKLE1","PRKACA","PRKCA","PRKCB","PRKCD","PRKCE","PRKDC","PRMT2","PRMT5","PRNP","PROX1","PRPF19","PRPF40A","PRR11","PRR19","PRR5","PSMA8","PSMC3IP","PSMD10","PSMD13","PSME1","PSME2","PSME3","PSMG2","PSRC1","PTCH1","PTEN","PTGS2","PTK6","PTP4A1","PTPA","PTPN11","PTPN3","PTPN6","PTPRC","PTPRK","PTTG1","PTTG2","PUM1","PUM2","PYHIN1","RAB11A","RAB11FIP3","RAB11FIP4","RAB35","RAB6C","RABGAP1","RACGAP1","RACK1","RAD1","RAD17","RAD18","RAD21","RAD21L1","RAD23A","RAD50","RAD51","RAD51AP1","RAD51B","RAD51C","RAD51D",
          "RAD54B","RAD54L","RAD9A","RAD9B","RAE1","RALA","RALB","RAN","RANBP1","RARA","RASA1","RASSF1","RASSF2","RASSF4","RB1","RB1CC1","RBBP4","RBBP8","RBL1","RBL2","RBM14","RBM38","RBM46","RBM7","RCBTB1","RCC1","RCC2","RDX","REC114","REC8","RECQL5","REEP3","REEP4","RFPL1","RFWD3","RGCC","RGS14","RGS2","RHEB","RHNO1","RHOA","RHOB","RHOC","RHOU","RIF1","RINT1","RIOK2","RIOK3","RIPOR2","RMDN1","RMI1","RMI2","RNASEH2B",
          "RNF103-CHMP3","RNF112","RNF167","RNF2","RNF212","RNF212B","RNF4","RNF8","ROCK1","ROCK2","ROPN1B","RPA1","RPA2","RPA3","RPA4","RPL10L","RPL23","RPL24","RPL26","RPRD1B","RPRM","RPS15A","RPS27L","RPS3","RPS6KA1","RPS6KA2","RPS6KA3","RPS6KB1","RPTOR","RRM1","RRM2","RRM2B","RRP8","RRS1","RSPH1","RTEL1","RTF2","RTKN","RTTN","RUNX3","RUVBL1","RUVBL2","RXFP3","SAPCD2","SASS6","SBDS","SCAND3","SCRIB","SDCBP","SDCCAG8","SDE2","SEH1L","SENP5","SENP6","SEPTIN1","SEPTIN10","SEPTIN11","SEPTIN12","SEPTIN14",
          "SEPTIN2","SEPTIN3","SEPTIN4","SEPTIN5","SEPTIN6","SEPTIN7","SEPTIN8","SEPTIN9","SERTAD1","SETD2","SETDB2","SETMAR","SFN","SFRP1","SGF29","SGO1","SGO2","SGSM3","SH2B1","SH3GLB1","SHB","SHCBP1L","SHOC1","SIAH1","SIAH2","SIK1","SIN3A","SIPA1","SIRT1","SIRT2","SIRT7","SIX3","SKA1","SKA2","SKA3","SKIL","SKP2","SLC16A1","SLC25A31","SLC25A5","SLC26A8","SLC39A5","SLC6A4","SLF1","SLF2","SLFN11","SLX4","
          SMARCA2","SMARCA4","SMARCA5","SMARCAD1","SMARCB1","SMARCC1","SMARCC2","SMARCD1","SMARCD2","SMARCD3","SMARCE1","SMC1A","SMC1B","SMC2","SMC3","SMC4","SMC5","SMC6","SMIM22","SMOC2","SMPD3","SND1","SNRK","SNX18","SNX33","SNX9","SOGA1","SON","SOX15","SOX2","SOX9","SPAG5","SPAG8","SPAST","SPATA17","SPATA22","SPC24","SPC25","SPDL1","SPDYA","SPDYC","SPDYE1","SPDYE10","SPDYE11","SPDYE12","SPDYE13","SPDYE14","SPDYE15","SPDYE16","SPDYE17","SPDYE18","SPDYE2","SPDYE21","SPDYE2B","SPDYE3","SPDYE4","SPDYE5","SPDYE6","SPDYE8","SPDYE9","SPECC1L-ADORA2A","SPECC1L","SPHK1","SPICE1","SPIN1","SPIN2A","SPIN2B","SPIRE1","SPIRE2","SPO11","SPOUT1","SPRY1","SPRY2","SPTBN1","SRA1","SRC","SRPK1","SRPK2","SSNA1","SSTR5","SSX2IP","STAG1","STAG2","STAG3","STAMBP","STARD9","STAT3","STAT5B","STEAP3","STIL","STK10","STK11","STK33","STMN1","STOX1","STRA8","STRADA","STRADB","STXBP4","SUGT1","SUN1","SUN2","SUSD2","SUV39H1","SUV39H2","SVIL","SYCE1","SYCE1L","SYCE2","SYCE3","SYCP1","SYCP2","SYCP2L","SYCP3","SYF2","TACC1","TACC2","TACC3","TADA2A","TADA3","TAF1","TAF10","TAF1L","TAF2","TAF6","TAL1","TAOK1","TAOK2","TAOK3","TARDBP","TAS1R2","TAS2R13","TBCD","TBCE","TBRG1","TBRG4","TBX1","TBX2","TBX20","TBX3","TCF3","TCIM","TDRD1","TDRD12","TDRD9","TDRKH","TELO2","TENT4A","TENT4B","TENT5B","TERB1","TERB2","TERF1","TERF2","TERT","TESMIN","TET2","TEX11","TEX12","TEX14","TEX15","TEX19","TFAP4","TFDP1","TFDP2",
          "TFDP3","TFPT","TGFA","TGFB1","TGFB2","TGFBR1","TGM1","THAP1","THAP5","THAP9","THOC1","THOC5","TICRR","TIMELESS","TIPIN","TIPRL","TK1","TLE6","TLK1","TLK2","TM4SF5","TMEM14B","TMEM250","TMEM67","TMEM8B","TMOD3","TMPRSS11A","TNF","TNKS","TOM1L1","TOM1L2","TOP1","TOP2A","TOP2B","TOP3A","TOP3B","TOPBP1","TP53","TP53BP1","TP53BP2","TP53I13","TP53INP1","TP63","TP73","TPD52L1","TPPP","TPR","TPRA1","TPX2","TRAPPC12","TREX1","TRIAP1","TRIM21","TRIM32","TRIM35","TRIM36","TRIM37","TRIM39","TRIM71","TRIM75","TRIOBP","TRIP13","TRNP1","TRRAP","TSC1","TSC2","TSC22D2","TSG101","TSPYL2","TSSK4","TTBK1","TTC19","TTC28","TTI1","TTI2","TTK","TTL","TTLL12","TTN","TTYH1","TUBA1A","TUBA1B","TUBA1C","TUBA3C","TUBA3D","TUBA3E","TUBA4A","TUBA4B","TUBA8","TUBAL3","TUBB","TUBB1","TUBB2A","TUBB2B","TUBB3","TUBB4A","TUBB4B","TUBB6","TUBB8","TUBB8B","TUBD1","TUBE1",
          "TUBG1","TUBG2","TUBGCP2","TUBGCP3","TUBGCP4","TUBGCP5","TUBGCP6","TUSC2","TXLNG","TXNIP","TXNL4A","TXNL4B","UBA3","UBB","UBD","UBE2A","UBE2B","UBE2C","UBE2E2","UBE2I","UBE2L3","UBE2S","UBR2","UBXN2B","UCHL5","UFL1","UHMK1","UHRF1","UHRF2","UIMC1","UNC119","UPF1","URGCP","USH1C","USP16","USP17L2","USP19","USP2","USP22","USP26","USP28","USP29","USP3","USP33","USP37","USP39","USP44","USP47","USP51","USP8","USP9X","UTP14C","UVRAG","UXT","VASH1","VCP","VPS4A","VPS4B","VPS72","VRK1","WAC","WAPL","WASHC5","WASL","WBP2NL","WDHD1","WDR12","WDR5","WDR6","WDR62","WDR76","WEE1","WEE2","WIZ","WNK1","WNT10B","WNT4","WNT5A","WRAP73","WRN","WTAP","XIAP","XPC","XPO1","XRCC2","XRCC3","YEATS2","YEATS4","YTHDC2","YTHDF2","YWHAE","YY1","YY1AP1","ZBED3","ZBTB17","ZBTB49","ZC3H12D","ZC3HC1","ZCWPW1","ZFP36L1","ZFP36L2","ZFP42","ZFYVE19","ZFYVE26","ZMPSTE24","ZMYND11","ZNF16","ZNF207","ZNF268","ZNF318","ZNF324","ZNF503","ZNF541","ZNF655","ZNF703","ZNF830","ZNRD2","ZNRF4","ZPR1","ZSCAN21","ZW10","ZWILCH","ZWINT","ZZZ3",
          "KI67","MKI67")




plt <- stats.time |> 
  dplyr::rename_with( ~ paste0(.x, "_time"), .cols=!matches("^protein_id$",perl = T)) |> 
  dplyr::full_join(stats.grade |> dplyr::rename_with( ~ paste0(.x, "_grade"), .cols=!matches("^protein_id$",perl = T)), by=c('protein_id'='protein_id')) |> 
  dplyr::full_join(stats.cgc |> dplyr::rename_with( ~ paste0(.x, "_cgc"), .cols=!matches("^protein_id$",perl = T)), by=c('protein_id'='protein_id')) |> 

  dplyr::full_join(stats.time.pat.corrected |> tibble::rownames_to_column("protein_id") |> dplyr::rename_with( ~ paste0(.x, "_time.cor"), .cols=!matches("^protein_id$",perl = T)), by=c('protein_id'='protein_id')) |> 
  dplyr::full_join(stats.grade.pat.corrected |> tibble::rownames_to_column("protein_id") |> dplyr::rename_with( ~ paste0(.x, "_grade.cor"), .cols=!matches("^protein_id$",perl = T)), by=c('protein_id'='protein_id')) |>  

  dplyr::mutate(cellcycling_go_ = protein_id %in% clng) |> 
  #dplyr::mutate(fn1 = grepl("FN1|FN|CIG|FINC|FIBRO", protein_id))
  dplyr::mutate(fn1 = grepl("FINC", protein_id))


#plt |> dplyr::filter(fn1) |> head()




plt <- plt |> dplyr::mutate(label = cellcycling_go_)
ggplot(plt, aes(x=t_time, y=t_grade, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  theme_nature



plt <- plt |> dplyr::mutate(label = cellcycling_go_)
ggplot(plt, aes(x=t_time, y=t_time.cor, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  theme_nature


plt <- plt |> dplyr::mutate(label = cellcycling_go_)
ggplot(plt, aes(x=t_grade, y=t_grade.cor, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  theme_nature



pltcor <- plt |> 
  dplyr::select(starts_with("t")) |> 
  dplyr::filter(!is.na(t_time)) |> 
  dplyr::filter(!is.na(t_grade)) |> 
  dplyr::filter(!is.na(t_cgc))|> 
  dplyr::filter(!is.na(t_time.cor)) |> 
  dplyr::filter(!is.na(t_grade.cor))
  #dplyr::filter(!is.na(t_cgc.cor))

corrplot::corrplot(cor(pltcor), order="hclust")


sum((plt |> dplyr::filter(!is.na(adj.P.Val_time)) |> dplyr::pull(adj.P.Val_time)) < 0.01)
sum((plt |> dplyr::filter(!is.na(adj.P.Val_time.cor)) |> dplyr::pull(adj.P.Val_time.cor)) < 0.01)
sum((plt |> dplyr::filter(!is.na(adj.P.Val_grade)) |> dplyr::pull(adj.P.Val_grade)) < 0.01)
sum((plt |> dplyr::filter(!is.na(adj.P.Val_grade.cor)) |> dplyr::pull(adj.P.Val_grade.cor)) < 0.01)
sum((plt |> dplyr::filter(!is.na(adj.P.Val_cgc)) |> dplyr::pull(adj.P.Val_cgc)) < 0.01)
sum((plt |> dplyr::filter(!is.na(adj.P.Val_cgc.cor)) |> dplyr::pull(adj.P.Val_cgc.cor)) < 0.01)




plt <- plt |> dplyr::mutate(label = protein_id %in% c("IDH1"))
p1 = ggplot(plt, aes(x=t, y=t.grade, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  ggrepel::geom_text_repel(data=subset(plt,  label == T), size=3, col="red") + 
  labs(x = "t-score CGC (continuous)", y="t-score WHO grade 2 - 3") + 
  theme_nature

plt <- plt |> dplyr::mutate(label = protein_id %in% c("IDH2"))
p2 = ggplot(plt, aes(x=t, y=t.grade, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  ggrepel::geom_text_repel(data=subset(plt,  label == T), size=3, col="red")  + 
  labs(x = "t-score CGC (continuous)", y="t-score WHO grade 2 - 3") + 
  theme_nature

plt <- plt |> dplyr::mutate(label = protein_id %in% c("GLUD1"))
p3 = ggplot(plt, aes(x=t, y=t.grade, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  ggrepel::geom_text_repel(data=subset(plt,  label == T), size=3, col="red") + 
  labs(x = "t-score CGC (continuous)", y="t-score WHO grade 2 - 3") +
  theme_nature


plt <- plt |> dplyr::mutate(label = protein_id %in% c("FINC"))
p4 = ggplot(plt, aes(x=t, y=t.grade, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  ggrepel::geom_text_repel(data=subset(plt,  label == T), size=3, col="red") + 
  labs(x = "t-score CGC (continuous)", y="t-score WHO grade 2 - 3") +
  theme_nature




p1 + p2 + p3 + p4

### corr t's ----


plt.c <- plt |> 
  dplyr::select(starts_with("t_")) |> 
  dplyr::filter(!is.na(t_time)) |> 
  dplyr::mutate(t_time = NULL ) |> 
  dplyr::mutate(t_grade = NULL) |> 
  cor()


corrplot::corrplot(plt.c, order="hclust") 





### FN1 ----



ggplot(plt, aes(x=t_cgc, y=t_time.cor, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_point(data=subset(plt, fn1 == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, fn1 == T), pch=19, cex=0.6, col="red") +
  ggrepel::geom_text_repel(data=subset(plt,  fn1 == T), size=3, col="red") + 
  labs(x = "t-score CGC (continuous)", y="t-score WHO grade 2 - 3") +
  theme_nature



#### top100 -----


top100 <- plt |> 
  dplyr::filter(!is.na(protein_id)) |> 
  dplyr::filter(adj.P.Val_cgc.cor < 0.01) |> 
  dplyr::arrange(abs(adj.P.Val_cgc.cor)) |> 
  dplyr::select(protein_id, adj.P.Val_cgc.cor) |> 
  head(n=100) |> 
  dplyr::pull(protein_id) |> 
  unique()


plt <- plt |> dplyr::mutate(label = protein_id %in% top100)
ggplot(plt, aes(x=t_time.cor, y=t_cgc.cor, label= protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.15) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  ggrepel::geom_text_repel(data=subset(plt,  label == T), size=3, col="red") + 
  labs(x = "t-score CGC (continuous)", y="t-score WHO grade 2 - 3") +
  theme_nature


#### enhanced volcano ----



EnhancedVolcano(plt,
                lab = plt$protein_id,
                FCcutoff = 0.5,
                x = 'logFC_cgc.cor',
                y = 'adj.P.Val_cgc.cor')




plt <- plt |> dplyr::mutate(significant = adj.P.Val_cgc.cor < 0.01 & abs(logFC_cgc.cor) > lfc_cut)
ggplot(plt, aes(x= `logFC_cgc.cor`, y= -log(`adj.P.Val_cgc.cor`), col=significant, label=protein_id)) +
  
  geom_hline(yintercept=0,          col="black", lwd=theme_nature_lwd, lty=1) +
  geom_hline(yintercept=-log(0.01), col="black", lwd=theme_nature_lwd, lty=2) +
  #geom_hline(yintercept=-log(0.05), col="gray", lwd=theme_nature_lwd) +
  
  geom_vline(xintercept=0,        col="black", lwd=theme_nature_lwd, lty=2) +
  geom_vline(xintercept=-lfc_cut, col="black", lwd=theme_nature_lwd, lty=2) +
  geom_vline(xintercept= lfc_cut, col="black", lwd=theme_nature_lwd, lty=2) +
  
  geom_point(pch=19, cex=0.5) +
  ggrepel::geom_text_repel(data=subset(plt, significant==T), size=theme_nature_size, col="black") +
  scale_color_manual(values=c(`TRUE`='red',`FALSE`='darkgray')) +
  theme_nature




lfc_cut <- log2(1.25)
plt <- plt |> dplyr::mutate(significant = adj.P.Val_cgc.cor < 0.01 & abs(logFC_cgc.cor) > lfc_cut)
ggplot(plt, aes(x= `logFC_cgc.cor`, y= -log(`adj.P.Val_cgc.cor`), col=significant, label=protein_id)) +
  
  geom_hline(yintercept=0,          col="black", lwd=theme_nature_lwd, lty=1) +
  geom_hline(yintercept=-log(0.01), col="black", lwd=theme_nature_lwd, lty=2) +
  #geom_hline(yintercept=-log(0.05), col="gray", lwd=theme_nature_lwd) +
  
  geom_vline(xintercept=0,        col="black", lwd=theme_nature_lwd, lty=2) +
  geom_vline(xintercept=-lfc_cut, col="black", lwd=theme_nature_lwd, lty=2) +
  geom_vline(xintercept= lfc_cut, col="black", lwd=theme_nature_lwd, lty=2) +
  
  geom_point(pch=19, cex=0.5) +
  ggrepel::geom_text_repel(data=subset(plt, protein_id %in% c("GLUD2", "IDH1","IDH2", "GLUL", "GLUD1", "ALDH18A1", "GLS", "PC", "GLS2", "GOT1L1", "GOT1", "CPS1", "GOT2",
                                                              "IDH3B",  "IDH3A", "IDH3G", "LDHA", "LDHB", "OGDH")), size=theme_nature_size, col="black") +
  scale_color_manual(values=c(`TRUE`='red',`FALSE`='darkgray')) +
  theme_nature





lfc_cut <- log2(1.25)
plt <- plt |> dplyr::mutate(significant = adj.P.Val_cgc.cor < 0.01 & abs(logFC_cgc.cor) > lfc_cut)
ggplot(plt, aes(x= `logFC_cgc.cor`, y= -log(`adj.P.Val_cgc.cor`), col=significant, label=protein_id)) +
  
  geom_hline(yintercept=0,          col="black", lwd=theme_nature_lwd, lty=1) +
  geom_hline(yintercept=-log(0.01), col="black", lwd=theme_nature_lwd, lty=2) +
  #geom_hline(yintercept=-log(0.05), col="gray", lwd=theme_nature_lwd) +
  
  geom_vline(xintercept=0,        col="black", lwd=theme_nature_lwd, lty=2) +
  geom_vline(xintercept=-lfc_cut, col="black", lwd=theme_nature_lwd, lty=2) +
  geom_vline(xintercept= lfc_cut, col="black", lwd=theme_nature_lwd, lty=2) +
  
  geom_point(pch=19, cex=0.5) +
  ggrepel::geom_text_repel(data=subset(plt, protein_id %in% c("DNMT1",  "DNMT3A", "APOBEC3C", "APOBEC3B", "APOBEC2"
                                                              # no TET's
                                                              # no TDG's
                                                              )), size=theme_nature_size, col="black") +
  scale_color_manual(values=c(`TRUE`='red',`FALSE`='darkgray')) +
  theme_nature




plt |> dplyr::filter(grepl("BER", protein_id)) |> dplyr::pull(protein_id)



plt |> dplyr::filter(significant) |> dplyr::arrange(adj.P.Val_cgc.cor) |>  dplyr::pull(protein_id)





#### PCNA / Ki67 (MKI67) ----




plt |> dplyr::filter(grepl("^HMGB", protein_id)) |> dplyr::pull(protein_id)



lfc_cut <- log2(1.25)
plt <- plt |> dplyr::mutate(significant = adj.P.Val_cgc.cor < 0.01 & abs(logFC_cgc.cor) > lfc_cut)
ggplot(plt, aes(x= `logFC_cgc.cor`, y= -log(`adj.P.Val_cgc.cor`), col=significant, label=protein_id)) +
  
  geom_hline(yintercept=0,          col="black", lwd=theme_nature_lwd, lty=1) +
  geom_hline(yintercept=-log(0.01), col="black", lwd=theme_nature_lwd, lty=2) +
  #geom_hline(yintercept=-log(0.05), col="gray", lwd=theme_nature_lwd) +
  
  geom_vline(xintercept=0,        col="black", lwd=theme_nature_lwd, lty=2) +
  geom_vline(xintercept=-lfc_cut, col="black", lwd=theme_nature_lwd, lty=2) +
  geom_vline(xintercept= lfc_cut, col="black", lwd=theme_nature_lwd, lty=2) +
  
  geom_point(pch=19, cex=0.5) +
  ggrepel::geom_text_repel(data=subset(plt, protein_id %in% c(
    "PCNA", "PNCA", "MKI67", "CDC6", "PLK1", "MCM3" ,  "MCM4" ,  "MCM5" ,  "MCM7"  , "MCM2" ,  "MCM6",
    
    "TMPO","DEK","RFC4",
    "HMGB3", "HMGB1", "HMGB2"
    
    #"H6PD", "H2AC25;H2AC6;H2AC8", "H2BC11;H2BC21", "H1-0", "H1-4", "H1-5", "H1-2", "H2BC5", "H4C16", "H3-3B;H3C12;H3C13", "H1-1", "H2AC19;H2AC20", "H3-7", "H2AZ2", "H2AC21", "H1-10", "H2BC1", "H1-6", "H2BC19P;H2BC20P"
    )), size=theme_nature_size, col="black") +
  scale_color_manual(values=c(`TRUE`='red',`FALSE`='darkgray')) +
  theme_nature



### recursive Cor N/A adj ----


lfc_cut <- 0.5

plt <- plt |>
  dplyr::mutate(significant = adj.P.Val_cgc < 0.01 & abs(logFC_cgc) > lfc_cut) 
table(plt$significant)


metadata <- glass_od.metadata.proteomics |> 
  dplyr::filter(is.na(proteomics_notes) | (proteomics_notes == "has replicate" & grepl("use this one first", proteomics_notes_zurich))) |> 
  dplyr::filter(is.na(resection_reason_excluded)) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  dplyr::filter(patient_study_name == "GLASS-OD") |> 
  dplyr::select(resection_id, starts_with("proteomics_")) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 140) 
    return(.)
  })()




data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id) |> 
  tibble::rownames_to_column("protein_id") |> 
  dplyr::filter(protein_id %in% (plt |> dplyr::filter(significant) |> dplyr::pull(protein_id))) |> 
  tibble::column_to_rownames("protein_id")




dim(data)



cordata <- metadata |> 
  dplyr::inner_join(
  glass_od.metadata.array_samples |> 
    filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
    dplyr::filter(resection_id %in% metadata$resection_id) |> 
    assertr::verify(!duplicated(resection_id)),
  by=c('resection_id'='resection_id')
) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 118) 
    return(.)
  })() |> 
  dplyr::select(resection_id, proteomics_id, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::left_join(
    as.data.frame(t(data)) |> tibble::rownames_to_column('proteomics_id'), by=c('proteomics_id'='proteomics_id')
  ) |> 
  tibble::column_to_rownames('resection_id') |> 
  dplyr::mutate(proteomics_id = NULL)


cordata <- cor(cordata, use="pairwise.complete.obs", method="pearson") |> 
  as.data.frame() |> 
  dplyr::select(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  tibble::rownames_to_column('protein_id') |> 
  dplyr::filter(protein_id != "array_A_IDH_HG__A_IDH_LG_lr__lasso_fit") |> 
  dplyr::rename(`cor CGC` = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)




labels <- data |>
  tibble::rownames_to_column("gene_id") |> 
  dplyr::select("gene_id") |> 
  dplyr::mutate(FN1 = gene_id == "FINC") |> 
  dplyr::mutate(cycling = gene_id %in% c(
    "PCNA", "PNCA", "MKI67", "CDC6", "PLK1", "MCM3" ,  "MCM4" ,  "MCM5" ,  "MCM7"  , "MCM2" ,  "MCM6",
    "TMPO","DEK","RFC4",
    "HMGB3", "HMGB1", "HMGB2"
  )) |> 
  dplyr::mutate(GO_0002526_acute_inflammatory_response = gene_id %in% c(
   "BTK","CD6","PTGER3","ITIH4","CREB3L3","TFRC","FCGR2B","PTGS2","GSTP1","B4GALT1","EIF2AK1","OSM","RHBDD3","PIK3CG","EPHB6","TFR2","LIPA","VNN1","IL4","IL1A","ACVR1","FN1","PARK7","ASH1L","PLA2G2D","F3","TNFSF4","CNR1","TNFSF11","IL1B","C3","FFAR2","CCR7","IL22","APOL2","LBP","EPO","ASS1","F12","SELENOS","PPARG","CRP","APCS","ALOX5AP","SAA2","IL6ST","EDNRB","IL6","IL1RN","PRCP","MYLK3","TNFRSF11A","REG3G","AHSG","OSMR","PTGES","SAA4","FCGR1A","ADAM8","IL6R","NLRP3","ADORA1","DNASE1L3","NPY5R","KLKB1","IL31RA","MBL2","SERPINF2","SCN11A","MRGPRX1","REG3A","CEBPB","SAA1","IL20RB","NLRP6","A2M","NUPR1","ANO6","CD163","CTNNBIP1","FCER1A","F2","FUT7","CXCR2","EXT1","F8","SIGIRR","FFAR3","PLSCR1","ZP3","SERPINA3","TRPV1","SERPINA1","SPN","ELANE","C2CD4A","FCGR3A","HLA-E","C2CD4B","IGHG1","DNASE1","ORM2","ORM1","TNF","UGT1A1","INS","SAA2-SAA4","HP"
  )) |>
  dplyr::mutate(GO_0062023_collagen_containing_extracellular_matrix = gene_id %in% c(
    "DCN","SEMA3B","MARCO","SERPINB1","MYOC","TIMP2","VCAN","TNC","USH2A","COL9A2","LTBP1","ELN","LAMC3","COL23A1","LAMA3","ENTPD2","CBLN4","ITIH4","ITIH1","F7","LAMC2","COL11A1","WNT8A","GPC1","CDON","NTN1","COL17A1","FGFR2","PKM","FGF10","SPP2","ADAM11","NTN4","LMAN1","DLG1","GPC4","FBLN1","COL5A3","HSP90AA1","COL4A4","IMPG2","COL19A1","COL16A1","FCN1","ACHE","ADAMTS2","MMP2","NID2","LTBP4","ICAM1","P3H2","LAMB4","LAMB1","APOH","ANGPT2","CCDC80","CMA1","COL9A3","TGFB2","AMELY","HNRNPM","LGALS1","TIMP3","PDGFB","CHADL","CTSG","COCH","MMP9","BMP7","CTSZ","COL20A1","ANGPT4","LAMA1","MXRA5","SRPX","F9","TIMP1","SRPX2","CBLN1","ZP2","CTSH","SFRP1","IL7","FGL1","CLC","APLP1","TGFB1","COMP","WNT2","PTPRZ1","PCOLCE","SERPINE1","PLOD3","RARRES2","AEBP1","MEGF9","OGN",
    "ASPN","ECM2","AMBP","CXCL12","ACTA2","KAZALD1","LGALS3BP","COL1A1","VTN","SOD3","CTSC","HPX","APOA4","APOC3","MDK","VWF","WNT5B","MGP","COL12A1","COL9A1","SMOC2","IMPG1","LAMA4","ERBIN","LOX","SPARC","THBS4","KNG1","HRG","WNT5A","COL7A1","LOXL3","EFEMP1","FN1","TNR","ANGPTL1","PRG4","NID1","ZP4","MFAP2","P3H1","F3","SERPINC1","CTSD","MMP8","APOA1","CCN2","LTBP2","TGFB3","TNN","TGFBI",
    "CLU","A1BG","FMOD","PLG","ANXA11","COL10A1","SLPI","MATN4","NDP","F13A1","SERPINB6","COL21A1","AMELX","CFP","PZP","FGL2","LOXL1","APOE","NCAN","PXDN","GDF15","COL5A1","LAMA5","FIBCD1","ANGPTL6","F12","LGALS3","PODNL1","RTBDN","MATN3","EMILIN2","SERPINF1","ITGB4","MATN2","BCAN","HAPLN2","APCS","ANGPTL3","POSTN","LOXL2","ADAMDEC1","WNT2B","COL4A2","ADAMTS8","ANXA1","CTSL","ADAM19","AGT","LAMC1","SERPINE2","ANGPTL2","CCN3","IGFBPL1","TINAG","SULF1","THBS1",
    "EMILIN1","LOXL4","ANXA7","CILP","SEMA7A","MMRN1","FRAS1","FBN2","COL2A1","INHBE","LUM","FBLN5","PCSK6","LMAN1L","HAPLN3","MFGE8","TGFB1I1","CDH13","COL6A1","COL6A2","ADAMTS10","FCN3","HSPG2","CCN1","TINAGL1","DPT","HDGF","HMCN1","ECM1","ADAMTSL4","ANXA9","S100A8","S100A7","FLG","LEFTY2","COL8A1","AHSG","SFRP2","HAPLN1","CASK","GPC3","HMCN2","SERPING1","SERPINH1","NCAM1","ZP1","FREM2","ITIH2","DST","SPARCL1","ABI3BP","ANGPT1","ADAMTS1","ADAMTS5","ADAMTS3","PRG3","ACAN","COLEC12","ADAMTS4","LAD1","C1QC","CSTB","FCN2","AZGP1","COL26A1","MATN1","SDC3","CTSS","S100A9","COL6A3","IGFBP7","FBLN2","ADAMTS9","CRELD1","PF4","CPA3","TGM4","CLEC3B","ANXA5","EDIL3","EGFLAM","RELL2","SHH","COL1A2","CTSB","SBSPON","CTHRC1","FREM1","VWA2","HTRA1","ADAMTS15","FBN1","SERPINB8","MFAP4","HSP90B1","SERPINB12","NAV2","GREM1","SERPINF2","KRT1","ANGPTL4",
    "SOST","LTBP3","TNXB","COL3A1","NPNT","LGALS9","COL4A3","THBS3","COL22A1","SDC2","MUC17","SERPINB9","CDH2","COL24A1","FGG","FGA","FGB","LGALS4","COL8A2","ANGPTL7","LAMB2","LRRC15","TPSAB1","EFEMP2","COL6A5","EGFL7","ADAMTS20","MMRN2","C1QA","DAG1","CSPG4","CTSF","PODN","FGFBP3","ZG16","NPPA","A2M","CLEC14A","CD151","CALR","GPC5","VWA1","SSC5D","F2","ADIPOQ","OTOL1","BGN","MXRA7","SNORC","ANXA2","COL18A1","FREM3","GPC6","EMILIN3","EFNA5","THBS2","PRG2","C17orf58","ANGPTL5","BCAM","COL4A1","HAPLN4","AMTN","THSD4","COL14A1","EYS","COL4A5","AGRN","PLSCR1","ZP3","IFNA2","SERPINA5","COL25A1","VWC2","PRELP","MMP23B","SERPINA3","S100A4","PRTN3","LAMA2","COL27A1","LAMB3","ANXA4","ANXA6","SERPINA1","TPSB2","COL13A1","SPN","ELANE","COL4A6","MFAP5","PSAP","S100A10","HRNR","S100A6","SMOC1","EGFL6","L1CAM","TGM2","COL11A2","COL5A2","COL15A1","PRSS1","VIT","DEFA1","COL6A6","COLQ","GPC2","ANG","COL28A1","ORM2","ORM1","DEFA1B","MARCOL","GH1","SPON1","ANXA8","RBP3","GDF10","MMP28"
  )) |> 
  dplyr::mutate(`GO_0072562_blood_microparticle` = gene_id %in% c(
    "CFH","SLC4A1","PON1","ITGA2B","CP","ITIH4","ITIH1","GRIPAP1","TFRC","CD5L","ACTB","AFM","PSMC5","TF","APOL1","TGFB1","AMBP","ENG","PFN1","LGALS3BP","VTN","HSPA8","HPX","APOA4","C9","KNG1","HRG","BCHE","FN1","CFHR3","SLC2A1","SERPINC1","APOA1","CLU","A1BG","PLG","C4BPA","DNPEP","F13A1","C3","HSPA2","PZP","INTS11","APOE","JCHAIN","APCS","AGT","CIB2","SDCBP","TMPRSS13","FCN3","OAZ3","ACTA1","EIF2A","AHSG","GC","MSN","STOM","GSN","SERPING1","ITIH2","C8A","APOA2","C1QC","ACTC1","FCN2","ACTG2","ALB","ANXA5","YWHAZ","ACSM1","SERPINF2","KRT1","ANGPTL4","FGG","FGA","FGB","HSPA6","A2M","C8G","ZBTB38","CPN2","F2","C1S","ACTG1","PROS1","KDM4D","POTEE","HBA2","ZNF177","SERPINA3","HBG2","POTEF","HSPA1B","HSPA1A","HSPA1L","PRSS1","HBA1","IGKV4-1","IGLV1-47","IGLV3-25","IGLV3-21","IGLC2","IGLC3","IGHA2","IGHG4","IGHG2","IGHA1","IGHG1","IGHM","IGHV3-7","IGHV3-13","IGHV3-23","CLIC1","HBE1","HBD","C4B","ORM2","ORM1","IGKV3-20","IGKV1D-33","IGKV1-17","IGKV3-11","IGKV1-33","IGKV1-39","IGKV2D-28","IGKV2-30","IGKV1-5","CFB","CFHR1","IGKV3-15","C4A","HBB","IGKV2D-40","HP","HPR","IGKV1D-12"
  )) |> 
  tibble::column_to_rownames("gene_id")


labels_continuous <- data |>
  tibble::rownames_to_column("gene_id") |> 
  dplyr::select("gene_id") |> 
  dplyr::mutate(fraction_NAs = rowSums(is.na(data))) |> 
  dplyr::left_join(cordata, by=c('gene_id'='protein_id')) |> 
  tibble::column_to_rownames("gene_id")




pplt <- data |>
  tibble::rownames_to_column('__hugo_symbol__') |>
  dplyr::filter(!duplicated(.data$`__hugo_symbol__`)) |>
  tibble::column_to_rownames('__hugo_symbol__')


pplt <- pplt |>
  base::as.matrix() |>
  base::t() |> 
  stats::cor(use="pairwise.complete.obs") # copes with N/A values



h <- stats::hclust(stats::as.dist(1 - stats::cor(pplt)), method = "ward.D2" ) # recursive cor-based cluastering !!!


o <- h$labels[h$order] |>
  base::rev()

#write.csv(data.frame(o), file="output/tables/proteomics_rcur.txt")




ph <- ggdendro::ggdendrogram(h, rotate = TRUE, theme_dendro = FALSE) +
  ggdendro::theme_dendro()


pplt <- pplt |>
  base::as.data.frame() |>
  dplyr::select(dplyr::all_of(o)) |>
  base::t() |>
  base::as.data.frame() |>
  dplyr::select(dplyr::all_of(o)) |>
  base::t() |>
  base::as.matrix()


o.join <- base::data.frame(name = o, i = 1:length(o))


plt.expanded2 <- reshape2::melt(pplt) |>
  dplyr::rename(y = .data$`Var1`) |>
  dplyr::rename(x = .data$`Var2`) |>
  dplyr::mutate(x = as.factor(.data$`x`)) |>
  dplyr::mutate(y = as.factor(.data$`y`)) |>
  dplyr::left_join(o.join |> dplyr::rename(x.order = .data$`i`), by = c("x" = "name")) |>
  dplyr::left_join(o.join |> dplyr::mutate(i = dplyr::n() - .data$i + 1) |> dplyr::rename(y.order = .data$i), by = c("y" = "name"))

base::rm(o.join)


font_scale = 4.5
legend_scale = 1.2


p1 <- ggplot2::ggplot(plt.expanded2, ggplot2::aes(
  x = .data$x.order,
  y = .data$y.order,
  radius = ((abs(.data$value) * 0.7) + 0.3) / 2 - 0.05,
  # [0.3 , 0.8] + 0.2 smoothened from lwd/border
  fill = .data$value,
  col = .data$value,
  label = .data$x
)
) +
  ggplot2::geom_tile(col = "gray", fill = "white", lwd = 0.15) +
  ggplot2::scale_fill_gradientn(colours = col2(200), na.value = "grey50", limits = c(-1, 1), guide = "none") + # guide = "colourbar",
  ggplot2::scale_color_gradientn(colours = col2(200), na.value = "grey50", limits = c(-1, 1), guide = "none") +
  recursiveCorPlot::geom_circle(radius.fixed = T) + # from THIS repo
  ggplot2::scale_x_discrete(labels = NULL, breaks = NULL) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.y = ggplot2::element_text(size = font_scale, angle = 0, hjust = 1, vjust = 0.5), # used to be [3,6] reduce font size here, should become argument
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5, color = "gray80"),
    text = ggplot2::element_text(size = 13),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank()
  ) +
  ggplot2::labs(y = NULL, x = NULL, main = NULL) +
  ggplot2::coord_fixed() +
  ggplot2::scale_y_continuous(name = NULL, breaks = base::length(o):1, labels = o)





pplt <- base::data.frame(gid = o, i = 1:length(o)) |>
  dplyr::left_join(labels |> tibble::rownames_to_column("gid"), by = c("gid" = "gid")) |>
  reshape2::melt(id.vars = c("gid", "i")) |>
  dplyr::mutate(variable = factor(.data$variable, levels = base::rev(base::colnames(labels))))


p2 <- ggplot2::ggplot(pplt, ggplot2::aes(
  x = .data$i,
  y = .data$variable,
  fill = .data$value,
  label = .data$gid
)) +
  ggplot2::geom_tile(col = "white", lwd = 0.15) +
  # scale_x_discrete(position = "bottom")  +
  ggplot2::scale_x_discrete(labels = NULL, breaks = NULL) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5),
    axis.text.y = ggplot2::element_text(size = font_scale, angle = 0, hjust = 1, vjust = 0.5),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank()
  ) +
  ggplot2::guides(fill = "none") +
  ggplot2::coord_fixed(ratio = legend_scale) + # used to be 2.75
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray98"))




pplt <- base::data.frame(gid = o, i = 1:length(o)) |>
  dplyr::left_join(labels_continuous |> tibble::rownames_to_column("gid"), by = c("gid" = "gid")) |>
  dplyr::mutate(fraction_NAs = 1.0 * fraction_NAs) |> 
  reshape2::melt(id.vars = c("gid", "i")) |>
  dplyr::mutate(variable = factor(.data$variable, levels = base::rev(base::colnames(labels_continuous))))

pplt <- rbind(
  pplt,
  pplt |> dplyr::mutate(value=0)
)



p3 <- ggplot2::ggplot(pplt, ggplot2::aes(
  x = reorder(gid, .data$i),
  y = .data$value,
  #fill = .data$value,
  label = .data$gid
)) +
  ggplot2::geom_line(col = "black", lwd = 0.15) +
  facet_grid(rows = vars(.data$variable), scales = "free") + 
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    axis.ticks.x=element_blank(), 
    axis.ticks.y=element_blank()
  ) +
  ggplot2::guides(fill = "none") +
  ggplot2::labs(x = NULL, y = NULL) 




patchwork::wrap_plots(
  A = p3,
  B = p2,
  C = p1
  #D = (ph + patchwork::plot_spacer()),
#   design = 'A
# B
# C'
  ) + 
  patchwork::plot_annotation(caption = 'caption') +
  patchwork::plot_layout(ncol=1)


ggsave("output/figures/vis_differential_proteomics_corrplot.pdf", width=8.4 , height=8.4*2, dpi=600)



#### GSEA ----


c("A","B;C")



write.table(
  plt |>
    dplyr::filter(!is.na(adj.P.Val_cgc)) |>
    dplyr::pull(protein_id) |> 
    stringr::str_split(";", n = Inf, simplify = FALSE) |> 
    unlist()
  ,
  "output/tables/GSEA_background_proteomics_cgc_cor.txt",  quote = F, row.names = F, col.names=F)



write.table(
  pplt |> 
    dplyr::filter(variable == 'cor CGC' & value != 0) |> 
    dplyr::filter(value > 0) |> 
    dplyr::pull(gid) |> 
    stringr::str_split(";", n = Inf, simplify = FALSE) |> 
    unlist()
  ,
  "output/tables/GSEA_proteins_up_proteomics_cgc_cor.txt",  quote = F, row.names = F, col.names=F)




write.table(
  pplt |> 
    dplyr::filter(variable == 'cor CGC' & value != 0) |> 
    dplyr::filter(value < 0) |> 
    dplyr::pull(gid) |> 
    stringr::str_split(";", n = Inf, simplify = FALSE) |> 
    unlist()
  ,
  "output/tables/GSEA_proteins_down_proteomics_cgc_cor.txt",  quote = F, row.names = F, col.names=F)





