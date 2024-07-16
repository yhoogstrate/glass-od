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



lfc_cut <- log2(1.25)



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








## WHO Grade ----
### regular ----


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
    assertthat::assert_that(nrow(.) == 119) 
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




### pat corrected ----




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
    assertthat::assert_that(nrow(.) == 119) 
    return(.)
  })() |> 
  
  dplyr::mutate(CGC = scale(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit))




data <- data.proteomics.glass_od |> 
  dplyr::select(metadata$proteomics_id)



stats.cgc.pat.corrected <- data |> 
  tibble::rownames_to_column("protein_id") |> 
  dplyr::select("protein_id") |> 
  dplyr::mutate(tmp = pbapply::pbapply(data, 1, pp, metadata$patient_id, metadata$CGC, "condition")) |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(adj.P.Val = NULL) |> 
  dplyr::mutate(adj.P.Val = p.adjust(P.Value, "fdr")) |> 
  tibble::column_to_rownames("protein_id")


sum((stats.cgc.pat.corrected |> dplyr::filter(!is.na(P.Value)) |> dplyr::pull(P.Value) ) < 0.01)
sum((stats.cgc.pat.corrected |> dplyr::filter(!is.na(adj.P.Val)) |> dplyr::pull(adj.P.Val) ) < 0.01)



rm(metadata, data)





## PC2 (methylation) ----


## plt tmp ----




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
          "TUBG1","TUBG2","TUBGCP2","TUBGCP3","TUBGCP4","TUBGCP5","TUBGCP6","TUSC2","TXLNG","TXNIP","TXNL4A","TXNL4B","UBA3","UBB","UBD","UBE2A","UBE2B","UBE2C","UBE2E2","UBE2I","UBE2L3","UBE2S","UBR2","UBXN2B","UCHL5","UFL1","UHMK1","UHRF1","UHRF2","UIMC1","UNC119","UPF1","URGCP","USH1C","USP16","USP17L2","USP19","USP2","USP22","USP26","USP28","USP29","USP3","USP33","USP37","USP39","USP44","USP47","USP51","USP8","USP9X","UTP14C","UVRAG","UXT","VASH1","VCP","VPS4A","VPS4B","VPS72","VRK1","WAC","WAPL","WASHC5","WASL","WBP2NL","WDHD1","WDR12","WDR5","WDR6","WDR62","WDR76","WEE1","WEE2","WIZ","WNK1","WNT10B","WNT4","WNT5A","WRAP73","WRN","WTAP","XIAP","XPC","XPO1","XRCC2","XRCC3","YEATS2","YEATS4","YTHDC2","YTHDF2","YWHAE","YY1","YY1AP1","ZBED3","ZBTB17","ZBTB49","ZC3H12D","ZC3HC1","ZCWPW1","ZFP36L1","ZFP36L2","ZFP42","ZFYVE19","ZFYVE26","ZMPSTE24","ZMYND11","ZNF16","ZNF207","ZNF268","ZNF318","ZNF324","ZNF503","ZNF541","ZNF655","ZNF703","ZNF830","ZNRD2","ZNRF4","ZPR1","ZSCAN21","ZW10","ZWILCH","ZWINT","ZZZ3")


dim(stats.time)

plt <- stats.time |>
  dplyr::rename_with( ~ paste0(.x, "_time"), .cols=!matches("^protein_id$",perl = T)) |> 
  dplyr::full_join(stats.grade |> dplyr::rename_with( ~ paste0(.x, "_grade"), .cols=!matches("^protein_id$",perl = T)), by=c('protein_id'='protein_id')) |> 
  dplyr::full_join(stats.cgc |> dplyr::rename_with( ~ paste0(.x, "_cgc"), .cols=!matches("^protein_id$",perl = T)), by=c('protein_id'='protein_id')) |> 

  dplyr::full_join(stats.time.pat.corrected |> tibble::rownames_to_column("protein_id") |> dplyr::rename_with( ~ paste0(.x, "_time.cor"), .cols=!matches("^protein_id$",perl = T)), by=c('protein_id'='protein_id')) |> 
  dplyr::full_join(stats.grade.pat.corrected |> tibble::rownames_to_column("protein_id") |> dplyr::rename_with( ~ paste0(.x, "_grade.cor"), .cols=!matches("^protein_id$",perl = T)), by=c('protein_id'='protein_id')) |> 
  dplyr::full_join(stats.cgc.pat.corrected |> tibble::rownames_to_column("protein_id") |> dplyr::rename_with( ~ paste0(.x, "_cgc.cor"), .cols=!matches("^protein_id$",perl = T)), by=c('protein_id'='protein_id')) |> 
  
  dplyr::mutate(cellcycling_go_ = protein_id %in% clng) 



background <- plt |>
  dplyr::filter(!is.na(adj.P.Val_cgc)) |>
  dplyr::pull(protein_id)
write.table(data.frame(background), "output/tables/GSEA_background_proteomics_cgc_cor.txt")


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
  dplyr::filter(!is.na(t_grade.cor)) |> 
  dplyr::filter(!is.na(t_cgc.cor))

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


plt <- plt |> dplyr::mutate(label = protein_id %in% c("FN1"))
p4 = ggplot(plt, aes(x=t, y=t.grade, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  ggrepel::geom_text_repel(data=subset(plt,  label == T), size=3, col="red") + 
  labs(x = "t-score CGC (continuous)", y="t-score WHO grade 2 - 3") +
  theme_nature




p1 + p2 + p3 + p4

#### FN1 ----


plt <- plt |> dplyr::mutate(label = protein_id %in% c("FN1"))
ggplot(plt, aes(x=t_cgc, y=t_cgc.cor, label=protein_id)) +
  geom_hline(yintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_vline(xintercept=0, col="darkgray", lwd=theme_nature_lwd) +
  geom_point(data=subset(plt, label == F), pch=19, cex=0.1, col="black", alpha=0.2) +
  geom_point(data=subset(plt, label == T), pch=19, cex=0.6, col="red") +
  ggrepel::geom_text_repel(data=subset(plt,  label == T), size=3, col="red") + 
  labs(x = "t-score CGC (continuous)", y="t-score WHO grade 2 - 3") +
  theme_nature



#### GSEA ----


plt |> 
  dplyr::filter(!is.na(protein_id)) |> 
  dplyr::filter(adj.P.Val_cgc.cor < 0.01) |> 
  dplyr::pull(protein_id) |> 
  unique()

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



#### recursive Cor N/A adj ----



plt <- plt |>
  dplyr::mutate(significant = adj.P.Val_cgc.cor < 0.01 & abs(logFC_cgc.cor) > lfc_cut)

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
    assertthat::assert_that(nrow(.) == 119) 
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
  dplyr::mutate(FN1 = gene_id == "FN1") |> 
  dplyr::mutate(cycling = gene_id %in% c(
    "PCNA", "PNCA", "MKI67", "CDC6", "PLK1", "MCM3" ,  "MCM4" ,  "MCM5" ,  "MCM7"  , "MCM2" ,  "MCM6",
    "TMPO","DEK","RFC4",
    "HMGB3", "HMGB1", "HMGB2"
  )) |> 
  dplyr::mutate(IDH = gene_id %in% c("GLUD2", "IDH1","IDH2", "GLUL", "GLUD1", "ALDH18A1", "GLS", "PC", "GLS2", "GOT1L1", "GOT1", "CPS1", "GOT2",
                                     "IDH3B",  "IDH3A", "IDH3G", "LDHA", "LDHB", "OGDH")) |> 
  dplyr::mutate(GO_0045202_synapse = gene_id %in% c(
    "SEMA3F","ICA1","ALS2","SARM1","PLXND1","ARHGAP33","PRKAR2B","TSPOAP1","ITGA3","ABCC8","CACNG3","AP2B1","TAC1","CACNA1G","RALA","USH1C","ARHGAP44","PAFAH1B1","BAIAP3","SYN1","CDKL5","ADAM22","SYPL1","MAPK8IP2","RPS20","SLC6A13","FYN","SLC6A7","APBA3","SYT7","GABRA3","CNTN1","ATP1A2","NCDN","SAMD4A","NRXN3","GABRA1","RNF10","ANK1","MUSK","ATP6V0A1","APBA2","RNF19A","CUL3","CYP46A1","SLC18A1","CASR","TRIO","VCAN","CDH1","RTN4R","CDH10","INPP4A","RAB27B","AP2S1","USH2A","ADRB1","ATP6V1H","DTNBP1","ANO2","GOPC","VAMP3","UTS2","SLC4A8","MAPK9","HOMER3","RAD51","PLEKHA5","PTPRN","KIF1B","CBLN4","CYFIP2","GDI2","CAMK2B","ATP2B4","PSD","YBX3","GNA15","RIMBP2","LZTS1",
    "GUCY1B1","SEZ6","GPC1","NGFR","LPAR2","SNCAIP","AP3D1","SNAP91","YBX1","ELAVL1","ATXN3","FGFR2","MYO9A","CACNB1","RHOA","PRKCZ","SYT1","ROGDI","ATP2B3","ROCK1","HDAC4","GRIPAP1","RASGRP2","SDK2","VPS35","KCNAB2","DRD4","NEDD4","GNB5","LRP6","PFN2","SPTB","AP3M2","CNGB1","CHAT","CAMK2A","CDC42","ATP2B1","NCK2","RPS6KA2","ATP6AP1","PRKACA","ADGRL1","ACTN1","FCGR2B","NDE1","CLCN4","ADAM11","PICALM","NSF","SNCB","KCNQ2","SEMA3A","ADD2","CACNG5","CACNG4","ACTB","PLD1","DLG1","RAB7A","GPC4","STXBP2","CTTNBP2","PAK3","ACTN2","CAPZB","ITGA8","AMPH","SYNJ2","GNB1","VDAC3","SLC1A3","OPHN1","SCGN","DNM2","RIMS1","STX7","RAB21","CHRNA3","KCNN2","PSEN1","DLGAP4",
    "IGSF9B","MEF2C","CADPS2","PCDHB4","DLG3","KCNK2","OPRK1","GSK3B","ERC1","P2RX5","NUFIP1","RPS5","CHMP2B","SSH1","EIF3I","CAD","SCAMP1","IGSF9","CTTN","EPS15","MTMR2","ACHE","MT3","ADD1","DNM1L","AURKA","GNA11","LZTS3","RPL6","P2RX7","TMEM230","RPLP0","RPH3A","NOS1","BLVRB","STRN4","RAB11FIP3","EFNB1","NRCAM","RAPGEF4","SLC17A6","RGS17","SNAP23","COMT","GABRP","WHRN","TDRD1","JAK2","ABL1","SYDE2","NRP1","PALMD","STX1B","EFNA2","CBARP","HNRNPM","PALM","SNAP29","CRKL","P2RX6","MAPK1","GRK3","EIF3L","PICK1","SEPTIN3","PACSIN2","NEFH","RASD2","SYNGR1","MYH9","EIF3D","PVALB","TRIM9","ATP6V1D","PSMA3","VTI1B","EIF5","SLC8A3","ZMYND8","RIMS4","DNAJC5","HRH3","PSMA7",
    "NTSR1","ARFGAP1","CHRNA4","EEF1A2","SNPH","MYLK2",
    "PDYN","PLCB4","PAK5","KIF3B","SNTA1","SLC32A1","MAP1LC3A","SYNDIG1","CELF4","USP14","PGRMC1","CHRDL1","GLRA2","SYP","KCND1","FMR1","RS1","RAB40AL","EEA1","CDK16","GABRE","PORCN","SRPX2","ZDHHC15","DRP2","HTR2A","ARHGEF7","CORO1A","MAPK3","CBLN1","VAC14","TSC2","EEF2K","DNAJA3","STX4","SLC6A2","TRIP4","AP3B2","RAB11A","HOMER2","DMXL2","SPG11","VPS18","BLOC1S6","ZDHHC2","FZD3","CALB1","PLAT","EIF3E","NEFM","RELB","LIN7B","SLC17A7","OLFM2","SYDE1","SLC1A6","RPS16","HNRNPUL1","RPS19","ICAM5","NAPA","RABAC1","ATP1A3","PTPRS","GRIN2D","SYNGR4","RAB3D","PPP2R1A","CACNG7","KCNN1","RAB3A","GSK3A","GRIK5","PTN","ATP6V0A4","CHN2","STX1A","CRHR2","NPTX2","EIF3B","PTPRZ1",
    "WASL","TFR2","EZH2","RHEB","MYL7","SLC1A1","DNM1","NCS1","MPDZ","APBA1","SH3GL2","ABHD17B","DVL1","DNMBP","EIF3A","MAPK8","CDH23","PPP3CB","CPEB3","NEURL1","RPL28","LGI1","GIT1","RPL19","FBXL20","RGS9",
    "SEPTIN4","P2RX1","RPS6KB1","PFN1","CHRNE","SLC6A4","SYNGR2","ASIC2","CNTNAP1","ABI3","PPP1R9B","MPP2","PRKAR1A","EFNB3","PMP22","UNC119","GABRA4","NMU","RPL34","WFS1","CLCN3","GLRB","RAPGEF2","CRYAB","HSPA8","P2RX3","NRXN2","APOA4","NECTIN1","SLC1A2","CALCA","RPS13","PPFIBP1","CORO1C","ASIC1","SYT10","LIN7A","PPM1H","SLC6A12","KCNA1","STX2","RAB5B","CHD4","GABRR2","OPRM1","TULP1","MAPK14","STK38","PHACTR1","MDGA1","CAP2","WASF1","ALDH5A1","SLC22A2","QKI","SLC29A1","LAMA4","ERBIN","PCDHB2","CDH9","PCDHB3","PCDHB5","PCDHB6","GRM6","GABRG2","DROSHA","CDH6","SLC12A7","PPP2CA","NR3C1","DPYSL3","CPEB4","HRH2","DBN1","KPNA1","UBE3A","WNT5A","FGF12",
    "GNAI2","RPL24","FXR1","SNX4","ATP6V1A","CSPG5","EIF4G1","ADAM23","OTOF","TANC1","SLC30A3","APC2","RPS15","NDUFS7","SPTBN1","RTN4","GLS","MOB4","SLC5A7","STRN","RAB3GAP1","DARS1","SLC1A4","AAK1","GRIN3B","ATP6V1B1","EPHA4","TNR","RPL22","STXBP3","PARK7","OPRD1","KCNC4","RAP1A","DLGAP3","FBXO2","DNAJC6","RIMS3","IGSF21","APH1A","SLC2A1","AKR1A1","PLPPR4","RPS6KA1","STX12","CHRNB4","SLC8A2","KPTN","RPS25","NRP2","FILIP1","CNR1","AKAP7","SLC16A7","GHRH","OLFM3","PCDH17","DENND1A","BCL11A","CRIPT","SLC17A5","KCNIP2","GRIA2","PCDHB10","PCDHB14","TNN","ARR3","PDZD11","IQSEC3","CLU","PTK2B","CHRNA2","PPP3CC","COPS5","CEP89","KCNJ8","PDZRN3","MAPK8IP1","HCRTR1","GHSR","SLITRK3",
    "HPCA","SV2C","PLG","NPY","CLTA","SSPN","GIPC1","NR4A1","DBH","PLP1","KCNJ2","SLC12A4","SLC12A5","PLCG1","ARFGEF2","STAU1","STX16","IQSEC2","VAMP7","ATP8A1","GRM4","PACSIN1",
    "NRN1","RAB17","PHAF1","RAP2A","EFNB2","SEPTIN6","GRIA3","RPL23","NAPB","FLRT3","LAMP5","EIF2S2","LRFN3","NR1D1","PRR12","PRKCG","ELK1","HNRNPH2","ZC4H2","CANX","RAP1B","PIN1","SYNGR3","PNKD","AKAP9","HIP1","LRFN1","DGCR8","YWHAH","ADORA2A","MPST","GALR3","RNF112","ATP6V1F","CDHR3","VGF","MKLN1","LRRC4","GAD1","ARHGAP22","TMOD2","PALLD","KCNC1","PLD2","ATP1B2","FXR2","NGDN","SNX6","CHRNA10","CDH15","ABHD17A","SYT5","APOE","RPL36","NCAN","KIF1A","SNX9","UNC13A","MAP1S","OLFM1","CALY","LAMA5","HIP1R","AKAP12","SYNE1","ARHGEF9","C1QL1","ATP6V1E1","PRR7","PPT1","SLC6A6","KCNC3","LRRC4B","PDLIM4","KIF3A","RPL27","NDFIP1","EXOC4","ANO1","PPFIA1","NINJ1","NPHP4","MAP1B","DNAJB1",
    "SLC6A11","PNISR","DLG4","FLOT2","SNAP25","PTPRA","BCAN","HAPLN2","SYT11","SYT4","WASF3","CHRM3","DCLK1","SPART","POSTN","EPHB2","RTN3","LARGE1","NTS","NUMB","EIF2S1","ARL8B","CNTN6","GNAT2","SYT6","SORT1","NGF","YWHAQ","RPS15A","HRH4","LRP4","MTNR1B",
    "DTNA","DAGLA","ARHGAP32","P2RX4","ADGRB3","HTR1B","SNX14","EPHA7","CAPRIN1","AGAP2","FAIM2","MICAL1","STX11","SEMA4F","DYSF","KCNMB4","SLC9A5","KCNK1","EGLN1","STX6","RGS8","CHRND","DOCK10","HTR2B","SERPINE2","PLXNC1","PCDH8","SCRN1","RAC1","DGKB","DBNL","ABHD17C","CIB2","MYCBPAP","SCN2A","VPS45","BIN1","GAD2","ABI1","TOR1A","STXBP1","ATP6V1G1","GABBR2","ARPC5L","RANBP6","PLAA","HCRTR2","TUBB2B","FLOT1","CPEB2","MYO7A","FCHSD2","SDCBP","FXYD6","UNC13C","ITPKA","PAK6","ADAM10","DTNB","MYOF","ATAD1","ASCC1","METTL5","CHRNA1","ABI2","SLC40A1","USP8","COPS4","HNRNPD","BLTP1","SEPTIN11","PPA2","PPP3CA","CLSTN3","VAMP1","PPFIA2","GIT2","NPFF","DIAPH3","CBLN3","RTN1","SLC38A6",
    "STON2","SLC12A6","HAPLN3","MCTP2","SH3GL3","PARN","CDH11","CDH13","PDPK1","ARRB2","MINK1","ARHGDIA","SLC16A3","RAB40B","EIF4A3","CBLN2","ERBB2","CACNA1A","SLC39A3","TPGS1","SH3GL1","APP","AKT1","SLC6A3","CACNG8","ACP4","RPS11","RPL13A","SYPL2",
    "SNX27","ANXA9","KCNH1","RAB13","SNAPIN","HCN3","SNAP47","ARF1","CNIH3","PSEN2","PPFIA4","SYT2","ARL8A","RHOB","HNRNPLL","RPS27A","VPS54","C1QL2","PKP4","SLC4A10","RAB5A","GRIP2","IQSEC1","RPL32","ARL6IP5","IGSF11","TMEM108","STXBP5L","DGKQ","SNCA","ANK2","RPS3A","GLRA3","HAPLN1","CRHBP","GABRA6","GABRB2","GLRA1","CPLX2","TENM2","KCNMB1","PSD2","LRRTM2","TPBG","GABRR1","SLC18B1","TIAM2","SLC22A3","SDK1","NLGN4X","SH3KBP1","CASK","OGT","HTR2C","ATP6V1B2","CHRNB3","CHRNA6","CRH","RPL7","MAL2","ARHGAP39","C9orf72","SIGMAR1","NTRK2","RPL7A","CACNA1B","HTR7","ADD3","EIF4EBP2","INA","RGS10","LIN7C","LRRC4C","SYT8","CAPN5","RPS3","DRD2","HTR3B","GRIK4","KIRREL3","SCN2B","PLCB3",
    "DOC2A","CNKSR2","CNTN5","ITGB1","PCDH15","CDH8","ADGRL3","LYPD1","ADRA2A","GPM6A","DLG2","ITPR1","GPR158","CACNA1C","KCNA6","ANK3","DLG5","EIF4E","PTPRO","EPS8","VTI1A","DRD3","GABRA2","TMEM163","GRID2","RIT2","SPOCK1","HOMER1","GRIA4","SPARCL1","FARP1","UTRN","GRM1","CNTNAP4","RAB3C","ASAP1","PTPRD","DGKE","GRAP","CEP112","TNIK","WNT3A","ENAH","SORBS2","WNT7A","RAB6B","EPHB1","PTPRN2","ATP6V1C1","PI4K2A","SLC16A1","GRIA1","ADCY8","GRIP1","KIF5A","PSD3","GNAQ","KCNMA1","TIAM1","SORCS3","RPL30","LRFN2","ZDHHC5","NPTN","LHFPL4","SST","GHRL","NMNAT2","ATP2B2","SLC6A1","FCHO2","SYN2","NRG1","LRP8","HTR5A","FAM81A","APPL1","DGKI","BRAF","ACAN",
    "CABP1","AP3S2","TPRG1L","WASF2","CLSTN2","EIF5B","KCNB1","PPP1R9A","HTR6","DMTN","MPZ","SYNJ1","SV2A","C1QC","IGF2BP1","GJD2","ACTC1","CTBP1","ATP6V0D1","RGS12","ABR","NPR2","KALRN","DIP2A","BRSK1","VPS11","CHRNB2","RUSC1","RPL8","AP2M1","GRIN2C","HCRT","ITGA5","SHANK1","TAMALIN","RPL26","TBC1D24","SHANK2","GNG3","RPS6KA4","ELAVL4","SLC1A7","SYNC","CAMK2N1","NTNG1","KCNJ9","NCSTN","TDRD5","MRPL55","DISC1","LRRTM1","KCNJ3","PDLIM5","FZD5","GABRG1","GABRB1","ATP1A1","ARPC2","ADORA1","PRKCI","CADPS","SYNPR","ABHD6","APBB2","UCN","GRIK3","CAMK2N2","SLC9B2","BSN","CAMKV","GRM2","PRSS12","GUCY1A1","NPY5R","ITGA2","HTR4","EGFLAM","CARTPT","SEPTIN8","GRIK2","STXBP5","IL31RA",
    "RPS14","HCN1","SLC29A4","FABP5","ADCY1","GPER1","CDK5","TMUB1","YWHAZ","BAALC","ALDH1A1","PTCHD1","ATP7A","NLGN4Y","VCP","SLITRK5","LRFN5","ARF6","DACT1","SLC18A2","TSC1","BORCS5","ZMYND19","NSMF","RAPSN","CPSF2","SLC6A5","C1QL3","CACNB2","KCNC2","GPR176","SVOP","RAB8B","GABRB3","STXBP4","APBB1","NETO1","PPFIBP2","RIC3","RPL27A","PRKCB","CPLX4","NDEL1","CHRFAM7A","HTR3A","AP1G1","CACNG2","STX3","TSC22D4","MAP1A","SAMD14","CRK","PRRT2","RAB8A","RPL13","CACNB3","TUBA1A","TTYH1","ATCAY","EEF2","RAB26","SLC3A2","CTNNB1","PNOC","RAB4A","KCNJ4","KIF5C","FAM107A","ARF4","SEPTIN2","LGI3","CHRM1","SEMA4C","SNTB2","NSG1","HTR1E","SLC35G2","MFF","GRM5","CPLX1","MECP2","ACTBL2","ROR2",
    "CPT1C","GSG1L","RAB3B","RGS14","IL1RAPL1","PTK2","KCNK9","PCBP1","DRD5","CHRNA5","RAC3","NLGN1","HNRNPF","P2RY1","CTNND2","SYAP1","FRMPD4","NLGN2","MYRIP","YWHAG","CHRNB1","USP50","GABARAP","SLC30A1","ZNF804A","ADORA2B","CDH2",
    "EMB","DLGAP1","SYT9","GPR37","RPS9","TANC2","ATP6V0E2","PRKCE","GRIK1","NETO2","NPTX1","KCNK3","CLCN5","TPPP","KCND3","PDE7B","LPAR3","BCL2L1","FGB","DSCAM","CLSTN1","GPHN","RPS21","PTEN","RPS7","PRNP","TLN2","SYNPO","RAB33B","GAP43","LAMB2","CALB2","SNTB1","RAB40A","PPP1CA","PURG","RPL38","NBEA","GRK2","RELA","ABLIM3","SYT12","SNCG","C1QB","C1QA","DAG1","LRFN4","EIF1AX","CNP","HAP1","SPTBN2","FBXO45","CHRNA9","ATP2A2","VWC2L","CNTNAP2","NPAS4","HRAS","FZD4","CNIH2","SLC22A1","CTBP2","DES","CADM2","CHRNA7","EIF3F","CLTB","PCSK1","MCTP1","CABP4","FOSL1","PRIMA1","BAIAP2","DOK7","SPHK1","BNIP3","LRRTM4","RIMS2","KCMF1","BDNF","CDK5R1","GRIN1","LY6H","SMCR8","DSCAML1",
    "RPS6KA3","KCNA3","KCNA2","PPFIA3","NECTIN3","HNRNPA0","KCNJ10","PCDHB9","AP3S1","RPS27","HTR3C","SLITRK1","HTR1A","UCN3","ERBB4","KCTD12","STX19","GAK","TUFM","HTR1F","EGR3","SLC17A8","SHARPIN","SLITRK4","HTR1D","AKAP5","NRXN1","SEPTIN1","EXOC3","LYNX1","TH","ZNRF2","PAK2","CHRM4","SCRIB","NQO1","RPH3AL","CHRM2","F2R","PENK","DDN","ACTL9","ADGRB1","SLC2A4","UNC5C","EXT1","ATP6AP2","KCNA4","GABRG3","CLN8","CACNB4","CAV3","PLCB1","SLC35D3","GRID1","RPS17","PLCXD3","RPL35A","GJC1","CADM1","AP2A2","SLC8A1","BEGAIN","GPC6","GABRR3","PMCH","GRIN2A","NOG","KCTD16","KCTD8","BTBD9","ACTG1","NIPSNAP1","CNTN2","KCNQ3","EFNA5","KCND2","PDE4B","SEPTIN5",
    "TMED9","DRD1","CHRM5","SORCS2","SEMA4B","ANKS1B","FLRT2","PURA","RAB11B","SCN10A","PRKN","STAC3","ROR1","SV2B","SHC4","NTF3","SYN3","NRG3","KCNQ5","ATP6V0C","LAMP1","SLITRK2","HTR3E","HTR3D","PIP5K1C","ZNRF1","GABRA5","PPP1CC","BACE1","RPS23","PCLO","RGS7BP","PDE2A","BCR","MAPT","ZDHHC17","P2RY4","ZACN","LYPD6","DCC","PCDHB13","MAGI2","KCNJ11","HAPLN4","ERC2","SLC18A3","GABRD","P2RX2","SHISA7","KLHL17","COL4A5","AGRN","PRKAR1B","POTEE","SRSF10","CLN3","VWC2","FZD9","SHISA6","RPL14","LRRK2","INSYN2A","IL1RAPL2","BTBD8","SH2D5","IL1RAP","PTPRT","GRM7","NLGN3","NTNG2","SLC35F1","SLC6A9","MME","LAMA2","POTEF","HRH1","TRAPPC4","TRPV1","NF1","DAPK1","CHRNG","POTEI","SCN8A",
    "FLNA","AP2A1","SLC6A17","PCBP2","SRC","ADGRA1","SYNGAP1","ADARB1","HTT","COL13A1","PCDHB11","SIPA1L1","RAB40C","RPS26","RPL37A","FCHSD1","RPL12","DNM3","MBP","DLGAP2","RPS4X","CD2AP","LPAR1","CACNA1E","RPL23A","ITSN2","CTNND1","ARC","SLC9A6","UNC13B","LRRTM3","GRIN3A","SCAMP5","GRM3","ARHGEF15","CD3E","NOS1AP","MAGEE1","DMD","PJA2","GDI1","MACO1","BMPR2","AGER","PRRT1","GABBR1","INSYN2B","TMEM240","PSENEN","INSYN1","ITSN1","COLQ","STK38L","GPC2","CPLX3","VDAC1","ATP6V1G2","NPS","CPEB1","HNRNPUL2","PHB2","GRID2IP","VAMP2","PPP3R1","PLXNA4","NPTXR","PPP2R2A","POTEJ","C4B","NTF4","ELFN1","RPS18","DYTN","RPS28","SHISA8","SHISA9","CORT","EIF6","STON1","LY6G6D","C4A","FER1L5",
    "SHANK3","LYN","OMP","ITGB3","FRRS1L","DYNLL2","SRGAP2","S1PR2","GABRQ","DOC2B","PCDHB16","GRIN2B","CYFIP1","SRCIN1","GPR179","NEFL","ADORA3","SLURP2","NPBWR1"
  )) |> 
  dplyr::mutate(GO_0062023_collagen_containing_ECM = gene_id %in% c(
    "DCN","SEMA3B","MARCO","SERPINB1","MYOC","TIMP2","VCAN","TNC","USH2A","COL9A2","LTBP1","ELN","LAMC3","COL23A1","LAMA3","ENTPD2","CBLN4","ITIH4","ITIH1","F7","LAMC2","COL11A1","WNT8A","GPC1","CDON","NTN1","COL17A1","FGFR2","PKM","FGF10","SPP2","ADAM11","NTN4","LMAN1","DLG1","GPC4","FBLN1","COL5A3","HSP90AA1","COL4A4","IMPG2","COL19A1","COL16A1","FCN1","ACHE","ADAMTS2","MMP2","NID2","LTBP4","ICAM1","P3H2","LAMB4","LAMB1","APOH","ANGPT2","CCDC80","CMA1","COL9A3","TGFB2","AMELY","HNRNPM","LGALS1","TIMP3","PDGFB","CHADL","CTSG","COCH","MMP9","BMP7","CTSZ","COL20A1","ANGPT4","LAMA1","MXRA5","SRPX","F9","TIMP1","SRPX2","CBLN1","ZP2","CTSH","SFRP1","IL7","FGL1","CLC","APLP1",
    "TGFB1","COMP","WNT2","PTPRZ1","PCOLCE","SERPINE1","PLOD3","RARRES2","AEBP1","MEGF9","OGN","ASPN","ECM2",
    "AMBP","CXCL12","ACTA2","KAZALD1","LGALS3BP","COL1A1","VTN","SOD3","CTSC","HPX","APOA4","APOC3","MDK","VWF","WNT5B","MGP","COL12A1","COL9A1","SMOC2","IMPG1","LAMA4","ERBIN","LOX","SPARC","THBS4","KNG1","HRG","WNT5A","COL7A1","LOXL3","EFEMP1","FN1","TNR","ANGPTL1","PRG4","NID1","ZP4","MFAP2","P3H1","F3","SERPINC1","CTSD","MMP8","APOA1","CCN2","LTBP2","TGFB3","TNN","TGFBI","CLU","A1BG","FMOD","PLG","ANXA11","COL10A1","SLPI","MATN4","NDP","F13A1","SERPINB6","COL21A1","AMELX","CFP","PZP","FGL2","LOXL1","APOE","NCAN","PXDN","GDF15","COL5A1","LAMA5","FIBCD1","ANGPTL6","F12","LGALS3","PODNL1","RTBDN","MATN3","EMILIN2","SERPINF1","ITGB4","MATN2","BCAN","HAPLN2","APCS","ANGPTL3",
    "POSTN","LOXL2","ADAMDEC1","WNT2B","COL4A2","ADAMTS8","ANXA1","CTSL","ADAM19","AGT","LAMC1",
    "SERPINE2","ANGPTL2","CCN3","IGFBPL1","TINAG","SULF1","THBS1","EMILIN1","LOXL4","ANXA7","CILP","SEMA7A","MMRN1","FRAS1","FBN2","COL2A1","INHBE","LUM","FBLN5","PCSK6","LMAN1L","HAPLN3","MFGE8","TGFB1I1","CDH13","COL6A1","COL6A2","ADAMTS10","FCN3","HSPG2","CCN1","TINAGL1","DPT","HDGF","HMCN1","ECM1","ADAMTSL4","ANXA9","S100A8","S100A7","FLG","LEFTY2","COL8A1","AHSG","SFRP2","HAPLN1","CASK","GPC3","HMCN2","SERPING1","SERPINH1","NCAM1","ZP1","FREM2","ITIH2","DST","SPARCL1","ABI3BP","ANGPT1","ADAMTS1","ADAMTS5","ADAMTS3","PRG3","ACAN","COLEC12","ADAMTS4","LAD1","C1QC","CSTB","FCN2","AZGP1","COL26A1","MATN1","SDC3","CTSS","S100A9","COL6A3","IGFBP7","FBLN2","ADAMTS9","CRELD1","PF4","CPA3","TGM4","CLEC3B","ANXA5","EDIL3","EGFLAM","RELL2","SHH","COL1A2","CTSB","SBSPON","CTHRC1","FREM1","VWA2","HTRA1","ADAMTS15","FBN1","SERPINB8","MFAP4","HSP90B1",
    "SERPINB12","NAV2","GREM1","SERPINF2","KRT1","ANGPTL4","SOST","LTBP3","TNXB","COL3A1","NPNT","LGALS9","COL4A3","THBS3","COL22A1","SDC2","MUC17","SERPINB9","CDH2","COL24A1","FGG","FGA","FGB","LGALS4","COL8A2","ANGPTL7","LAMB2","LRRC15","TPSAB1","EFEMP2","COL6A5","EGFL7","ADAMTS20","MMRN2","C1QA","DAG1","CSPG4","CTSF","PODN","FGFBP3","ZG16","NPPA","A2M","CLEC14A","CD151","CALR","GPC5","VWA1","SSC5D","F2","ADIPOQ","OTOL1","BGN","MXRA7","SNORC","ANXA2","COL18A1","FREM3","GPC6","EMILIN3","EFNA5","THBS2","PRG2","C17orf58","ANGPTL5","BCAM","COL4A1","HAPLN4","AMTN","THSD4","COL14A1","EYS","COL4A5","AGRN","PLSCR1","ZP3","IFNA2","SERPINA5","COL25A1","VWC2","PRELP","MMP23B","SERPINA3","S100A4","PRTN3","LAMA2","COL27A1","LAMB3","ANXA4","ANXA6","SERPINA1","TPSB2","COL13A1","SPN","ELANE","COL4A6","MFAP5","PSAP","S100A10","HRNR","S100A6","SMOC1","EGFL6","L1CAM","TGM2","COL11A2","COL5A2","COL15A1","PRSS1","VIT","DEFA1","COL6A6","COLQ","GPC2","ANG","COL28A1","ORM2","ORM1","DEFA1B","MARCOL","GH1","SPON1","ANXA8","RBP3","GDF10","MMP28"
  )) |> 
  #dplyr::mutate(GO_0006955_immune_reponse = gene_id %in% c(
  #  "TSPAN6","FGR","CFH","SARM1","CD38","DHX33","SPPL2B","CREBBP","ITGAL","LAMP2","TMEM98","MAP3K14","CX3CL1","CCL26","MARK4","NOS2","CD79B","PGLYRP1","ST3GAL1","IL32","TENM1","MASP2","PRSS3","CD4","BTK","HFE","FYN","TYROBP","ERCC1","CD22","LTF","NR1H4","ALOX5","KDM5D","NUB1","POLR3B","CD6","WAS","MATR3","SLC11A1","MARCO","CD74","C8B","NLRP2","BIRC3","NR1H3","VIM","FAS","RNASET2","SLAMF7","BTN3A1","PRKCH","IFNGR1","TNFRSF1B","POU2F2","HMGB3","GRN","GAB2","STAP1","C6","PARP3","LCP2","C12orf4","TNFRSF17","VAMP3","FOXP3","TNIP3","MAPK9","BCAR1","CYBA","PHPT1","LY75","EIF2AK2","PUM2","RC3H2","TRAF3IP2","PRDM1","ATG5","POLR3E","WNK1","SPHK2","CASP8","SBNO2","IPO5","KARS1","MAP2K4","PRKCQ","CD84","SPI1","SP100","TP53BP1","PRKCZ","HDAC4","IFI35","TFE3","SIRT2","BCL3","RORA","TGFBR3","NEDD4","RAB27A","UFD1","EIF2B3","CAMK2A","CDC42","FAM3A","TFRC","FCGR2B","PVR","ALPK1","CD5L","ST6GAL1","TBX21","DPP8","TXK","SEMA3C","GRAMD4","FOSL2","MKRN2","MLH1","UNG","TRAF4","PAG1","RAP1GAP","STXBP2","MAP2K7","NFKB2","IL4R","PAK3","APBB1IP","LRCH4","CST7","LAMP3","UBE2K","EDN1","ITM2A","ITCH","TOLLIP","CDH17","THOC1","XRCC5","RAPGEF3","CEACAM1","DDX1","STX7","RIF1","MID2","PSEN1","HSP90AA1","CFHR2","CXCL2","TCF7","ARG2","MEF2C","PTPRC","FYB1","OPRK1","CYLD","CD59","FCN1","WDFY1","LAT2","BAX","DEFB127","FKBP1A","MAVS","P2RX7","OAS1","LAG3","RGS1","ICAM1","IRAK3","LYZ","MUL1","SPG21","CD209","MCOLN1","NLRC4","IL5RA","RAPGEF4","TF","NLRP1","ESR1","CD200","PPP2R3C","RNF31","TYRO3","SNAP23","EZR","UNC13D","TGFB2","MSH2","MAP3K1","BLNK","TREM2","CRISP3","SRPK1","PGC","SIRT1","JAK2","IL12RB1","ABL1","TMED1","MARCHF2","MADCAM1","CRKL","OSM","MAPK1","MFNG","PLA2G3","LGALS1","XBP1","GTPBP1",
  #  "TTLL1","APOBEC3H","TTLL12","APOL1","CSF2RB","IL2RB","EP300","POLR3H","CTSG","GZMB","CNIH1","CHGA","NFKBIA","CD40","SLA2","NFATC2","PTK6","HCK","SAMHD1","JAG1","BPI","CST9L","WFDC2","USP14","SMAD7","RNF125","RIOK3","TLR8","XIAP","BMX","PQBP1","GATA1","CD40LG","TNFSF13B","ACP5","RGCC","ACOD1","CORO1A","N4BP1","CCL22","CCL17","FOXF1","SLC7A5","MEFV","DNAJA3","PYCARD","STX4","IL21R","CSK","CTSH","CD276","TRIM35","RIPK2","NBN","UBE2W","IKBKB","IL7","GSDMD","FGL1","KCNN4","RELB","FCGRT","FCER2","CLEC4M","IL4I1","LILRA1","IL27RA","AKAP8","CLC","GPI","EBI3","PRKD2","CACTIN","TGFB1","CD79A","RPS19","NOP53","NKG7","CD33","TYK2","CDC37","PTPRS","ZNRF4","CARD8","LILRB5","DNASE2","JAK3","PIK3R2","HAMP","PIK3CG","ITGB8","GSDME","ZC3HAV1","CAV1","LFNG","TAX1BP1","NOD1","EPHB6","CCL24","PHF14","RARRES2","AHR","LIMK1","TRIM14","C5","AMBP","TNFSF8","RIGI","RAPGEF1","PTGDS","SHB","EXOSC3","GATA3","RAB11FIP2","CXCL12","MAPK8","PLEKHA1","VSIR","PPP3CB","BMPR1A","LIPA","GBF1","CUEDC2","CSF3","RPS6KB1","C1QBP","CCL7","CCL2","CCL8","CCL1","DHX58","NAGLU","DUSP3","YWHAE","VTN","TNFAIP1","SUPT6H","TMEM33","PF4V1","NFKB1","MAPK10","IL2","DHX15","CLNK","NSD2","BST1","CTSC","CRTAM","HSPA8","LPXN","DTX4","UNC93B1","KMT5B","HPX","TRIM3","PANX1","APOA4","BIRC2","CBL","NECTIN1","SLC15A3","MDK","CD81","TCIRG1","POU2AF1","IL23A","CYP27B1","ITFG2","LTBR","OAS3","OAS2","EIF2B1","ENDOU","IFNG","GAPDH","PTPN6","HCFC2","CLEC4A","AICDA","KLRB1","BTN3A3","FRK","ULBP1","MAPK14","IL17A","IL17F","RNF8","FBXO9","CD83","BACH2","GPLD1","VNN1","TRIM38","CCR6","TFEB","DUSP22","BTN2A1","LY86","ERBIN","C7","CERT1","ITK","IL12B","BTNL8","POLR3G","PDE4D","IL4","IL5","PPP2CA","TRIM23","C9","CSNK1A1","HRH2","HRG","BCL6","CD86","WNT5A","CBLB","HHLA2","SNX4","SCAP","MAPKAPK3","FOXP1","IL1A","CCL20","ZAP70","ACTR3","EIF2B4","ITGB6","ITGA4","IFIH1","RTN4","LOXL3","REG1A","STAT1","GNLY","IL1R2","IL1R1","IL1RL2","IL1RL1","IL18R1","IL18RAP","SOS1","KYNU","NFE2L2","MSH6","STXBP3","PARK7","OPRD1","RNF19B","TRIM62","SFPQ","ARHGEF2","MAD2L2","PRG4","NCF2","PLA2G4A","CFHR3","PLEKHM2","CD58","CD2","ZP4","SLAMF1","UAP1","PLA2G2D","GBP3","GBP1","RAB29","CD160","CR2","CD46","MPL","PRDX1","PIK3R3","FASLG","TNFSF4","COLEC11","STK11","KMT2A","APOA1","CNR1","TNFAIP3","MYB","ARG1","VAMP8","CCDC92","TRIM32","PPP6C","SLC46A2","NR4A3","TGFB3","EIF2B2","IFIT3","NKX2-3","IFIT2","CRHR1","CD274","IFNA6","IFNA8","TASL","TNFSF18","GPR31","TNFSF11","ELF1","CLU","PTK2B","EPX","AKAP1","TRIM25","TRIM6","ADCY7","KCNJ8","PSPC1","CD80","KHDRBS1","CCRL2","CCR2","
  #  TNFSF10","PIK3CA","CXCR4","SASH3","LAX1","CD244","LY9","SHLD2","PMS2","TWIST1","PKN1","ADGRE5","OPTN","NCKAP1L","LRP1","IL13RA2","NMI","EXOSC9","C4BPA","C4BPB","PI3","SLPI","WFDC3","SEMG2","PLCG1","ZNFX1","SEMG1","PCK1","ZBP1","VAMP7","IL9R","NAGK","CEACAM8","BTN2A2","WRNIP1","BTN1A1","H2BC11","TREM1","RAB17","CXCL6","EREG","TRIM51","IRF1","KIR2DL1","SRMS","IL1B","IL37","TNFSF9","MED1","CD70","C3","GPR108","TNFSF14","DEFB126","RBCK1","S1PR4","CITED1","ROMO1","BPIFB1","FFAR2","HCST","CCR7","NR1D1","IRF3","IFI6","CFP","ZBTB1","TRAF2","MASP1","PLA2G5","ADGRE2","GFER","TICAM1","IL17B","FGL2","VPREB3","RFPL1","RFPL2","RFPL3","IGLL1","RAC2","LIF","APOBEC3A","APOBEC3F","IRF5","EIF2AK4","ACKR4","PSMA1","PLD2","CD68","RIPK3","RAB2B","FOXJ1","SEC14L1","RHBDF2","LBP","CRACR2A","ECSIT","NECTIN2","APOE","BST2","FCHO1","PXDN","TRPM4","ASS1","TRIM28","SMPDL3B","THEMIS2","SHFL","ULBP2","ULBP3","
  #  LILRB2","DEFB118","CCL25","F12","PRR7","IDO1","TRAF3","RFTN1","NR1H2","ANKHD1","NDFIP1","NINJ1","RARA","USP29","SELENOS","GCH1","LGALS3","RFX1","TRIM21","RAF1","PPARG","TRIM5","TRIM22","RAP1GAP2","JCHAIN","ANKRD17","KDM6B","CLEC10A","GPS2","POLR3F","CRP","APCS","CTNNBL1","SERINC3","AP3B1","CHIT1","EPHB2","PRAM1","KMT5C","SFTPD","SPINK5","SWAP70","MEN1","ADAMDEC1","CD180","IRAK2","ARL8B","VAV3","PTPN22","NOTCH2","VTCN1","SPIRE1","RSAD2","IL6ST","CFHR4","CFHR5","GRP","IL2RA","DOCK2","KLRD1","KLRC1","PUM1","C5AR2","ETS1","UBQLN1","ANXA1","CTSL","HAVCR2","OASL","TRAFD1","CD36","SRPK2","SYNCRIP","AKIRIN2","MAP3K7","TESPA1","CD164","NMBR","TEC","DYSF","USP15","DHX9","RC3H1","DOCK10","EIF4E2","USP44","APPL2","FLNB","LCP1","IL6","DBNL","MYO1G","CALCOCO2","TANK","BLK","EPRS1","HLX","IL10","IL36G","IL1RN","IL36A","IL36RN","IL36B","IL1F10","DAB2IP","STXBP1","TLR4","USP20","CTSV","ENPP2","IL33","CCL21","SIT1","IFNA21","CNPY3","IRF4","FLOT1","ATAT1","TLR2","IL18BP","DDX60","TRIM29","CASP1","THBS1","SMAD6","IFI44L","IFI44","ACTR2","CH25H","STAT4","SENP7","PARP9","SPPL2A","SEMA7A","UBL7","HERC5","IL21","CXCL9","G3BP2","ANXA3","CASP6","LEF1","RNF185","C1RL","KLRG1","CD27","PIANP","SLC15A4","TARBP2","ITGB7","RAB20","GPR65","TCF12","PSTPIP1","IGF1R","PML","FURIN","IGSF6","CHST4","NLRC5","CMTM3","IRF8","PDPK1","SKAP1","SLC39A6","GATA6","ARRB2","PIK3R5","TP53","TRIM65","SECTM1","TNFRSF11A","VAV1","IFITM3","IFNAR1","APP","AKT1","IL19","GPR32","SIGLEC10","KLK3","GARIN5A","RPL13A","FCN3","BCL10","ATP1B1","XCL1","XCL2","FCGR2A","RORC","ECM1","CTSK","SUSD4","DUSP10","ADAM15","S100A8","S100A7","LYST","PARP1","REG3G","TRIM43B","TRIM43","GPR17","HSPD1","ACKR3","ACKR2","NFKBIZ","EIF2B5","LYAR","SNCA","TIFA","MARCHF1","SKP2","PLK2","GZMA","PIK3R1","TSLP","ATG12","CXCL14","IL9","FBXO38","TNIP1","G3BP1","N4BP3","CPLX2","TRIM7","TRIM41","TNFRSF21","IRAK1BP1","TRIM4","NONO","IL2RG","DOCK11","MFHAS1","EBAG9","IFNA5","IFNA16","C9orf72","IFNK","GSN","LCN2","NOTCH1","POLR3A","DGKZ","SERPING1","PTPRJ","ENDOD1","PAK1","RPS3","DRD2","TKFC","MS4A2","CDC42EP2","FAU","MMP3","TRIM48","FCGR1A","TIRAP","CD226","IL18","NEK7","FER","VPS26B","ADAM8","ADAM17","TMEM45B","CCL28","CYSLTR2","EPG5","DCLRE1C","CAMK4","IFIT5","MR1","BANK1","BMP6","CD96","CYRIB","FBXL2","CD8A","PTPRD","CEBPG","MCOLN2","SDHAF4","THY1","OTULIN","TOMM70","ANGPT1","ENPP3","TRIM11","GBP5","LY96","RABGEF1","PLCL2","SAMSN1","PIK3AP1","VSIG4","RAET1L","BATF","CXCL13","RPL30","HK1","PRG3","UBE2L6","ZDHHC5","MS4A1","ITGAD","C8A","LRP8","PAXIP1","IL34","KIT","APPL1","MX1","BRAF","SPPL3","TNFRSF14","NCK1","COLEC12","CD1D","CD1A","CD1C","CD1B","CD1E","HMHB1","NCF1","SLAMF8","RNF166","PLA2G2F","EDA","FCER1G","APOA2","IFNAR2","IFNGR2","C1QC","MORC3","PADI4","PSMB4","C1R","PGLYRP3","SPON2","ZDHHC1","ZYX","TNFRSF13C","UBASH3A","ICOSLG","AIRE","ITGB2","VAV2","ADAMTS13","FCN2","ZDHHC12","NLRP4","CD3G","CXCR5","ZBTB7B","NLRX1","ADAR","IL6R","CCR5","FCRL3","AZGP1","CYP11B1","LRRC14","SH2B2","PGLYRP2","ALOX15","TREML1","CXCL16","SCIMP","TNFSF13","POLR3K","CLPB","JAK1","LAPTM5","IL23R","GBP2","GBP4","GFI1","NLRP3","TRIM58","SLAMF6","FCRLB","FCGR3B","MAPKAPK2","FCMR","PIGR","FCAMR","REL","SANBR","NEURL3","CTSS","TNFAIP8L2","PGLYRP4","S100A9","S100A12","FZD5","CLDN1","DCST1","INAVA","SLC15A2","CXCR1","EOMES","AZI2","FCRL4","TRAT1","MNDA","PYHIN1","IFI16","AIM2","CTLA4","ICOS","CD200R1","PTX3","DNASE1L3","RBM47","CXCL3","CXCL5","PPBP","PF4","CXCL1","CCR1","DTX3L","ZC3H12A","PRKCD","RNF168","OTOP1","ERMAP","CAMP","MST1R","HMGB2","NPY5R","IL15","OTUD4","F2RL1","WDR41","SERINC5","ERAP1","ERAP2","TLR3","IL3","CSF2","LEAP2","CGAS","IL31RA","RAET1E","SLC30A8","DEFA5","DEFA4","DEFA6","DEFB1","GPER1","SPAG11B","YWHAZ","GEM","SYK","NFIL3","SVEP1","CYBB","ATP7A","STOML2","MARCHF8","INPPL1","MBL2","SLC18A2","TSC1","HPRT1","DDX21","RNASE7","BTNL9","METTL3","IFI27","JAM3","IL25","STXBP4","C2","SMPD1","RAG1","PLD4","PRKCB","CLEC4E","CLEC4D","B2M","AP1G1","STAT6","SMAD3","PDIA3",
  #  "PHB1","SP2","GPRC5B","CRK","NOD2","DUS2","CD3D","TRIM68","CACNB3","AXL","NFKBID","LAIR1","KIR3DL1","DAPK3","TMIGD2","KLK5","KRT1","OTUB1","CD300A","NLRC3","FTH1","FADD","MAP4K2","RBPJ","PTGDR","BMI1","CX3CR1","SCN11A","TAP1","SCNN1B","POLR3D","COL3A1","STAT3","WFDC13","IL7R","WFDC12","IL12A","TNIP2","BTNL3","INPP5D","TRIM49","LGALS9","KLK7","IL13","GCSAML","CXCL10","CXCL11","UMOD","GP2","SIN3A","RNASE2","RNASE3","PTK2","PTAFR","RNASE6","CXCL8","GPR183","VPREB1","APLF","TRIM56","ITGAM","PYDC1","ALCAM","USP38","USP50","STX8","ADORA2B","CD14","SERPINB9","STAT2","RNF34","KIF5B","TMEM43","PLA2G1B","OSCAR","FPR2","FPR1","PRKCE","SOCS5","MUC7","TMEM126A","TRIM8","JUNB","APLN","PTGER4","FGA","FGB","PIK3CD","DEFB4A","LGALS4","BCL2","IFNB1","C3AR1","PRNP","AQP4","SHLD1","REG3A","REG1B","CD8B","CCL11","MALT1","ISG20","CXCR6","CEBPB","AZU1","CLEC7A","IL16","CTSW","RASGRP1","KLHL6","THEMIS","CCL19","RAB43","MYD88","SLC22A13","KAT5","RELA",
  #  "PARP14","CYSLTR1","GPR151","CST9","C1QA","RNASE8","PPP1R14B","MST1","XCR1","CCR9","STAT5B","CD7","NRROS","TLR10","TLR1","TLR6","EXO1","LIG4","GCSAM","PELI3","IL20RB","CMKLR1","LEP","HRAS","ADGRE1","PDE12","NLRP6","DHX36","RAG2","TRAF6","WFDC5","BANF1","PTPN2","CLCF1","FOSL1","GAPT","A2M","ZNF683","YES1","ATAD5","DEFB104A","PNMA1","TCIM","C8G","SMCR8","DEFB104B","IFNW1","RPS6KA3","DEFB4B","PAWR","CD19","IL17RA","DEFB114","CHID1","GRB2","UBE2N","PLEC","SPAG11A","CD28","DEFB125","CYBC1","AURKB","TRIML2","PTPN11","HLA-DQB1","GATA2","CIITA","ZFPM1","CDC42EP4","LACC1","FCER1A","APOBEC3B","CCR8","SSC5D","WFDC11","WFDC9","MYL11","F2","PRKRA","WFDC10A","PAK2","DEFB124","DEFB123","DEFB119","FUT7","H2BC4","PRF1","CXCR2","DEFB112","NQO1","CCL13","RNF135","SETD2","TNFSF15","ADGRB1","TRIM49B","UBA7","EXT1","SHMT2","NLRP10","C1S","CACNB4","IFNL1","NPLOC4","FES","CLEC4G","CSF1R","LCK","CEP63","WFDC10B","CADM1","SPNS2","PTGDR2","GBP6","MX2","KLRC4","HMCES","CCR3","PSG9","MRGPRX2","IFNL2","ASCL2","TBK1","CCR4","SH2D1A","ACTG1","TRIML1","IRAK1","DEFB108B","PRKD1","CSF1","COLEC10","PLA2G6","LRRC19","CCR10","BPIFC","SOCS3","STING1","PDE4B","H2BC21","USP18","TMEM106A","IFNE","IFITM2","ATP6V0A2","IFNLR1","IRF7","BRCC3","YTHDF3","IFIT1","NLRP9","IFITM1","LAMP1","FFAR3","DEFB128","CD300LF","POLR3C","DEFB131A","ZNRF1","BPIFB3","BTLA","FCAR","DEFB132","BTN3A2","DEFB105A","DEFB107A","DEFB106A","DEFB105B","PRG2","BCR","IFNA10","CXCR3","LILRB4","TNFRSF4","KRT16","HEXIM1","DEFB106B","LILRA5","FPR3","TLR5","ISG15","SEMA4D","CARD9","FYB2","DMBT1","C17orf99","TREML4","NCR3LG1","TUBB4B","PLSCR1","ZP3","GPR15LG","IFNA2","PDCD1","CLEC2A","IGHV1OR15-9","ZDHHC11","CNR2","KIR2DL4","LITAF","S100A13","TRIM64B","S100A14","HMGB1","NCR1","IL1RAP","HLA-DRB1","SEMA4A","TUBB","CD55","PTPN1","PRTN3","XRCC6","TLR7","DAPK1","HLA-DQA1","CD47","ILRUN","ADA","ARID5A","SLC39A10","CASP4","IFNL3","PCBP2","SRC","IL27","PELI1","SLC22A5","ADARB1","C5AR1","MAP3K5","SPN","GALP","GZMM","ELANE","ENPP1","MPEG1","PDCD1LG2","NMB","CR1L","CFD","ATAD3A","MYO1C","H2BC12","IFNA1","PLCG2","MBP","IRAK4","ENTPD7","SIRPA","DEFB107B","CLEC4C","BPIFA1","CARD11","YTHDF2","HLA-DRB5","SH2D1B","CALM1","DLL1","MSRB1","GPATCH3","CNOT7","MTOR","PNP","CD247","SUCNR1","HMGN2","UBE2J1","CD3E","RPL39","TBKBP1","PJA2","CR1","RAET1G","FCGR3A","PLPP4","LIME1","DEFB110","GIGYF2","ZDHHC18","HLA-DOA","HLA-DMA","TAP2","HLA-DRA","BTNL2","AGER","HSPA1B","HSPA1A","TRIM49C","TRIM64","BAG6","AIF1","NCR3","LST1","NFKBIL1","MICB","HLA-C","DEFB121","DHX16","LILRB3","GNL1","HLA-E","TRIM15","TRIM10","TRIM40","TRIM31","HLA-G","HLA-F","MOG","TRIM27","CD177","TRIM13","PKHD1L1","LGR4","PSMB10","IPO7","CFI","SAMD9","KRT6A","PLPP6","KLRC2","KLRC3","CLEC6A","DEFB134","DEFB135","DEFB136","DEFA1","SERPINB4","HLA-A","IGKJ1","IGKV4-1","IGKV5-2","IGKV6-21","IGKV2D-26","IGKV3D-20","IGKV6D-41",
  #  "IGKV3D-11","IGKV1D-42","IGLV4-69","IGLV8-61","IGLV4-60","IGLV6-57","IGLV11-55","IGLV10-54","IGLV5-52","IGLV1-51","IGLV1-50","IGLV5-48","IGLV1-47","IGLV7-46","IGLV5-45","IGLV1-44","IGLV7-43","IGLV1-40","IGLV5-37","IGLV1-36","IGLV2-33","IGLV3-32","IGLV3-27","IGLV3-25","IGLV2-23","IGLV3-22","IGLV3-21","IGLV3-19","IGLV2-18","IGLV3-16","IGLV2-14","IGLV3-12","IGLV2-11","IGLV3-10","IGLV3-9","IGLV4-3","IGLV3-1","IGLJ1","IGLC2","IGLC3","TRGJ1","TRGV11","TRGV10","TRGV9","TRGV8","TRGV5","TRGV4","TRGV3","TRGV1","TRBV6-1","TRBV7-1","TRBV4-1","TRBV6-4","TRBV7-3","TRBV5-3","TRBV9","TRBV10-1","TRBV11-1","TRBV6-5","TRBV6-6","TRBV5-5","TRBV7-6","TRBV5-6","TRBV5-7","TRBV5-1","TRBV4-2","TRBV19","TRBV20-1","TRBV23-1","TRBV24-1","TRBV27","TRBV28","TRBJ2-1","TRBJ2-3","TRBJ2-7","TRAV2","TRAV3","TRAV4","TRAV5","TRAV6","TRAV7","TRAV8-1","TRAV9-1","TRAV10","TRAV12-1","TRAV8-2","TRAV8-3","TRAV13-1","TRAV12-2","TRAV8-4","TRAV13-2","TRAV14DV4","TRAV9-2","TRAV12-3","TRAV8-6","TRAV16","TRAV17","TRAV18","TRAV19","TRAV20","TRAV21","TRAV22","TRAV23DV6","TRDV1","TRAV24","TRAV25","TRAV26-1","TRAV27","TRAV29DV5","TRAV26-2","TRAV34","TRAV36DV7","TRAV38-1","TRAV38-2DV8","TRAV39","TRAV40","TRAV41","TRDV2","TRDJ1","TRAJ42","TRAJ31","TRAJ3","IGHA2","IGHE","IGHG4","IGHG2","IGHA1","IGHG1","IGHM","IGHJ1","IGHV6-1","IGHV1-2","IGHV1-3","IGHV2-5","IGHV3-7","IGHV3-11","IGHV3-13","IGHV3-15","IGHV3-16","IGHV1-18","IGHV3-20","IGHV3-21","IGHV3-23","IGHV1-24","IGHV2-26","IGHV4-28","IGHV3-33","IGHV4-34","IGHV3-35","IGHV3-38","IGHV4-39","IGHV1-45","IGHV1-46","IGHV3-48","IGHV3-49","IGHV5-51","IGHV3-53","IGHV1-58","IGHV4-61","IGHV3-66","IGHV1-69","IGHV2-70D","IGHV3-73","IGHV7-81","DENND1B","CRIP1","TRIM59","CHUK","PVRIG","SIPA1","GBP7","LAT","TREX1","KLRK1","EMP2","UBD","LTB4R","DNASE1","CCL27","IFNA7","ANG","SCART1","TRIM77","DEFB113","TRIM64C","GPR33","PHB2","NLRP2B","DDX3X","DEFB108C","DEFB116","VAMP2","TRIM51G","HMSD","IGLV9-49","TRIM49D1","USP17L2","EXOSC6","RFPL4A","IGHV3-64","HLA-DPB1","TRDD1","IGKV3D-15","CPTP","IGHV4-59","C4B","IGHV3-74","IGKV6D-21","SLC26A6","IGHV3-72","DEFB131B","TRBV2","LTA","LTB","IFNA14","IGKV3D-7","RFPL4AL1","TRBV10-2","TRBV5-4","IGKV1OR2-108","HLA-DPA1","IGHV4-31","IGHV3-43","HLA-DQB2","TNF","TRBV29-1","GPX1","TRGV2","IGHV3OR16-10","TRIM49D2","IFNA13","IGKV3OR2-268","TRIM26","H2BC12L","HLA-B","IFNA17","NFAM1","IGHD1-1","IFNA4","TRBV30","HLA-DQA2","IRGM","TRBV3-1","RBM14","IGKV2D-30","TNFSF12","APOBEC3G","TLR9","IGKV1D-8","DEFA3","IGKV1-6","IGKV1-37","IGKV3-20","LILRA4","IGKV1D-33","LILRA2","DEFA1B","IGKV1-17","TNFRSF13B","IGKV1-8","IGKV1-16","MIF","HLA-DOB","TDGF1","IGKV1D-16","CRCP","IGKV2-24","IGKV3-11","IGKV2D-24","IGKV1-9","SPRR2A","IGKV1-33","IGKV1-39","IGKV2D-28","HLA-DMB","IGKV1D-43","IGKV1D-17","IGKV3-7","IGKV2-30","IGKV2D-29","IGKV1-12","STMP1","TICAM2","IGKV1-5","IL10RB","CFB","TTC4","KIR2DL3","APOBEC3D","IGKV2-28","CFHR1","IGKV3-15","LILRA6","APOBEC3C","IGKV1-27","C4A","ZCCHC3","TARM1","TNFSF12-TNFSF13","NAIP","IGKV1D-37","IGKV2D-40","TMED7-TICAM2","IGKV1D-39","TRBV6-7","SHLD3","TRBV7-7","TRBV7-4","TRBV6-8","PYDC2","PRKDC","LYN","CD8B2","INS","IGLL5","TRAV1-1","RAB44","TRIL","TRAV1-2","TRDV3","CLEC12B","KLRF2","LSM14A","CLEC5A","TRIM34","RNASE4","TRAV30","IGHV4OR15-8","IGHV2OR16-5","IGHV3OR15-7","IGHV3OR16-17","MRC1","MMP12","IKBKE","OTUD7B","IKBKG","IGHV3OR16-12","IGHV1OR15-1","IGHV3-30","IGHV3OR16-8","IGHV3OR16-13","CCL5","MILR1","CD24","H2BC8","USP27X","IGKV2-40","H2BC6","IGHV2-70","CCL23","TRBV12-3","CCL16","TRBV12-5","TRBV16","CCL4","CCL18","CCL15-CCL14","CCL15","TRBV14","TRBV10-3","PRSS2","ORAI1","CCL4L2","CCL3L3","PIK3R6","TRBV13","CCL14","TRBV18","IGKV1D-13","TRBV11-3","RAB7B","IGHV4-4","TRBV12-4","H2BC7","IGHV1OR21-1","CCL3","TRBV17","TRBV7-9","IGLV2-8","H2BC10","IGKV1D-12","IGHV1-69D","TRBJ1-4","IGHV1-69-2","IGHV7-4-1","TRBJ1-3","TRBJ1-5","TRBJ1-1","TRBJ1-2","TRBD1","TRBV25-1","IGHV3-64D","IGHV5-10-1","TRBJ1-6","TRBV7-2","TRBV6-2","PYDC5"
  #)) |> 
  #dplyr::mutate(GO_0002526_acute_inflammatory_response = gene_id %in% c(
  #  "BTK","CD6","PTGER3","ITIH4","CREB3L3","TFRC","FCGR2B","PTGS2","GSTP1","B4GALT1","EIF2AK1","OSM","RHBDD3","PIK3CG","EPHB6","TFR2","LIPA","VNN1","IL4","IL1A","ACVR1","FN1","PARK7","ASH1L","PLA2G2D","F3","TNFSF4","CNR1","TNFSF11","IL1B","C3","FFAR2","CCR7","IL22","APOL2","LBP","EPO","ASS1","F12","SELENOS","PPARG","CRP","APCS","ALOX5AP","SAA2","IL6ST","EDNRB","IL6","IL1RN","PRCP","MYLK3","TNFRSF11A","REG3G","AHSG","OSMR","PTGES","SAA4","FCGR1A","ADAM8","IL6R","NLRP3","ADORA1","DNASE1L3","NPY5R","KLKB1","IL31RA","MBL2","SERPINF2","SCN11A","MRGPRX1","REG3A","CEBPB","SAA1","IL20RB","NLRP6","A2M","NUPR1","ANO6","CD163","CTNNBIP1","FCER1A","F2","FUT7","CXCR2","EXT1","F8","SIGIRR","FFAR3","PLSCR1","ZP3","SERPINA3","TRPV1","SERPINA1","SPN","ELANE","C2CD4A","FCGR3A","HLA-E","C2CD4B","IGHG1","DNASE1","ORM2","ORM1","TNF","UGT1A1","INS","SAA2-SAA4","HP" 
  #)) |> 
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
  geom_circle(radius.fixed = T) + # from THIS repo
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


ggsave("output/figures/vis_differential_proteomics_corrplot.png", width=8.4, height=8.4*2, dpi=600)




