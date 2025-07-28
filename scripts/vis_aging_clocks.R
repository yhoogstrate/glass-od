#!/usr/bin/env R


# load ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# if(!exists('data.mvalues.probes')) {
#   source('scripts/load_mvalues_hq_samples.R')
# }
# 


library(ggplot2)


# corr + forest mainly PCs ----
## load data ----


plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::mutate(`age_at_diagnosis_days` = as.Date(patient_diagnosis_date) - as.Date(patient_birth_date)) |> 
  #dplyr::filter(!is.na(age_at_diagnosis_days)) |> 
  #dplyr::mutate(age_at_diagnosis_days = age_at_diagnosis_days - min(na.omit(age_at_diagnosis_days))) |> 
  #dplyr::mutate(age_at_diagnosis_days = as.numeric(age_at_diagnosis_days)) |> 
  #dplyr::mutate(time_tissue_in_ffpe =  ifelse(isolation_material == "ffpe", time_between_resection_and_array, 0)) |> 
  
  dplyr::rename_with(~ gsub("array_","",.x)) |>
  dplyr::rename_with(~ gsub("dnaMethyAge__","dnaMethyAge: ",.x)) |>
  
  dplyr::select(resection_id,
                
                PC1, 
                PC2, 
                PC3, 
                PC4, 
                PC5, 
                PC6, 

                `dnaMethyAge: PCHorvathS2018`,
                `dnaMethyAge: PCHannumG2013`,
                percentage.detP.signi,
                `qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`,
                GLASS_NL_g2_g3_sig, A_IDH_HG__A_IDH_LG_lr__lasso_fit,
                
  ) |> 
  tibble::column_to_rownames('resection_id') |>  
  
  dplyr::mutate(`percentage.detP.signi` = log(percentage.detP.signi)) |> 
  dplyr::rename(`log(% det-P)` = percentage.detP.signi) |> 
  
  dplyr::rename(`GLASS-NL p-r signature` = GLASS_NL_g2_g3_sig) |> 
  dplyr::rename(`CGC[Ac]` = `A_IDH_HG__A_IDH_LG_lr__lasso_fit`) |> 
  dplyr::rename(`QC: Spec I GT MM 6` = `qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`)




## corrplot ----


# <3
ggcorrplot(plt, abs=T) +
  labs(caption=paste0("n=",nrow(plt)," samples"),
       subtitle=format_subtitle("correlation PCs"))

ggsave("output/figures/vis_aging_clocks__PC_clustering.pdf", height = 1.9, width=4)




## stats grade ----


data <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_first_G2_and_last_G3(154) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"_remainder")))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = factor(ifelse(resection_tumor_grade == 2, "Grade2", "Grade3"), levels=c("Grade2", "Grade3"))) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 154) 
    return(.)
  })() |> 
  
  dplyr::select(resection_id, gr.status, patient) |> 
  dplyr::left_join(
    plt |> 
      tibble::rownames_to_column('resection_id'),
    by=c('resection_id'='resection_id')
  ) |> 
  tibble::column_to_rownames('resection_id')


n <- nrow(data)


design <- model.matrix(~~patient + gr.status, data=data)
fit <- limma::lmFit(data |> dplyr::mutate(patient = NULL) |> dplyr::mutate(gr.status = NULL) |> t() |>  as.data.frame()
                    , design)
fit <- limma::eBayes(fit, trend=T)
stats <- limma::topTable(fit,
                         n=nrow(data),
                         coef="gr.statusGrade3",
                         sort.by = "none",
                         confint=T,
                         adjust.method="fdr") |> 
  tibble::rownames_to_column('covar') 


order = data.frame(covar = c(
  "PC6",
  "PC4",
  "PC5",
  "log(% det-P)",
  "PC1",
  "QC: Spec I GT MM 6"    ,
  "GLASS-NL p-r signature" ,
  "PC2"     ,
  "CGC[Ac]"  ,
  "dnaMethyAge: PCHannumG2013" ,
  "dnaMethyAge: PCHorvathS2018",
  "PC3" 
)) |> 
  dplyr::mutate(rank = 1:dplyr::n())


stats <- stats |> 
  dplyr::left_join(order, by=c('covar'='covar'), suffix=c('',''))




stats.all <- rbind(
  stats |>
    dplyr::mutate(logFC = 0) |> 
    dplyr::mutate(t = 0) |> 
    dplyr::mutate(type = "zero")
  ,
  stats |>
    dplyr::mutate(type = "datapoint")
  ,
  stats |>
    dplyr::mutate(logFC = CI.L) |> 
    dplyr::mutate(type = "CI")
  ,
  stats |>
    dplyr::mutate(logFC = CI.R) |> 
    dplyr::mutate(type = "CI")
)




ggplot(stats.all, aes(x = t, y=reorder(covar, -rank), label=paste0("q=",format.pval(adj.P.Val,nsmall=3, digits = 1)), group=covar)) +
  geom_line(data=subset(stats.all, type %in% c("datapoint", "zero")), lwd=theme_nature_lwd) + 
  #geom_line(data=subset(stats.all, type %in% c("CI")),lwd=theme_nature_lwd*3) + 
  geom_point(data=subset(stats.all, type == "datapoint"), size=theme_nature_size/3) +
  geom_text(data=subset(stats.all, type %in% c("datapoint")), size=theme_nature_size, family=theme_nature_font_family) +
  #ggrepel::geom_text_repel(
  #  data=subset(stats.all, type %in% c("datapoint")), size=theme_nature_size, nudge_y = 0, family=theme_nature_font_family) + 
  labs(x="t WHO grade 2 - 3",
       y = NULL,
       caption = paste0("First Grade 2 and last Grade 3: n=",n, " samples")) +
  theme_nature


ggsave("output/figures/vis_aging_clocks__PC_clustering_limma.pdf", width=3.5, height=2.028)





# corr + forest with QC & clocks ----
## corrplot ----


plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::mutate(`age_at_diagnosis_days` = as.Date(patient_diagnosis_date) - as.Date(patient_birth_date)) |> 
  dplyr::filter(!is.na(age_at_diagnosis_days)) |> 
  dplyr::mutate(age_at_diagnosis_days = age_at_diagnosis_days - min(na.omit(age_at_diagnosis_days))) |> 
  dplyr::mutate(age_at_diagnosis_days = as.numeric(age_at_diagnosis_days)) |> 
  dplyr::mutate(time_tissue_in_ffpe =  ifelse(isolation_material == "ffpe", time_between_resection_and_array, 0)) |> 
  dplyr::select(resection_id,
                
                array_epiTOC2_tnsc, array_epiTOC2_hypoSC, contains("array_dnaMethyAge"), array_RepliTali,
                array_percentage.detP.signi, array_PC1, `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`,
                time_tissue_in_ffpe,
                array_GLASS_NL_g2_g3_sig,array_PC2, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit,
                
                array_PC3,
                age_at_diagnosis_days,
                array_PC4,
                array_PC5,
                array_PC6,
                
                array_median.overall.methylation
  ) |> 
  tibble::column_to_rownames('resection_id') |> 
  dplyr::filter(!is.na(time_tissue_in_ffpe)) |> 
  
  dplyr::mutate(array_GLASS_OD_g2_g3_sig2 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig2 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_g2_g3_sig3 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig3 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_g2_g3_sig4 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig4 = NULL) |> 
  
  
  dplyr::mutate(array_PC3 = -1 * array_PC3) |> 
  dplyr::mutate(array_percentage.detP.signi = log(array_percentage.detP.signi)) |> 
  
  dplyr::mutate(`CGC[Ac]` = -1 * array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit = NULL) |> 
  dplyr::mutate(`-1 * array_dnaMethyAge__ZhangY2017` = -1 * array_dnaMethyAge__ZhangY2017 , array_dnaMethyAge__ZhangY2017 = NULL) |> 
  dplyr::mutate(`-1 * array_epiTOC2_hypoSC` = -1 * array_epiTOC2_hypoSC, array_epiTOC2_hypoSC = NULL) |>  # seems at inversed scale
  dplyr::mutate(`-1 * array_dnaMethyAge__LuA2019` = -1 * array_dnaMethyAge__LuA2019, array_dnaMethyAge__LuA2019 = NULL) |>  # seems at inversed scale
  dplyr::mutate(`-1 * QC: SPECIFICITY I GT MM 6` = -1 * `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`, `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`=NULL) |> 
  dplyr::select(-array_dnaMethyAge__YangZ2016) |>  # epiTOC1, near identical to epiTOC2
  dplyr::select(-array_dnaMethyAge__PCHorvathS2013) |>  # very similar to its 2018 equivalent
  dplyr::mutate(`array_qc_BISULFITE_CONVERSION_I_Beta_I-3_Beta_larger_0.12_0.18` = NULL) |> # contains N/A's
  dplyr::mutate(`array_qc_BISULFITE_CONVERSION_I_Beta_I-6_Beta_larger_0.2_0.3` = NULL) |> # contains N/A's
  dplyr::rename(`% detP significant probes (log)` = array_percentage.detP.signi) |> 
  
  dplyr::mutate(`-1 * CGC[Ac]` = -1 * `CGC[Ac]`, `CGC[Ac]`= NULL) |> 
  dplyr::mutate(`array_dnaMethyAge__ZhangY2017` = -1 * `-1 * array_dnaMethyAge__ZhangY2017`, `-1 * array_dnaMethyAge__ZhangY2017`=NULL) |>
  dplyr::mutate(`-1 * array_GLASS_NL_g2_g3_sig` = -1 * array_GLASS_NL_g2_g3_sig, array_GLASS_NL_g2_g3_sig = NULL) |>
  
  dplyr::mutate(`-1 * array_median.overall.methylation` = -1 * array_median.overall.methylation) |> 
  dplyr::mutate(array_median.overall.methylation = NULL) |> 

  
  dplyr::mutate(`-1 * array_PC2` = -1 * array_PC2, array_PC2 = NULL)


colnames(plt) <- gsub("array_","",colnames(plt))



p1 = ggcorrplot(plt)
p1
ggsave(plot = p1, "output/figures/vis_aging_clocks__ggcorrplot.pdf", width= 8.5 * 0.975 / 2 , height = 4)




order <- c(
  "PC5",
  "PC4",
  "PC6",
  
  
  "-1 * PC2",
  "-1 * GLASS_NL_g2_g3_sig",
  "-1 * CGC[Ac]",
  "-1 * median.overall.methylation",
  
  "dnaMethyAge__ZhangY2017",
  
  "dnaMethyAge__LuA2023p3",
  "dnaMethyAge__LuA2023p2",

  "dnaMethyAge__PCPhenoAge",
  "dnaMethyAge__PCHorvathS2018",
  "dnaMethyAge__PCHannumG2013",
  "PC3",
  "dnaMethyAge__ShirebyG2020",
  "dnaMethyAge__McEwenL2019",
  "dnaMethyAge__LuA2023p1",
  "-1 * dnaMethyAge__LuA2019",
  "dnaMethyAge__Cortex_common",
  "dnaMethyAge__CBL_common",
  "RepliTali",
  "dnaMethyAge__CBL_specific",
  "dnaMethyAge__HorvathS2018",
  "dnaMethyAge__ZhangQ2019",
  "dnaMethyAge__HannumG2013",
  "epiTOC2_tnsc",
  "dnaMethyAge__LevineM2018",
  "age_at_diagnosis_days",
  
  "-1 * epiTOC2_hypoSC",
  "PC1",
  "% detP significant probes (log)",
  "-1 * QC: SPECIFICITY I GT MM 6",
  "time_tissue_in_ffpe"
)





## forests ----
### grade ----


tmp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_first_G2_and_last_G3(154) |> 
  dplyr::mutate(`age_at_diagnosis_days` = as.Date(patient_diagnosis_date) - as.Date(patient_birth_date)) |> 
  dplyr::filter(!is.na(age_at_diagnosis_days)) |> 
  dplyr::mutate(age_at_diagnosis_days = age_at_diagnosis_days - min(na.omit(age_at_diagnosis_days))) |> 
  dplyr::mutate(age_at_diagnosis_days = as.numeric(age_at_diagnosis_days)) |> 
  dplyr::mutate(time_tissue_in_ffpe =  ifelse(isolation_material == "ffpe", time_between_resection_and_array, 0)) |> 
  dplyr::select(resection_id,
                resection_tumor_grade,
                
                array_median.overall.methylation,
                
                array_epiTOC2_tnsc, array_epiTOC2_hypoSC, contains("array_dnaMethyAge"), array_RepliTali,
                array_percentage.detP.signi, array_PC1, `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`,
                time_tissue_in_ffpe,
                array_GLASS_NL_g2_g3_sig, array_PC2, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit,
                
                array_PC3,
                age_at_diagnosis_days,
                
                array_PC4,
                array_PC5,
                array_PC6
  ) |> 
  tibble::column_to_rownames('resection_id') |> 
  dplyr::filter(!is.na(time_tissue_in_ffpe)) |> 
  
  dplyr::mutate(array_GLASS_OD_g2_g3_sig2 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig2 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_g2_g3_sig3 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig3 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_g2_g3_sig4 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig4 = NULL) |> 
  
  dplyr::mutate(`-1 * array_median.overall.methylation` = -1 * array_median.overall.methylation) |> 
  dplyr::mutate(array_median.overall.methylation = NULL) |> 
  
  
  dplyr::mutate(array_PC3 = -1 * array_PC3) |> 
  dplyr::mutate(array_percentage.detP.signi = log(array_percentage.detP.signi)) |> 
  
  dplyr::mutate(`CGC[Ac]` = -1 * array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit = NULL) |> 
  dplyr::mutate(`-1 * array_dnaMethyAge__ZhangY2017` = -1 * array_dnaMethyAge__ZhangY2017 , array_dnaMethyAge__ZhangY2017 = NULL) |> 
  dplyr::mutate(`-1 * array_epiTOC2_hypoSC` = -1 * array_epiTOC2_hypoSC, array_epiTOC2_hypoSC = NULL) |>  # seems at inversed scale
  dplyr::mutate(`-1 * array_dnaMethyAge__LuA2019` = -1 * array_dnaMethyAge__LuA2019, array_dnaMethyAge__LuA2019 = NULL) |>  # seems at inversed scale
  dplyr::mutate(`-1 * QC: SPECIFICITY I GT MM 6` = -1 * `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`, `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`=NULL) |> 
  dplyr::select(-array_dnaMethyAge__YangZ2016) |>  # epiTOC1, near identical to epiTOC2
  dplyr::select(-array_dnaMethyAge__PCHorvathS2013) |>  # very similar to its 2018 equivalent
  dplyr::mutate(`array_qc_BISULFITE_CONVERSION_I_Beta_I-3_Beta_larger_0.12_0.18` = NULL) |> # contains N/A's
  dplyr::mutate(`array_qc_BISULFITE_CONVERSION_I_Beta_I-6_Beta_larger_0.2_0.3` = NULL) |> # contains N/A's
  dplyr::rename(`% detP significant probes (log)` = array_percentage.detP.signi) |> 
  
  dplyr::mutate(resection_tumor_grade = factor(paste0("Grade",resection_tumor_grade), levels=c("Grade2","Grade3"))) |> 
  
  dplyr::mutate(`-1 * CGC[Ac]` = -1 * `CGC[Ac]`, `CGC[Ac]`= NULL) |> 
  dplyr::mutate(`array_dnaMethyAge__ZhangY2017` = -1 * `-1 * array_dnaMethyAge__ZhangY2017`, `-1 * array_dnaMethyAge__ZhangY2017`=NULL) |>
  dplyr::mutate(`-1 * array_GLASS_NL_g2_g3_sig` = -1 * array_GLASS_NL_g2_g3_sig, array_GLASS_NL_g2_g3_sig = NULL) |>
  dplyr::mutate(`-1 * array_PC2` = -1 * array_PC2, array_PC2 = NULL)

colnames(tmp) <- gsub("array_","",colnames(tmp))


data <- tmp |>
  dplyr::mutate(`resection_tumor_grade` = NULL) |> 
  dplyr::mutate(array_PC1 = NULL) |> 
  dplyr::mutate(PC1 = NULL) |> 
  as.data.frame() |> 
  scale(center=T, scale=T)  |> 
  t() |>
  as.data.frame() 


design <- model.matrix(~factor(resection_tumor_grade), data=tmp) # no correcion, you can' time in ffpe for quality
#design <- model.matrix(~age_at_diagnosis_days + factor(resection_tumor_grade), data=tmp)
#design <- model.matrix(~factor(resection_tumor_grade), data=tmp)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit, trend=T)
stats <- limma::topTable(fit,
                         n=nrow(data),
                         coef="factor(resection_tumor_grade)Grade3",
                         sort.by = "none",
                         confint=T,
                         adjust.method="fdr") |> 
  tibble::rownames_to_column('covar') 


stats = rbind(stats,
              
              data.frame(covar=c('array_PC1'), logFC=0, CI.L=0, CI.R=0, AveExpr=0,t=0,P.Value=NA, adj.P.Val=NA,B=NA)
              
) |> 
  dplyr::mutate(covar =  gsub("array_", "", covar))


plt <- data.frame(covar = order) |>
  dplyr::mutate(y = dplyr::n():1) |> 
  assertr::verify(covar %in% stats$covar) |> 
  assertr::verify(stats$covar %in% covar) |> 
  dplyr::left_join(stats, by=c('covar'='covar')) 


plt <- rbind(
  plt |>
    dplyr::mutate(logFC = 0) |> 
    dplyr::mutate(t = 0) |> 
    dplyr::mutate(type = "zero")
  ,
  plt |>
    dplyr::mutate(type = "datapoint")
  ,
  plt |>
    dplyr::mutate(logFC = CI.L) |> 
    dplyr::mutate(type = "CI")
  ,
  plt |>
    dplyr::mutate(logFC = CI.R) |> 
    dplyr::mutate(type = "CI")
) |> 
  dplyr::mutate(stars = gtools::stars.pval(adj.P.Val))


#order


ggplot(plt, aes(x = t, y=reorder(covar, y), label=stars, group=covar)) +
  geom_line(data=subset(plt, type %in% c("datapoint", "zero")),lwd=theme_nature_lwd) + 
  #geom_line(data=subset(plt, type %in% c("CI")),lwd=theme_nature_lwd*3) + 
  geom_point(data=subset(plt, type == "datapoint"), size=theme_nature_size/3) +
  geom_text(data=subset(plt, type %in% c("datapoint")), x = 7.5, size=theme_nature_size, family=theme_nature_font_family) + 
  xlim(-5,7.7) + 
  labs(y = NULL, x="t statistic WHO Grade 2 vs. Grade 3") +
  theme_nature


ggsave("output/figures/vis_aging_clocks__forest__grade.pdf", width=2.75,height = 2.665)


# https://support.bioconductor.org/p/37524/
# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/forest.html
# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/metareg.html


### prim-rec ----


tmp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(179) |> 
  dplyr::mutate(`age_at_diagnosis_days` = as.Date(patient_diagnosis_date) - as.Date(patient_birth_date)) |> 
  dplyr::filter(!is.na(age_at_diagnosis_days)) |> 
  dplyr::mutate(age_at_diagnosis_days = age_at_diagnosis_days - min(na.omit(age_at_diagnosis_days))) |> 
  dplyr::mutate(age_at_diagnosis_days = as.numeric(age_at_diagnosis_days)) |> 
  dplyr::mutate(time_tissue_in_ffpe =  ifelse(isolation_material == "ffpe", time_between_resection_and_array, 0)) |> 
  dplyr::select(resection_id,
                resection_number,
                
                array_median.overall.methylation,
                
                array_epiTOC2_tnsc, array_epiTOC2_hypoSC, contains("array_dnaMethyAge"), array_RepliTali,
                array_percentage.detP.signi, array_PC1, `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`,
                time_tissue_in_ffpe,
                array_GLASS_NL_g2_g3_sig, array_PC2, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit,
                
                array_PC3,
                age_at_diagnosis_days,
                
                array_PC4,
                array_PC5,
                array_PC6
  ) |> 
  tibble::column_to_rownames('resection_id') |> 
  dplyr::filter(!is.na(time_tissue_in_ffpe)) |> 
  
  dplyr::mutate(array_GLASS_OD_g2_g3_sig2 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig2 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_g2_g3_sig3 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig3 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_g2_g3_sig4 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig4 = NULL) |> 
  
  
  dplyr::mutate(`-1 * array_median.overall.methylation` = -1 * array_median.overall.methylation) |> 
  dplyr::mutate(array_median.overall.methylation = NULL) |> 
  
  dplyr::mutate(array_PC3 = -1 * array_PC3) |> 
  dplyr::mutate(array_percentage.detP.signi = log(array_percentage.detP.signi)) |> 
  
  dplyr::mutate(`CGC[Ac]` = -1 * array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit = NULL) |> 
  dplyr::mutate(`-1 * array_dnaMethyAge__ZhangY2017` = -1 * array_dnaMethyAge__ZhangY2017 , array_dnaMethyAge__ZhangY2017 = NULL) |> 
  dplyr::mutate(`-1 * array_epiTOC2_hypoSC` = -1 * array_epiTOC2_hypoSC, array_epiTOC2_hypoSC = NULL) |>  # seems at inversed scale
  dplyr::mutate(`-1 * array_dnaMethyAge__LuA2019` = -1 * array_dnaMethyAge__LuA2019, array_dnaMethyAge__LuA2019 = NULL) |>  # seems at inversed scale
  dplyr::mutate(`-1 * QC: SPECIFICITY I GT MM 6` = -1 * `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`, `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`=NULL) |> 
  dplyr::select(-array_dnaMethyAge__YangZ2016) |>  # epiTOC1, near identical to epiTOC2
  dplyr::select(-array_dnaMethyAge__PCHorvathS2013) |>  # very similar to its 2018 equivalent
  dplyr::mutate(`array_qc_BISULFITE_CONVERSION_I_Beta_I-3_Beta_larger_0.12_0.18` = NULL) |> # contains N/A's
  dplyr::mutate(`array_qc_BISULFITE_CONVERSION_I_Beta_I-6_Beta_larger_0.2_0.3` = NULL) |> # contains N/A's
  dplyr::rename(`% detP significant probes (log)` = array_percentage.detP.signi) |> 
  
  dplyr::mutate(resection_prim_rec = factor(ifelse(resection_number == 1,"primary", "recurrent"), levels=c("primary","recurrent"))) |> 
  dplyr::mutate(resection_number = NULL) |> 
  
  dplyr::mutate(`-1 * CGC[Ac]` = -1 * `CGC[Ac]`, `CGC[Ac]`= NULL) |> 
  dplyr::mutate(`array_dnaMethyAge__ZhangY2017` = -1 * `-1 * array_dnaMethyAge__ZhangY2017`, `-1 * array_dnaMethyAge__ZhangY2017`=NULL) |>
  dplyr::mutate(`-1 * array_GLASS_NL_g2_g3_sig` = -1 * array_GLASS_NL_g2_g3_sig, array_GLASS_NL_g2_g3_sig = NULL) |>
  dplyr::mutate(`-1 * array_PC2` = -1 * array_PC2, array_PC2 = NULL)

colnames(tmp) <- gsub("array_","",colnames(tmp))


data <- tmp |>
  dplyr::mutate(`resection_prim_rec` = NULL) |> 
  dplyr::mutate(array_PC1 = NULL) |> 
  dplyr::mutate(PC1 = NULL) |> 
  as.data.frame() |> 
  scale(center=T, scale=T)  |> 
  t() |>
  as.data.frame() 


design <- model.matrix(~ factor(resection_prim_rec), data=tmp)
#design <- model.matrix(~age_at_diagnosis_days + factor(resection_tumor_grade), data=tmp)
#design <- model.matrix(~factor(resection_tumor_grade), data=tmp)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit, trend=T)
stats <- limma::topTable(fit,
                         n=nrow(data),
                         coef="factor(resection_prim_rec)recurrent",
                         sort.by = "none",
                         confint=T,
                         adjust.method="fdr") |> 
  tibble::rownames_to_column('covar') 


stats = rbind(stats,
              
              data.frame(covar=c('array_PC1'), logFC=0, CI.L=0, CI.R=0, AveExpr=0,t=0,P.Value=NA, adj.P.Val=NA,B=NA)
              
) |> 
  dplyr::mutate(covar =  gsub("array_", "", covar))


plt <- data.frame(covar = order) |>
  dplyr::mutate(y = dplyr::n():1) |> 
  assertr::verify(covar %in% stats$covar) |> 
  assertr::verify(stats$covar %in% covar) |> 
  dplyr::left_join(stats, by=c('covar'='covar')) 


plt <- rbind(
  plt |>
    dplyr::mutate(logFC = 0) |> 
    dplyr::mutate(t = 0) |> 
    dplyr::mutate(type = "zero")
  ,
  plt |>
    dplyr::mutate(type = "datapoint")
  ,
  plt |>
    dplyr::mutate(logFC = CI.L) |> 
    dplyr::mutate(type = "CI")
  ,
  plt |>
    dplyr::mutate(logFC = CI.R) |> 
    dplyr::mutate(type = "CI")
) |> 
  dplyr::mutate(stars = gtools::stars.pval(adj.P.Val))


#order


ggplot(plt, aes(x = t, y=reorder(covar, y), label=stars, group=covar)) +
  geom_line(data=subset(plt, type %in% c("datapoint", "zero")),lwd=theme_nature_lwd) + 
  #geom_line(data=subset(plt, type %in% c("CI")),lwd=theme_nature_lwd*3) + 
  geom_point(data=subset(plt, type == "datapoint"), size=theme_nature_size/3) +
  geom_text(data=subset(plt, type %in% c("datapoint")), x = 7.5, size=theme_nature_size, family=theme_nature_font_family) + 
  xlim(-5, 7.7) + 
  labs(y = NULL, x="t statistic primary vs. recurrent") +
  theme_nature


ggsave("output/figures/vis_aging_clocks__forest__prim_rec.pdf", width=2.75,height = 2.665)


# https://support.bioconductor.org/p/37524/
# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/forest.html
# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/metareg.html




# DMP time_in_tissue_ffpe

# plot order: 
plt.order <- data.frame(name = c(
  "dnaMethyAge__ZhangY2017",
  "CGC[Ac]",
  "GLASS_NL_g2_g3_sig",
  "PC2",
  "time_tissue_in_ffpe",
  "QC: SPECIFICITY I GT MM 6",
  "detP significant probes (log)",
  "PC1",
  "epiTOC2_hypoSC",
  "dnaMethyAge__LuA2023p2",
  "dnaMethyAge__LuA2023p3",
  "dnaMethyAge__LuA2019",
  "dnaMethyAge__CBL_specific",
  "RepliTali",
  "epiTOC2_tnsc",
  "dnaMethyAge__HannumG2013",
  "dnaMethyAge__CBL_common",
  "dnaMethyAge__Cortex_common",
  "dnaMethyAge__McEwenL2019",
  "dnaMethyAge__ZhangQ2019",
  "dnaMethyAge__HorvathS2018",
  "dnaMethyAge__LevineM2018",
  "dnaMethyAge__LuA2023p1",
  "PC3",
  "dnaMethyAge__ShirebyG2020",
  "dnaMethyAge__PCHorvathS2018",
  "dnaMethyAge__PCHannumG2013",
  "dnaMethyAge__PCPhenoAge"
)) |> 
  dplyr::mutate(order = 1:dplyr::n()) |> 
  dplyr::mutate(colname = dplyr::recode(name, 
    `PC1` = "DMP__PCs__pp_nc__PC1_adj.P.Val",
    `PC2` = "DMP__PCs__pp_nc__PC2_adj.P.Val",
    `PC3` = "DMP__PCs__pp_nc__PC3_adj.P.Val",
    `detP significant probes (log)` = "DMP__pct_detP_signi__pp_nc__adj.P.Val",
    `CGC[Ac]` = "DMP__AcCGAP__pp_nc__adj.P.Val",
    `GLASS_NL_g2_g3_sig` = "DMP__GLASS_NL_signature__pp_nc__adj.P.Val",
    `RepliTali` = "DMP__RepliTali__up_nc__adj.P.Val",
    
    `QC: SPECIFICITY I GT MM 6` = "DMP__mnp_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5__pp_nc__adj.P.Val",
    `epiTOC2_hypoSC` = "DMP__epiTOC2_hypoSC__up_nc__adj.P.Val",
    `epiTOC2_tnsc` = "DMP__epiTOC2_tnsc__up_nc__adj.P.Val",
    
    `dnaMethyAge__CBL_common` = "DMP__dnaMethyAge__CBL_common__up_nc__adj.P.Val",
    `dnaMethyAge__CBL_specific` = "DMP__dnaMethyAge__CBL_specific__up_nc__adj.P.Val",
    `dnaMethyAge__Cortex_common` = "DMP__dnaMethyAge__Cortex_common__up_nc__adj.P.Val",
    `dnaMethyAge__HannumG2013` = "DMP__dnaMethyAge__HannumG2013__up_nc__adj.P.Val",
    `dnaMethyAge__HorvathS2018` = "DMP__dnaMethyAge__HorvathS2018__up_nc__adj.P.Val",
    `dnaMethyAge__LevineM2018` = "DMP__dnaMethyAge__LevineM2018__up_nc__adj.P.Val",
    `dnaMethyAge__LuA2019` = "DMP__dnaMethyAge__LuA2019__up_nc__adj.P.Val",
    `dnaMethyAge__LuA2023p1` = "DMP__dnaMethyAge__LuA2023p1__up_nc__adj.P.Val",
    `dnaMethyAge__LuA2023p2` = "DMP__dnaMethyAge__LuA2023p2__up_nc__adj.P.Val",
    `dnaMethyAge__LuA2023p3` = "DMP__dnaMethyAge__LuA2023p3__up_nc__adj.P.Val",
    `dnaMethyAge__McEwenL2019` = "DMP__dnaMethyAge__McEwenL2019__up_nc__adj.P.Val",
    `dnaMethyAge__PCHannumG2013` = "DMP__dnaMethyAge__PCHannumG2013__up_nc__adj.P.Val",
    `dnaMethyAge__PCHorvathS2013` = "DMP__dnaMethyAge__PCHorvathS2013__up_nc__adj.P.Val",
    `dnaMethyAge__PCHorvathS2018` = "DMP__dnaMethyAge__PCHorvathS2018__up_nc__adj.P.Val",
    `dnaMethyAge__PCPhenoAge` = "DMP__dnaMethyAge__PCPhenoAge__up_nc__adj.P.Val",
    `dnaMethyAge__ShirebyG2020` = "DMP__dnaMethyAge__ShirebyG2020__up_nc__adj.P.Val",
    `dnaMethyAge__YangZ2016` = "DMP__dnaMethyAge__YangZ2016__up_nc__adj.P.Val",
    `dnaMethyAge__ZhangQ2019` = "DMP__dnaMethyAge__ZhangQ2019__up_nc__adj.P.Val",
    `dnaMethyAge__ZhangY2017` = "DMP__dnaMethyAge__ZhangY2017__up_nc__adj.P.Val",
    `time_tissue_in_ffpe` = "DMP__FFPE_decay_time__pp_nc__adj.P.Val"
                            ))




plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  tibble::tibble() |> 
  dplyr::select(
    `DMP__PCs__pp_nc__PC1_adj.P.Val`,
    `DMP__PCs__pp_nc__PC2_adj.P.Val`,
    `DMP__PCs__pp_nc__PC3_adj.P.Val`,
    `DMP__pct_detP_signi__pp_nc__adj.P.Val`,
    `DMP__AcCGAP__pp_nc__adj.P.Val`,
    `DMP__GLASS_NL_signature__pp_nc__adj.P.Val`,
    `DMP__RepliTali__up_nc__adj.P.Val`,
    
    `DMP__mnp_qc_SPECIFICITY_I_GT_Mismatch_6_PM_Red_smaller_NA_0_5__pp_nc__adj.P.Val`,
    `DMP__epiTOC2_hypoSC__up_nc__adj.P.Val`,
    `DMP__epiTOC2_tnsc__up_nc__adj.P.Val`,
    
    `DMP__dnaMethyAge__CBL_common__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__CBL_specific__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__Cortex_common__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__HannumG2013__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__HorvathS2018__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__LevineM2018__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__LuA2019__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__LuA2023p1__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__LuA2023p2__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__LuA2023p3__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__McEwenL2019__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__PCHannumG2013__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__PCHorvathS2018__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__PCPhenoAge__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__ShirebyG2020__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__ZhangQ2019__up_nc__adj.P.Val`,
    `DMP__dnaMethyAge__ZhangY2017__up_nc__adj.P.Val`,
    `DMP__FFPE_decay_time__pp_nc__adj.P.Val`
   ) |> 
  tidyr::pivot_longer(cols = dplyr::everything(), values_to = "adj.P.Val", names_to = "colname") |> 
  dplyr::mutate(one.min.adj.P.val = 1 - adj.P.Val) |> 
  dplyr::left_join(plt.order, by=c('colname'='colname')) |> 
  dplyr::rename(order_y = order) |> 
  dplyr::group_by(colname) |> 
  dplyr::mutate(order_x = rank(adj.P.Val)) |> 
  dplyr::ungroup() |> 
  dplyr::arrange(order_y, -order_x) |> 
  dplyr::mutate(facet = factor(name, levels = plt.order |> dplyr::arrange(order) |> dplyr::pull(name))) |> 
  dplyr::filter(1:dplyr::n() %% 450 == 1)



plt2 <- rbind(plt |> 
               dplyr::mutate(group = paste0(order_x,colname)),
             plt  |> 
               dplyr::mutate(group = paste0(order_x,colname)) |>
               dplyr::mutate(adj.P.Val=0,
                             one.min.adj.P.val=0))



ggplot(plt2, aes(x=order_x, y=one.min.adj.P.val, group=group, col = adj.P.Val < 0.01)) +
  facet_grid(rows = vars(facet)) + # ,  levels=plt.order |> dplyr::arrange(order) |> dplyr::pull(name)
  geom_line(lwd=theme_nature_lwd) +
  theme_nature +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(strip.text.y.right = element_text(angle = 0)) +
  coord_cartesian(ylim = c(0,1))




# generic statistics ----

## prim - rec ----
tmp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(179) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"remainder")))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 
  
  dplyr::mutate(`age_at_diagnosis_days` = as.Date(patient_diagnosis_date) - as.Date(patient_birth_date)) |> 
  dplyr::filter(!is.na(age_at_diagnosis_days)) |> 
  dplyr::mutate(age_at_diagnosis_days = age_at_diagnosis_days - min(na.omit(age_at_diagnosis_days))) |> 
  dplyr::mutate(age_at_diagnosis_days = as.numeric(age_at_diagnosis_days)) |> 
  dplyr::mutate(time_tissue_in_ffpe =  ifelse(isolation_material == "ffpe", time_between_resection_and_array, 0)) |> 
  dplyr::select(resection_id,
                resection_tumor_grade,
                patient,
                pr.status,
                array_PC1,
                array_dnaMethyAge__PCHorvathS2018
  ) |> 
  tibble::column_to_rownames('resection_id') |> 
  dplyr::rename(PCHorvathS2018 = array_dnaMethyAge__PCHorvathS2018)


data <- tmp |>
  dplyr::select(PCHorvathS2018) |> 
  t() |> 
  as.data.frame()


stopifnot(rownames(tmp) == colnames(data))


design <- model.matrix(~array_PC1 +# factor(patient) + 
                          factor(patient) + 
                           factor(pr.status), data=tmp)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit, trend=T)
stats <- limma::topTable(fit,
                         n=nrow(data),
                         coef="factor(pr.status)recurrence",
                         sort.by = "none",
                         confint=T,
                         adjust.method="fdr") |> 
  tibble::rownames_to_column('covar') 
stats





## WHO grade ----


tmp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(180) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"remainder")))) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2", "Grade3")) |> 
  
  dplyr::mutate(`age_at_diagnosis_days` = as.Date(patient_diagnosis_date) - as.Date(patient_birth_date)) |> 
  dplyr::filter(!is.na(age_at_diagnosis_days)) |> 
  dplyr::mutate(age_at_diagnosis_days = age_at_diagnosis_days - min(na.omit(age_at_diagnosis_days))) |> 
  dplyr::mutate(age_at_diagnosis_days = as.numeric(age_at_diagnosis_days)) |> 
  dplyr::mutate(time_tissue_in_ffpe =  ifelse(isolation_material == "ffpe", time_between_resection_and_array, 0)) |> 
  dplyr::select(resection_id,
                resection_tumor_grade,
                patient,
                gr.status,
                array_PC1,
                array_PC2,
                array_PC3,
                array_dnaMethyAge__PCHorvathS2018,
                array_A_IDH_HG__A_IDH_LG_lr__lasso_fit
  ) |> 
  tibble::column_to_rownames('resection_id') |> 
  dplyr::rename(PCHorvathS2018 = array_dnaMethyAge__PCHorvathS2018)


data <- tmp |>
  dplyr::select(PCHorvathS2018) |> 
  #dplyr::select(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  #dplyr::select(array_PC3) |> 
  t() |> 
  as.data.frame()


stopifnot(rownames(tmp) == colnames(data))


design <- model.matrix(~array_PC1 + factor(patient) + factor(gr.status), data=tmp)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit, trend=T)
stats <- limma::topTable(fit,
                         n=nrow(data),
                         coef="factor(gr.status)Grade3",
                         sort.by = "none",
                         confint=T,
                         adjust.method="fdr") |> 
  tibble::rownames_to_column('covar') 
stats



