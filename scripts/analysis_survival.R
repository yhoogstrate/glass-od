#!/usr/bin/env R


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


## ggsurvplot export bugfix - https://github.com/kassambara/survminer/issues/152#issuecomment-938941051
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}


# set data ----
## from primary ----


stats.primary <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::mutate(cgc = cut(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, breaks=2)) |> 
  dplyr::filter(resection_number == 1) |> 
  assertr::verify(!duplicated(patient_id))


## from last recurrence ----


stats.last_rec <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(resection_number > 1) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::filter(resection_number == max(resection_number)) |> 
  dplyr::ungroup() |> 
  assertr::verify(!duplicated(patient_id))





# CGC Ac low high ----
## primary ----



stats.primary <- stats.primary |> 
  dplyr::mutate(cgc = cut(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, breaks=2)) |> 
  dplyr::mutate(cgc = ifelse(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit < -2.35, "CGC low", "CGC high")) |> 
  dplyr::filter(!is.na(time_between_resection_and_last_event) & !is.na(patient_last_follow_up_event))


table(stats.primary$cgc)
table(stats.primary$resection_tumor_grade)



#   dplyr::select(resection_id, cgc,  time_between_resection_and_last_event, patient_last_follow_up_event)
# 



s1 <- survfit(Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ cgc, data = stats.primary)
s2 <- survfit(Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ resection_tumor_grade, data = stats.primary)


#ggsurvfit(s1) + add_confidence_interval()
survminer::ggsurvplot(s1, data = stats.primary, pval = TRUE, risk.table=T, tables.y.text = FALSE)
survminer::ggsurvplot(s2, data = stats.primary, pval = TRUE, risk.table=T, tables.y.text = FALSE)



## last rec ----


stats.last_rec <- stats.last_rec |> 
  dplyr::mutate(cgc = cut(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, breaks=2)) |> 
  dplyr::mutate(cgc = ifelse(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit < -4.9, "CGC[Astro] low", "CGC[Astro] high")) |> 
  dplyr::filter(!is.na(time_between_resection_and_last_event) & !is.na(patient_last_follow_up_event))


table(stats.last_rec$cgc)
table(stats.last_rec$resection_tumor_grade)





s1 <- survfit(Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ cgc, data = stats.last_rec)
s2 <- survfit(Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ resection_tumor_grade, data = stats.last_rec)


p2 = survminer::ggsurvplot(s2, data = stats.last_rec, pval = TRUE, risk.table=T, tables.y.text = FALSE, xlab="Time (days) last recurrence -> death")
p1 = survminer::ggsurvplot(s1, data = stats.last_rec, pval = TRUE, risk.table=T, tables.y.text = FALSE, xlab="Time (days) last recurrence -> death", newpage = FALSE)


grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}
ggsave("output/figures/analysis_survival_cgc_at_recurrence.pdf", width=8.3 / 2 * 0.8,height=2.6, scale=2, plot=p1)


# GMS Ac low high ----
## primary ----

## last rec ----


stats.last_rec <- stats.last_rec |> 
  dplyr::mutate(gms = ifelse(array_GLASS_NL_g2_g3_sig < 1.7175, "GMS[GLASS-NL] low", "GMS[GLASS-NL] high")) |> 
  dplyr::filter(!is.na(time_between_resection_and_last_event) & !is.na(patient_last_follow_up_event))


table(stats.last_rec$gms)
table(stats.last_rec$resection_tumor_grade)





s1 <- survfit(Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ gms, data = stats.last_rec)
s2 <- survfit(Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ resection_tumor_grade, data = stats.last_rec)


p2 = survminer::ggsurvplot(s2, data = stats.last_rec, pval = TRUE, risk.table=T, tables.y.text = FALSE, xlab="Time (days) last recurrence -> death", newpage = FALSE)
p1 = survminer::ggsurvplot(s1, data = stats.last_rec, pval = TRUE, risk.table=T, tables.y.text = FALSE, xlab="Time (days) last recurrence -> death", newpage = FALSE)


grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}
ggsave("output/figures/analysis_survival_gms_at_recurrence.pdf", width=8.3 / 2 * 0.8,height=2.6, scale=2, plot=p1)



