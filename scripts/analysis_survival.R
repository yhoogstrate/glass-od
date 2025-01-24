#!/usr/bin/env R


library(survival)
library(survminer)


source('scripts/load_themes.R')


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


## ggsurvplot export bugfix - https://github.com/kassambara/survminer/issues/152#issuecomment-938941051
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}


# set data ----
## glass-od ----
### from primary ----


glass_od.stats.primary <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::mutate(CGC_scaled = scale(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, center=T)) |> 
  dplyr::filter(resection_number == 1) |> 
  assertr::verify(!duplicated(patient_id)) |> 
  
  (function(.) { print(dim(.)); assertthat::assert_that(nrow(.) == 76); return(.) })() |> 
  dplyr::filter(!is.na(time_between_resection_and_last_event)) |>  # excl none
  (function(.) { print(dim(.)); assertthat::assert_that(nrow(.) == 76); return(.) })() |> 
  assertr::verify(!is.na(time_between_resection_and_last_event) & !is.na(patient_last_follow_up_event)) |> 
  
  dplyr::mutate(resection_tumor_grade = paste0("G",resection_tumor_grade)) |> 
  dplyr::mutate(resection_tumor_grade = factor(resection_tumor_grade, levels=c("G2","G3"))) |> 
  
  dplyr::mutate(sex = factor(patient_sex, levels=c("female", "male")))




### from last recurrence ----



glass_od.stats.last_rec <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::mutate(CGC_scaled = scale(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, center=T)) |> 
  dplyr::filter(resection_number > 1) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::filter(resection_number == max(resection_number)) |> 
  dplyr::ungroup() |> 
  assertr::verify(!duplicated(patient_id)) |>
  
  (function(.) { print(dim(.)); assertthat::assert_that(nrow(.) == 103); return(.) })() |> 
  dplyr::filter(!is.na(time_between_resection_and_last_event)) |>  # excl 0091(-R2)
  (function(.) { print(dim(.)); assertthat::assert_that(nrow(.) == 102); return(.) })() |> 
  assertr::verify(!is.na(time_between_resection_and_last_event) & !is.na(patient_last_follow_up_event)) |> 
  
  dplyr::mutate(resection_tumor_grade = paste0("G",resection_tumor_grade)) |> 
  dplyr::mutate(resection_tumor_grade = factor(resection_tumor_grade, levels=c("G2","G3"))) |> 
  
  dplyr::mutate(`CDKN2ABHD` = as.factor(`CDKN2A/B deletion status last recurrence` == "HD")) |>  # no spaces or slashes allowed in surv packages >.<
  dplyr::mutate(`CDKN2ABiWT` = as.factor(`CDKN2A/B deletion status last recurrence` == "wt")) |> 
  dplyr::mutate(`CDKN2ABnSD_arm` = as.factor(`CDKN2A/B deletion status last recurrence` == "SD, arm")) |> 
  dplyr::mutate(`CDKN2ABnSD_arm_or_chr` = as.factor(`CDKN2A/B deletion status last recurrence` %in% c("SD, arm", "SD, chr"))) |> 
  dplyr::mutate(`CDKN2ABnHD_or_focal` = as.factor(`CDKN2A/B deletion status last recurrence` %in% c("HD", "SD, focal"))) |> 
  
  assertr::verify(!is.na(time_between_resection_and_last_event) & !is.na(patient_last_follow_up_event)) |> 
  
  dplyr::mutate(sex = factor(patient_sex, levels=c("female", "male")))



## validation ----
### from primary ----


# glass_od.stats.primary <- glass_od.metadata.array_samples |> 
#   filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
#   dplyr::mutate(cgc = cut(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, breaks=2)) |> 
#   dplyr::filter(resection_number == 1) |> 
#   assertr::verify(!duplicated(patient_id))


### from last recurrence ----


# glass_od.stats.last_rec <- glass_od.metadata.array_samples |> 
#   filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
#   dplyr::filter(resection_number > 1) |> 
#   dplyr::group_by(patient_id) |> 
#   dplyr::filter(resection_number == max(resection_number)) |> 
#   dplyr::ungroup() |> 
#   assertr::verify(!duplicated(patient_id))








# full split sweep ----
### primary ----


tmp.stats.primary <- glass_od.stats.primary |> 
  dplyr::arrange(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(i = 1:dplyr::n())


tmp.primary.CGC <- data.frame(
  k = 1:nrow(tmp.stats.primary),
  
  survdiff_MH_p = rep(1.0, nrow(tmp.stats.primary)),
  survdiff_PP_p = rep(1.0, nrow(tmp.stats.primary)),
  
  logrank_Z = rep(0.0, nrow(tmp.stats.primary)),
  logrank_p = rep(1.0, nrow(tmp.stats.primary)),
  
  coxph_logtest_p = rep(1.0, nrow(tmp.stats.primary)),
  coxph_waldtest_p = rep(1.0, nrow(tmp.stats.primary)),
  
  name = rep("patient", nrow(tmp.stats.primary)),
  
  censored = tmp.stats.primary$patient_last_follow_up_event == 0
)



tmp.s.primary <- survival::Surv(tmp.stats.primary$time_between_resection_and_last_event,
                                tmp.stats.primary$patient_last_follow_up_event)



for(split in ((2:nrow(glass_od.stats.primary))-0.5) ) {
  
  tmp.stats.primary <- tmp.stats.primary |> 
    dplyr::mutate(cgc = ifelse(i < split, "CGC low", "CGC high")) |> 
    dplyr::mutate(cgc = factor(cgc, levels=c("CGC low", "CGC high")))
  
  
  tmp.1 <- coin::logrank_test(tmp.s.primary ~ cgc, data=tmp.stats.primary)
  tmp.2 <- survival::survdiff(tmp.s.primary ~ cgc, data=tmp.stats.primary, rho=0)
  tmp.3 <- survival::coxph(   tmp.s.primary ~ cgc, data=tmp.stats.primary)
  tmp.4 <- survival::survdiff(tmp.s.primary ~ cgc, data=tmp.stats.primary, rho=1)
  
  
  df <- data.frame(
    
    k = split ,
    
    survdiff_MH_p = tmp.2 |> purrr::pluck('pvalue'),
    survdiff_PP_p = tmp.4 |> purrr::pluck('pvalue'),
    
    logrank_Z = tmp.1@statistic@teststatistic,
    logrank_p = coin::pvalue(tmp.1),
    
    coxph_logtest_p = summary(tmp.3)$logtest['pvalue'],
    coxph_waldtest_p = summary(tmp.3)$waldtest['pvalue'],
    
    name="CGC",
    
    censored = F
  )
  
  tmp.primary.CGC <- rbind(tmp.primary.CGC, df)
  rm(sd, tmp.1, tmp.2, tmp.3, tmp.4)
}



tmp.1_who_grade <- coin::logrank_test(tmp.s.primary ~ resection_tumor_grade, data=tmp.stats.primary)
tmp.2_who_grade <- survival::survdiff(tmp.s.primary ~ resection_tumor_grade, data=tmp.stats.primary, rh=0)
tmp.3_who_grade <- survival::coxph(   tmp.s.primary ~ resection_tumor_grade, data=tmp.stats.primary)
tmp.4_who_grade <- survival::survdiff(tmp.s.primary ~ resection_tumor_grade, data=tmp.stats.primary, rh=1)


rm(tmp.stats.last_rec)


plt.primary <- rbind(tmp.primary.CGC,
                     
                     data.frame(
                       k=sum(glass_od.stats.primary$resection_tumor_grade == "G2"),
                       
                       logrank_Z = tmp.1_who_grade@statistic@teststatistic,
                       logrank_p = coin::pvalue(tmp.1_who_grade),
                       
                       survdiff_MH_p = tmp.2_who_grade |> purrr::pluck('pvalue'),
                       survdiff_PP_p = tmp.2_who_grade |> purrr::pluck('pvalue'),
                       
                       coxph_logtest_p = summary(tmp.3_who_grade)$logtest['pvalue'],
                       coxph_waldtest_p = summary(tmp.3_who_grade)$waldtest['pvalue'],
                       
                       name = "WHO Grade",
                       censored=F)
)


rm(tmp.primary.CGC)
rm(tmp.1_who_grade, tmp.2_who_grade, tmp.3_who_grade, tmp.4_who_grade)





plt.facet <- plt.primary |> 
  dplyr::mutate(logrank_Z = NULL) |> 
  tidyr::pivot_longer(cols=c(survdiff_MH_p, survdiff_PP_p ,  logrank_p, coxph_logtest_p, coxph_waldtest_p ), names_to="pval_method", values_to="pvalue")




ggplot(plt.facet, aes(x=k, y=-log10(pvalue), col=name, group=name, pch=censored)) +
  facet_grid(cols = vars(pval_method)) +
  geom_hline(yintercept=-log10(0.01), col="red", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.05), col="gray60", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(1), col="black", lwd=theme_nature_lwd) +
  geom_line(data=subset(plt.facet, name != "patient"), lwd=theme_nature_lwd, alpha=0.5) +
  geom_point(data=subset(plt.facet, name == "patient"),y = -0.15, size=theme_nature_size/3) + 
  geom_point(data=subset(plt.facet, name != "patient"),size=theme_nature_size/3) + 
  theme_nature + 
  scale_y_continuous(breaks = -log10(c(0, 0.01, 0.05, 1)),
                     labels =        c(0, 0.01, 0.05, 1)) +
  scale_shape_manual(values=c(20, 3)) +
  labs(y = "p-value logrank test", x=paste0("# primary samples (of ",nrow(tmp.stats.primary),") stratified low grade"))







#### sandbox ----



# tmp.primary.GMS <- data.frame()
# tmp.stats.primary <- stats.primary |> 
#   dplyr::arrange(-array_GLASS_NL_g2_g3_sig) |> 
#   dplyr::mutate(i = 1:dplyr::n())
# 
# for(split in ((2:nrow(stats.primary))-0.5) ) {
#   tmp.stats.primary <- tmp.stats.primary |> 
#     dplyr::mutate(gms = ifelse(i < split, "GMS lg", "GMS hg"))
#   
#   tmp.primary.GMS <- rbind(tmp.primary.GMS, data.frame(
#     k = split - 0.5,
#     pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ gms, tmp.stats.primary) |> purrr::pluck('pvalue'),
#     name="GMS"
#   ))
# }
# rm(tmp.stats.primary)



# 
# tmp.primary.PC2 <- data.frame()
# tmp.stats.primary <- stats.primary |> 
#   dplyr::arrange(-array_PC2) |> 
#   dplyr::mutate(i = 1:dplyr::n())
# 
# for(split in ((2:nrow(stats.primary))-0.5) ) {
#   tmp.stats.primary <- tmp.stats.primary |> 
#     dplyr::mutate(pc2 = ifelse(i < split, "PC2 lg", "PC2 hg"))
#   
#   tmp.primary.PC2 <- rbind(tmp.primary.PC2, data.frame(
#     k = split - 0.5,
#     pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ pc2, tmp.stats.primary) |> purrr::pluck('pvalue'),
#     name="PC2"
#   ))
# }
# rm(tmp.stats.primary)
# 
# 
# 
# 
# tmp.primary.PC3 <- data.frame()
# tmp.stats.primary <- stats.primary |> 
#   dplyr::arrange(-array_PC3) |> 
#   dplyr::mutate(i = 1:dplyr::n())
# 
# for(split in ((2:nrow(stats.primary))-0.5) ) {
#   tmp.stats.primary <- tmp.stats.primary |> 
#     dplyr::mutate(pc3 = ifelse(i < split, "PC3 lg", "PC3 hg"))
#   
#   tmp.primary.PC3 <- rbind(tmp.primary.PC3, data.frame(
#     k = split - 0.5,
#     pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ pc3, tmp.stats.primary) |> purrr::pluck('pvalue'),
#     name="PC3"
#   ))
# }
# rm(tmp.stats.primary)



### last_rec ----



tmp.stats.last_rec <- glass_od.stats.last_rec |> 
  dplyr::arrange(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(i = 1:dplyr::n())


tmp.last_rec.CGC <- data.frame(
  k = 1:nrow(tmp.stats.last_rec),
  
  survdiff_MH_p = rep(1.0, nrow(tmp.stats.last_rec)),
  survdiff_PP_p = rep(1.0, nrow(tmp.stats.last_rec)),
  
  logrank_Z = rep(0.0, nrow(tmp.stats.last_rec)),
  logrank_p = rep(1.0, nrow(tmp.stats.last_rec)),
  
  coxph_logtest_p = rep(1.0, nrow(tmp.stats.last_rec)),
  coxph_waldtest_p = rep(1.0, nrow(tmp.stats.last_rec)),
  
  name = rep("patient", nrow(tmp.stats.last_rec)),
  
  censored = tmp.stats.last_rec$patient_last_follow_up_event == 0
)



tmp.s.last_rec <- survival::Surv(tmp.stats.last_rec$time_between_resection_and_last_event,
                                 tmp.stats.last_rec$patient_last_follow_up_event)



for(split in ((2:nrow(glass_od.stats.last_rec)) - 0.5)) {

  tmp.stats.last_rec <- tmp.stats.last_rec |> 
    dplyr::mutate(cgc = factor(ifelse(i < split, "CGC low", "CGC high"), levels=c("CGC low", "CGC high")))
  
  tmp.1 <- coin::logrank_test(tmp.s.last_rec ~ cgc, data=tmp.stats.last_rec, 
                              #ties="Hothorn-Lausen"
                              #type="Gehan-Breslow",
                              type="Fleming-Harrington"
                              )
  tmp.2 <- survival::survdiff(tmp.s.last_rec ~ cgc, data=tmp.stats.last_rec, rho=0) # fits ggsurvplot https://www.datacamp.com/tutorial/survival-analysis-R
  tmp.3 <- survival::coxph(   tmp.s.last_rec ~ cgc, data=tmp.stats.last_rec)
  tmp.4 <- survival::survdiff(tmp.s.last_rec ~ cgc, data=tmp.stats.last_rec, rho=1)
  
  
  df <- data.frame(
    
    k = split ,
    
    survdiff_MH_p = tmp.2 |> purrr::pluck('pvalue'),
    survdiff_PP_p = tmp.4 |> purrr::pluck('pvalue'),
    
    logrank_Z = tmp.1@statistic@teststatistic,
    logrank_p = coin::pvalue(tmp.1),
    
    coxph_logtest_p = summary(tmp.3)$logtest['pvalue'],
    coxph_waldtest_p = summary(tmp.3)$waldtest['pvalue'],
    
    name="CGC",
    
    censored = F
  )
  
  tmp.last_rec.CGC <- rbind(tmp.last_rec.CGC, df)
  rm(df, tmp.1, tmp.2, tmp.3, tmp.4)
}



tmp.1_who_grade <- coin::logrank_test(tmp.s.last_rec ~ resection_tumor_grade, data=tmp.stats.last_rec)
tmp.2_who_grade <- survival::survdiff(tmp.s.last_rec ~ resection_tumor_grade, data=tmp.stats.last_rec, rh=0)
tmp.3_who_grade <- survival::coxph(   tmp.s.last_rec ~ resection_tumor_grade, data=tmp.stats.last_rec)
tmp.4_who_grade <- survival::survdiff(tmp.s.last_rec ~ resection_tumor_grade, data=tmp.stats.last_rec, rh=1)

tmp.1_CDKN2ABHD <- coin::logrank_test(tmp.s.last_rec ~ `CDKN2ABHD`, data=tmp.stats.last_rec)
tmp.2_CDKN2ABHD <- survival::survdiff(tmp.s.last_rec ~ `CDKN2ABHD`, data=tmp.stats.last_rec, rho=0)
tmp.3_CDKN2ABHD <- survival::coxph(   tmp.s.last_rec ~ `CDKN2ABHD`, data=tmp.stats.last_rec)
tmp.4_CDKN2ABHD <- survival::survdiff(tmp.s.last_rec ~ `CDKN2ABHD`, data=tmp.stats.last_rec, rho=1)

tmp.1_CDKN2ABiWT <- coin::logrank_test(tmp.s.last_rec ~ `CDKN2ABiWT`, data=tmp.stats.last_rec)
tmp.2_CDKN2ABiWT <- survival::survdiff(tmp.s.last_rec ~ `CDKN2ABiWT`, data=tmp.stats.last_rec, rho=0)
tmp.3_CDKN2ABiWT <- survival::coxph(   tmp.s.last_rec ~ `CDKN2ABiWT`, data=tmp.stats.last_rec)
tmp.4_CDKN2ABiWT <- survival::survdiff(tmp.s.last_rec ~ `CDKN2ABiWT`, data=tmp.stats.last_rec, rho=1)

tmp.1_CDKN2ABnSD_arm <- coin::logrank_test(tmp.s.last_rec ~ `CDKN2ABnSD_arm`, data=tmp.stats.last_rec)
tmp.2_CDKN2ABnSD_arm <- survival::survdiff(tmp.s.last_rec ~ `CDKN2ABnSD_arm`, data=tmp.stats.last_rec, rho=0)
tmp.3_CDKN2ABnSD_arm <- survival::coxph(   tmp.s.last_rec ~ `CDKN2ABnSD_arm`, data=tmp.stats.last_rec)
tmp.4_CDKN2ABnSD_arm <- survival::survdiff(tmp.s.last_rec ~ `CDKN2ABnSD_arm`, data=tmp.stats.last_rec, rho=1)

tmp.1_CDKN2ABnSD_arm_or_chr <- coin::logrank_test(tmp.s.last_rec ~ `CDKN2ABnSD_arm_or_chr`, data=tmp.stats.last_rec)
tmp.2_CDKN2ABnSD_arm_or_chr <- survival::survdiff(tmp.s.last_rec ~ `CDKN2ABnSD_arm_or_chr`, data=tmp.stats.last_rec, rho=0)
tmp.3_CDKN2ABnSD_arm_or_chr <- survival::coxph(   tmp.s.last_rec ~ `CDKN2ABnSD_arm_or_chr`, data=tmp.stats.last_rec)
tmp.4_CDKN2ABnSD_arm_or_chr <- survival::survdiff(tmp.s.last_rec ~ `CDKN2ABnSD_arm_or_chr`, data=tmp.stats.last_rec, rho=1)

tmp.1_CDKN2ABnHD_or_focal <- coin::logrank_test(tmp.s.last_rec ~ `CDKN2ABnHD_or_focal`, data=tmp.stats.last_rec)
tmp.2_CDKN2ABnHD_or_focal <- survival::survdiff(tmp.s.last_rec ~ `CDKN2ABnHD_or_focal`, data=tmp.stats.last_rec, rho=0)
tmp.3_CDKN2ABnHD_or_focal <- survival::coxph(   tmp.s.last_rec ~ `CDKN2ABnHD_or_focal`, data=tmp.stats.last_rec)
tmp.4_CDKN2ABnHD_or_focal <- survival::survdiff(tmp.s.last_rec ~ `CDKN2ABnHD_or_focal`, data=tmp.stats.last_rec, rho=1)


rm(tmp.stats.last_rec)


plt.last_rec = tmp.last_rec.CGC


plt.last_rec <- rbind(tmp.last_rec.CGC,
                      
                      data.frame(
                        k=sum(glass_od.stats.last_rec$resection_tumor_grade == "G2"),

                        logrank_Z = tmp.1_who_grade@statistic@teststatistic,
                        logrank_p = coin::pvalue(tmp.1_who_grade),

                        survdiff_MH_p = tmp.2_who_grade |> purrr::pluck('pvalue'),
                        survdiff_PP_p = tmp.2_who_grade |> purrr::pluck('pvalue'),

                        coxph_logtest_p = summary(tmp.3_who_grade)$logtest['pvalue'],
                        coxph_waldtest_p = summary(tmp.3_who_grade)$waldtest['pvalue'],

                        name = "WHO Grade",
                        censored=F)
                      ,
  
                      data.frame(
                        k=sum(glass_od.stats.last_rec$`CDKN2ABHD` == "FALSE"),

                        logrank_Z = tmp.1_CDKN2ABHD@statistic@teststatistic,
                        logrank_p = coin::pvalue(tmp.1_CDKN2ABHD),

                        survdiff_MH_p = tmp.2_CDKN2ABHD |> purrr::pluck('pvalue'),
                        survdiff_PP_p = tmp.2_CDKN2ABHD |> purrr::pluck('pvalue'),
                        
                        coxph_logtest_p = summary(tmp.3_CDKN2ABHD)$logtest['pvalue'],
                        coxph_waldtest_p = summary(tmp.3_CDKN2ABHD)$waldtest['pvalue'],
                        
                        name = "CDKN2A/B HD",
                        censored=F)
                      ,

                      data.frame(
                        k=sum(glass_od.stats.last_rec$`CDKN2ABiWT` == "TRUE"),

                        logrank_Z = tmp.1_CDKN2ABiWT@statistic@teststatistic,
                        logrank_p = coin::pvalue(tmp.1_CDKN2ABiWT),

                        survdiff_MH_p = tmp.2_CDKN2ABiWT |> purrr::pluck('pvalue'),
                        survdiff_PP_p = tmp.2_CDKN2ABiWT |> purrr::pluck('pvalue'),

                        coxph_logtest_p = summary(tmp.3_CDKN2ABiWT)$logtest['pvalue'],
                        coxph_waldtest_p = summary(tmp.3_CDKN2ABiWT)$waldtest['pvalue'],
                        
                        name = "CDKN2A/B wt",
                        censored=F)
                      ,

                      data.frame(
                        k=sum(glass_od.stats.last_rec$`CDKN2ABnSD_arm` == "FALSE"),
                        
                        logrank_Z = tmp.1_CDKN2ABnSD_arm@statistic@teststatistic,
                        logrank_p = coin::pvalue(tmp.1_CDKN2ABnSD_arm),

                        survdiff_MH_p = tmp.2_CDKN2ABnSD_arm |> purrr::pluck('pvalue'),
                        survdiff_PP_p = tmp.2_CDKN2ABnSD_arm |> purrr::pluck('pvalue'),

                        coxph_logtest_p = summary(tmp.3_CDKN2ABnSD_arm)$logtest['pvalue'],
                        coxph_waldtest_p = summary(tmp.3_CDKN2ABnSD_arm)$waldtest['pvalue'],
                        
                        name = "CDKN2A/B SD (arm)",
                        censored=F)
                      ,

                      data.frame(
                        k=sum(glass_od.stats.last_rec$`CDKN2ABnSD_arm_or_chr` == "FALSE"),
                        
                        survdiff_MH_p = tmp.2_CDKN2ABnSD_arm_or_chr |> purrr::pluck('pvalue'),
                        survdiff_PP_p = tmp.2_CDKN2ABnSD_arm_or_chr |> purrr::pluck('pvalue'),
                        
                        logrank_Z = tmp.1_CDKN2ABnSD_arm_or_chr@statistic@teststatistic,
                        logrank_p = coin::pvalue(tmp.1_CDKN2ABnSD_arm_or_chr),
                        
                        coxph_logtest_p = summary(tmp.3_CDKN2ABnSD_arm_or_chr)$logtest['pvalue'],
                        coxph_waldtest_p = summary(tmp.3_CDKN2ABnSD_arm_or_chr)$waldtest['pvalue'],
                        
                        name = "CDKN2A/B SD (arm | chr)",
                        censored=F)
                      ,

                      data.frame(
                        k=sum(glass_od.stats.last_rec$`CDKN2ABnHD_or_focal` == "FALSE"),

                        survdiff_MH_p = tmp.2_CDKN2ABnHD_or_focal |> purrr::pluck('pvalue'),
                        survdiff_PP_p = tmp.2_CDKN2ABnHD_or_focal |> purrr::pluck('pvalue'),
                        
                        logrank_Z = tmp.1_CDKN2ABnHD_or_focal@statistic@teststatistic,
                        logrank_p = coin::pvalue(tmp.1_CDKN2ABnHD_or_focal),

                        coxph_logtest_p = summary(tmp.3_CDKN2ABnHD_or_focal)$logtest['pvalue'],
                        coxph_waldtest_p = summary(tmp.3_CDKN2ABnHD_or_focal)$waldtest['pvalue'],

                        name = "CDKN2A/B HD or focal",
                        censored=F)

                      )



rm(tmp.1_who_grade, tmp.2_who_grade, tmp.3_who_grade, tmp.4_who_grade)
rm(tmp.1_CDKN2ABHD, tmp.2_CDKN2ABHD, tmp.3_CDKN2ABHD, tmp.4_CDKN2ABHD)
rm(tmp.1_CDKN2ABiWT, tmp.2_CDKN2ABiWT, tmp.3_CDKN2ABiWT, tmp.4_CDKN2ABiWT)
rm(tmp.1_CDKN2ABnSD_arm, tmp.2_CDKN2ABnSD_arm, tmp.3_CDKN2ABnSD_arm, tmp.4_CDKN2ABnSD_arm)
rm(tmp.1_CDKN2ABnSD_arm_or_chr, tmp.2_CDKN2ABnSD_arm_or_chr, tmp.3_CDKN2ABnSD_arm_or_chr, tmp.4_CDKN2ABnSD_arm_or_chr)
rm(tmp.1_CDKN2ABnHD_or_focal, tmp.2_CDKN2ABnHD_or_focal, tmp.3_CDKN2ABnHD_or_focal, tmp.4_CDKN2ABnHD_or_focal)
rm(tmp.last_rec.CGC) # ,tmp.last_rec.GMS , tmp.last_rec.PC2, tmp.last_rec.PC3



plt.facet <- plt.last_rec |> 
  dplyr::mutate(logrank_Z = NULL) |> 
  tidyr::pivot_longer(cols=c(survdiff_MH_p, survdiff_PP_p ,  logrank_p, coxph_logtest_p, coxph_waldtest_p ), names_to="pval_method", values_to="pvalue") |> 
  dplyr::filter(pval_method != "coxph_waldtest_p")



ggplot(plt.facet, aes(x=k, y=-log10(pvalue), col=name, group=name, pch=censored)) +
  facet_grid(cols = vars(pval_method)) +
  geom_hline(yintercept=-log10(0.01), col="red", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.05), col="gray60", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(1), col="black", lwd=theme_nature_lwd) +
  geom_line(data=subset(plt.facet, name != "patient"), lwd=theme_nature_lwd, alpha=0.5) +
  geom_point(data=subset(plt.facet, name == "patient"),y = -0.15, size=theme_nature_size/3) + 
  geom_point(data=subset(plt.facet, name != "patient"),size=theme_nature_size/3) + 
  theme_nature + 
  scale_y_continuous(breaks = -log10(c(0, 0.01, 0.05, 1)),
                     labels =        c(0, 0.01, 0.05, 1)) +
  scale_shape_manual(values=c(20, 3)) +
  labs(y = "p-value logrank test", x=paste0("# last_rec samples (of ",nrow(glass_od.stats.last_rec),") stratified low grade"))





ggsave("output/figures/vis_analysis_survival__last_rec__split_sweep.pdf", width = 8.5 * 0.95 / 2, height = 2.5)


#### sandbox ----




# tmp.last_rec.GMS <- data.frame()
# tmp.stats.last_rec <- stats.last_rec |> 
#   dplyr::arrange(-array_GLASS_NL_g2_g3_sig) |> 
#   dplyr::mutate(i = 1:dplyr::n())
# 
# for(split in ((2:nrow(stats.last_rec))-0.5) ) {
#   tmp.stats.last_rec <- tmp.stats.last_rec |> 
#     dplyr::mutate(gms = ifelse(i < split, "GMS lg", "GMS hg"))
#   
#   tmp.last_rec.GMS <- rbind(tmp.last_rec.GMS, data.frame(
#     k = split - 0.5,
#     pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ gms, tmp.stats.last_rec) |> purrr::pluck('pvalue'),
#     name="GMS"
#   ))
# }
# rm(tmp.stats.last_rec)




# tmp.last_rec.PC2 <- data.frame()
# tmp.stats.last_rec <- stats.last_rec |> 
#   dplyr::arrange(-array_PC2) |> 
#   dplyr::mutate(i = 1:dplyr::n())
# 
# for(split in ((2:nrow(stats.last_rec))-0.5) ) {
#   tmp.stats.last_rec <- tmp.stats.last_rec |> 
#     dplyr::mutate(pc2 = ifelse(i < split, "PC2 lg", "PC2 hg"))
#   
#   tmp.last_rec.PC2 <- rbind(tmp.last_rec.PC2, data.frame(
#     k = split - 0.5,
#     pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ pc2, tmp.stats.last_rec) |> purrr::pluck('pvalue'),
#     name="PC2"
#   ))
# }
# rm(tmp.stats.last_rec)
# 
# 
# 
# 
# tmp.last_rec.PC3 <- data.frame()
# tmp.stats.last_rec <- stats.last_rec |> 
#   dplyr::arrange(-array_PC3) |> 
#   dplyr::mutate(i = 1:dplyr::n())
# 
# for(split in ((2:nrow(stats.last_rec))-0.5) ) {
#   tmp.stats.last_rec <- tmp.stats.last_rec |> 
#     dplyr::mutate(pc3 = ifelse(i < split, "PC3 lg", "PC3 hg"))
#   
#   tmp.last_rec.PC3 <- rbind(tmp.last_rec.PC3, data.frame(
#     k = split - 0.5,
#     pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ pc3, tmp.stats.last_rec) |> purrr::pluck('pvalue'),
#     name="PC3"
#   ))
# }
# rm(tmp.stats.last_rec)


### prim + last_rec merged ----


plt.facet <- rbind(
  plt.primary |> 
    dplyr::mutate(sample_selection = paste0("A. n=",(plt.primary |> dplyr::filter(name=="patient") |> nrow())," primary resections"))
  ,
  plt.last_rec |> 
    dplyr::mutate(sample_selection = paste0("B. n=",(plt.last_rec |> dplyr::filter(name=="patient") |> nrow())," last avaialable recurrent resection"))
  ) |> 
  dplyr::mutate(sample_selection = factor(sample_selection)) |> 
  dplyr::select(k, sample_selection, survdiff_MH_p, name, censored) |> 
  dplyr::filter(name != "CDKN2A/B SD (arm)") |> 
  dplyr::filter(name != "CDKN2A/B SD (arm | chr)") |> 
  dplyr::filter(name != "CDKN2A/B wt") |> 
  dplyr::filter(name != "CDKN2A/B SD (arm)") |> 
  dplyr::filter(name != "CDKN2A/B HD or focal")



ggplot(plt.facet, aes(x=k, y=-log10(survdiff_MH_p), col=name, group=name, pch=censored)) +
  facet_grid(cols = vars(sample_selection), scales="free_x", space="free_x") +
  geom_hline(yintercept=-log10(0.0000000001), col="gray80", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.01), col="red", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.05), col="gray60", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(1), col="black", lwd=theme_nature_lwd) +
  geom_line(data=subset(plt.facet, name != "patient"), lwd=theme_nature_lwd, alpha=0.5) +
  geom_point(data=subset(plt.facet, name == "patient" & censored),y = -0.25, size=theme_nature_size/8) + 
  geom_point(data=subset(plt.facet, name == "patient" & !censored),y = -0.25, size=theme_nature_size/3 * 0.75, alpha=0.5) + 
  geom_point(data=subset(plt.facet, name != "patient"),size=theme_nature_size/3) + 
  theme_nature + 
  scale_y_continuous(breaks = -log10(c(0,0.0000000001, 0.01, 0.05, 1)),
                     labels =        c(0,0.0000000001, 0.01, 0.05, 1)) +
  scale_shape_manual(values=c(16, 3)) +
  labs(y = "p-value logrank test survdiff", x=paste0("resection sample ordered by CGC"))




ggsave("output/figures/analysis_survival_sweep__prim__last_rec.pdf", width=8.5 * 0.975 * 0.9* 1.02,height=2.9 * 0.85)







### plot ----


plt <- rbind(
  plt.primary |> dplyr::mutate(facet = paste0("a. primaries (n=",nrow(stats.primary),")")),
  plt.last_rec |> dplyr::mutate(facet = paste0("b. last recurrences (n=",nrow(stats.last_rec),")"))
)


ggplot(plt, aes(x=k, y=-log10(pval), col=name, group=name)) +
  #facet_wrap(~facet, space="free_x") +
  facet_grid(cols = vars(facet), scales="free_x", space="free_x") +
  geom_line(lwd=theme_nature_lwd) +
  geom_point(size=theme_nature_size/3) + 
  geom_hline(yintercept=-log10(0.01), col="red", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.05), col="gray60", lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.5), col="gray60", lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(1), col="gray60", lwd=theme_nature_lwd) +
  scale_y_continuous(breaks = -log10(c(0,0.01,0.05,0.5,1)),
                     labels =        c(0,0.01,0.05,0.5,1)) +
  theme_nature + 
  labs(y = "p-value logrank test", x=paste0("P-value for cut-off samples stratified as low grade"))



ggsave("output/figures/analysis_survival_sweep.pdf", width=8.5 * 0.975,height=3.2)


# MV + TRT ----
## prim ----


tmp <- glass_od.stats.primary |> 
  dplyr::mutate(`Chemotherapy before` = factor(ifelse(resection_treatment_status_chemo, "Yes", "No"), levels=c("No", "Yes"))) |> 
  dplyr::mutate(`Radiotherapy before` = factor(ifelse(resection_treatment_status_radio, "Yes", "No"), levels=c("No", "Yes"))) |> 
  dplyr::mutate(`CGC` = CGC_scaled) |> 
  dplyr::rename(Sex = sex) |> 
  dplyr::rename(`WHO Grade` = resection_tumor_grade) |> 
  dplyr::mutate(placeholder_2factor = factor(rep("a",nrow(glass_od.stats.primary)),levels=c("a","b"))) |> 
  dplyr::mutate(placeholder_continuous = rep(1.0,nrow(glass_od.stats.primary)))



### all ----

svvl.prim <- survival::Surv(tmp$time_between_resection_and_last_event,
                            tmp$patient_last_follow_up_event)

fit.cox <- survival::coxph(svvl.prim ~
                             `CGC`+
                             `WHO Grade` + 
                             #`CDKN2A/B HD` +
                             Sex +
                             `Chemotherapy before` +
                             `Radiotherapy before`
                           ,
                           data = tmp)

survminer::ggforest(fit.cox, data = as.data.frame(tmp))


ggsave("output/figures/analysis_survival__forest__prim__all.pdf", width=8.5 * 0.975 * 0.67,height=3.2)




### CGC only ----


svvl.prim <- survival::Surv(tmp$time_between_resection_and_last_event,
                            tmp$patient_last_follow_up_event)

fit.cox <- survival::coxph(svvl.prim ~
                             `CGC`+
                             #`WHO Grade` + 
                             placeholder_2factor +
                             Sex +
                             `Chemotherapy before` +
                             `Radiotherapy before`
                           ,
                           data = tmp)

survminer::ggforest(fit.cox, data = as.data.frame(tmp))


ggsave("output/figures/analysis_survival__forest__prim__CGC-only.pdf", width=8.5 * 0.975 * 0.67,height=3.2)




### WHO only ----


svvl.prim <- survival::Surv(tmp$time_between_resection_and_last_event,
                            tmp$patient_last_follow_up_event)

fit.cox <- survival::coxph(svvl.prim ~
                             #`CGC` +
                             placeholder_continuous +
                             `WHO Grade` + 
                             Sex +
                             `Chemotherapy before` +
                             `Radiotherapy before`
                           ,
                           data = tmp)

survminer::ggforest(fit.cox, data = as.data.frame(tmp))


ggsave("output/figures/analysis_survival__forest__prim__who-only.pdf", width=8.5 * 0.975 * 0.67,height=3.2)




## last_rec ----


tmp <- glass_od.stats.last_rec |> 
  dplyr::mutate(`Chemotherapy before` = factor(ifelse(resection_treatment_status_chemo, "Yes", "No"), levels=c("No", "Yes"))) |> 
  dplyr::mutate(`Radiotherapy before` = factor(ifelse(resection_treatment_status_radio, "Yes", "No"), levels=c("No", "Yes"))) |> 
  dplyr::mutate(`CGC` = CGC_scaled) |> 
  dplyr::rename(Sex = sex) |> 
  dplyr::rename(`WHO Grade` = resection_tumor_grade) |> 
  dplyr::mutate(`CDKN2AB HD` = ifelse(`CDKN2A/B deletion status last recurrence` == "HD","Yes","No")) |> 
  dplyr::mutate(placeholder_2factor = factor(rep("a",nrow(glass_od.stats.last_rec)),levels=c("a","b"))) |> 
  dplyr::mutate(placeholder_continuous = rep(1.0,nrow(glass_od.stats.last_rec)))




### all ----


svvl.last_rec <- survival::Surv(tmp$time_between_resection_and_last_event,
                                tmp$patient_last_follow_up_event)

fit.cox <- survival::coxph(svvl.last_rec ~
                             `CGC` + 
                             Sex +
                             `WHO Grade` +
                             `CDKN2AB HD` +
                             `Chemotherapy before` +
                             `Radiotherapy before`
                           ,
                           data = as.data.frame(tmp))

survminer::ggforest(fit.cox, data = as.data.frame(tmp)) # https://stackoverflow.com/questions/75152731/r-ggforest-error-in-match-namesclabs-namesxi-names-do-not-match-previou

ggsave("output/figures/analysis_survival__forest__last_rec__all.pdf", width=8.5 * 0.975 * 0.67,height=3.2)



### CGC only ----


svvl.last_rec <- survival::Surv(tmp$time_between_resection_and_last_event,
                                tmp$patient_last_follow_up_event)

fit.cox <- survival::coxph(svvl.last_rec ~
                             `CGC` + 
                             Sex +
                             #`WHO Grade` +
                             `placeholder_2factor` + 
                           `CDKN2AB HD` +
                             `Chemotherapy before` +
                             `Radiotherapy before`
                           ,
                           data = as.data.frame(tmp))

survminer::ggforest(fit.cox, data = as.data.frame(tmp)) # https://stackoverflow.com/questions/75152731/r-ggforest-error-in-match-namesclabs-namesxi-names-do-not-match-previou

ggsave("output/figures/analysis_survival__forest__last_rec__CGC-only.pdf", width=8.5 * 0.975 * 0.67,height=3.2)



### WHO only ----


svvl.last_rec <- survival::Surv(tmp$time_between_resection_and_last_event,
                                tmp$patient_last_follow_up_event)

fit.cox <- survival::coxph(svvl.last_rec ~
                             #`CGC` + 
                             placeholder_continuous +
                             Sex +
                             `WHO Grade` +
                           `CDKN2AB HD` +
                             `Chemotherapy before` +
                             `Radiotherapy before`
                           ,
                           data = as.data.frame(tmp))

survminer::ggforest(fit.cox, data = as.data.frame(tmp)) # https://stackoverflow.com/questions/75152731/r-ggforest-error-in-match-namesclabs-namesxi-names-do-not-match-previou

ggsave("output/figures/analysis_survival__forest__last_rec__who-only.pdf", width=8.5 * 0.975 * 0.67,height=3.2)






tmp |> 
  dplyr::select(
    `WHO Grade`,
    #`Chemotherapy before`,
    `Radiotherapy before`
  ) |> 
  table()





### all ----


svvl.last_rec <- survival::Surv(tmp$time_between_resection_and_last_event,
                                tmp$patient_last_follow_up_event)

fit.cox <- survival::coxph(svvl.last_rec ~
                             `CGC` + 
                             Sex +
                             `WHO Grade` +
                             `CDKN2AB_HD` +
                             `Chemotherapy before` +
                             `Radiotherapy before`
                           ,
                           data = as.data.frame(tmp))

survminer::ggforest(fit.cox, data = as.data.frame(tmp)) # https://stackoverflow.com/questions/75152731/r-ggforest-error-in-match-namesclabs-namesxi-names-do-not-match-previou

ggsave("output/figures/analysis_survival__forest__last_rec__all.pdf", width=8.5 * 0.975 * 0.67,height=3.2)






# full sweep intersect Ki67 ----



tmp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::mutate(CGC_scaled = scale(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, center=T)) |> 
  dplyr::filter(resection_number > 1) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::filter(resection_number == max(resection_number)) |> 
  dplyr::ungroup() |> 
  assertr::verify(!duplicated(patient_id)) |>
  
  (function(.) { print(dim(.)); assertthat::assert_that(nrow(.) == 103); return(.) })() |> 
  dplyr::filter(!is.na(time_between_resection_and_last_event)) |>  # excl 0091(-R2)
  (function(.) { print(dim(.)); assertthat::assert_that(nrow(.) == 102); return(.) })() |> 
  assertr::verify(!is.na(time_between_resection_and_last_event) & !is.na(patient_last_follow_up_event)) |> 
  
  dplyr::mutate(resection_tumor_grade = paste0("G",resection_tumor_grade)) |> 
  dplyr::mutate(resection_tumor_grade = factor(resection_tumor_grade, levels=c("G2","G3"))) |> 
  
  dplyr::mutate(`CDKN2ABHD` = as.factor(`CDKN2A/B deletion status last recurrence` == "HD")) |>  # no spaces or slashes allowed in surv packages >.<
  dplyr::mutate(`CDKN2ABiWT` = as.factor(`CDKN2A/B deletion status last recurrence` == "wt")) |> 
  dplyr::mutate(`CDKN2ABnSD_arm` = as.factor(`CDKN2A/B deletion status last recurrence` == "SD, arm")) |> 
  dplyr::mutate(`CDKN2ABnSD_arm_or_chr` = as.factor(`CDKN2A/B deletion status last recurrence` %in% c("SD, arm", "SD, chr"))) |> 
  dplyr::mutate(`CDKN2ABnHD_or_focal` = as.factor(`CDKN2A/B deletion status last recurrence` %in% c("HD", "SD, focal"))) |> 
  
  assertr::verify(!is.na(time_between_resection_and_last_event) & !is.na(patient_last_follow_up_event)) |> 
  
  dplyr::mutate(sex = factor(patient_sex, levels=c("female", "male"))) |> 
  
  dplyr::filter(!is.na(staining_KI67_filename)) |> 
  dplyr::filter(!is.na(staining_KI67_Num.KI67pos)) |>  # one sample has N/A as file was too poor for analysis?
  
  as.data.frame() |> 
  dplyr::select(resection_id, isolation_id, patient_id, 
                
                time_between_resection_and_last_event, patient_last_follow_up_event,
                
                array_A_IDH_HG__A_IDH_LG_lr__lasso_fit,
                
                staining_KI67_filename,
                
                staining_KI67_lr_pos_neg_cells, staining_KI67_pos_per_area_um2, staining_KI67_pos_neg_cell_density
                
                )


#ki67 = tmp |> dplyr::select(resection_id, contains("i67")) 
#rank(ki67$staining_KI67_pos_per_detected_cells) == rank(ki67$staining_KI67_lr_pos_neg_cells)
#rank(ki67$staining_KI67_pos_per_area_um2)
##rank(ki67$staining_KI67_pos_neg_cell_density)


# 
# order_staining_KI67_lr_pos_neg_cells <- tmp |> 
#   dplyr::arrange(staining_KI67_lr_pos_neg_cells) |> 
#   dplyr::mutate(i = 1:dplyr::n())

order_staining_KI67_pos_per_area_um2 <- tmp |> 
  dplyr::arrange(staining_KI67_pos_per_area_um2) |> 
  dplyr::mutate(i = 1:dplyr::n())

order_cgc <- tmp |> 
  dplyr::arrange(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(i = 1:dplyr::n())



# svvl.s1  <- survival::Surv(order_staining_KI67_lr_pos_neg_cells$time_between_resection_and_last_event,
#                            order_staining_KI67_lr_pos_neg_cells$patient_last_follow_up_event)

svvl.s2  <- survival::Surv(order_staining_KI67_pos_per_area_um2$time_between_resection_and_last_event,
                           order_staining_KI67_pos_per_area_um2$patient_last_follow_up_event)

svvl.cgc <- survival::Surv(order_cgc$time_between_resection_and_last_event,
                           order_cgc$patient_last_follow_up_event)




df <- data.frame()


for(split in ((2:nrow(tmp)) - 0.5)) {
  
  # tmp.s1 <- order_staining_KI67_lr_pos_neg_cells |> 
  #   dplyr::mutate(ki67 = factor(ifelse(i < split, "Ki67 low", "Ki67 high"), levels=c("Ki67 low", "Ki67 high")))

  tmp.s2 <- order_staining_KI67_pos_per_area_um2 |> 
    dplyr::mutate(ki67 = factor(ifelse(i < split, "Ki67 low", "Ki67 high"), levels=c("Ki67 low", "Ki67 high")))

  tmp.cgc <- order_cgc |> 
    dplyr::mutate(cgc = factor(ifelse(i < split, "CGC low", "CGC high"), levels=c("CGC low", "CGC high")))


  #sv_diff.s1  <- survival::survdiff(svvl.s1  ~ ki67, data=tmp.s1)
  sv_diff.s2  <- survival::survdiff(svvl.s2  ~ ki67, data=tmp.s2)
  sv_diff.cgc <- survival::survdiff(svvl.cgc ~ cgc,  data=tmp.cgc)
  
  # 
  # df.s1 <- data.frame(
  #   k = split ,
  #   survdiff_p = sv_diff.s1 |> purrr::pluck('pvalue'),
  #   name="Ki67 log(pos/neg)",
  #   censored = F
  # )

  df.s2 <- data.frame(
    k = split ,
    survdiff_p = sv_diff.s2 |> purrr::pluck('pvalue'),
    name="Ki67 pos / um2",
    censored = F
  )
  
  df.cgc <- data.frame(
    k = split ,
    survdiff_p = sv_diff.cgc |> purrr::pluck('pvalue'),
    name="CGC",
    censored = F
  )
  
  df <- rbind(df, df.s2, df.cgc)
  
  #rm(tmp.s1, sv_diff.s1, df.s1)
  rm(tmp.s2, sv_diff.s2, df.s2)
  rm(tmp.cgc, sv_diff.cgc, df.cgc)

}





ggplot(df, aes(x=k, y=-log10(survdiff_p), col=name, group=name)) +
  #facet_grid(cols = vars(pval_method)) +
  geom_hline(yintercept=-log10(0.01), col="red", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.05), col="gray60", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(1), col="black", lwd=theme_nature_lwd) +
  geom_line(data=subset(df, name != "patient"), lwd=theme_nature_lwd, alpha=0.5) +
  geom_point(data=subset(df, name == "patient"),y = -0.15, size=theme_nature_size/3) + 
  geom_point(data=subset(df, name != "patient"),size=theme_nature_size/3) + 
  theme_nature + 
  scale_y_continuous(breaks = -log10(c(0, 0.01, 0.05, 1)),
                     labels =        c(0, 0.01, 0.05, 1)) +
  #scale_shape_manual(values=c(20, 3)) +
  labs(y = "p-value logrank test", x=paste0("# last_rec samples with Ki-67 staining (of ",nrow(order_cgc),") stratified low grade / low Ki67"))



ggsave("output/figures/analysis_survival__Ki67_sweep.pdf", width = 11 * 1/3 * 0.85, height = 2.36)




# mair berghoff ----


tmp.mair_berghoff.CGC <- data.frame()
tmp.stats.mair_berghoff <- stats.mair_berghoff |> 
  dplyr::arrange(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(i = 1:dplyr::n())

for(split in ((2:nrow(stats.mair_berghoff))-0.5) ) {
  tmp.stats.mair_berghoff <- tmp.stats.mair_berghoff |> 
    dplyr::mutate(cgc = ifelse(i < split, "CGC low", "CGC high"))
  
  tmp.mair_berghoff.CGC <- rbind(tmp.mair_berghoff.CGC, data.frame(
    k = split - 0.5,
    pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ cgc, tmp.stats.mair_berghoff) |> purrr::pluck('pvalue'),
    name="CGC"
  ))
}
rm(tmp.stats.mair_berghoff)



tmp.mair_berghoff.GMS <- data.frame()
tmp.stats.mair_berghoff <- stats.mair_berghoff |> 
  dplyr::arrange(-array_GLASS_NL_g2_g3_sig) |> 
  dplyr::mutate(i = 1:dplyr::n())

for(split in ((2:nrow(stats.mair_berghoff))-0.5) ) {
  tmp.stats.mair_berghoff <- tmp.stats.mair_berghoff |> 
    dplyr::mutate(gms = ifelse(i < split, "GMS lg", "GMS hg"))
  
  tmp.mair_berghoff.GMS <- rbind(tmp.mair_berghoff.GMS, data.frame(
    k = split - 0.5,
    pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ gms, tmp.stats.mair_berghoff) |> purrr::pluck('pvalue'),
    name="GMS"
  ))
}
rm(tmp.stats.mair_berghoff)




tmp.mair_berghoff.PC2 <- data.frame()
tmp.stats.mair_berghoff <- stats.mair_berghoff |> 
  dplyr::arrange(-array_PC2) |> 
  dplyr::mutate(i = 1:dplyr::n())

for(split in ((2:nrow(stats.mair_berghoff))-0.5) ) {
  tmp.stats.mair_berghoff <- tmp.stats.mair_berghoff |> 
    dplyr::mutate(pc2 = ifelse(i < split, "PC2 lg", "PC2 hg"))
  
  tmp.mair_berghoff.PC2 <- rbind(tmp.mair_berghoff.PC2, data.frame(
    k = split - 0.5,
    pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ pc2, tmp.stats.mair_berghoff) |> purrr::pluck('pvalue'),
    name="PC2"
  ))
}
rm(tmp.stats.mair_berghoff)




tmp.mair_berghoff.PC3 <- data.frame()
tmp.stats.mair_berghoff <- stats.mair_berghoff |> 
  dplyr::arrange(-array_PC3) |> 
  dplyr::mutate(i = 1:dplyr::n())

for(split in ((2:nrow(stats.mair_berghoff))-0.5) ) {
  tmp.stats.mair_berghoff <- tmp.stats.mair_berghoff |> 
    dplyr::mutate(pc3 = ifelse(i < split, "PC3 lg", "PC3 hg"))
  
  tmp.mair_berghoff.PC3 <- rbind(tmp.mair_berghoff.PC3, data.frame(
    k = split - 0.5,
    pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ pc3, tmp.stats.mair_berghoff) |> purrr::pluck('pvalue'),
    name="PC3"
  ))
}
rm(tmp.stats.mair_berghoff)




tmp.stats.mair_berghoff <- stats.mair_berghoff |> 
  dplyr::mutate(resection_tumor_grade = ifelse(resection_tumor_grade == 2, "G2", "G3"))



plt.mair_berghoff <- rbind(tmp.mair_berghoff.CGC,
                           tmp.mair_berghoff.GMS,
                           #tmp.mair_berghoff.PC2,
                           #tmp.mair_berghoff.PC3,
                           
                           data.frame(
                             k=sum(tmp.stats.mair_berghoff$resection_tumor_grade == "G2"),
                             pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ resection_tumor_grade, tmp.stats.mair_berghoff) |> purrr::pluck('pvalue'),
                             name = "WHO Grade"))

rm(tmp.stats.mair_berghoff, tmp.mair_berghoff.CGC,tmp.mair_berghoff.GMS)



ggplot(plt.mair_berghoff, aes(x=k, y=-log10(pval), col=name, group=name)) +
  geom_line(lwd=theme_nature_lwd) +
  geom_point() + 
  geom_hline(yintercept=-log10(0.01), col="red", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.05), col="gray60", lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.5), col="gray60", lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(1), col="gray60", lwd=theme_nature_lwd) +
  scale_y_continuous(breaks = -log10(c(0,0.01,0.05,0.5,1)),
                     labels =        c(0,0.01,0.05,0.5,1)) +
  theme_nature + 
  labs(y = "p-value logrank test", x=paste0("# mair_berghoff samples (of ",nrow(stats.mair_berghoff),") stratified low grade"))







s1 <- survfit(Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ gms, data = stats.last_rec)
s2 <- survfit(Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ resection_tumor_grade, data = stats.last_rec)


p2 = survminer::ggsurvplot(s2, data = stats.last_rec, pval = TRUE, risk.table=T, tables.y.text = FALSE, xlab="Time (days) last recurrence -> death", newpage = FALSE)
p1 = survminer::ggsurvplot(s1, data = stats.last_rec, pval = TRUE, risk.table=T, tables.y.text = FALSE, xlab="Time (days) last recurrence -> death", newpage = FALSE)


grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}
ggsave("output/figures/analysis_survival_gms_at_recurrence.pdf", width=8.3 / 2 * 0.8,height=2.6, scale=2, plot=p1)



