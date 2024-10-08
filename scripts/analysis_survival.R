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



## mair ----


stats.mair_berghoff <- glass_od.metadata.array_samples |> 
  filter_OD_validation_idats(CONST_N_OD_VALIDATION_INCLUDED_SAMPLES) |> 
  dplyr::filter(grepl("35998208", array_notes)) |> 
  assertr::verify(!duplicated(patient_id)) |> 
  dplyr::filter(!is.na(patient_last_follow_up_date)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 21)
    return(.)
  })()
 


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


stats.primary <- stats.primary |> 
  dplyr::mutate(cgc = cut(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, breaks=2)) |> 
  dplyr::mutate(cgc = ifelse(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit < -4.9, "CGC[Astro] low", "CGC[Astro] high")) |> 
  dplyr::filter(!is.na(time_between_resection_and_last_event) & !is.na(patient_last_follow_up_event))


table(stats.primary$cgc)
table(stats.primary$resection_tumor_grade)





s1 <- survfit(Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ cgc, data = stats.primary)
s2 <- survfit(Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ resection_tumor_grade, data = stats.primary)


p2 = survminer::ggsurvplot(s2, data = stats.primary, pval = TRUE, risk.table=T, tables.y.text = FALSE, xlab="Time (days) last recurrence -> death")
p1 = survminer::ggsurvplot(s1, data = stats.primary, pval = TRUE, risk.table=T, tables.y.text = FALSE, xlab="Time (days) last recurrence -> death", newpage = FALSE)


grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}
ggsave("output/figures/analysis_survival_cgc_at_recurrence.pdf", width=8.3 / 2 * 0.8,height=2.6, scale=2, plot=p1)




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





## full split sweep ----
### primary ----


tmp.primary.CGC <- data.frame()
tmp.stats.primary <- stats.primary |> 
  dplyr::arrange(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(i = 1:dplyr::n())

for(split in ((2:nrow(stats.primary))-0.5) ) {
  tmp.stats.primary <- tmp.stats.primary |> 
    dplyr::mutate(cgc = ifelse(i < split, "CGC low", "CGC high"))
  
  tmp.primary.CGC <- rbind(tmp.primary.CGC, data.frame(
    k = split - 0.5,
    pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ cgc, tmp.stats.primary) |> purrr::pluck('pvalue'),
    name="CGC"
  ))
}
rm(tmp.stats.primary)



tmp.primary.GMS <- data.frame()
tmp.stats.primary <- stats.primary |> 
  dplyr::arrange(-array_GLASS_NL_g2_g3_sig) |> 
  dplyr::mutate(i = 1:dplyr::n())

for(split in ((2:nrow(stats.primary))-0.5) ) {
  tmp.stats.primary <- tmp.stats.primary |> 
    dplyr::mutate(gms = ifelse(i < split, "GMS lg", "GMS hg"))
  
  tmp.primary.GMS <- rbind(tmp.primary.GMS, data.frame(
    k = split - 0.5,
    pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ gms, tmp.stats.primary) |> purrr::pluck('pvalue'),
    name="GMS"
  ))
}
rm(tmp.stats.primary)




tmp.primary.PC2 <- data.frame()
tmp.stats.primary <- stats.primary |> 
  dplyr::arrange(-array_PC2) |> 
  dplyr::mutate(i = 1:dplyr::n())

for(split in ((2:nrow(stats.primary))-0.5) ) {
  tmp.stats.primary <- tmp.stats.primary |> 
    dplyr::mutate(pc2 = ifelse(i < split, "PC2 lg", "PC2 hg"))
  
  tmp.primary.PC2 <- rbind(tmp.primary.PC2, data.frame(
    k = split - 0.5,
    pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ pc2, tmp.stats.primary) |> purrr::pluck('pvalue'),
    name="PC2"
  ))
}
rm(tmp.stats.primary)




tmp.primary.PC3 <- data.frame()
tmp.stats.primary <- stats.primary |> 
  dplyr::arrange(-array_PC3) |> 
  dplyr::mutate(i = 1:dplyr::n())

for(split in ((2:nrow(stats.primary))-0.5) ) {
  tmp.stats.primary <- tmp.stats.primary |> 
    dplyr::mutate(pc3 = ifelse(i < split, "PC3 lg", "PC3 hg"))
  
  tmp.primary.PC3 <- rbind(tmp.primary.PC3, data.frame(
    k = split - 0.5,
    pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ pc3, tmp.stats.primary) |> purrr::pluck('pvalue'),
    name="PC3"
  ))
}
rm(tmp.stats.primary)




tmp.stats.primary <- stats.primary |> 
  dplyr::mutate(resection_tumor_grade = ifelse(resection_tumor_grade == 2, "G2", "G3"))




plt.primary <- rbind(tmp.primary.CGC,
                     tmp.primary.GMS,
                     tmp.primary.PC2,
                     tmp.primary.PC3,
                     
                     data.frame(
                       k=sum(tmp.stats.primary$resection_tumor_grade == "G2"),
                       pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ resection_tumor_grade, tmp.stats.primary) |> purrr::pluck('pvalue'),
                       name = "WHO Grade"))

rm(tmp.stats.primary, tmp.primary.CGC,tmp.primary.GMS, tmp.primary.PC2, tmp.primary.PC3)




p1 <- ggplot(plt.primary, aes(x=k, y=-log10(pval), col=name, group=name)) +
  geom_line(lwd=theme_nature_lwd) +
  geom_point() + 
  geom_hline(yintercept=-log10(0.01), col="red", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.05), col="gray60", lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.5), col="gray60", lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(1), col="gray60", lwd=theme_nature_lwd) +
  scale_y_continuous(breaks = -log10(c(0,0.01,0.05,0.5,1)),
                     labels =        c(0,0.01,0.05,0.5,1)) +
  theme_nature + 
  labs(y = "p-value logrank test", x=paste0("# primary samples (of ",nrow(stats.primary),") stratified low grade"))
p1




### last_rec ----



tmp.last_rec.CGC <- data.frame()
tmp.stats.last_rec <- stats.last_rec |> 
  dplyr::arrange(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(i = 1:dplyr::n())

for(split in ((2:nrow(stats.last_rec))-0.5) ) {
  tmp.stats.last_rec <- tmp.stats.last_rec |> 
    dplyr::mutate(cgc = ifelse(i < split, "CGC low", "CGC high"))
  
  tmp.last_rec.CGC <- rbind(tmp.last_rec.CGC, data.frame(
    k = split - 0.5,
    pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ cgc, tmp.stats.last_rec) |> purrr::pluck('pvalue'),
    name="CGC"
  ))
}
rm(tmp.stats.last_rec)



tmp.last_rec.GMS <- data.frame()
tmp.stats.last_rec <- stats.last_rec |> 
  dplyr::arrange(-array_GLASS_NL_g2_g3_sig) |> 
  dplyr::mutate(i = 1:dplyr::n())

for(split in ((2:nrow(stats.last_rec))-0.5) ) {
  tmp.stats.last_rec <- tmp.stats.last_rec |> 
    dplyr::mutate(gms = ifelse(i < split, "GMS lg", "GMS hg"))
  
  tmp.last_rec.GMS <- rbind(tmp.last_rec.GMS, data.frame(
    k = split - 0.5,
    pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ gms, tmp.stats.last_rec) |> purrr::pluck('pvalue'),
    name="GMS"
  ))
}
rm(tmp.stats.last_rec)




tmp.last_rec.PC2 <- data.frame()
tmp.stats.last_rec <- stats.last_rec |> 
  dplyr::arrange(-array_PC2) |> 
  dplyr::mutate(i = 1:dplyr::n())

for(split in ((2:nrow(stats.last_rec))-0.5) ) {
  tmp.stats.last_rec <- tmp.stats.last_rec |> 
    dplyr::mutate(pc2 = ifelse(i < split, "PC2 lg", "PC2 hg"))
  
  tmp.last_rec.PC2 <- rbind(tmp.last_rec.PC2, data.frame(
    k = split - 0.5,
    pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ pc2, tmp.stats.last_rec) |> purrr::pluck('pvalue'),
    name="PC2"
  ))
}
rm(tmp.stats.last_rec)




tmp.last_rec.PC3 <- data.frame()
tmp.stats.last_rec <- stats.last_rec |> 
  dplyr::arrange(-array_PC3) |> 
  dplyr::mutate(i = 1:dplyr::n())

for(split in ((2:nrow(stats.last_rec))-0.5) ) {
  tmp.stats.last_rec <- tmp.stats.last_rec |> 
    dplyr::mutate(pc3 = ifelse(i < split, "PC3 lg", "PC3 hg"))
  
  tmp.last_rec.PC3 <- rbind(tmp.last_rec.PC3, data.frame(
    k = split - 0.5,
    pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ pc3, tmp.stats.last_rec) |> purrr::pluck('pvalue'),
    name="PC3"
  ))
}
rm(tmp.stats.last_rec)




tmp.stats.last_rec <- stats.last_rec |> 
  dplyr::mutate(resection_tumor_grade = ifelse(resection_tumor_grade == 2, "G2", "G3"))



plt.last_rec <- rbind(tmp.last_rec.CGC,
                      tmp.last_rec.GMS,
                      tmp.last_rec.PC2,
                      tmp.last_rec.PC3,
                      
                      data.frame(
                        k=sum(tmp.stats.last_rec$resection_tumor_grade == "G2"),
                        pval = survival::survdiff(survival::Surv(time_between_resection_and_last_event, patient_last_follow_up_event) ~ resection_tumor_grade, tmp.stats.last_rec) |> purrr::pluck('pvalue'),
                        name = "WHO Grade"))

rm(tmp.stats.last_rec, tmp.last_rec.CGC,tmp.last_rec.GMS, tmp.last_rec.PC2, tmp.last_rec.PC3)



p2 <- ggplot(plt.last_rec, aes(x=k, y=-log10(pval), col=name, group=name)) +
  geom_hline(yintercept=-log10(0.01), col="red", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.05), col="gray60", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.5), col="gray60", lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(1), col="black", lwd=theme_nature_lwd) +
  geom_line(lwd=theme_nature_lwd, alpha=0.5) +
  geom_point(size=theme_nature_size/3) + 
  scale_y_continuous(breaks = -log10(c(0,0.01,0.05,0.5,1)),
                     labels =        c(0,0.01,0.05,0.5,1)) +
  theme_nature + 
  labs(y = "p-value logrank test", x=paste0("# last_rec samples (of ",nrow(stats.last_rec),") stratified low grade"))
p2


ggsave("output/figures/vis_analysis_survival__last_rec__split_sweep.pdf", width = 8.5*0.95/2, height = 2.5)




# simplified version for pres
ggplot(plt.last_rec |> 
         dplyr::filter(name %in% c("PC3", "GMS") == F)
         , aes(x=k, y=-log10(pval), col=name, group=name)) +
  geom_hline(yintercept=-log10(0.01), col="red", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.05), col="gray60", lty=2, lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(0.5), col="gray60", lwd=theme_nature_lwd) +
  geom_hline(yintercept=-log10(1), col="black", lwd=theme_nature_lwd) +
  geom_line(lwd=theme_nature_lwd, alpha=0.5) +
  geom_point(size=theme_nature_size/3) + 
  scale_y_continuous(breaks = -log10(c(0,0.01,0.05,0.5,1)),
                     labels =        c(0,0.01,0.05,0.5,1)) +
  theme_nature + 
  labs(y = "p-value logrank test", x=paste0("# last_rec samples (of ",nrow(stats.last_rec),") stratified low grade"))



ggsave("output/figures/vis_analysis_survival__last_rec__split_sweep_simplified.pdf", width = 8.5*0.95/2, height = 2.5)




### mair berghoff ----


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



