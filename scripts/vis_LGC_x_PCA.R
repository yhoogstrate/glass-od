#!/usr/bin/env R

# load ----


if(!exists('glass_od.metadata.idats')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}




# GLASS-OD / OD ----


metadata <- glass_od.metadata.idats |> 
  filter_GLASS_OD_idats(163)




## PCA ----


plt <- metadata


plt.split <- rbind(
  plt |>
    dplyr::mutate(col = as.factor(paste0("Grade ",resection_tumor_grade))) |> 
    dplyr::mutate(facet = "Histological grade")
  ,
  plt |> 
    dplyr::mutate(col = ifelse(mnp_predictBrain_v2.0.1_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", mnp_predictBrain_v2.0.1_cal_class)) |> 
    dplyr::mutate(facet = "MNP CNS Classifier 114b/2.0.1")  ,
  plt |> 
    dplyr::mutate(col = ifelse(mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", mnp_predictBrain_v12.8_cal_class)) |> 
    dplyr::mutate(facet = "MNP CNS Classifier 12.8"),
  
  plt |> 
    dplyr::mutate(col = ifelse(isolation_person_name == "USA / Duke", "Batch [US]", "Batch [EU]")) |> 
    dplyr::mutate(facet = "Batch")
)


ggplot(plt.split, aes(x=A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=-PC2, col=col)) + 
  facet_grid(cols = vars(facet), scales = "free", space="free") +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  geom_point() +
  theme_bw() + 
  scale_color_manual(values=c(
    "Grade 2"="blue","O_IDH"="blue",    "Batch [EU]"="darkgreen",
    "OLIGOSARC_IDH" = "orange",
    "NA" = "gray",
    "other" = "purple",
    "Grade 3"="red","A_IDH_HG"="red",     "Batch [US]"="brown"
  ))



ggplot(plt, aes(x=A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=-median.overall.methylation)) + 
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  geom_point() +
  theme_bw() 



ggplot(plt, aes(x=PC2, y=median.overall.methylation, col=isolation_person_name)) + 
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), cor.coef.name ="rho") +
  geom_point() +
  theme_bw() 


ggplot(plt, aes(x=PC2, y=median.glass_nl_supervised.methylation, col=isolation_person_name)) + 
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  geom_point() +
  theme_bw() 




## stats ----


stats <- metadata |> 
  dplyr::left_join(data.all.pca.obj.data, by=c('sentrix_id'='sentrix_id'), suffix=c('','')) |> 
  dplyr::filter(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(resection_tumor_grade__hg = ifelse(resection_tumor_grade == 3, 1 , 0)) |> 
  dplyr::mutate(batch_us = ifelse(isolation_person_name == "USA / Duke", 1, 0))  |> 
  dplyr::mutate(PC2inv = -PC2)



### t-test LR ----

t.test(
  stats |>
    dplyr::filter(resection_tumor_grade == 2) |> 
    dplyr::pull(A_IDH_HG__A_IDH_LG_lr__lasso_fit),
  stats |>
    dplyr::filter(resection_tumor_grade == 3) |> 
    dplyr::pull(A_IDH_HG__A_IDH_LG_lr__lasso_fit)
)

### t-test PC2 ----


t.test(
  stats |>
    dplyr::filter(resection_tumor_grade == 2) |> 
    dplyr::pull(PC2),
  stats |>
    dplyr::filter(resection_tumor_grade == 3) |> 
    dplyr::pull(PC2)
)

### logistic LR ----


model <- glm(resection_tumor_grade__hg ~ A_IDH_HG__A_IDH_LG_lr__lasso_fit, data = stats,family = binomial)
summary(model)
#sjPlot::plot_model(model)

pval <- as.data.frame(summary(model)$coefficients)$`Pr(>|z|)`[2]

#Data frame with hp in ascending order
Predicted_data <- data.frame(A_IDH_HG__A_IDH_LG_lr__lasso_fit=seq(min(stats$A_IDH_HG__A_IDH_LG_lr__lasso_fit), max(stats$A_IDH_HG__A_IDH_LG_lr__lasso_fit), len=500))

# Fill predicted values using regression model
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")

# Plot Predicted data and original data points
#plot((resection_tumor_grade__hg) ~ A_IDH_HG__A_IDH_LG_lr__lasso_fit, data=stats)
#lines((resection_tumor_grade__hg) ~ A_IDH_HG__A_IDH_LG_lr__lasso_fit, Predicted_data, lwd=2, col="darkgreen")


plt.logit <- rbind(
  stats |> 
    dplyr::select( resection_tumor_grade__hg, A_IDH_HG__A_IDH_LG_lr__lasso_fit, mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data"),
  Predicted_data |> 
    dplyr::mutate(mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2)


p1 <- ggplot(plt.logit, aes(x = A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=`resection tumor grade`, col= mnp_predictBrain_v12.8_cal_class)) +
  geom_point(data = plt.logit |> dplyr::filter(type == "data"),pch="|",size=5) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), col="darkgreen") +
  theme_bw() + 
  annotate("text", x = 7.5, y = 2.5, label = paste0("p = ",format.pval(pval))) + 
  theme(legend.position = 'bottom') + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col="mnp v12.8")




### logistic PC2 ----


model <- glm(resection_tumor_grade__hg ~ PC2inv, data = stats,family = binomial)
summary(model)
#sjPlot::plot_model(model)

pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "PC2inv") |> 
  dplyr::pull(`Pr(>|z|)`)


#Data frame with hp in ascending order
Predicted_data <- data.frame(PC2inv=seq(min(stats$PC2inv), max(stats$PC2inv), len=500))

# Fill predicted values using regression model
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")

# Plot Predicted data and original data points
#plot((resection_tumor_grade__hg) ~ PC2inv, data=stats)
#lines((resection_tumor_grade__hg) ~ PC2inv, Predicted_data, lwd=2, col="darkgreen")

plt.logit <- rbind(
  stats |> 
    dplyr::select( resection_tumor_grade__hg, PC2inv, mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data"),
  Predicted_data |> 
    dplyr::mutate(mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2)


p2 <- ggplot(plt.logit, aes(x = PC2inv, y=`resection tumor grade`, col= mnp_predictBrain_v12.8_cal_class)) +
  geom_point(data = plt.logit |> dplyr::filter(type == "data"),pch="|",size=5) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), col="darkgreen") +
  theme_bw() + 
  annotate("text", x = 750, y = 2.5, label = paste0("p = ",format.pval(pval))) + 
  theme(legend.position = 'bottom') + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col="mnp v12.8")



### logistic median.overall.methylation ----


model <- glm(resection_tumor_grade__hg ~  median.overall.methylation, data = stats,family = binomial)
summary(model)
#sjPlot::plot_model(model)

pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "median.overall.methylation") |> 
  dplyr::pull(`Pr(>|z|)`)

Predicted_data <- data.frame(median.overall.methylation=seq(min(stats$median.overall.methylation), max(stats$median.overall.methylation), len=500)) |> 
  dplyr::mutate(batch_us = 0.5)

Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  stats |> 
    dplyr::select(resection_tumor_grade__hg, median.overall.methylation, mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(batch_us = 0.5)
  ,
  Predicted_data |> 
    dplyr::mutate(mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2)


p3 <- ggplot(plt.logit, aes(x = median.overall.methylation, y=`resection tumor grade`, col= mnp_predictBrain_v12.8_cal_class)) +
  geom_point(data = plt.logit |> dplyr::filter(type == "data"),pch="|",size=5) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), col="darkgreen") +
  theme_bw() + 
  annotate("text", x = 2.5, y = 2.5, label = paste0("p = ",format.pval(pval))) + 
  theme(legend.position = 'bottom') + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col="mnp v12.8")



p1 + p2 + p3


### logistic median.glass_nl_supervised.methylation ----


model <- glm(resection_tumor_grade__hg ~  median.glass_nl_supervised.methylation, data = stats,family = binomial)
summary(model)
#sjPlot::plot_model(model)

pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "median.glass_nl_supervised.methylation") |> 
  dplyr::pull(`Pr(>|z|)`)

Predicted_data <- data.frame(median.glass_nl_supervised.methylation=seq(min(stats$median.glass_nl_supervised.methylation), max(stats$median.glass_nl_supervised.methylation), len=500)) |> 
  dplyr::mutate(batch_us = 0.5)

Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  stats |> 
    dplyr::select(resection_tumor_grade__hg, median.glass_nl_supervised.methylation, mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(batch_us = 0.5)
  ,
  Predicted_data |> 
    dplyr::mutate(mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2)


p4 <- ggplot(plt.logit, aes(x = median.glass_nl_supervised.methylation, y=`resection tumor grade`, col= mnp_predictBrain_v12.8_cal_class)) +
  geom_point(data = plt.logit |> dplyr::filter(type == "data"),pch="|",size=5) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), col="darkgreen") +
  theme_bw() + 
  annotate("text", x = 2.5, y = 2.5, label = paste0("p = ",format.pval(pval))) + 
  theme(legend.position = 'bottom') + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col="mnp v12.8")



p1 + p2 + p3 + p4


# confirm - methylation lowers
# ggplot(stats,aes(
#   x = resection_tumor_grade,
#   y = median.glass_nl_supervised.methylation
# )) +
#   ggbeeswarm::geom_quasirandom()
# 


### multi LM ----




model <- lm(resection_tumor_grade__hg ~ PC2 + A_IDH_HG__A_IDH_LG_lr__lasso_fit, data = stats)
model.g <- glm(resection_tumor_grade__hg ~ PC2 + A_IDH_HG__A_IDH_LG_lr__lasso_fit, data = stats,family = binomial)
summary(model)
summary(model.g)
sjPlot::plot_model(model)


model2 <- lm(A_IDH_HG__A_IDH_LG_lr__lasso_fit ~ resection_tumor_grade__hg + batch_us, data = stats)
summary(model2)
sjPlot::plot_model(model2)

model3 <- lm(-PC2 ~ resection_tumor_grade__hg + batch_us, data = stats)
summary(model3)
sjPlot::plot_model(model3)


model3 <- lm(PC2 ~ resection_tumor_grade__hg, data = stats)
summary(model3)



# G-SAM ----



metadata <- gsam.metadata.idats |> 
  filter_GSAM_idats(73) |> 
  dplyr::mutate(resection__rec = ifelse(resection == "R1",0,1))


data <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata$sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (694299))
    return(.)
  })() |> 
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |> # this synthax, oh my
  dplyr::arrange(mad) |> 
  dplyr::mutate(mad = NULL)


## logistic LR & recurrence ----



ggplot(metadata, aes(x=resection, y=A_IDH_HG__A_IDH_LG_lr, group=patient)) +
  ggbeeswarm::geom_quasirandom() +
  geom_line(col="black", alpha=0.25, lwd=0.25) +
  theme_bw()


model <- glm(resection__rec ~  A_IDH_HG__A_IDH_LG_lr, 
             data = metadata,
             family = binomial)
summary(model)


pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "A_IDH_HG__A_IDH_LG_lr") |> 
  dplyr::pull(`Pr(>|z|)`)


Predicted_data <- data.frame(A_IDH_HG__A_IDH_LG_lr=modelr::seq_range(metadata$A_IDH_HG__A_IDH_LG_lr, 500))
Predicted_data$resection__rec = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  metadata |> 
    dplyr::select(resection__rec, A_IDH_HG__A_IDH_LG_lr, mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data")
  ,
  Predicted_data |> 
    dplyr::mutate(mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
)


ggplot(plt.logit, aes(x = A_IDH_HG__A_IDH_LG_lr, y=`resection__rec`, col= mnp_predictBrain_v12.8_cal_class)) +
  geom_point(data = plt.logit |> dplyr::filter(type == "data"),pch="|",size=5) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), col="darkgreen") +
  theme_bw() + 
  annotate("text", x = 2.5, y = 0.5, label = paste0("p = ",format.pval(pval))) + 
  theme(legend.position = 'bottom') + 
  scale_y_continuous(breaks = c(0,1)) + 
  labs(col="mnp v12.8", y = "Resection (0 primary, 1 recurrent)")



## survival ----
### R1 ----


ggplot(
  metadata |> 
    dplyr::filter(resection == "R1"),
  aes(y= survivalDays , x=A_IDH_HG__A_IDH_LG_lr, shape=status)) + 
  geom_point()



### R2 ----

#plot(sort(stats$A_IDH_HG__A_IDH_LG_lr__lasso_fit)) + abline(h=4.65)


stats <-  metadata |> 
  dplyr::filter(resection == "R2") |> 
  dplyr::mutate(A_IDH_HG__A_IDH_LG__cut = cut(A_IDH_HG__A_IDH_LG_lr__lasso_fit, breaks=2)) |> 
  dplyr::mutate(LR_status_median = ifelse(A_IDH_HG__A_IDH_LG_lr__lasso_fit > 6.1, "LR high", "LR low")) |> 
  dplyr::mutate(LR_status_L1 = ifelse(A_IDH_HG__A_IDH_LG_lr__lasso_fit > 4.7, "LR high", "LR low")) |> 
  dplyr::mutate(LR_status_L2 = ifelse(A_IDH_HG__A_IDH_LG_lr__lasso_fit > 4.1, "LR high", "LR low"))


ggplot(stats, aes(y= survivalFromSecondSurgeryDays , x=A_IDH_HG__A_IDH_LG_lr, shape=status)) + 
  geom_point()


surv_object <- survival::Surv(time = stats$survivalFromSecondSurgeryDays, event=stats$os.event)
#surv_object <- survival::Surv(time = stats$survivalDays, event=stats$os.event)


fit1 <- survival::survfit(surv_object ~  A_IDH_HG__A_IDH_LG__cut , data = stats)
print(survminer::surv_pvalue(fit1)$pval)

survminer::ggsurvplot(fit1, data = stats, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                      xlab="Survival from recurrence")


fit2 <- survival::survfit(surv_object ~  LR_status_median , data = stats)
print(survminer::surv_pvalue(fit2)$pval)

survminer::ggsurvplot(fit2, data = stats, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                      xlab="Survival from recurrence")

fit3 <- survival::survfit(surv_object ~  LR_status_L1 , data = stats)
print(survminer::surv_pvalue(fit3)$pval)

survminer::ggsurvplot(fit3, data = stats, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                      xlab="Survival from recurrence")

fit4 <- survival::survfit(surv_object ~  LR_status_L2 , data = stats)
print(survminer::surv_pvalue(fit4)$pval)

survminer::ggsurvplot(fit4, data = stats, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                      xlab="Survival from recurrence")




fit_l <- survival::coxph(surv_object ~  A_IDH_HG__A_IDH_LG_lr__lasso_fit , data = stats)
summary(fit_l)




# GLASS-NL / AC ----


metadata <- glass_nl.metadata.idats |> 
  dplyr::filter(!qc.pca.detP.outlier) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (218))
    return(.)
  })()


## PCA ----


masked <- readRDS('cache/mvalues.HQ_samples.detP_mask.Rds') |> 
  is.na() |> 
  rowSums() |> 
  as.data.frame() |> 
  dplyr::rename(count = 1) |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(count > 0) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (66116))
    return(.)
  })() |> 
  dplyr::pull('probe_id')




data <- readRDS('cache/mvalues.HQ_samples.Rds') 

data <- data |> 
  dplyr::select(metadata$sentrix_id)

data <- data |> 
  tibble::rownames_to_column('probe_id')

data <- data |> 
  dplyr::filter(probe_id %in% (masked) == F) |> 
  tibble::column_to_rownames('probe_id')

dim(data)


data <- data |> 
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |> # this synthax, oh my
  dplyr::arrange(mad) |> 
  dplyr::mutate(mad = NULL)


data.all.pca.obj <- data |> 
  #dplyr::slice_head(n=100000) |> 
  t() |> 
  prcomp()

data.all.pca <- data.all.pca.obj |> 
  purrr::pluck('x') |> 
  as.data.frame(stringsAsFactors=F) |> 
  dplyr::select(paste0("PC",1:40)) |> 
  tibble::rownames_to_column('sentrix_id')


plt <- metadata |> 
  dplyr::left_join(data.all.pca, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))


plt <- plt |> 
  dplyr::mutate(mnp_predictBrain_v12.8_cal_A_IDH_LG_odds = (mnp_predictBrain_v12.8_cal_A_IDH_LG/100) / (1 - (mnp_predictBrain_v12.8_cal_A_IDH_LG/100))) |> 
  dplyr::mutate(mnp_predictBrain_v12.8_cal_A_IDH_HG_odds = (mnp_predictBrain_v12.8_cal_A_IDH_HG/100) / (1 - (mnp_predictBrain_v12.8_cal_A_IDH_HG/100))) |> 
  dplyr::mutate(A_IDH_HG__A_IDH_LG_lr = log(mnp_predictBrain_v12.8_cal_A_IDH_HG_odds / mnp_predictBrain_v12.8_cal_A_IDH_LG_odds))


ggplot(plt, aes(x=A_IDH_HG__A_IDH_LG_lr, y=PC1)) +
  geom_point()




# detP mask filter
# order on MAD
# pick top 200.000
# prcomp


# plot(plt$PC2, plt$median.overall.methylation)
# plot(plt$PC2, plt$A_IDH_HG__A_IDH_LG_lr)
# plot(plt$PC2, plt$A_IDH_HG__A_IDH_LG_lr__lasso_fit)
# plot(plt$median.overall.methylation, plt$A_IDH_HG__A_IDH_LG_lr__lasso_fit)
# plt$A_IDH_HG__O_IDH_lr, plt$A_IDH_HG__A_IDH_LG_lr)


# c <- plt |>
#   dplyr::mutate(mnp_predictBrain_v12.8_cal_A_IDH_HG = -log(mnp_predictBrain_v12.8_cal_A_IDH_HG)) |> 
#   dplyr::select(
#   PC2,
#   A_IDH_HG__A_IDH_LG_lr,
#   A_IDH_HG__A_IDH_LG_lr__lasso_fit,
#   A_IDH_HGoligsarc__O_IDH_lr,
#   median.overall.methylation,
#   median.glass_nl_supervised.methylation,
#   A_IDH_HG__O_IDH_lr
#   #mnp_predictBrain_v12.8_cal_A_IDH_HG
# ) |> 
#   as.matrix()
# 
# 
corrplot::corrplot(abs(cor(c, method="spearman")), order="hclust")
corrplot::corrplot(abs(cor(c, method="kendall")), order="hclust")
# 
# 
# plot(plt$A_IDH_HGoligsarc__O_IDH_lr, plt$PC2)
# plot(plt$A_IDH_HG__O_IDH_lr, plt$PC2)
# 
# plot(plt$A_IDH_HGoligsarc__O_IDH_lr, plt$PC2)
# plot(plt$PC2)
# 
