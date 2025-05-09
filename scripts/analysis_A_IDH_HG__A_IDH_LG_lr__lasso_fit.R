#!/usr/bin/env R

# load data ----


library(ggplot2)


## m-values ----


if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}


## metadata ----


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}

if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}

if(!exists('gsam.metadata.array_samples')) {
  source('scripts/load_G-SAM_metadata.R')
}




metadata.glass_od <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES)

data.glass_od <- data.mvalues.hq_samples |>
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.glass_od$array_sentrix_id)  |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_GLASS_OD_INCLUDED_SAMPLES)
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()




metadata.tcga_lgg <- tcga_lgg_od.metadata.array_samples
data.tcga_lgg <- readRDS(file=paste0("cache/mvalues.tcga-lgg.Rds")) |> 
  tibble::column_to_rownames('probe_id')




metadata.od_validation <- glass_od.metadata.array_samples |> 
  filter_OD_validation_idats(CONST_N_OD_VALIDATION_INCLUDED_SAMPLES)

data.od_validation <- data.mvalues.hq_samples |>
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.od_validation$array_sentrix_id)  |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_OD_VALIDATION_INCLUDED_SAMPLES) 
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()





metadata.glass_nl <- glass_nl.metadata.array_samples |> 
  filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES) |> 
  dplyr::mutate(i = 1:dplyr::n()) |> 
  dplyr::mutate(slice = i %% 10, i = NULL)


data.glass_nl <- data.mvalues.hq_samples |>
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.glass_nl$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_GLASS_NL_INCLUDED_SAMPLES) 
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()






metadata.gsam <- gsam.metadata.array_samples |> 
  filter_GSAM_idats(CONST_N_GSAM_INCLUDED_SAMPLES)

data.gsam <- data.mvalues.hq_samples |>
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.gsam$array_sentrix_id)  |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_GSAM_INCLUDED_SAMPLES)
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP)
    return(.)
  })()




metadata.catnon <- glass_od.metadata.array_samples |> 
  dplyr::filter(patient_study_name == "CATNON")

data.catnon <- data.mvalues.hq_samples |>
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.catnon$array_sentrix_id)  |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_CATNON_ALL_SAMPLES)
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()






## add PCA to metadata ----


data.glass_nl.pca.obj <- data.glass_nl |> 
  t() |> 
  prcomp()
data.glass_nl.pca <- data.glass_nl.pca.obj |> 
  purrr::pluck('x') |> 
  t() |> 
  as.data.frame(stringsAsFactors=F)



# model 1: purely on m-values ----
#' https://www.statology.org/lasso-regression-in-r/
#' PCA based works a little better and considerably faster

## train parameters ----


set.seed(123456)
cv_model_probe_based <- glmnet::cv.glmnet(
  data.glass_nl |>
    dplyr::select(metadata.glass_nl$array_sentrix_id) |> 
    t(),
  metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, alpha = 1, relax=F)

saveRDS(cv_model_probe_based, file="cache/analysis_A_IDH_HG__A_IDH_LG_lr__lasso_fit__probe_based__train_paramters.Rds")

# @todo export & ggplot
pdf("output/figures/analysis_A_IDH_HG__A_IDH_LG_l__prediction_err.pdf", width = (8.5/3) * 0.95, height = 2.80)
plot(cv_model_probe_based)
dev.off()




## test predicted vs real ----

cv_model_probe_based <- readRDS(file="cache/analysis_A_IDH_HG__A_IDH_LG_lr__lasso_fit__probe_based__train_paramters.Rds")

out.probe_based <- data.frame()
for(cv in 0:9) {
  
  metadata.train <- metadata.glass_nl |> dplyr::filter(slice != cv)
  metadata.test <- metadata.glass_nl |> dplyr::filter(slice == cv)
  
  set.seed(123456)
  lm.probe_based <- glmnet::glmnet(
    data.glass_nl |>
      dplyr::select(metadata.train$array_sentrix_id) |> 
      t(),
    metadata.train$array_A_IDH_HG__A_IDH_LG_lr_v12.8,
    alpha = 1, # 1 = lasso; 0 = ridge
    lambda = cv_model_probe_based$lambda.min # 0.01
    )
  
  
  metadata.test$LGC.lasso.probe_based <- predict(lm.probe_based, data.glass_nl |>
                                                   dplyr::select(metadata.test$array_sentrix_id) |> 
                                                   t())|> 
    as.numeric()
  
  out.probe_based <- rbind(
    metadata.test |>  
      dplyr::select(array_sentrix_id, LGC.lasso.probe_based, array_A_IDH_HG__A_IDH_LG_lr_v12.8, array_mnp_predictBrain_v12.8_cal_class)
    , out.probe_based)
}
rm(cv)


#sqrt(sum((out.probe_based$LGC.lasso.probe_based - out.probe_based$A_IDH_HG__A_IDH_LG_lr)^2) / length(out.probe_based$A_IDH_HG__A_IDH_LG_lr))



### export 10xCV data points ----


tmp <- out.probe_based |> 
  dplyr::rename(array_A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV = LGC.lasso.probe_based) |> 
  dplyr::select(array_sentrix_id, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV)

saveRDS(tmp, file="cache/analysis_A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV.Rds")



## build model on full data ----


lm.full.probe_based <- glmnet::glmnet(data.glass_nl |> 
                                      dplyr::select(metadata.glass_nl$array_sentrix_id) |>
                                      t(),
                                    metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8,
                                    alpha = 1,
                                    lambda = cv_model_probe_based$lambda.min)




### export final AC-trained classifier ----


saveRDS(rownames(data.glass_nl), file="cache/LGC_predictor_probe_based_ordered_probes.Rds")
saveRDS(lm.full.probe_based, file="cache/LGC_predictor_probe_based_lm.Rds")




# probes <- names(lm.full.probe_based$beta[which(lm.full.probe_based$beta != 0),1])


## build model on full data [450k probes] ----


tmp <- data.glass_nl |> 
  dplyr::select(metadata.glass_nl$array_sentrix_id) |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% rownames(data.tcga_lgg)) |> 
  tibble::column_to_rownames('probe_id')


cv_model_probe_based_450k <- glmnet::cv.glmnet(
  tmp |> t(),
  metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, alpha = 1, relax=F)


lm.full.probe_based.450k <- glmnet::glmnet(tmp |> t(),
                                           metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8,
                                           alpha = 1,
                                           lambda = cv_model_probe_based_450k$lambda.min)



saveRDS(rownames(tmp), file="cache/LGC_predictor_probe_based_ordered_probes_450k.Rds")
saveRDS(lm.full.probe_based.450k, file="cache/LGC_predictor_probe_based_lm_450k.Rds")


rm(tmp)




# model 2: on PCA ----

## train parameters ----

rm(cv_model.pca_based, lm)
metadata.test$LGC.lasso.pca_based <- NULL

set.seed(123456)
cv_model.pca_based <- glmnet::cv.glmnet(
  data.glass_nl.pca |>  dplyr::select(metadata.glass_nl$sentrix_id) |> 
    t(),
                              metadata.glass_nl$A_IDH_HG__A_IDH_LG_lr, alpha = 1, relax=F)
plot(cv_model.pca_based) 



## simple cv ----


rm(out.pca_based)
out.pca_based <- data.frame()
for(cv in 0:9) {
  
  metadata.train <- metadata.glass_nl |> dplyr::filter(slice != cv)
  metadata.test <- metadata.glass_nl |> dplyr::filter(slice == cv)
  
  set.seed(123456)
  lm.pca_based <- glmnet::glmnet(
    data.glass_nl.pca |>
      dplyr::select(metadata.train$sentrix_id) |> 
      t(),
              metadata.train$A_IDH_HG__A_IDH_LG_lr, 
    alpha = 1,
    lambda = cv_model.pca_based$lambda.min)
  
  
  metadata.test$LGC.lasso.pca_based <- predict(lm.pca_based, 
                                               data.glass_nl.pca |> 
                                                 dplyr::select(metadata.test$sentrix_id) |> t()
                                               ) |> 
    as.numeric()
  sqrt(sum((metadata.test$LGC.lasso.pca_based - metadata.test$A_IDH_HG__A_IDH_LG_lr)^2) / length(metadata.test$A_IDH_HG__A_IDH_LG_lr))
  
  
  out.pca_based <- rbind(
    metadata.test |>  dplyr::select(sentrix_id, LGC.lasso.pca_based, A_IDH_HG__A_IDH_LG_lr, mnp_predictBrain_v12.8_cal_class)
    , out.pca_based
  )
}



out.pca_based <- out.pca_based |> 
  dplyr::mutate(mnp_predictBrain_v12.8_cal_class = ifelse(mnp_predictBrain_v12.8_cal_class %in% c("A_IDH","A_IDH_LG","A_IDH_HG") == F, "other", mnp_predictBrain_v12.8_cal_class))


ggplot(out.pca_based, aes(
    y=LGC.lasso.pca_based,
    x=A_IDH_HG__A_IDH_LG_lr,
    col=mnp_predictBrain_v12.8_cal_class
    #label=sentrix_id
  )) +
    geom_point(size=1.5) +
    ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1") +
    labs(subtitle = "Linear PCA-based classifier predicting log(IDH_HG / IDH[_LG]) ratio [10x CV]")+
    theme_bw() +
    theme(legend.position = 'bottom')
  
  



## combined ----


plt <- out.probe_based |> 
  dplyr::left_join(out.pca_based |> dplyr::select(sentrix_id, `LGC.lasso.pca_based`),
                   by=c('sentrix_id'='sentrix_id'))


ggplot(plt, aes(x=LGC.lasso.probe_based, y=LGC.lasso.pca_based )) +
  geom_point()


## build model on full data ----


lm.full.pca_based <- glmnet::glmnet(data.glass_nl.pca |> 
                            dplyr::select(metadata.glass_nl$sentrix_id) |>
                            t(),
                     metadata.glass_nl$A_IDH_HG__A_IDH_LG_lr,
                     alpha = 1, lambda = cv_model.pca_based$lambda.min)

saveRDS(data.glass_nl.pca.obj, file="cache/LGC_predictor_PCA_based_prcomp.Rds")
saveRDS(lm.full.pca_based, file="cache/LGC_predictor_PCA_based_lm.Rds")




# apply to GLASS-OD & G-SAM ----


classifier <- readRDS("cache/LGC_predictor_probe_based_lm.Rds")


data.1 <- data.glass_od |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% rownames(classifier$beta)) |> 
  tibble::column_to_rownames('probe_id') |> 
  t() |>
  as.data.frame() |> 
  dplyr::select(readRDS("cache/LGC_predictor_probe_based_ordered_probes.Rds")) |> 
  as.matrix()


data.2 <- data.od_validation  |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% rownames(classifier$beta)) |> 
  tibble::column_to_rownames('probe_id') |> 
  t() |>
  as.data.frame() |> 
  dplyr::select(readRDS("cache/LGC_predictor_probe_based_ordered_probes.Rds")) |> 
  as.matrix()


data.3 <- data.gsam |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% rownames(classifier$beta)) |> 
  tibble::column_to_rownames('probe_id') |> 
  t() |>
  as.data.frame() |> 
  dplyr::select(readRDS("cache/LGC_predictor_probe_based_ordered_probes.Rds")) |> 
  as.matrix()



data.4 <- data.catnon |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% rownames(classifier$beta)) |> 
  tibble::column_to_rownames('probe_id') |> 
  t() |>
  as.data.frame() |> 
  dplyr::select(readRDS("cache/LGC_predictor_probe_based_ordered_probes.Rds")) |> 
  as.matrix()


stopifnot(colnames(data.1) == colnames(data.2))
stopifnot(colnames(data.1) == colnames(data.3))
stopifnot(colnames(data.1) == colnames(data.4))
data <- rbind(data.1, data.2, data.3, data.4)



stopifnot(rownames(classifier$beta) == colnames(data))
stopifnot(readRDS("cache/LGC_predictor_probe_based_ordered_probes.Rds") == colnames(data))


p <- predict(classifier, data) |> 
  as.data.frame() |> 
  dplyr::rename(`array_A_IDH_HG__A_IDH_LG_lr__lasso_fit` = 1) |> 
  tibble::rownames_to_column('array_sentrix_id')


saveRDS(p, file="cache/analysis_A_IDH_HG__A_IDH_LG_lr__lasso_fit.Rds")



rm(p, data, data.1, data.2, data.3, data.4)



# apply 450k model to TCGA-LGG oligo samples ----


classifier <- readRDS("cache/LGC_predictor_probe_based_lm_450k.Rds")

data <- data.tcga_lgg |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% rownames(classifier$beta)) |> 
  tibble::column_to_rownames('probe_id') |> 
  t() |>
  as.data.frame() |> 
  dplyr::select(readRDS("cache/LGC_predictor_probe_based_ordered_probes_450k.Rds")) |> 
  as.matrix()


p <- predict(classifier, data) |> 
  as.data.frame() |> 
  dplyr::rename(`array_A_IDH_HG__A_IDH_LG_lr__lasso_fit_450k` = 1) |> 
  tibble::rownames_to_column('array_sentrix_id')


saveRDS(p, file="cache/analysis_A_IDH_HG__A_IDH_LG_lr__lasso_fit__TCGA-LGG_450k.Rds")


rm(classifier, data, p)
rm(cv_model_probe_based_450k)
rm(lm.full.probe_based.450k)
rm(metadata.tcga_lgg)
rm(data.tcga_lgg)
