#!/usr/bin/env R

# make linear model to predict A_IDH ~ A_IDH_HG

if(!exists('glass_nl.data.mvalues')) {
  source('scripts/load_mvalues.R')
}



metadata.all <- glass_nl.metadata.resections |> 
  dplyr::filter(!is.na(IDH_HG_IDH_ratio)) |> 
  dplyr::mutate(i = 1:n()) |> 
  dplyr::mutate(m = i %% 10)


data.all <- glass_nl.data.mvalues |>
  dplyr::select(metadata.all$methylation.sid) |>
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |> # this synthax, oh my
  dplyr::arrange(mad) |> 
  dplyr::mutate(mad = NULL)

data.all.pca.obj <- data.all |> 
  t() |> 
  prcomp()
data.all.pca <- data.all.pca.obj |> 
  purrr::pluck('x') |> 
  t() |> 
  as.data.frame(stringsAsFactors=F)



# https://www.statology.org/lasso-regression-in-r/




# purely on m-values ----
#' PCA based works a little better and considerably faster



set.seed(123456)
cv_model_probe_based <- glmnet::cv.glmnet(t(data.all |>  dplyr::select(metadata.all$methylation.sid)),
                                          metadata.all$IDH_HG_IDH_ratio, alpha = 1, relax=F)
plot(cv_model_probe_based) 

## simple cv ----


out.probe_based <- data.frame()
for(cv in 0:9) {
  
  metadata.train <- metadata.all |> dplyr::filter(m != cv)
  metadata.test <- metadata.all |> dplyr::filter(m == cv)
  
  set.seed(123456)
  lm.probe_based <- glmnet::glmnet(
    data.all |>
      dplyr::select(metadata.train$methylation.sid) |> 
      t(),
    metadata.train$IDH_HG_IDH_ratio,
    alpha = 1, # 1 = lasso; 0 = ridge
    lambda = cv_model_probe_based$lambda.min)
  
  
  metadata.test$LGC.lasso.probe_based <- predict(lm.probe_based, data.all |>
                                                   dplyr::select(metadata.test$methylation.sid) |> 
                                                   t())|> 
    as.numeric()
  
  out.probe_based <- rbind(
    metadata.test |>  
      dplyr::select(methylation.sid, LGC.lasso.probe_based, IDH_HG_IDH_ratio, methylation.sub.diagnosis)
    , out.probe_based)
}
rm(cv)

sqrt(sum((out.probe_based$LGC.lasso.probe_based - out.probe_based$IDH_HG_IDH_ratio)^2) / length(out.probe_based$IDH_HG_IDH_ratio))


out.probe_based <- out.probe_based |> 
  dplyr::mutate(methylation.sub.diagnosis = ifelse(methylation.sub.diagnosis %in% c("A_IDH","A_IDH_LG","A_IDH_HG") == F, "other", methylation.sub.diagnosis))


ggplot(out.probe_based, aes(
  y=LGC.lasso.probe_based, 
  x=IDH_HG_IDH_ratio, 
  col=methylation.sub.diagnosis
)) +
  geom_point() +
  labs(subtitle = "Linear probe-based classifier predicting log(IDH_HG / IDH[_LG]) ratio [10x CV]")


rm(out.probe_based)




## export ----




lm.full.probe_based <- glmnet::glmnet(data.all |> 
                            dplyr::select(metadata.all$methylation.sid) |>
                            t(),
                          metadata.all$IDH_HG_IDH_ratio,
                          alpha = 1,
                          lambda = cv_model_probe_based$lambda.min)

#saveRDS(rownames(data.all), file="cache/LGC_predictor_probe_based_ordered_probes.Rds")
#saveRDS(lm.full.probe_based, file="cache/LGC_predictor_probe_based_lm.Rds")
#lm.full.probe_based <- readRDS("cache/LGC_predictor_probe_based_lm.Rds")

# probes <- names(lm.full.probe_based$beta[which(lm.full.probe_based$beta != 0),1])


# on PCA ----

rm(cv_model.pca_based, lm)
metadata.test$LGC.lasso.pca_based <- NULL

set.seed(123456)
cv_model.pca_based <- glmnet::cv.glmnet(
  data.all.pca |>  dplyr::select(metadata.all$methylation.sid) |> 
    t(),
                              metadata.all$IDH_HG_IDH_ratio, alpha = 1, relax=F)
plot(cv_model.pca_based) 



## simple cv ----


rm(out.pca_based)
out.pca_based <- data.frame()
for(cv in 0:9) {
  
  metadata.train <- metadata.all |> dplyr::filter(m != cv)
  metadata.test <- metadata.all |> dplyr::filter(m == cv)
  
  set.seed(123456)
  lm.pca_based <- glmnet::glmnet(
    data.all.pca |>
      dplyr::select(metadata.train$methylation.sid) |> 
      t(),
              metadata.train$IDH_HG_IDH_ratio, 
    alpha = 1,
    lambda = cv_model.pca_based$lambda.min)
  
  
  metadata.test$LGC.lasso.pca_based <- predict(lm.pca_based, 
                                               data.all.pca |> 
                                                 dplyr::select(metadata.test$methylation.sid) |> t()
                                               ) |> 
    as.numeric()
  sqrt(sum((metadata.test$LGC.lasso.pca_based - metadata.test$IDH_HG_IDH_ratio)^2) / length(metadata.test$IDH_HG_IDH_ratio))
  
  
  out.pca_based <- rbind(
    metadata.test |>  dplyr::select(methylation.sid, LGC.lasso.pca_based, IDH_HG_IDH_ratio, methylation.sub.diagnosis)
    , out.pca_based
  )
}



out.pca_based <- out.pca_based |> 
  dplyr::mutate(methylation.sub.diagnosis = ifelse(methylation.sub.diagnosis %in% c("A_IDH","A_IDH_LG","A_IDH_HG") == F, "other", methylation.sub.diagnosis))


ggplot(out.pca_based, aes(
  y=LGC.lasso.pca_based, 
  x=IDH_HG_IDH_ratio, 
  col=methylation.sub.diagnosis
)) +
  geom_point() +
  labs(subtitle = "Linear PCA-based classifier predicting log(IDH_HG / IDH[_LG]) ratio [10x CV]")



## combined ----


plt <- out.probe_based |> 
  dplyr::left_join(out.pca_based |> dplyr::select(methylation.sid, `LGC.lasso.pca_based`),
                   by=c('methylation.sid'='methylation.sid'))


ggplot(plt, aes(x=LGC.lasso.probe_based, y=LGC.lasso.pca_based )) +
  geom_point()


## export to reproduce ----


lm.full.pca_based <- glmnet::glmnet(data.all.pca |> 
                            dplyr::select(metadata.all$methylation.sid) |>
                            t(),
                     metadata.all$IDH_HG_IDH_ratio,
                     alpha = 1, lambda = cv_model$lambda.min)

saveRDS(data.all.pca.obj, file="cache/LGC_predictor_PCA_based_prcomp.Rds")
saveRDS(lm.full.pca_based, file="cache/LGC_predictor_PCA_based_lm.Rds")


# apply to OD samples ----

classifier <- readRDS("cache/LGC_predictor_probe_based_lm.Rds")
data <- glass_od.data.mvalues |> 
  tibble::rownames_to_column('probeID') |> 
  dplyr::filter(probeID %in% rownames(classifier$beta)) |> 
  tibble::column_to_rownames('probeID') |> 
  t() |>
  as.data.frame() |> 
  dplyr::select(readRDS("cache/LGC_predictor_probe_based_ordered_probes.Rds")) |> 
  as.matrix()


stopifnot(rownames(classifier$beta) == colnames(data))
stopifnot(readRDS("cache/LGC_predictor_probe_based_ordered_probes.Rds") == colnames(data))

p = predict(classifier, data) |> 
  as.data.frame() |> 
  dplyr::rename(`A_IDH_A_IDH_HG_lasso_fit` = 1) |> 
  tibble::rownames_to_column('sentrix_id')

plt <- glass_od.metadata.idats |> 
  dplyr::left_join(p, by=c('sentrix_id'='sentrix_id'),suffix=c('','')) |> 
  dplyr::left_join(plt.pca, by=c('sentrix_id'='sentrix_id'),suffix=c('','')) |> 
  
  dplyr::mutate(label = paste0(as.character(plt$resection_id),""))



ggplot(plt, aes(x=A_IDH_A_IDH_HG_lasso_fit, y=resection_tumor_grade)) + 
  geom_point()


#ggplot(plt, aes(x=A_IDH_A_IDH_HG_lasso_fit, y=PC3, col=as.factor(predictBrain_12.5_cal_class), label=label)) + 
ggplot(plt, aes(x=A_IDH_A_IDH_HG_lasso_fit, y=PC3, col=qc.pca.outlier, label=label)) + 
  geom_point() +
  labs(y = "Generic PCA / PC3 in all Oligo samples")
  #ggrepel::geom_text_repel(col="blue", size=2.5 , segment.size=0.35, show.legend  = FALSE) 


ggplot(plt |>  dplyr::filter(qc.pca.outlier == F), aes(x=A_IDH_A_IDH_HG_lasso_fit, y=PC3, col=as.factor(predictBrain_12.5_cal_class), label=label)) + 
  ggpubr::stat_cor(method = "spearman", aes(label = ..r.label..), col="1") +
  geom_point() +
  labs(y = "Generic PCA / PC3 in all Oligo samples") + 
  theme_bw()


ggplot(plt |>  dplyr::filter(qc.pca.outlier == F), aes(x=A_IDH_A_IDH_HG_lasso_fit, y=PC3, col=as.factor(resection_tumor_grade))) + 
  ggpubr::stat_cor(method = "spearman", aes(label = ..r.label..), col="1") +
  geom_point() +
  labs(y = "Generic PCA / PC3 in all Oligo samples") + 
  theme_bw()


ggplot(plt, aes(x=A_IDH_A_IDH_HG_lasso_fit, y=PC4)) + 
  geom_point()


