#!/usr/bin/env R


# load data ----


library(ggplot2)


## m-values ----


if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}


## metadata ----


if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}

# if(!exists('glass_od.metadata.array_samples')) {
#   source('scripts/load_GLASS-OD_metadata.R')
# }
# 
# if(!exists('gsam.metadata.array_samples')) {
#   source('scripts/load_G-SAM_metadata.R')
# }




metadata.glass_nl <- glass_nl.metadata.array_samples |> 
  filter_GLASS_NL_idats(218) |> 
  dplyr::filter(!is.na(rnaseq_cell.cycling.signature)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 156) # + 10x astro
    return(.)
  })() |> 
  dplyr::mutate(i = 1:dplyr::n()) |> 
  dplyr::mutate(slice = i %% 10, i = NULL)
  #dplyr::select(!starts_with("array_PC.GLASS_OD_NL") & !starts_with("array_PC") & !starts_with("array_mnp_predictBrain_v12.8_cal_"))


data.glass_nl <- data.mvalues.hq_samples |>
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.glass_nl$array_sentrix_id)




## EM iteration 1 ----


metadata.glass_nl <- metadata.glass_nl |> 
  dplyr::select(!starts_with("rnaseq_cell.cycling.signature_it"))



### params ----


set.seed(123456)
cv_model_probe_based <- glmnet::cv.glmnet(
  data.glass_nl |>
    dplyr::select(metadata.glass_nl$array_sentrix_id) |> 
    t(),
  metadata.glass_nl$rnaseq_cell.cycling.signature, alpha = 1, relax=F)


# @todo export & ggplot
plot(cv_model_probe_based)



### prediction ----


tmp <- data.frame()

for(cv in 0:9) {
  
  metadata.train <- metadata.glass_nl |> dplyr::filter(slice != cv)
  metadata.test <- metadata.glass_nl |> dplyr::filter(slice == cv)
  
  set.seed(123456)
  lm.probe_based <- glmnet::glmnet(
    data.glass_nl |>
      dplyr::select(metadata.train$array_sentrix_id) |> 
      t(),
    metadata.train$rnaseq_cell.cycling.signature,
    alpha = 1, # 1 = lasso; 0 = ridge
    lambda = 0.01 # cv_model_probe_based$lambda.min
  )
  
  
  metadata.test$rnaseq_cell.cycling.signature_it1 <- predict(lm.probe_based, data.glass_nl |>
                                                               dplyr::select(metadata.test$array_sentrix_id) |> 
                                                               t())|> 
    as.numeric()
  
  tmp <- rbind(tmp, 
    metadata.test |>  
      dplyr::select(array_sentrix_id,
                    rnaseq_cell.cycling.signature_it1
      ))
}
rm(cv)


metadata.glass_nl <- metadata.glass_nl |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))



## EM iteration 2 ----

### params ----

set.seed(123456)
cv_model_probe_based <- glmnet::cv.glmnet(
  data.glass_nl |>
    dplyr::select(metadata.glass_nl$array_sentrix_id) |> 
    t(),
  metadata.glass_nl$rnaseq_cell.cycling.signature_it1, alpha = 1, relax=F)


# @todo export & ggplot
plot(cv_model_probe_based)



### prediction ----


tmp <- data.frame()
for(cv in 0:9) {
  
  metadata.train <- metadata.glass_nl |> dplyr::filter(slice != cv)
  metadata.test <- metadata.glass_nl |> dplyr::filter(slice == cv)
  
  set.seed(123456)
  lm.probe_based <- glmnet::glmnet(
    data.glass_nl |>
      dplyr::select(metadata.train$array_sentrix_id) |> 
      t(),
    metadata.train$rnaseq_cell.cycling.signature_it1,
    alpha = 1, # 1 = lasso; 0 = ridge
    lambda = 0.01 # cv_model_probe_based$lambda.min
  )
  
  
  metadata.test$rnaseq_cell.cycling.signature_it2 <- predict(lm.probe_based, data.glass_nl |>
                                                               dplyr::select(metadata.test$array_sentrix_id) |> 
                                                               t())|> 
    as.numeric()
  
  tmp <- rbind(
    metadata.test |>  
      dplyr::select(array_sentrix_id,
                    rnaseq_cell.cycling.signature_it2
      ),
    tmp)
}
rm(cv)


metadata.glass_nl <- metadata.glass_nl |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))


plot(metadata.glass_nl$rnaseq_cell.cycling.signature, metadata.glass_nl$rnaseq_cell.cycling.signature_it1)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it1, metadata.glass_nl$rnaseq_cell.cycling.signature_it2)

plot(metadata.glass_nl$rnaseq_cell.cycling.signature, metadata.glass_nl$rnaseq_cell.cycling.signature_it2)


metadata.glass_nl |> dplyr::select(starts_with("rnaseq_cell.cycling.signature")) |> head(n=1)


## EM iteration 3 ----

### params ----


set.seed(123456)
cv_model_probe_based <- glmnet::cv.glmnet(
  data.glass_nl |>
    dplyr::select(metadata.glass_nl$array_sentrix_id) |> 
    t(),
  metadata.glass_nl$rnaseq_cell.cycling.signature_it2, alpha = 1, relax=F)


# @todo export & ggplot
plot(cv_model_probe_based)



### prediction ----


tmp <- data.frame()
for(cv in 0:9) {
  
  metadata.train <- metadata.glass_nl |> dplyr::filter(slice != cv)
  metadata.test <- metadata.glass_nl |> dplyr::filter(slice == cv)
  
  set.seed(123456)
  lm.probe_based <- glmnet::glmnet(
    data.glass_nl |>
      dplyr::select(metadata.train$array_sentrix_id) |> 
      t(),
    metadata.train$rnaseq_cell.cycling.signature_it2,
    alpha = 1, # 1 = lasso; 0 = ridge
    lambda = 0.01 # cv_model_probe_based$lambda.min
  )
  
  
  metadata.test$rnaseq_cell.cycling.signature_it3 <- predict(lm.probe_based, data.glass_nl |>
                                                               dplyr::select(metadata.test$array_sentrix_id) |> 
                                                               t())|> 
    as.numeric()
  
  tmp <- rbind(
    metadata.test |>  
      dplyr::select(array_sentrix_id,
                    rnaseq_cell.cycling.signature_it3
      ),
    tmp)
}
rm(cv)



metadata.glass_nl <- metadata.glass_nl |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))


plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it2, metadata.glass_nl$rnaseq_cell.cycling.signature_it3)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it1, metadata.glass_nl$rnaseq_cell.cycling.signature_it2)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature, metadata.glass_nl$rnaseq_cell.cycling.signature_it1)



metadata.glass_nl |> dplyr::select(starts_with("rnaseq_cell.cycling.signature")) |> head(n=1)



## EM iteration 4 ----

### params ----


set.seed(123456)
cv_model_probe_based <- glmnet::cv.glmnet(
  data.glass_nl |>
    dplyr::select(metadata.glass_nl$array_sentrix_id) |> 
    t(),
  metadata.glass_nl$rnaseq_cell.cycling.signature_it3, alpha = 1, relax=F)


# @todo export & ggplot
plot(cv_model_probe_based)



### prediction ----


tmp <- data.frame()
for(cv in 0:9) {
  
  metadata.train <- metadata.glass_nl |> dplyr::filter(slice != cv)
  metadata.test <- metadata.glass_nl |> dplyr::filter(slice == cv)
  
  set.seed(123456)
  lm.probe_based <- glmnet::glmnet(
    data.glass_nl |>
      dplyr::select(metadata.train$array_sentrix_id) |> 
      t(),
    metadata.train$rnaseq_cell.cycling.signature_it3,
    alpha = 1, # 1 = lasso; 0 = ridge
    lambda = 0.01 # cv_model_probe_based$lambda.min
  )
  
  
  metadata.test$rnaseq_cell.cycling.signature_it4 <- predict(lm.probe_based, data.glass_nl |>
                                                               dplyr::select(metadata.test$array_sentrix_id) |> 
                                                               t())|> 
    as.numeric()
  
  tmp <- rbind(
    metadata.test |>  
      dplyr::select(array_sentrix_id,
                    rnaseq_cell.cycling.signature_it4
      ),
    tmp)
}
rm(cv)



metadata.glass_nl <- metadata.glass_nl |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))


plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it3, metadata.glass_nl$rnaseq_cell.cycling.signature_it4)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it2, metadata.glass_nl$rnaseq_cell.cycling.signature_it3)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it1, metadata.glass_nl$rnaseq_cell.cycling.signature_it2)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature, metadata.glass_nl$rnaseq_cell.cycling.signature_it1)



metadata.glass_nl |> dplyr::select(starts_with("rnaseq_cell.cycling.signature")) |> head(n=1)




## EM iteration 5 ----

### params ----


set.seed(123456)
cv_model_probe_based <- glmnet::cv.glmnet(
  data.glass_nl |>
    dplyr::select(metadata.glass_nl$array_sentrix_id) |> 
    t(),
  metadata.glass_nl$rnaseq_cell.cycling.signature_it4, alpha = 1, relax=F)


# @todo export & ggplot
plot(cv_model_probe_based)



### prediction ----


tmp <- data.frame()
for(cv in 0:9) {
  
  metadata.train <- metadata.glass_nl |> dplyr::filter(slice != cv)
  metadata.test <- metadata.glass_nl |> dplyr::filter(slice == cv)
  
  set.seed(123456)
  lm.probe_based <- glmnet::glmnet(
    data.glass_nl |>
      dplyr::select(metadata.train$array_sentrix_id) |> 
      t(),
    metadata.train$rnaseq_cell.cycling.signature_it4,
    alpha = 1, # 1 = lasso; 0 = ridge
    lambda = 0.01 # cv_model_probe_based$lambda.min
  )
  
  
  metadata.test$rnaseq_cell.cycling.signature_it5 <- predict(lm.probe_based, data.glass_nl |>
                                                               dplyr::select(metadata.test$array_sentrix_id) |> 
                                                               t())|> 
    as.numeric()
  
  tmp <- rbind(
    metadata.test |>  
      dplyr::select(array_sentrix_id,
                    rnaseq_cell.cycling.signature_it5
      ),
    tmp)
}
rm(cv)



metadata.glass_nl <- metadata.glass_nl |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))


plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it4, metadata.glass_nl$rnaseq_cell.cycling.signature_it5)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it3, metadata.glass_nl$rnaseq_cell.cycling.signature_it4)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it2, metadata.glass_nl$rnaseq_cell.cycling.signature_it3)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it1, metadata.glass_nl$rnaseq_cell.cycling.signature_it2)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature, metadata.glass_nl$rnaseq_cell.cycling.signature_it1)



metadata.glass_nl |> dplyr::select(starts_with("rnaseq_cell.cycling.signature")) |> head(n=1)


## EM iteration 6 ----

### params ----


set.seed(123456)
cv_model_probe_based <- glmnet::cv.glmnet(
  data.glass_nl |>
    dplyr::select(metadata.glass_nl$array_sentrix_id) |> 
    t(),
  metadata.glass_nl$rnaseq_cell.cycling.signature_it5, alpha = 1, relax=F)


# @todo export & ggplot
plot(cv_model_probe_based)



### prediction ----


tmp <- data.frame()
for(cv in 0:9) {
  
  metadata.train <- metadata.glass_nl |> dplyr::filter(slice != cv)
  metadata.test <- metadata.glass_nl |> dplyr::filter(slice == cv)
  
  set.seed(123456)
  lm.probe_based <- glmnet::glmnet(
    data.glass_nl |>
      dplyr::select(metadata.train$array_sentrix_id) |> 
      t(),
    metadata.train$rnaseq_cell.cycling.signature_it5,
    alpha = 1, # 1 = lasso; 0 = ridge
    lambda = 0.01 # cv_model_probe_based$lambda.min
  )
  
  
  metadata.test$rnaseq_cell.cycling.signature_it6 <- predict(lm.probe_based, data.glass_nl |>
                                                               dplyr::select(metadata.test$array_sentrix_id) |> 
                                                               t())|> 
    as.numeric()
  
  tmp <- rbind(
    metadata.test |>  
      dplyr::select(array_sentrix_id,
                    rnaseq_cell.cycling.signature_it6
      ),
    tmp)
}
rm(cv)



metadata.glass_nl <- metadata.glass_nl |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))


plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it5, metadata.glass_nl$rnaseq_cell.cycling.signature_it6)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it4, metadata.glass_nl$rnaseq_cell.cycling.signature_it5)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it3, metadata.glass_nl$rnaseq_cell.cycling.signature_it4)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it2, metadata.glass_nl$rnaseq_cell.cycling.signature_it3)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it1, metadata.glass_nl$rnaseq_cell.cycling.signature_it2)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature, metadata.glass_nl$rnaseq_cell.cycling.signature_it1)



metadata.glass_nl |> dplyr::select(starts_with("rnaseq_cell.cycling.signature")) |> head(n=1)





# plots ----


plot(metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, metadata.glass_nl$rnaseq_cell.cycling.signature)
plot(metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, metadata.glass_nl$rnaseq_cell.cycling.signature_it1)
plot(metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, metadata.glass_nl$rnaseq_cell.cycling.signature_it2)
plot(metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, metadata.glass_nl$rnaseq_cell.cycling.signature_it3)
plot(metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, metadata.glass_nl$rnaseq_cell.cycling.signature_it4)
plot(metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, metadata.glass_nl$rnaseq_cell.cycling.signature_it5)
plot(metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, metadata.glass_nl$rnaseq_cell.cycling.signature_it6)

plot(metadata.glass_nl$rnaseq_cell.cycling.signature, metadata.glass_nl$array_PC3)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it1, metadata.glass_nl$array_PC3)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it2, metadata.glass_nl$array_PC3)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it3, metadata.glass_nl$array_PC3)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it4, metadata.glass_nl$array_PC3)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it5, metadata.glass_nl$array_PC3)
plot(metadata.glass_nl$rnaseq_cell.cycling.signature_it6, metadata.glass_nl$array_PC3)


ggplot(metadata.glass_nl, aes(x=array_A_IDH_HG__A_IDH_LG_lr_v12.8, y=array_PC3, col=array_mnp_predictBrain_v12.8_cal_class)) +
  geom_point()


ggplot(metadata.glass_nl, aes(x=rnaseq_cell.cycling.signature_it6, y=array_A_IDH_HG__A_IDH_LG_lr_v12.8, col=array_mnp_predictBrain_v12.8_cal_class)) +
  geom_point()


ggplot(metadata.glass_nl, aes(x=rnaseq_cell.cycling.signature_it6, y=metadata.glass_nl$array_median.glass_nl_supervised.methylation, col=array_mnp_predictBrain_v12.8_cal_class)) +
  geom_point()





