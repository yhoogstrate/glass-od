#!/usr/bin/env R


# load data ----


library(ggplot2)



lin_repredict <- function(data, vector, output_name) {
  stopifnot(rownames(data) == names(vector))
  
  metadata <- data.frame(vector) |> 
    tibble::rownames_to_column('array_sentrix_id') |> 
    dplyr::mutate(i = 1:dplyr::n()) |> 
    dplyr::mutate(slice = i %% 10, i = NULL)
  
  
  set.seed(123456)
  cv_model <- glmnet::cv.glmnet(data, metadata$vector, alpha = 1, relax=F)
  
  
  tmp <- data.frame()
  for(cv in 0:9) {
    metadata.train <- metadata |> dplyr::filter(slice != cv)
    metadata.test <- metadata |> dplyr::filter(slice == cv)
    
    set.seed(123456)
    lm.probe_based <- glmnet::glmnet(data |>t() |> as.data.frame() |> dplyr::select(metadata.train$array_sentrix_id) |> t()      ,
      metadata.train |> dplyr::pull(vector, name="array_sentrix_id"),
      alpha = 1, # 1 = lasso; 0 = ridge
      lambda = cv_model$lambda.min
    )
    
    
    out.df <- predict(lm.probe_based, data |>t() |> as.data.frame() |> dplyr::select(metadata.test$array_sentrix_id) |> t()) |> 
      as.data.frame() |> 
      dplyr::rename("out" = 1) |> 
      tibble::rownames_to_column("array_sentrix_id")
    
    tmp <- rbind(tmp, out.df)
  }

  # reorder
  tmp <- data.frame(array_sentrix_id = rownames(data)) |> 
    dplyr::left_join(tmp, by=c('array_sentrix_id' = 'array_sentrix_id'))

  return(tmp)
}




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



## Cell cycling signature ----


tmp.cc.1 = lin_repredict(data, metadata.glass_nl |> dplyr::pull(rnaseq_cell.cycling.signature, name="array_sentrix_id"))
tmp.cc.2 = lin_repredict(data, tmp.cc.1 |> dplyr::pull(out, name="array_sentrix_id"))
tmp.cc.3 = lin_repredict(data, tmp.cc.2 |> dplyr::pull(out, name="array_sentrix_id"))
tmp.cc.4 = lin_repredict(data, tmp.cc.3 |> dplyr::pull(out, name="array_sentrix_id"))
tmp.cc.5 = lin_repredict(data, tmp.cc.4 |> dplyr::pull(out, name="array_sentrix_id"))
tmp.cc.6 = lin_repredict(data, tmp.cc.5 |> dplyr::pull(out, name="array_sentrix_id"))



### export ----

exp_Cell_Cycling <- tmp.cc.1 |> dplyr::rename(rnaseq_Cell_Cycling.signature_EM_it1 = out) |> 
  dplyr::left_join(
    tmp.cc.2 |> dplyr::rename(rnaseq_Cell_Cycling.signature_EM_it2 = out),
    by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')
  ) |> 
  dplyr::left_join(
    tmp.cc.3 |> dplyr::rename(rnaseq_Cell_Cycling.signature_EM_it3 = out),
    by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')
  ) |> 
  dplyr::left_join(
    tmp.cc.4 |> dplyr::rename(rnaseq_Cell_Cycling.signature_EM_it4 = out),
    by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')
  ) |> 
  dplyr::left_join(
    tmp.cc.5 |> dplyr::rename(rnaseq_Cell_Cycling.signature_EM_it5 = out),
    by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')
  ) |> 
  dplyr::left_join(
    tmp.cc.6 |> dplyr::rename(rnaseq_Cell_Cycling.signature_EM_it6 = out),
    by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')
  )

saveRDS(exp_Cell_Cycling, file="cache/analysis_DNAm_cell_cycling_signature_EM_CC.Rds")



## ECM signature ----



tmp.ecm.1 = lin_repredict(data, metadata.glass_nl |> dplyr::pull(rnaseq_ECM.signature, name="array_sentrix_id"))
tmp.ecm.2 = lin_repredict(data, tmp.ecm.1 |> dplyr::pull(out, name="array_sentrix_id"))
tmp.ecm.3 = lin_repredict(data, tmp.ecm.2 |> dplyr::pull(out, name="array_sentrix_id"))
tmp.ecm.4 = lin_repredict(data, tmp.ecm.3 |> dplyr::pull(out, name="array_sentrix_id"))
tmp.ecm.5 = lin_repredict(data, tmp.ecm.4 |> dplyr::pull(out, name="array_sentrix_id"))
tmp.ecm.6 = lin_repredict(data, tmp.ecm.5 |> dplyr::pull(out, name="array_sentrix_id"))



### export ----


exp_ecm <- tmp.ecm.1 |> dplyr::rename(rnaseq_ECM.signature_EM_it1 = out) |> 
  dplyr::left_join(
    tmp.ecm.2 |> dplyr::rename(rnaseq_ECM.signature_EM_it2 = out),
    by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')
  ) |> 
  dplyr::left_join(
    tmp.ecm.3 |> dplyr::rename(rnaseq_ECM.signature_EM_it3 = out),
    by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')
  ) |> 
  dplyr::left_join(
    tmp.ecm.4 |> dplyr::rename(rnaseq_ECM.signature_EM_it4 = out),
    by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')
  ) |> 
  dplyr::left_join(
    tmp.ecm.5 |> dplyr::rename(rnaseq_ECM.signature_EM_it5 = out),
    by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')
  ) |> 
  dplyr::left_join(
    tmp.ecm.6 |> dplyr::rename(rnaseq_ECM.signature_EM_it6 = out),
    by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('','')
  )



saveRDS(exp_Cell_Cycling, file="cache/analysis_DNAm_cell_cycling_signature_EM_ECM.Rds")




# sandbox plots ----




plot(metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, metadata.glass_nl$rnaseq_cell.cycling.signature)
plot(metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, metadata.glass_nl$rnaseq_cell.cycling.signature_it1)
plot(metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, metadata.glass_nl$rnaseq_cell.cycling.signature_it2)
plot(metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, metadata.glass_nl$rnaseq_cell.cycling.signature_it3)
plot(metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, metadata.glass_nl$rnaseq_cell.cycling.signature_it4)
plot(metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, metadata.glass_nl$rnaseq_cell.cycling.signature_it5)
plot(metadata.glass_nl$array_A_IDH_HG__A_IDH_LG_lr_v12.8, metadata.glass_nl$rnaseq_cell.cycling.signature_it6)


plot(metadata.glass_nl$array_median.glass_nl_supervised.methylation, metadata.glass_nl$rnaseq_cell.cycling.signature)
plot(metadata.glass_nl$array_median.glass_nl_supervised.methylation, metadata.glass_nl$rnaseq_cell.cycling.signature_it1)
plot(metadata.glass_nl$array_median.glass_nl_supervised.methylation, metadata.glass_nl$rnaseq_cell.cycling.signature_it2)
plot(metadata.glass_nl$array_median.glass_nl_supervised.methylation, metadata.glass_nl$rnaseq_cell.cycling.signature_it3)
plot(metadata.glass_nl$array_median.glass_nl_supervised.methylation, metadata.glass_nl$rnaseq_cell.cycling.signature_it4)
plot(metadata.glass_nl$array_median.glass_nl_supervised.methylation, metadata.glass_nl$rnaseq_cell.cycling.signature_it5)
plot(metadata.glass_nl$array_median.glass_nl_supervised.methylation, metadata.glass_nl$rnaseq_cell.cycling.signature_it6)



plot(metadata.glass_nl$rnaseq_cell.cycling.signature, metadata.glass_nl$array_PC3)

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


## convergence plot CC ----


metadata.glass_nl <- metadata.glass_nl |> 
  dplyr::left_join(exp_Cell_Cycling, by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('',''))


plt <- rbind(
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_cell.cycling.signature) |> 
    dplyr::rename(y = rnaseq_Cell_Cycling.signature_EM_it1) |> 
    dplyr::mutate(facet = "Iteration 1") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it1) |> 
    dplyr::rename(y = rnaseq_Cell_Cycling.signature_EM_it2) |> 
    dplyr::mutate(facet = "Iteration 2") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it2) |> 
    dplyr::rename(y = rnaseq_Cell_Cycling.signature_EM_it3) |> 
    dplyr::mutate(facet = "Iteration 3") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it3) |> 
    dplyr::rename(y = rnaseq_Cell_Cycling.signature_EM_it4) |> 
    dplyr::mutate(facet = "Iteration 4") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it4) |> 
    dplyr::rename(y = rnaseq_Cell_Cycling.signature_EM_it5) |> 
    dplyr::mutate(facet = "Iteration 5") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it5) |> 
    dplyr::rename(y = rnaseq_Cell_Cycling.signature_EM_it6) |> 
    dplyr::mutate(facet = "Iteration 6") |> 
    dplyr::select(x, y, facet)
)



ggplot(plt, aes(x=x, y=y)) +
  facet_wrap(~facet, scales="free_x", ncol=22) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho", size=theme_cellpress_size) +
  geom_vline(xintercept = -20, lwd = theme_cellpress_lwd) +
  geom_point() +
  labs(x = "cell cycling signature (input)", y="cell cycling signature (predicted)") +
  theme_cellpress




## convergence plot ECM ----


metadata.glass_nl <- metadata.glass_nl |> 
  dplyr::left_join(exp_ecm, by=c('array_sentrix_id'='array_sentrix_id'),suffix=c('',''))


plt <- rbind(
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature) |> 
    dplyr::rename(y = rnaseq_ECM.signature_EM_it1) |> 
    dplyr::mutate(facet = "Iteration 1") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it1) |> 
    dplyr::rename(y = rnaseq_ECM.signature_EM_it2) |> 
    dplyr::mutate(facet = "Iteration 2") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it2) |> 
    dplyr::rename(y = rnaseq_ECM.signature_EM_it3) |> 
    dplyr::mutate(facet = "Iteration 3") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it3) |> 
    dplyr::rename(y = rnaseq_ECM.signature_EM_it4) |> 
    dplyr::mutate(facet = "Iteration 4") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it4) |> 
    dplyr::rename(y = rnaseq_ECM.signature_EM_it5) |> 
    dplyr::mutate(facet = "Iteration 5") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it5) |> 
    dplyr::rename(y = rnaseq_ECM.signature_EM_it6) |> 
    dplyr::mutate(facet = "Iteration 6") |> 
    dplyr::select(x, y, facet)
)



ggplot(plt, aes(x=x, y=y)) +
  facet_wrap(~facet, scales="free_x", ncol=22) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho", size=theme_cellpress_size) +
  geom_vline(xintercept = -20, lwd = theme_cellpress_lwd) +
  geom_point() +
  labs(x = "ECM signature (input)", y="ECM signature (predicted)") +
  theme_cellpress



## ECM x CC ----


plt <- rbind(
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_cell.cycling.signature) |> 
    dplyr::rename(y = rnaseq_ECM.signature) |> 
    dplyr::mutate(facet = "RNA data") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it1) |> 
    dplyr::rename(y = rnaseq_ECM.signature_EM_it1) |> 
    dplyr::mutate(facet = "Iteration 1") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it2) |> 
    dplyr::rename(y = rnaseq_ECM.signature_EM_it2) |> 
    dplyr::mutate(facet = "Iteration 2") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it3) |> 
    dplyr::rename(y = rnaseq_ECM.signature_EM_it3) |> 
    dplyr::mutate(facet = "Iteration 3") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it4) |> 
    dplyr::rename(y = rnaseq_ECM.signature_EM_it4) |> 
    dplyr::mutate(facet = "Iteration 4") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it5) |> 
    dplyr::rename(y = rnaseq_ECM.signature_EM_it5) |> 
    dplyr::mutate(facet = "Iteration 5") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it6) |> 
    dplyr::rename(y = rnaseq_ECM.signature_EM_it6) |> 
    dplyr::mutate(facet = "Iteration 6") |> 
    dplyr::select(x, y, facet)
)



ggplot(plt, aes(x=x, y=y)) +
  facet_wrap(~facet, scales="free_x", ncol=22) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho", size=theme_cellpress_size) +
  geom_vline(xintercept = -20, lwd = theme_cellpress_lwd) +
  geom_point() +
  labs(x = "Cell Cycling signature", y="ECM signature") +
  theme_cellpress




## ECM x demeth wies ----


plt <- rbind(
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature) |> 
    dplyr::rename(y = array_median.glass_nl_supervised.methylation) |> 
    dplyr::mutate(facet = "a RNA data") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it1) |> 
    dplyr::rename(y = array_median.glass_nl_supervised.methylation) |> 
    dplyr::mutate(facet = "Iteration 1") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it2) |> 
    dplyr::rename(y = array_median.glass_nl_supervised.methylation) |> 
    dplyr::mutate(facet = "Iteration 2") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it3) |> 
    dplyr::rename(y = array_median.glass_nl_supervised.methylation) |> 
    dplyr::mutate(facet = "Iteration 3") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it4) |> 
    dplyr::rename(y = array_median.glass_nl_supervised.methylation) |> 
    dplyr::mutate(facet = "Iteration 4") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it5) |> 
    dplyr::rename(y = array_median.glass_nl_supervised.methylation) |> 
    dplyr::mutate(facet = "Iteration 5") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it6) |> 
    dplyr::rename(y = array_median.glass_nl_supervised.methylation) |> 
    dplyr::mutate(facet = "Iteration 6") |> 
    dplyr::select(x, y, facet)
)



ggplot(plt, aes(x=x, y=y)) +
  facet_wrap(~facet, scales="free_x", ncol=22) +
  geom_vline(xintercept = -20, lwd = theme_cellpress_lwd) +
  geom_point() +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="red", cor.coef.name ="rho", size=theme_cellpress_size) +
  labs(x = "ECM signature [RNA & EM classified]", y="Wies demeth signature") +
  theme_cellpress






## CC x demeth wies ----


plt <- rbind(
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_cell.cycling.signature) |> 
    dplyr::rename(y = array_median.glass_nl_supervised.methylation) |> 
    dplyr::mutate(facet = "a RNA data") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it1) |> 
    dplyr::rename(y = array_median.glass_nl_supervised.methylation) |> 
    dplyr::mutate(facet = "Iteration 1") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it2) |> 
    dplyr::rename(y = array_median.glass_nl_supervised.methylation) |> 
    dplyr::mutate(facet = "Iteration 2") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it3) |> 
    dplyr::rename(y = array_median.glass_nl_supervised.methylation) |> 
    dplyr::mutate(facet = "Iteration 3") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it4) |> 
    dplyr::rename(y = array_median.glass_nl_supervised.methylation) |> 
    dplyr::mutate(facet = "Iteration 4") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it5) |> 
    dplyr::rename(y = array_median.glass_nl_supervised.methylation) |> 
    dplyr::mutate(facet = "Iteration 5") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it6) |> 
    dplyr::rename(y = array_median.glass_nl_supervised.methylation) |> 
    dplyr::mutate(facet = "Iteration 6") |> 
    dplyr::select(x, y, facet)
)



ggplot(plt, aes(x=x, y=y)) +
  facet_wrap(~facet, scales="free_x", ncol=22) +
  geom_vline(xintercept = -20, lwd = theme_cellpress_lwd) +
  geom_point() +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), cor.coef.name ="rho", size=theme_cellpress_size, col="red") +
  labs(x = "Cell Cycling signature", y="ECM signature") +
  theme_cellpress





## CC x CGC ----


plt <- rbind(
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_cell.cycling.signature) |> 
    dplyr::rename(y = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(facet = "a RNA data") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it1) |> 
    dplyr::rename(y = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(facet = "Iteration 1") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it2) |> 
    dplyr::rename(y = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(facet = "Iteration 2") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it3) |> 
    dplyr::rename(y = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(facet = "Iteration 3") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it4) |> 
    dplyr::rename(y = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(facet = "Iteration 4") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it5) |> 
    dplyr::rename(y = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(facet = "Iteration 5") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_Cell_Cycling.signature_EM_it6) |> 
    dplyr::rename(y = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(facet = "Iteration 6") |> 
    dplyr::select(x, y, facet)
)



ggplot(plt, aes(x=x, y=y)) +
  facet_wrap(~facet, scales="free_x", ncol=22) +
  geom_vline(xintercept = -20, lwd = theme_cellpress_lwd) +
  geom_point() +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), cor.coef.name ="rho", size=theme_cellpress_size, col="red") +
  labs(x = "Cell Cycling signature", y="ECM signature") +
  theme_cellpress




## ECM x CGC ----



plt <- rbind(
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature) |> 
    dplyr::rename(y = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(facet = "a RNA data") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it1) |> 
    dplyr::rename(y = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(facet = "Iteration 1") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it2) |> 
    dplyr::rename(y = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(facet = "Iteration 2") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it3) |> 
    dplyr::rename(y = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(facet = "Iteration 3") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it4) |> 
    dplyr::rename(y = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(facet = "Iteration 4") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it5) |> 
    dplyr::rename(y = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(facet = "Iteration 5") |> 
    dplyr::select(x, y, facet),
  metadata.glass_nl |> 
    dplyr::rename(x = rnaseq_ECM.signature_EM_it6) |> 
    dplyr::rename(y = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(facet = "Iteration 6") |> 
    dplyr::select(x, y, facet)
)



ggplot(plt, aes(x=x, y=y)) +
  facet_wrap(~facet, scales="free_x", ncol=22) +
  geom_vline(xintercept = -13, lwd = theme_cellpress_lwd) +
  geom_point() +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="red", cor.coef.name ="rho", size=theme_cellpress_size) +
  labs(x = "ECM signature [RNA & EM classified]", y="Wies demeth signature") +
  theme_cellpress








