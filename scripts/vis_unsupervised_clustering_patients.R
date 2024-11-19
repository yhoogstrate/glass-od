#!/usr/bin/env R

# load data ----


library(ggplot2)

source('scripts/load_constants.R')
source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# plt ----


## PC hierarchical ----


library(ComplexHeatmap)

plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::mutate(heidi_label = dplyr::case_when(
    array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") ~ array_mnp_predictBrain_v12.8_cal_class,
    T ~ "other"
  ))




mat <- plt |> 
  dplyr::select(resection_id, contains("array_PC")) |> 
  dplyr::mutate(array_PC1 = NULL) |> 
  tibble::column_to_rownames('resection_id')  |> 
  #scale() |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(plt$resection_id) |> 
  head(n=19) |> 
  tibble::rownames_to_column('rownames') |> 
  dplyr::mutate(rownames = gsub("^array_","", rownames)) |> 
  tibble::column_to_rownames('rownames')




heidi <- HeatmapAnnotation(`1. MNP` = (plt$heidi_label),
                           `2. Radio Therapy` = ifelse(plt$resection_treatment_status_radio, 'Yes', 'No'),
                           `3. Chemo` = ifelse(plt$resection_treatment_status_chemo, 'Yes', 'No'),
                           `4. WHO grade` = paste0("Grade ", plt$resection_tumor_grade),
                           `5. resection` = ifelse(plt$resection_number == 1, "Primary", "Recurrent"),
                           `6. log(frac. detP)` = log(plt$array_percentage.detP.signi),
                           `7. PC1` = plt$array_PC1 ,
                           `8. PCHorvathS2018` = plt$array_dnaMethyAge__PCHorvathS2018,
                           `9. CGC[ac]` = plt$array_A_IDH_HG__A_IDH_LG_lr__lasso_fit
                           ,
                           
                           col  = list(
                             `1. MNP` = palette_mnp_12.8_7,
                             `2. Radio Therapy` = palette_yes_no_1,
                             `3. Chemo` = palette_yes_no_1,
                             `4. WHO grade` = palette_g2_g3,
                             `5. resection` = palette_p_r,
                             `6. log(frac. detP)` = circlize::colorRamp2(seq(min(log(plt$array_percentage.detP.signi)),max(log(plt$array_percentage.detP.signi)),length.out=100), col6(100)) ,
                             `7. PC1` = circlize::colorRamp2(seq(min(plt$array_PC1),max(plt$array_PC1),length.out=100), col6(100)),
                             `8. PCHorvathS2018` = circlize::colorRamp2(seq(max(plt$array_dnaMethyAge__PCHorvathS2018),min(plt$array_dnaMethyAge__PCHorvathS2018),length.out=100), col6(100)) ,
                             `9. CGC[ac]` = circlize::colorRamp2(seq(
                              
                               0.8 * max(abs(plt$array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)),
                              -0.8 * max(abs(plt$array_A_IDH_HG__A_IDH_LG_lr__lasso_fit)),
                              
                              length.out=100), col6(100))
                            )
)
Heatmap(mat, name = "mat", 
        top_annotation = heidi,
        cluster_rows = F,
        clustering_method_columns = "ward.D2",
        
        col = circlize::colorRamp2(seq(-410, 410, length.out=100), col6(100))
)



## MDS ----


mat.mds_PC2_PC20 <- plt |> 
  dplyr::select(array_sentrix_id, contains("array_PC")) |> 
  dplyr::mutate(array_PC1 = NULL) |> 
  #dplyr::mutate(array_PC2 = NULL) |> 
  tibble::column_to_rownames('array_sentrix_id')  |> 
  #scale() |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(plt$array_sentrix_id) |> 
  head(n=19) |> 
  tibble::rownames_to_column('rownames') |> 
  dplyr::mutate(rownames = gsub("^array_","", rownames)) |> 
  tibble::column_to_rownames('rownames')

mat.mds_PC3_PC20 <- plt |> 
  dplyr::select(array_sentrix_id, contains("array_PC")) |> 
  dplyr::mutate(array_PC1 = NULL) |> 
  dplyr::mutate(array_PC2 = NULL) |> 
  tibble::column_to_rownames('array_sentrix_id')  |> 
  #scale() |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(plt$array_sentrix_id) |> 
  head(n=19) |> 
  tibble::rownames_to_column('rownames') |> 
  dplyr::mutate(rownames = gsub("^array_","", rownames)) |> 
  tibble::column_to_rownames('rownames')


fit_PC2_PC20 <- cmdscale(dist(t(mat.mds_PC2_PC20)), eig=TRUE, k=2)
fit_PC3_PC20 <- cmdscale(dist(t(mat.mds_PC3_PC20)), eig=TRUE, k=2)



plt.mds_PC2_PC20 <- fit_PC2_PC20$points |> 
  as.data.frame() |> 
  dplyr::rename('MDS: dim 1' = 1, 'MDS: dim 2' = 2) |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::left_join(plt, by=c('array_sentrix_id'='array_sentrix_id'))



plt.mds_PC3_PC20 <- fit_PC3_PC20$points |> 
  as.data.frame() |> 
  dplyr::rename('MDS: dim 1' = 1, 'MDS: dim 2' = 2) |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::left_join(plt, by=c('array_sentrix_id'='array_sentrix_id'))


plt.mds <- rbind(
  plt.mds_PC2_PC20 |> 
    dplyr::mutate(`used components` = 'Included principal components: PC2 - PC20')
  ,
  plt.mds_PC3_PC20 |> 
    dplyr::mutate(`used components` = 'Included principal components: PC3 - PC20')
)



ggplot(plt.mds, aes(x=`MDS: dim 1`, y=`MDS: dim 2`, col=heidi_label)) +
  facet_wrap(~`used components`) +
  geom_point(size=theme_nature_size/3) +
  scale_color_manual(values=palette_mnp_12.8_7) +
  theme_nature +
  labs(subtitle=format_subtitle("Unsupervised MDS GLASS-OD"))

ggsave("output/figures/vis_unsupervised_clustering__MDS_on_PCs.pdf", width=3.8,height=2.5)



