#!/usr/bin/env R




source('scripts/load_constants.R')
source('scripts/load_functions.R')

source('scripts/load_palette.R')
source('scripts/load_themes.R')



library(ggplot2)



# load all M-values ----



if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}



# load GLASS-OD ----


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}



metadata.glass_od <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES)  # only include FFPE's - FF are different regardless
  


# load GLASS-NL ----


if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}


metadata.glass_nl <- glass_nl.metadata.array_samples |> 
  filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES)
  


# plt UMAP ----




metadata <- rbind(
  metadata.glass_od |> 
    dplyr::rename(`CGC` = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
    dplyr::select(array_sentrix_id, patient_sex, resection_id, array_mnp_predictBrain_v12.8_cal_class, array_RFpurity_estimate , CGC) |> 
    dplyr::mutate(disease = "Oligodendroglioma"),
  metadata.glass_nl |> 
    dplyr::rename(CGC = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
    dplyr::mutate(resection_id = Sample_Name) |>
    dplyr::rename(array_RFpurity_estimate = array_RFpurity.estimate) |> 
    dplyr::mutate(patient_sex = ifelse(Sample_Sex == "M", "male", "female")) |> 
    dplyr::select(array_sentrix_id, patient_sex, resection_id, array_mnp_predictBrain_v12.8_cal_class, array_RFpurity_estimate , CGC) |> 
    dplyr::mutate(disease = "Astrocytoma")
) |> 
  dplyr::mutate(col2 = dplyr::case_when(
    array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG", "OLIGOSARC_IDH") ~ "high grade classes",
    array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_LG", "O_IDH") ~ "low grade classes",
    T ~ "other"
  ))



cpgs_x_y <- metadata.cg_probes.epic |> 
  dplyr::filter(CHR_hg38 %in% c("chrM" ,"chrX", "chrY") | is.na(CHR_hg38)) |> 
  dplyr::pull('probe_id')



plt <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::filter(array_sentrix_id %in% cpgs_x_y == F) |> 
  tibble::column_to_rownames('array_sentrix_id') |> 
  dplyr::select(metadata$array_sentrix_id)


plt$mad_value <- pbapply::pbapply(plt, 1, mad) # sparser mem than rowwise?


plt <- plt |> 
  dplyr::arrange(desc(mad_value)) |> 
  dplyr::mutate(mad_value = NULL) |> 
  dplyr::slice_head(n = 7500)




umap_result <- plt |> 
  t() |> 
  uwot::umap(n_neighbors = 10)



# Convert to data frame for plotting
umap_df <- as.data.frame(umap_result) |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::rename(UMAP1 = V1, UMAP2 = V2) |> 
  dplyr::left_join(metadata,by=c('array_sentrix_id'='array_sentrix_id')) |> 
  dplyr::mutate(col = paste0(disease , ": ", col2))  |> 
  transform(
    RFpurity_group = ifelse(
      array_RFpurity_estimate > median(array_RFpurity_estimate, na.rm = TRUE),
      "High",
      "Low"
    ))  |> 
  transform(
    CGC_group = ifelse(
      CGC > median(CGC, na.rm = TRUE),
      "High",
      "Low"
    ))




# Plot with ggplot2
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, col = col, label=resection_id)) + # color by purity and CGC?
  geom_point(size = theme_nature_size/3) +
  theme_nature +
  labs(x = "UMAP 1", y = "UMAP 2", col="") +
  scale_color_manual(values=c(
    #`Astrocytoma: high grade classes` = as.character(palette_mnp_12.8_6["A_IDH_HG"]),
    #`Astrocytoma: low grade classes` = as.character(palette_mnp_12.8_6["A_IDH_LG"]),
    
    #`Oligodendroglioma: high grade classes` = as.character(palette_mnp_12.8_6["OLIGOSARC_IDH"]),
    #`Oligodendroglioma: low grade classes` = as.character(palette_mnp_12.8_6["O_IDH"]),
    
    `Astrocytoma: high grade classes` = 'darkblue',
    `Astrocytoma: low grade classes` =  'lightblue',
    
    `Oligodendroglioma: high grade classes` = 'darkgreen',
    `Oligodendroglioma: low grade classes` = 'lightgreen',
    
    
    `Astrocytoma: other` = "red",
    `Oligodendroglioma: other` = "orange"
  )) +
  coord_fixed(ratio = 1) # keep natural umap aspect ratio


ggsave("output/figures/vis_UMAP_oligo_astro.pdf", width = 3.5, height = 3.5)





# Plot with ggplot2
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, col = RFpurity_group, label=resection_id)) + # color by purity and CGC?
  geom_point(size = theme_nature_size/3) +
  theme_nature +
  labs(x = "UMAP 1", y = "UMAP 2", col="") +
  scale_color_manual(values=c(
    `High` = 'red',
    `Low` = 'darkgreen'
  )) +
  coord_fixed(ratio = 1) # keep natural umap aspect ratio



# Plot with ggplot2
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, col = CGC_group, label=resection_id)) + # color by purity and CGC?
  geom_point(size = theme_nature_size/3) +
  theme_nature +
  labs(x = "UMAP 1", y = "UMAP 2", col="") +
  scale_color_manual(values=c(
    `High` = 'red',
    `Low` = 'darkgreen'
  )) +
  coord_fixed(ratio = 1) # keep natural umap aspect ratio





