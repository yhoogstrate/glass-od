#!/usr/bin/env R


source('scripts/load_constants.R')
source('scripts/load_functions.R')



if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# GEO ----
## GLASS-OD ----


tmp.1 <-  glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) 

tmp.2 <- glass_od.metadata.array_samples |> 
  filter_OD_validation_idats(CONST_N_OD_VALIDATION_INCLUDED_SAMPLES) |> 
  dplyr::filter(patient_center_name == "Rotterdam")


main <- glass_od.metadata.array_samples |> 
  dplyr::filter(array_sentrix_id %in% c(c(tmp.1$array_sentrix_id, tmp.2$array_sentrix_id))) |> 
  dplyr::mutate(`Sample name` = resection_id) |> 
  dplyr::mutate(`title` = resection_id) |> 
  dplyr::mutate(`organism` = 'Homo sapiens') |> 
  dplyr::mutate(`*idat file 1` = base::basename(array_channel_green)) |> 
  dplyr::mutate(`*idat file 2` = base::basename(array_channel_red)) |> 
  dplyr::mutate(`Sex` = patient_sex) |> 
  dplyr::mutate(set = ifelse(array_sentrix_id %in% tmp.1$array_sentrix_id, "GLASS-OD", "GLASS-OD validation set")) |> 
  dplyr::mutate(`disease state` = ifelse(patient_suspected_noncodel == F, paste0("Oligodendroglioma, WHO grade ",resection_tumor_grade , ", resection #",resection_number), "-")) |> 
  dplyr::mutate(Tissue = "tumor genomic DNA") |> 
  dplyr::mutate(`**cell line`= "") |> 
  dplyr::mutate(`**cell type`= "") |> 
  dplyr::mutate(`genotype`= "") |> 
  dplyr::mutate(treatment = "") |> 
  dplyr::mutate(molecule = "genomic DNA") |> 
  dplyr::mutate(label = "Cy3 and Cy5") |> 
  dplyr::mutate(description = paste0(Tissue, " from patient ",patient_id,": ", `disease state`)) |> 
  dplyr::mutate(`BeadChip or GEO platform id` = 'GPL21145') |> 
  
  dplyr::select(`Sample name`, title, organism, `*idat file 1`, `*idat file 2`, Sex, `disease state`,
                
                Tissue, `**cell line`, `**cell type`, genotype, treatment,
                
                
                molecule, label, description, `BeadChip or GEO platform id`, set) |> 
  dplyr::arrange(set, title)






write.table(main, file="/tmp/geo_main_dataset.txt", sep="\t", quote=F, row.names=F)





cp <- glass_od.metadata.array_samples |> 
  dplyr::filter(array_sentrix_id %in% c(c(tmp.1$array_sentrix_id, tmp.2$array_sentrix_id))) |>
  #dplyr::filter(array_sentrix_id %in% c(c(tmp.2$array_sentrix_id))) |>
  dplyr::mutate(`*idat file 1` = base::basename(array_channel_green)) |> 
  dplyr::mutate(`*idat file 2` = base::basename(array_channel_red)) |> 
  
  dplyr::mutate(array_channel_green = paste0("/home/youri/projects/glass-od/", array_channel_green)) |> 
  dplyr::mutate(array_channel_red = paste0("/home/youri/projects/glass-od/", array_channel_red)) |> 
  
  dplyr::select(array_channel_green, `*idat file 1`,  array_channel_red, `*idat file 2`)


write.table(cp, file="/tmp/cp.txt", sep="\t", quote=T, row.names=F)






if(!exists('data.mvalues.probes')) {
  source('scripts/load_mvalues_hq_samples.R')
}




tmp <- glass_od.metadata.array_samples |> 
  dplyr::filter(array_sentrix_id %in% c(c(tmp.1$array_sentrix_id, tmp.2$array_sentrix_id))) |>
  dplyr::select(resection_id, array_sentrix_id) |> 
  dplyr::arrange(resection_id)


mvalues_to_table <- read.delim("cache/GPL21145-48548.txt", comment.char = "#", stringsAsFactors = F) |> # preserve order and probes from this file
  dplyr::select(ID) |> 
  dplyr::rename(ID_REF = ID) |> 
  dplyr::left_join(
    data.mvalues.hq_samples |> 
      dplyr::select(tmp$array_sentrix_id) |> 
      tibble::rownames_to_column('probe_id'),
    by=c("ID_REF" = 'probe_id'), suffix=c('','')
  ) |> 
  tibble::column_to_rownames('ID_REF')



test <- mvalues_to_table |> 
  dplyr::rename(!!! tibble::deframe(tmp)) |> 
  dplyr::rename_with( ~ paste0(.x))

stopifnot(na.omit(sum(mvalues_to_table != test) == 0))

test <- test |> 
  tibble::rowid_to_column('ID_REF')


write.table(test, file="geo/GLASS-OD_processed_M-values.txt", sep="\t", quote=T, row.names=F)


