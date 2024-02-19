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




if(!exists('data.intensities.probes')) {
  source('scripts/load_intensities_hq_samples.R')
}




if(!exists('data.beta.values.hq_samples')) {
  source('scripts/load_beta.values_hq_samples.R')
}



# primary - recurrence ----


metadata <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(177) |> 
  
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  #dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"remainder")))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 
  dplyr::select(array_sentrix_id, pr.status) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (177)) 
    return(.)
  })()



## total ----


median.combined.intensity.primary <- data.intensities.combined.hq_samples |> 
  dplyr::select(
    metadata |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.combined.intensity.primary = 1) |> 
  tibble::rownames_to_column('probe_id') 


median.combined.intensity.recurrent <- data.intensities.combined.hq_samples |> 
  dplyr::select(
    metadata |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.combined.intensity.recurrent = 1) |> 
  tibble::rownames_to_column('probe_id') 




## methylated ----


median.methylated.intensity.primary <- data.intensities.methylated.hq_samples |> 
  dplyr::select(
    metadata |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.methylated.intensity.primary = 1) |> 
  tibble::rownames_to_column('probe_id') 


median.methylated.intensity.recurrent <- data.intensities.methylated.hq_samples |> 
  dplyr::select(
    metadata |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.methylated.intensity.recurrent = 1) |> 
  tibble::rownames_to_column('probe_id') 



## unmethylated ----


median.unmethylated.intensity.primary <- data.intensities.unmethylated.hq_samples |> 
  dplyr::select(
    metadata |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.unmethylated.intensity.primary = 1) |> 
  tibble::rownames_to_column('probe_id') 


median.unmethylated.intensity.recurrent <- data.intensities.unmethylated.hq_samples |> 
  dplyr::select(
    metadata |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.unmethylated.intensity.recurrent = 1) |> 
  tibble::rownames_to_column('probe_id') 



## beta ----


median.beta.primary <- data.beta.values.hq_samples |> 
  dplyr::select(
    metadata |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.beta.primary = 1) |> 
  tibble::rownames_to_column('probe_id') 


median.beta.recurrent <- data.beta.values.hq_samples |> 
  dplyr::select(
    metadata |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.beta.recurrent = 1) |> 
  tibble::rownames_to_column('probe_id') 



## export ----


tmp <- median.beta.primary |> 
  dplyr::left_join(median.beta.recurrent, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(median.combined.intensity.primary, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(median.combined.intensity.recurrent, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(median.methylated.intensity.primary, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(median.methylated.intensity.recurrent, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(median.unmethylated.intensity.primary, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(median.unmethylated.intensity.recurrent, by=c('probe_id'='probe_id'))


head(tmp)





rm(median.beta.primary)
rm(median.beta.recurrent)
rm(median.combined.intensity.primary)
rm(median.combined.intensity.recurrent)
rm(median.methylated.intensity.primary)
rm(median.methylated.intensity.recurrent)
rm(median.unmethylated.intensity.primary)
rm(median.unmethylated.intensity.recurrent)


saveRDS(tmp, file="cache/analysis_probe_median_statistics.Rds")


rm(data.intensities.combined.hq_samples)
rm(data.intensities.methylated.hq_samples)
rm(data.intensities.unmethylated.hq_samples)
rm(data.beta.values.hq_samples)


gc()




# sandbox ----



plt <- rbind(
  data.frame(median = methylated.primary.medians) |> 
    dplyr::mutate(signal = "methylated") |> 
    dplyr::mutate(resection = "primary") |> 
    dplyr::arrange(median) |> 
    dplyr::mutate(x = 1:dplyr::n()) 
  ,
  data.frame(median = methylated.recurrent.medians) |> 
    dplyr::mutate(signal = "methylated") |> 
    dplyr::mutate(resection = "recurrent") |> 
    dplyr::arrange(median) |> 
    dplyr::mutate(x = 1:dplyr::n())
  ,
  data.frame(median = unmethylated.primary.medians) |> 
    dplyr::mutate(signal = "unmethylated") |> 
    dplyr::mutate(resection = "primary") |> 
    dplyr::arrange(median) |> 
    
    dplyr::mutate(x = 1:dplyr::n()) 
    #head(n=10000)
  ,
  data.frame(median = unmethylated.recurrent.medians) |> 
    dplyr::mutate(signal = "unmethylated") |> 
    dplyr::mutate(resection = "recurrent") |> 
    dplyr::arrange(median) |> 
    
    dplyr::mutate(x = 1:dplyr::n())
    #head(n=10000)
)


ggplot(plt, aes(x=x,y=median, col=signal, group=paste0(resection,": ",signal), lty = resection)) +
  #facet_wrap(~resection) +
  geom_line()


