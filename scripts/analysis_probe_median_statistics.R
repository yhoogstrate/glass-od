#!/usr/bin/env R

# load data ----


library(ggplot2)

source('scripts/load_constants.R')
source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')



## load clinical data ----


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}




### primary - recurrence ----


metadata.glass_od <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(180) |> 
  
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary", "recurrence"),levels=c("primary", "recurrence"))) |> 
  dplyr::select(array_sentrix_id, pr.status) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 180) 
    return(.)
  })()


metadata.glass_nl <- glass_nl.metadata.array_samples |> 
  filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(178) |> 
  
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary", "recurrence"),levels=c("primary", "recurrence"))) |> 
  dplyr::select(array_sentrix_id, pr.status) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 178) 
    return(.)
  })()



## total ----
### load ----


if(!exists('data.intensities.combined.hq_samples')) {
  source('scripts/load_intensities_hq_samples.R')
}


### GLASS-OD ----


median.combined.intensity.primary.glass_od <- data.intensities.combined.hq_samples |> 
  dplyr::select(
    metadata.glass_od |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.combined.intensity.primary = 1) |> 
  tibble::rownames_to_column('probe_id') 


median.combined.intensity.recurrent.glass_od <- data.intensities.combined.hq_samples |> 
  dplyr::select(
    metadata.glass_od |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.combined.intensity.recurrent = 1) |> 
  tibble::rownames_to_column('probe_id') 


### GLASS-NL ----


median.combined.intensity.primary.glass_nl <- data.intensities.combined.hq_samples |> 
  dplyr::select(
    metadata.glass_nl |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.combined.intensity.primary = 1) |> 
  tibble::rownames_to_column('probe_id') 


median.combined.intensity.recurrent.glass_nl <- data.intensities.combined.hq_samples |> 
  dplyr::select(
    metadata.glass_nl |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.combined.intensity.recurrent = 1) |> 
  tibble::rownames_to_column('probe_id') 


### clean-up ----



rm(data.intensities.combined.hq_samples)
gc()




## methylated ----

### load ----


if(!exists('data.intensities.methylated.hq_samples')) {
  source('scripts/load_intensities_hq_samples.R')
}


### GLASS-OD ----


median.methylated.intensity.primary.glass_od <- data.intensities.methylated.hq_samples |> 
  dplyr::select(
    metadata.glass_od |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.methylated.intensity.primary = 1) |> 
  tibble::rownames_to_column('probe_id') 


median.methylated.intensity.recurrent.glass_od <- data.intensities.methylated.hq_samples |> 
  dplyr::select(
    metadata.glass_od |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.methylated.intensity.recurrent = 1) |> 
  tibble::rownames_to_column('probe_id') 


### GLASS-NL ----


median.methylated.intensity.primary.glass_nl <- data.intensities.methylated.hq_samples |> 
  dplyr::select(
    metadata.glass_nl |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.methylated.intensity.primary = 1) |> 
  tibble::rownames_to_column('probe_id') 


median.methylated.intensity.recurrent.glass_nl <- data.intensities.methylated.hq_samples |> 
  dplyr::select(
    metadata.glass_nl |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.methylated.intensity.recurrent = 1) |> 
  tibble::rownames_to_column('probe_id') 



### clean-up ----



rm(data.intensities.methylated.hq_samples)
gc()



## unmethylated ----



### GLASS-OD ----


median.unmethylated.intensity.primary.glass_od <- data.intensities.unmethylated.hq_samples |> 
  dplyr::select(
    metadata.glass_od |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.unmethylated.intensity.primary = 1) |> 
  tibble::rownames_to_column('probe_id') 


median.unmethylated.intensity.recurrent.glass_od <- data.intensities.unmethylated.hq_samples |> 
  dplyr::select(
    metadata.glass_od |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.unmethylated.intensity.recurrent = 1) |> 
  tibble::rownames_to_column('probe_id') 




### GLASS-NL ----


median.unmethylated.intensity.primary.glass_nl <- data.intensities.unmethylated.hq_samples |> 
  dplyr::select(
    metadata.glass_nl |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.unmethylated.intensity.primary = 1) |> 
  tibble::rownames_to_column('probe_id') 


median.unmethylated.intensity.recurrent.glass_nl <- data.intensities.unmethylated.hq_samples |> 
  dplyr::select(
    metadata.glass_nl |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.unmethylated.intensity.recurrent = 1) |> 
  tibble::rownames_to_column('probe_id') 



### clean-up ----


rm(data.intensities.unmethylated.hq_samples)
gc()



## m-values ----

### load ----


if(!exists('data.intensities.unmethylated.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}



### GLASS-OD ----


median.mvalue.primary.glass_od <- data.mvalues.hq_samples |> 
  dplyr::select(
    metadata.glass_od |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.mvalue.primary = 1) |> 
  tibble::rownames_to_column('probe_id') 


median.mvalue.recurrent.glass_od <- data.mvalues.hq_samples |> 
  dplyr::select(
    metadata.glass_od |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.mvalue.recurrent = 1) |> 
  tibble::rownames_to_column('probe_id') 




### GLASS-NL ----



median.mvalue.primary.glass_nl <- data.mvalues.hq_samples |> 
  dplyr::select(
    metadata.glass_nl |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.mvalue.primary = 1) |> 
  tibble::rownames_to_column('probe_id') 


median.mvalue.recurrent.glass_nl <- data.mvalues.hq_samples |> 
  dplyr::select(
    metadata.glass_nl |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.mvalue.recurrent = 1) |> 
  tibble::rownames_to_column('probe_id') 



### clean-up ----


rm(data.mvalues.hq_samples)
gc()




## beta values ----

### load ----

if(!exists('data.beta.values.hq_samples')) {
  source('scripts/load_beta.values_hq_samples.R')
}


### GLASS-OD ----


median.beta.primary.glass_od <- data.beta.values.hq_samples |> 
  dplyr::select(
    metadata.glass_od |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.beta.primary = 1) |> 
  tibble::rownames_to_column('probe_id') 


median.beta.recurrent.glass_od <- data.beta.values.hq_samples |> 
  dplyr::select(
    metadata.glass_od |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.beta.recurrent = 1) |> 
  tibble::rownames_to_column('probe_id') 



### GLASS-NL ----


median.beta.primary.glass_nl <- data.beta.values.hq_samples |> 
  dplyr::select(
    metadata.glass_nl |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.beta.primary = 1) |> 
  tibble::rownames_to_column('probe_id') 


median.beta.recurrent.glass_nl <- data.beta.values.hq_samples |> 
  dplyr::select(
    metadata.glass_nl |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T) |> 
  data.frame() |> 
  dplyr::rename(median.beta.recurrent = 1) |> 
  tibble::rownames_to_column('probe_id') 




# stopifnot(
#   glass_nl.metadata.array_samples |> 
#     filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES) |> 
#     dplyr::filter(resection_number == 1) |> 
#     dplyr::pull(array_sentrix_id)
#   %in% colnames(data.beta.values.hq_samples))
# 
# 
# median.combined.intensity.primary.AC <- data.beta.values.hq_samples |>
#   dplyr::select(
#     glass_nl.metadata.array_samples |> 
#       filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES) |> 
#       dplyr::filter(resection_number == 1) |> 
#       dplyr::pull(array_sentrix_id)
#     ) |> 
#   as.matrix() |> 
#   MatrixGenerics::rowMedians(useNames = T) |> 
#   data.frame() |> 
#   dplyr::rename(median.combined.intensity.primary.AC = 1) |> 
#   tibble::rownames_to_column('probe_id') 



### clean-up ----


rm(data.beta.values.hq_samples)
gc()



# export ----


tmp <- median.beta.primary |> 
  dplyr::left_join(median.beta.recurrent, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(median.mvalue.primary, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(median.mvalue.recurrent, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(median.combined.intensity.primary, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(median.combined.intensity.recurrent, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(median.methylated.intensity.primary, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(median.methylated.intensity.recurrent, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(median.unmethylated.intensity.primary, by=c('probe_id'='probe_id')) |> 
  dplyr::left_join(median.unmethylated.intensity.recurrent, by=c('probe_id'='probe_id')) |> 
  
  dplyr::left_join(median.combined.intensity.primary.AC, by=c('probe_id'='probe_id')) |> 


head(tmp)





rm(median.beta.primary)
rm(median.beta.recurrent)
rm(median.mvalue.primary)
rm(median.mvalue.recurrent)
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


# 
# tmp <- readRDS("cache/analysis_probe_median_statistics.Rds") |> 
#   dplyr::left_join(median.combined.intensity.primary.AC, by=c('probe_id'='probe_id')) |> 
#   assertr::verify(!is.na(median.combined.intensity.primary.AC))
# 
# saveRDS(tmp, file="cache/analysis_probe_median_statistics.Rds")
# 



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



## signal reduction plot ----


# grep getMeth
# grep getUnmeth


meth <- data.intensities.methylated.hq_samples |> 
  head()
umeth <- data.intensities.unmethylated.hq_samples |> 
  head()


2^meth # how can these not be integers??
2^umeth # how can these not be integers??


meth |> 
  (function(.) { return(.^2) })()


### x-check methylated ----

a = readRDS("cache/intensities_m_hq/203293640061_R08C01.Rds") |> 
  dplyr::mutate(int = 2^`203293640061_R08C01`)


#head(a)
# A tibble: 6 × 2
#probe_id   `203293640061_R08C01`
#<chr>                      <dbl>
#  1 cg18478105                  7.30
#2 cg09835024                  7.92
#3 cg14361672                 12.7 
#4 cg01763666                 11.6 
#5 cg02115394                  7.98
#6 cg25813447                  6.74


# A tibble: 6 × 3
#probe_id   `203293640061_R08C01`   int
#<chr>                      <dbl> <dbl>
#  1 cg18478105                  7.30  157.
#2 cg09835024                  7.92  243.


b = readRDS("cache/intensities_m_hq/intensities_m_hq.Rds") |>
  dplyr::select( `203293640061_R08C01`) |>
  head(n=15)

# 203293640061_R08C01
# cg18478105            7.296209
# cg09835024            7.923823
# cg14361672           12.654370
# cg01763666           11.586341
# cg02115394            7.984861
# cg25813447            6.740520
# 
# 

a$`203293640061_R08C01`[1:15]
b$`203293640061_R08C01`[1:15]



readRDS("cache/intensities_m_hq/intensities_m_hq.Rds") |>
  dplyr::select( `203293640061_R08C01`) |>
  head() |> 
  (function(.) { return(.^2) })()



### x-check unmethylated ----


readRDS("cache/intensities_um_hq/intensities_um_hq.Rds") |>
  dplyr::select( `203293640061_R08C01`) |>
  head()



## bal rat ----




metadata <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES)





tot <- data.intensities.combined.hq_samples |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  (function(.) { return(2^.) })() |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::left_join(
    metadata.cg_probes.epic |> 
      dplyr::select(probe_id, probe_type_orientation, MASK_general), # probe_type, probe_type_orientation
    by=c('probe_id'='probe_id')
  ) |> 
  dplyr::filter(MASK_general == F & !is.na(probe_type_orientation)) |> 
  dplyr::mutate(MASK_general = NULL) |> 
  tibble::column_to_rownames('probe_id')
tot <- tot |> 
  dplyr::group_by(probe_type_orientation) |> 
  dplyr::summarise_all(sum) |> 
  tidyr::pivot_longer(cols = -c(probe_type_orientation))



m = data.intensities.methylated.hq_samples |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  (function(.) { return(2^.) })() |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::left_join(
    metadata.cg_probes.epic |> 
      dplyr::select(probe_id, probe_type_orientation, MASK_general), # probe_type, probe_type_orientation
    by=c('probe_id'='probe_id')
  ) |> 
  dplyr::filter(MASK_general == F & !is.na(probe_type_orientation)) |> 
  dplyr::mutate(MASK_general = NULL) |> 
  tibble::column_to_rownames('probe_id')
m <- m |> 
  dplyr::group_by(probe_type_orientation) |> 
  dplyr::summarise_all(sum) |> 
  tidyr::pivot_longer(cols = -c(probe_type_orientation))



um = data.intensities.unmethylated.hq_samples |> 
  dplyr::select(metadata$array_sentrix_id) |> 
  (function(.) { return(2^.) })() |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::left_join(
    metadata.cg_probes.epic |> 
      dplyr::select(probe_id, probe_type_orientation, MASK_general), # probe_type, probe_type_orientation
    by=c('probe_id'='probe_id')
  ) |> 
  dplyr::filter(MASK_general == F & !is.na(probe_type_orientation)) |> 
  dplyr::mutate(MASK_general = NULL) |> 
  tibble::column_to_rownames('probe_id')
um <- um |> 
  dplyr::group_by(probe_type_orientation) |> 
  dplyr::summarise_all(sum) |> 
  tidyr::pivot_longer(cols = -c(probe_type_orientation))



plt.tot <- tot |> 
  dplyr::left_join(metadata |> dplyr::select(array_sentrix_id, array_percentage.detP.signi, array_PC1), by=c('name'='array_sentrix_id'))

ggplot(plt.tot, aes(x=log(array_percentage.detP.signi), y=log(value), col=probe_type_orientation)) +
  geom_point() + 
  theme_bw()


plt.m <- m |> 
  dplyr::left_join(metadata |> dplyr::select(array_sentrix_id, array_percentage.detP.signi, array_PC1), by=c('name'='array_sentrix_id'))

ggplot(plt.m, aes(x=log(array_percentage.detP.signi), y=log(value), col=probe_type_orientation)) +
  geom_point() + 
  theme_bw()



plt.um <- um |> 
  dplyr::left_join(metadata |> dplyr::select(array_sentrix_id, array_percentage.detP.signi, array_PC1), by=c('name'='array_sentrix_id'))

ggplot(plt.um, aes(x=log(array_percentage.detP.signi), y=log(value), col=probe_type_orientation)) +
  geom_point() + 
  theme_bw()






plt <- metadata |> 
  dplyr::left_join(
    data.frame(meth = m) |> 
      tibble::rownames_to_column('array_sentrix_id'),
    by=c('array_sentrix_id'='array_sentrix_id')
  ) |> 
  dplyr::left_join(
    data.frame(unmeth = um) |> 
      tibble::rownames_to_column('array_sentrix_id'),
    by=c('array_sentrix_id'='array_sentrix_id')
  ) |> 
  dplyr::left_join(
    data.frame(total = tot) |> 
      tibble::rownames_to_column('array_sentrix_id'),
    by=c('array_sentrix_id'='array_sentrix_id')
  ) |> 
  dplyr::mutate(lr1 = log(meth / unmeth)) |> 
  dplyr::mutate(lr2 = log(log(meth) / log(unmeth)))



ggplot(plt, aes(x=array_PC1, y = log(meth))) +
  geom_point()

ggplot(plt, aes(x=array_PC1, y = log(unmeth))) +
  geom_point()

ggplot(plt, aes(x=array_PC1, y = lr1)) +
  geom_point()

ggplot(plt, aes(x=array_PC1, y = lr2)) +
  geom_point()


plt2 <- rbind(plt |> 
                dplyr::mutate(signal = meth) |> 
                dplyr::mutate(type = "meth")
              ,
              plt |> 
                dplyr::mutate(signal = unmeth) |> 
                dplyr::mutate(type = "unmeth")
              ,
              plt |> 
                dplyr::mutate(signal = total) |> 
                dplyr::mutate(type = "total"))


#ggplot(plt2, aes(x=array_PC1, y=log(signal), col=type)) +
ggplot(plt2, aes(x=log(array_percentage.detP.signi), y=log(signal), col=type)) +
  geom_point() +
  ggpubr::stat_cor(method = "pearson", aes(label = after_stat(r.label)))  +
  theme_minimal()




ggplot(plt, aes(x=array_PC1, y = log(total))) +
  geom_point()

ggplot(plt, aes(x=array_PC1, y = lr)) +
  geom_point()


ggplot(plt, aes(x=log(array_percentage.detP.signi), y = meth)) +
  geom_point()

ggplot(plt, aes(x=log(array_percentage.detP.signi), y = total)) +
  geom_point()


corrplot::corrplot(cor(  plt |> dplyr::select(meth, unmeth, total) ))

