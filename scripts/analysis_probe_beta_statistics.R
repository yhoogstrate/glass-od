#!/usr/bin/env R

# load ----


library(ggplot2)

source('scripts/load_constants.R')
source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}



# if(!exists('data.mvalues.probes')) {
#   source('scripts/load_mvalues_hq_samples.R')
# }


if(!exists('data.intensities.probes')) {
  source('scripts/load_intensities_hq_samples.R')
}


# median beta value ----
# median beta value at primary ----
# median beta value at recurrence ----



# median intensity value ----
# median intensity at primary ----
# median intensity at recurrence ----



# per probe per channel median/mean intensity, for primary and recurrence separate ----



metadata <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  filter_primaries_and_last_recurrences(178) |> 
  
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  
  #dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"remainder")))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 
  dplyr::select(array_sentrix_id, pr.status) |> 
  
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (178)) 
    return(.)
  })()


methylated.primary.medians <- data.intensities.methylated.hq_samples |> 
  dplyr::select(
    metadata |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T)


methylated.recurrent.medians <- data.intensities.methylated.hq_samples |> 
  dplyr::select(
    metadata |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T)





unmethylated.primary.medians <- data.intensities.unmethylated.hq_samples |> 
  dplyr::select(
    metadata |> dplyr::filter(pr.status == "primary") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T)


unmethylated.recurrent.medians <- data.intensities.unmethylated.hq_samples |> 
  dplyr::select(
    metadata |> dplyr::filter(pr.status == "recurrence") |>  dplyr::pull(array_sentrix_id)
  ) |> 
  as.matrix() |> 
  MatrixGenerics::rowMedians(useNames = T)



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


