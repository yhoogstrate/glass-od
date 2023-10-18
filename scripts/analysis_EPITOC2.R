#!/usr/bin/env R

# load m-values ----


source('scripts/load_beta.values_hq_samples.R')



# load EPITOC2 probes ----
#' https://github.com/perishky/meffonym/blob/master/inst/epitoc2/epitoc2.r

load('data/epiTOC2/dataETOC2.Rd')
source('data/epiTOC2/epiTOC2.R')



### run epitoc ----


dat <- data.beta.values.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.beta.values.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 693017)
    assertthat::assert_that(ncol(.) == 510)
    return(.)
  })()


epiTOC2.data.out <- epiTOC2(as.matrix(dat))
stopifnot(names(epiTOC2.data.out$tnsc) == names(epiTOC2.data.out$tnsc2))
stopifnot(names(epiTOC2.data.out$tnsc) == names(epiTOC2.data.out$pcgtAge))
stopifnot(names(epiTOC2.data.out$tnsc) == names(epiTOC2.data.out$hypoSC))


out <- data.frame(
  tnsc = epiTOC2.data.out$tnsc,
  tnsc2 = epiTOC2.data.out$tnsc2,
  pcgtAge = epiTOC2.data.out$pcgtAge,
  hypoSC = epiTOC2.data.out$hypoSC
) |> 
  dplyr::rename_with( ~ paste0("array_epiTOC2_", .x)) |> 
  tibble::rownames_to_column('array_sentrix_id')



saveRDS(out, file="cache/analysis_EPITOC2.Rds")



rm(epiTOC2.data.out, out)
gc()



# dnaMethyAge package ----

# devtools::install_github("yiluyucheng/dnaMethyAge")

dnaMethyAge::availableClock()
library(dnaMethyAge)
data(list='PC-clocks', envir=environment())

age_HannumG2013    <- dnaMethyAge::methyAge(dat, clock='HannumG2013')    |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__HannumG2013 = mAge)
age_HorvathS2013   <- dnaMethyAge::methyAge(dat, clock='HorvathS2013')   |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__HorvathS2013 = mAge)
age_LevineM2018    <- dnaMethyAge::methyAge(dat, clock='LevineM2018')    |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__LevineM2018 = mAge)
age_ZhangQ2019     <- dnaMethyAge::methyAge(dat, clock='ZhangQ2019')     |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__ZhangQ2019 = mAge)
age_ShirebyG2020   <- dnaMethyAge::methyAge(dat, clock='ShirebyG2020')   |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__ShirebyG2020 = mAge)
age_ZhangY2017     <- dnaMethyAge::methyAge(dat, clock='ZhangY2017')     |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__ZhangY2017 = mAge)
age_LuA2019        <- dnaMethyAge::methyAge(dat, clock='LuA2019')        |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__LuA2019 = mAge)
age_HorvathS2018   <- dnaMethyAge::methyAge(dat, clock='HorvathS2018')   |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__HorvathS2018 = mAge)
age_McEwenL2019    <- dnaMethyAge::methyAge(dat, clock='McEwenL2019')    |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__McEwenL2019 = mAge)
age_CBL_specific   <- dnaMethyAge::methyAge(dat, clock='CBL_specific')   |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__CBL_specific = mAge)
age_PCHorvathS2013 <- dnaMethyAge::methyAge(dat, clock='PCHorvathS2013') |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__PCHorvathS2013  = mAge)
age_PCHannumG2013  <- dnaMethyAge::methyAge(dat, clock='PCHannumG2013')  |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__PCHannumG2013 = mAge)
age_PCHorvathS2018 <- dnaMethyAge::methyAge(dat, clock='PCHorvathS2018') |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__PCHorvathS2018 = mAge)
age_PCPhenoAge     <- dnaMethyAge::methyAge(dat, clock='PCPhenoAge')     |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__PCPhenoAge = mAge)
age_CBL_common     <- dnaMethyAge::methyAge(dat, clock='CBL_common')     |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__CBL_common = mAge)
age_Cortex_common  <- dnaMethyAge::methyAge(dat, clock='Cortex_common')  |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__Cortex_common = mAge)
age_epiTOC2        <- dnaMethyAge::methyAge(dat, clock='epiTOC2')        |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__epiTOC2 = mAge)
age_LuA2023p1      <- dnaMethyAge::methyAge(dat, clock='LuA2023p1')      |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__LuA2023p1 = mAge)
age_LuA2023p2      <- dnaMethyAge::methyAge(dat, clock='LuA2023p2')      |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__LuA2023p2 = mAge)
age_LuA2023p3      <- dnaMethyAge::methyAge(dat, clock='LuA2023p3')      |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__LuA2023p3 = mAge)
age_YangZ2016      <- dnaMethyAge::methyAge(dat, clock='YangZ2016')      |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__YangZ2016 = mAge)

# The following fail/errr out
#age_DunedinPACE    <- dnaMethyAge::methyAge(dat, clock='DunedinPACE') |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__DunedinPACE = mAge)
#age_PCGrimAge      <- dnaMethyAge::methyAge(dat, clock='PCGrimAge') |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__PCGrimAge = mAge)
#age_BernabeuE2023c <- dnaMethyAge::methyAge(dat, clock='BernabeuE2023c') |> assertr::verify(is.numeric(mAge) & !is.na(mAge)) |> dplyr::rename(dnaMethyAge__BernabeuE2023c  = mAge)



out <- age_HannumG2013 |> 
  dplyr::left_join(age_HorvathS2013   , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_LevineM2018    , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_ZhangQ2019     , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_ShirebyG2020   , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_YangZ2016      , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_ZhangY2017     , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_LuA2019        , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_HorvathS2018   , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_McEwenL2019    , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_CBL_specific   , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_PCHorvathS2013 , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_PCHannumG2013  , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_PCHorvathS2018 , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_PCPhenoAge     , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_CBL_common     , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_Cortex_common  , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_epiTOC2        , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_LuA2023p1      , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_LuA2023p2      , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  dplyr::left_join(age_LuA2023p3      , by=c('Sample'='Sample'), suffix=c('',''))  |> 
  tibble::column_to_rownames('Sample')

#dplyr::left_join(age_DunedinPACE      , by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))  |> 
#dplyr::left_join(age_PCGrimAge        , by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))  |> 
#dplyr::left_join(age_BernabeuE2023c   , by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))  |> 

corrplot::corrplot(cor(out, method = "spearman"))


# LuA2019 seems inversed, ZhangY2017 seems different signal
saveRDS(out, file="cache/analysis_dnaMethyAge.Rds")


rm(out)




