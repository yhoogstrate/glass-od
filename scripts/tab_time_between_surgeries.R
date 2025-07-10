#!/usr/bin/env R


tmp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::select(patient_id, resection_id, resection_date, resection_number) |> 
  dplyr::arrange(resection_id) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(pat_count = dplyr::n()) |> 
  dplyr::ungroup() |> 
  dplyr::filter(pat_count > 1)


a = tmp |> 
  dplyr::mutate(date = as.Date(resection_date)) |> # ensure 'date' is in Date format
  dplyr::group_by(patient_id) |> 
  dplyr::summarise(timespan_days = as.numeric(max(date) - min(date))) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(timespan_months = timespan_days / CONST_DAYS_PER_MONTH) 

a |> View()

a$timespan_days |>  median() / CONST_DAYS_PER_MONTH

