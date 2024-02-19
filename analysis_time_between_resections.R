#!/usr/bin/env R

if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# analysis ----

#' Patient 77 had some inconsistent survival or follow-up data, not allowing to deduce everything. The 
#' Block nrs however indicate that one operation was in 2003 and the other in 2009.


tmp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::filter(patient_id != "0077")


dd <- c()
for(patient in unique(tmp$patient_id)) {
  print(patient)
  
  res <- tmp |> 
    dplyr::filter(patient_id == patient) |> 
    dplyr::select(array_sentrix_id, isolation_id, resection_number, resection_date) |> 
    dplyr::arrange(resection_number)
  
  if(nrow(res) > 1) {
    dates1 <- res$resection_date[1:(nrow(res)-1)]
    dates2 <- res$resection_date[2:nrow(res)]
  }
  
  delta.dates <- (dates2 - dates1)
  dd <- c(dd, delta.dates)
  for(d in delta.dates) {
    if(d < (365/2)) {
      print("******")
    }
    print(d)
  }
  print("")
}


