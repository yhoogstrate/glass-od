#!/usr/bin/env/R

filter_GLASS_OD_idats <- function(metadata, nrow.check = 0, exclude.suspected.noncodels = T) {
  out <- metadata |> 
    dplyr::filter(is.na(patient_reason_excluded)) |> 
    dplyr::filter(is.na(resection_reason_excluded)) |> 
    dplyr::filter(is.na(isolation_reason_excluded)) |> 
    dplyr::filter(is.na(array_reason_excluded)) |> # formally there are a few replicates that are o.k. but of which the other is better
    dplyr::filter(patient_study_name == "GLASS-OD") |> # oligosarcoma's from CATNON excl
    assertr::verify(!is.na(array_qc.pca.detP.outlier)) |> 
    dplyr::filter(array_qc.pca.detP.outlier == F) |> 
    #assertr::verify(!duplicated(resection_id)) |> 
    assertr::verify(array_sentrix_id != "204808700074_R04C01")
  
  if(exclude.suspected.noncodels == T) {
    print(dim(out))
    
    out <- out |> 
      dplyr::filter(patient_suspected_noncodel == F)
    
    print(dim(out))
  }
  
  
  if(nrow.check > 0) {
    out <- out |> 
        (function(.) {
          print(dim(.))
          assertthat::assert_that(nrow(.) == nrow.check)
          return(.)
        })()
  }
  
  return (out)
}


