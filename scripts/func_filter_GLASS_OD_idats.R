#!/usr/bin/env/R

filter_GLASS_OD_idats <- function(metadata, nrow.check = 0) {
  out <- metadata |> 
    dplyr::filter(is.na(reason_excluded_patient)) |> 
    dplyr::filter(is.na(reason_excluded_resection)) |> 
    dplyr::filter(is.na(reason_excluded_resection_isolation)) |> 
    dplyr::filter(is.na(reason_excluded_array_sample)) |> # formally there are a few replicates that are o.k. but of which the other is better
    
    dplyr::filter(study_name == "GLASS-OD") |> # oligosarcoma's from CATNON excl
    
    assertr::verify(!is.na(qc.pca.outlier)) |> 
    dplyr::filter(qc.pca.outlier == F) |> 
    
    assertr::verify(!duplicated(resection_id)) |> 
    
    assertr::verify(sentrix_id != "204808700074_R04C01")
  
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

