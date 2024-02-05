#!/usr/bin/env/R




filter_OD_validation_idats <- function(metadata, nrow.check = 0, exclude.suspected.noncodels = T) {
  out <- metadata |> 
    dplyr::filter(is.na(patient_reason_excluded)) |> 
    dplyr::filter(is.na(resection_reason_excluded)) |> 
    dplyr::filter(is.na(isolation_reason_excluded)) |> 
    dplyr::filter(is.na(array_reason_excluded)) |> # formally there are a few replicates that are o.k. but of which the other is better
    dplyr::filter(patient_study_name == "OD-validation") |> # oligosarcoma's incl
    dplyr::filter(arraychip_version == "EPICv1") |> 
    assertr::verify(!is.na(array_qc.pca.detP.outlier)) |> 
    dplyr::filter(array_qc.pca.detP.outlier == F)
  
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


