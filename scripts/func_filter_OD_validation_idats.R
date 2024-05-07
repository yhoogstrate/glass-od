#!/usr/bin/env/R




filter_OD_validation_idats <- function(metadata, nrow.check = 0, exclude.suspected.noncodels = T) {
  out <- metadata |> 
    dplyr::filter(patient_study_name == "OD-validation") |> # oligosarcoma's incl
    (function(.) { print(dim(.)) ; return(.) })() |> 
    dplyr::filter(arraychip_version == "EPICv1") |> 
    (function(.) { print(dim(.)) ; return(.) })() |> 
    dplyr::filter(is.na(patient_reason_excluded)) |> 
    (function(.) { print(dim(.)) ; return(.) })() |> 
    dplyr::filter(is.na(resection_reason_excluded)) |> 
    (function(.) { print(dim(.)) ; return(.) })() |> 
    dplyr::filter(is.na(isolation_reason_excluded)) |> 
    (function(.) { print(dim(.)) ; return(.) })() |> 
    dplyr::filter(is.na(array_reason_excluded)) |>
    (function(.) { print(dim(.)) ; return(.) })() |> 
    assertr::verify(!is.na(array_qc.pca.detP.outlier)) |> 
    dplyr::filter(array_qc.pca.detP.outlier == F) |> 
    (function(.) { print(dim(.)) ; return(.) })() 
    
  if(exclude.suspected.noncodels == T) {
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


