#!/usr/bin/env/R

filter_GSAM_idats <- function(metadata, nrow.check = 0) {
  out <- metadata |> 
    assertr::verify(!is.na(array_sentrix_id)) |> 
    assertr::verify(!is.na(array_qc.pca.detP.outlier)) |> 
    dplyr::filter(array_qc.pca.detP.outlier == F) |> 
    dplyr::filter(IDH == F) |> 
    
    assertr::verify(!duplicated(Sample_Name)) |> 
    assertr::verify(array_sentrix_id != "206467010123_R02C01") |> 
    assertr::verify(array_sentrix_id != "206467010123_R03C01")
  
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

