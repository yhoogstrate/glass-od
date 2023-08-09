#!/usr/bin/env/R

filter_GLASS_NL_idats <- function(metadata, nrow.check = 0) {
  out <- metadata |> 
    dplyr::filter(!is.na(sentrix_id)) |> 
    assertr::verify(!is.na(qc.pca.detP.outlier)) |> 
    dplyr::filter(qc.pca.detP.outlier == F) |> 
    #dplyr::filter(qc.pca.pc3purity.outlier == F) |> 
    assertr::verify(!duplicated(sentrix_id))
  
  
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


