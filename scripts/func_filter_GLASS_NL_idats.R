#!/usr/bin/env/R

filter_GLASS_NL_idats <- function(metadata, nrow.check = 0) {
  out <- metadata |> 
    
    dplyr::filter(!is.na(array_sentrix_id)) |> 
    (function(.) { print(dim(.)); return(.) })() |> 
    
    assertr::verify(!is.na(array_qc.pca.detP.outlier)) |> 
    dplyr::filter(array_qc.pca.detP.outlier == F) |> 
    (function(.) { print(dim(.)); return(.) })() |> 
    
    dplyr::filter(qc.pca.pc3purity.outlier == F) |> 
    (function(.) { print(dim(.)); return(.) })() |> 
    
    assertr::verify(!duplicated(array_sentrix_id))
  
  
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

