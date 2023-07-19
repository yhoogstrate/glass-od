#!/usr/bin/env R


query_mnp_12.8_CNVP_segment_csv <- function(path, n, sentrix_check = NULL) {
  tmp <- list.files(path = path, pattern = "*.seg", recursive = TRUE) |> 
    data.frame(heidelberg_cnvp_segments = _) |> 
    dplyr::mutate(heidelberg_cnvp_segments = paste0(path, heidelberg_cnvp_segments)) |> 
    assertr::verify(file.exists(heidelberg_cnvp_segments)) |> 
    dplyr::mutate(heidelberg_cnvp_version = gsub("^.+/cnvp_([^/]+)/.+$","\\1", heidelberg_cnvp_segments)) |> 
    dplyr::mutate(sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", heidelberg_cnvp_segments)) |> 
    assertr::verify(!is.na(sentrix_id)) |> 
    assertr::verify(!duplicated(sentrix_id)) |> # only one version per sample needed
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == n)
      return(.)
    })() 

  if(!is.null(sentrix_check)) {
    stopifnot(tmp$sentrix_id %in% sentrix_check)
      
  }
  
  return(tmp)
}
