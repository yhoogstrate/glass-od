#!/usr/bin/env R


query_mnp_12.8_CNVP_segment_csv <- function(path, n, sentrix_check = NULL, prefix="") {
  tmp <- list.files(path = path, pattern = "*.seg", recursive = TRUE) |> 
    data.frame(CNVP_segments = _) |> 
    dplyr::mutate(CNVP_segments = paste0(path, CNVP_segments)) |> 
    assertr::verify(file.exists(CNVP_segments)) |> 
    dplyr::mutate(CNVP_version = gsub("^.+/cnvp_([^/]+)/.+$","\\1", CNVP_segments)) |> 
    dplyr::mutate(array_sentrix_id = gsub("^.+([0-9]{12}_[A-Z][0-9]+[A-Z][0-9]+).+$","\\1", CNVP_segments)) |> 
    assertr::verify(!is.na(array_sentrix_id)) |> 
    assertr::verify(!duplicated(array_sentrix_id)) |> # only one version per sample needed
    dplyr::rename_with( ~ paste0(prefix, .x), .cols=!matches("^array_sentrix_id$",perl = T)) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == n)
      return(.)
    })()

  if(!is.null(sentrix_check)) {
    stopifnot(tmp$array_sentrix_id %in% sentrix_check)
      
  }
  
  return(tmp)
}
