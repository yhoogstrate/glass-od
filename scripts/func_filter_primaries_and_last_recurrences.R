#!/usr/bin/env R 



filter_primaries_and_last_recurrences <- function(metadata, nrow.check = 0) {
  out <- metadata |> 
    dplyr::group_by(patient_id) |> 
    dplyr::mutate(last_recurrence = resection_number != 1 & resection_number == max(resection_number)) |>
    dplyr::ungroup() |>
    dplyr::filter(resection_number == 1 | last_recurrence)
  
  
  if(nrow.check > 0) {
    out <- out |> 
      (function(.) {
        print(dim(.))
        assertthat::assert_that(nrow(.) == nrow.check)
        return(.)
      })()
  }
  
  
  return(out)
}



