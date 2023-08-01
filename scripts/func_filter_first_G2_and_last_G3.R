#!/usr/bin/env R


filter_first_G2_and_last_G3 <- function(metadata, nrow.check = 0) {
  out <- metadata |> 
    dplyr::group_by(patient_id) |> 
    dplyr::mutate(first_G2 = 
                    !is.na(resection_tumor_grade) & 
                    resection_tumor_grade == 2 &
                    resection_number == min(resection_number)) |>
    dplyr::mutate(last_G3 = 
                    !is.na(resection_tumor_grade) & 
                    resection_tumor_grade == 3 &
                    resection_number == max(resection_number)) |>
    dplyr::ungroup() |>
    dplyr::filter(first_G2 | last_G3)
  
  
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



