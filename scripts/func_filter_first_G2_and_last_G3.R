#!/usr/bin/env R


filter_first_G2_and_last_G3 <- function(metadata, nrow.check = 0) {
  out <- metadata
  
  if('resection_tumor_grade' %in% colnames(metadata)) {
    out <- out |> 
      dplyr::mutate(LG_HG_status = dplyr::case_when(
        resection_tumor_grade %in% c(2) ~ "LG",
        resection_tumor_grade %in% c(3) ~ "HG",
        T ~ as.character(NA),
        ))
  } else if('WHO_Classification2021' %in% colnames(metadata)) {
    out <- out |> 
      dplyr::mutate(LG_HG_status = dplyr::case_when(
        WHO_Classification2021 %in% c("Astrocytoma, IDH-mutant, WHO grade 2", "Astrocytoma, IDH-mutant, WHO grade 3") ~ "LG",
        WHO_Classification2021 %in% c("Astrocytoma, IDH-mutant, WHO grade 4") ~ "HG",
        T ~ as.character(NA),
      ))
  }
  
  out <- out |> 
    dplyr::group_by(patient_id, LG_HG_status) |> 
    dplyr::mutate(first_LG = 
                    !is.na(LG_HG_status) & 
                    LG_HG_status == "LG" &
                    resection_number == min(resection_number)) |>
    dplyr::mutate(last_HG = 
                    !is.na(LG_HG_status) & 
                    LG_HG_status == "HG" &
                    resection_number == max(resection_number)) |>
    dplyr::ungroup() |>
    dplyr::filter(first_LG | last_HG)
  
  
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



