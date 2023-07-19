
parse_mnp_RsGender_csv <- function(fn, prefix) {
  
  a <- read.csv(fn,header=T) |> 
    dplyr::mutate(idat = NULL, array = NULL) |>  # also returns status if it is inconfident 
    dplyr::mutate(predicted = ifelse(predicted == F, "F", predicted)) |> 
    dplyr::rename_with( ~ paste0(prefix, .x)) 
  
  return(a)
}

