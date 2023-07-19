
parse_mnp_predictMGMT_csv <- function(fn, prefix) {
  
  a <- read.csv(fn,header=T) |> 
    dplyr::mutate(X=NULL, Status = NULL) |>  # also returns status if it is inconfident 
    dplyr::mutate(status = dplyr::case_when(
      Estimated < Cutoff & CI_Lower < Cutoff & CI_Upper < Cutoff ~ "unmethylated",
      Estimated > Cutoff & CI_Lower > Cutoff & CI_Upper > Cutoff ~ "methylated",
      T ~ as.character(NA)
    )) |> 
    dplyr::rename_with( ~ paste0(prefix, .x)) 
  
  return(a)

}

