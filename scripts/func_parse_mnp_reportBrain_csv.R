
parse_mnp_reportBrain_csv <- function(fn, prefix) {
  
  a <- read.csv(fn) |> 
    tibble::column_to_rownames('X') |> 
    dplyr::rename(pval = 1) |> 
    dplyr::arrange(-pval)
    #dplyr::mutate(pval = pval * 100) 
  
  top <- a |> 
    tibble::rownames_to_column('class') |> 
    dplyr::slice_head(n=1) |> 
    dplyr::pull(class)
  
  a <- a |> 
    t() |> 
    as.data.frame() |> 
    dplyr::mutate(class = top) |> 
    dplyr::rename_with( ~ paste0(prefix,"cal_", .x)) 
  
  return(a)
}

