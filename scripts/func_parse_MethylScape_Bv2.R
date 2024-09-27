#!/usr/bin/env R

#fn = "data/GLASS_OD/DNA Methylation - EPIC arrays - MethylScape Bethesda classifier/v2.0/0017-R3_Report-CNS-Bv2_203989100149_R05C01.html"

parse_MethylScape_Bv2 <- function(fn) {
  if(length(fn) > 1) {
    out <-  do.call(rbind, pbapply::pblapply(fn, parse_MethylScape_Bv2))
  } else {
    
    dat <- xml2::read_html(fn) |> 
      rvest::html_elements("div #bethesda-classifier-v2")
    
    sentrix_id <- dat |> 
      rvest::html_element("table") |> 
      rvest::html_element("thead") |> 
      rvest::html_elements("th") |> 
      purrr::pluck(2) |> 
      rvest::html_text2()
    
    class <- dat |> 
      rvest::html_elements("table") |> 
      purrr::pluck(2) |> 
      rvest::html_element("tbody") |> 
      rvest::html_element("tr") |> 
      rvest::html_elements("td")  |> 
      purrr::pluck(2) |> 
      rvest::html_text2()
    
    class_score <- dat |> 
      rvest::html_element("table") |> 
      rvest::html_element("tbody") |> 
      rvest::html_elements("tr") |> 
      purrr::pluck(5) |> 
      rvest::html_elements("td") |> 
      purrr::pluck(2) |> 
      rvest::html_text2() |> 
      as.numeric()
    
    out <- data.frame(sentrix_id = sentrix_id,
                      class = class,
                      class_score = class_score)
  }
  
  return(out)
}


# html_structure()


