#!/usr/bin/env R

parse_sentrix_id <- function(idat) {
  tmp <- illuminaio::readIDAT(idat)
  
  return(paste0(tmp$Barcode, "_", tmp$Unknowns$MostlyA))
}

