#!/usr/bin/env R

split_match <- function(cc, gid) {
  params <- unlist(stringr::str_split(cc, ";"))
  return(sum(gid %in% params) > 0)
}

