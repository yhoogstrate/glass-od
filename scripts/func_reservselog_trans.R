#!/usr/bin/env R

# https://stackoverflow.com/questions/11053899/how-to-get-a-reversed-log10-scale-in-ggplot2

#library(scales)

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv, 
            scales::log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

