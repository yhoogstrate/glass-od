#!/usr/bin/env R

load_masks_hq <- function() {
  return( readRDS("cache/masks_hq/masks_hq.Rds") |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(ncol(.) == CONST_N_SAMPLES)
      assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED)
      return(.)
    })() )
}
