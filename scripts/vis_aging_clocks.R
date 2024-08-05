#!/usr/bin/env R


# DMP time_in_tissue_ffpe

data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  tibble::tibble() |> 
  dplyr::select(
    DMP__PCs__pp_nc__PC1_t,
    DMP__PCs__pp_nc__PC2_t,
    DMP__PCs__pp_nc__PC3_t,
    
    DMP__pct_detP_signi__pp_nc__t,
    DMP__AcCGAP__pp_nc__t,
    DMP__GLASS_NL_signature__pp_nc__t,
    
    DMP__RepliTali__up_nc__t,
    
    DMP__FFPE_decay_time__pp_nc__t,

    contains("dnaMe")
      ) |> 
  colnames()


