#!/usr/bin/env R

# per entity statistical power ----


# DMP time_in_tissue_ffpe

# plot order: 
data.frame(name = c(
  "dnaMethyAge__ZhangY2017",
  "CGC[Ac]",
  "GLASS_NL_g2_g3_sig",
  "PC2",
  "time_tissue_in_ffpe",
  "QC: SPECIFICITY I GT MM 6",
  "detP significant probes (log)",
  "PC1",
  "epiTOC2_hypoSC",
  "dnaMethyAge__LuA2023p2",
  "dnaMethyAge__LuA2023p3",
  "dnaMethyAge__LuA2019",
  "dnaMethyAge__CBL_specific",
  "RepliTali",
  "epiTOC2_tnsc",
  "dnaMethyAge__HannumG2013",
  "dnaMethyAge__CBL_common",
  "dnaMethyAge__Cortex_common",
  "dnaMethyAge__McEwenL2019",
  "dnaMethyAge__ZhangQ2019",
  "dnaMethyAge__HorvathS2018",
  "dnaMethyAge__LevineM2018",
  "dnaMethyAge__LuA2023p1",
  "PC3",
  "dnaMethyAge__ShirebyG2020",
  "dnaMethyAge__PCHorvathS2018",
  "dnaMethyAge__PCHannumG2013",
  "dnaMethyAge__PCPhenoAge"
)) |> 
  dplyr::mutate(order = 1:dplyr::n()) |> 
  dplyr::mutate(colname = dplyr::recode(name, 
    `PC1` = "DMP__PCs__pp_nc__PC1_t",
    `PC2` = "DMP__PCs__pp_nc__PC2_t",
    `PC3` = "DMP__PCs__pp_nc__PC3_t",
    `detP significant probes (log)` = "DMP__pct_detP_signi__pp_nc__t",
    `CGC[Ac]` = "DMP__AcCGAP__pp_nc__t",
    `GLASS_NL_g2_g3_sig` = "DMP__GLASS_NL_signature__pp_nc__t",
    `RepliTali` = "DMP__RepliTali__up_nc__t",
    
    `QC: SPECIFICITY I GT MM 6` = "array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5",
    `epiTOC2_hypoSC` = "DMP__epiTOC2_hypoSC__up_nc__t",
    `epiTOC2_tnsc` = "DMP__epiTOC2_tnsc__up_nc__t",
    
    `dnaMethyAge__CBL_common` = "DMP__dnaMethyAge__CBL_common__up_nc__t",
    `dnaMethyAge__CBL_specific` = "DMP__dnaMethyAge__CBL_specific__up_nc__t",
    `dnaMethyAge__Cortex_common` = "DMP__dnaMethyAge__Cortex_common__up_nc__t",
    `dnaMethyAge__HannumG2013` = "DMP__dnaMethyAge__HannumG2013__up_nc__t",
    `dnaMethyAge__HorvathS2018` = "DMP__dnaMethyAge__HorvathS2018__up_nc__t",
    `dnaMethyAge__LevineM2018` = "DMP__dnaMethyAge__LevineM2018__up_nc__t",
    `dnaMethyAge__LuA2019` = "DMP__dnaMethyAge__LuA2019__up_nc__t",
    `dnaMethyAge__LuA2023p1` = "DMP__dnaMethyAge__LuA2023p1__up_nc__t",
    `dnaMethyAge__LuA2023p2` = "DMP__dnaMethyAge__LuA2023p2__up_nc__t",
    `dnaMethyAge__LuA2023p3` = "DMP__dnaMethyAge__LuA2023p3__up_nc__t",
    `dnaMethyAge__McEwenL2019` = "DMP__dnaMethyAge__McEwenL2019__up_nc__t",
    `dnaMethyAge__PCHannumG2013` = "DMP__dnaMethyAge__PCHannumG2013__up_nc__t",
    `dnaMethyAge__PCHorvathS2013` = "DMP__dnaMethyAge__PCHorvathS2013__up_nc__t",
    `dnaMethyAge__PCHorvathS2018` = "DMP__dnaMethyAge__PCHorvathS2018__up_nc__t",
    `dnaMethyAge__PCPhenoAge` = "DMP__dnaMethyAge__PCPhenoAge__up_nc__t",
    `dnaMethyAge__ShirebyG2020` = "DMP__dnaMethyAge__ShirebyG2020__up_nc__t",
    `dnaMethyAge__YangZ2016` = "DMP__dnaMethyAge__YangZ2016__up_nc__t",
    `dnaMethyAge__ZhangQ2019` = "DMP__dnaMethyAge__ZhangQ2019__up_nc__t",
    `dnaMethyAge__ZhangY2017` = "DMP__dnaMethyAge__ZhangY2017__up_nc__t",
    `time_tissue_in_ffpe` = "DMP__FFPE_decay_time__pp_nc__t"
                            ))




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
    
    DMP__epiTOC2_hypoSC__up_nc__t,
    DMP__epiTOC2_tnsc__up_nc__t,
    
    contains("dnaMe") |
      contains("SPECIFICITY")
      ) |> 
  colnames()


