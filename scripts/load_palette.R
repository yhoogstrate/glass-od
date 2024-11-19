#!/usr/bin/env R



if(!exists('mixcol')) {
  source('scripts/load_functions.R')
}



# center with gray - original also used in corrplot and some HI-publications
col5 <- grDevices::colorRampPalette(c(       mixcol("#4DCB12", mixcol("#FFFFFF","black",0.125),0.10),
                                             mixcol("#FFFFFF","black",0.125),
                                             mixcol("#DF1BA7", "white",0.025)
))



# center with darker gray - original also used in corrplot and some HI-publications
col4 <- grDevices::colorRampPalette(c(       "#67001F", 
                                             mixcol("#B2182B","black",0.025*2),
                                             mixcol("#D6604D","black",0.050*2),
                                             mixcol("#F4A582","black",0.075*2),
                                             mixcol("#FDDBC7","black",0.100*2),
                                             
                                             mixcol("#FFFFFF","black",0.125*2),
                                             
                                             mixcol("#D1E5F0","black",0.100*2),
                                             mixcol("#92C5DE","black",0.075*2),
                                             mixcol("#4393C3","black",0.050*2),
                                             mixcol("#2166AC","black",0.025*2),
                                             
                                             "#053061"
))




# center with gray - original also used in corrplot and some HI-publications
col3 <- grDevices::colorRampPalette(c(       "#67001F", 
                                             mixcol("#B2182B","black",0.025),
                                             mixcol("#D6604D","black",0.050),
                                             mixcol("#F4A582","black",0.075),
                                             mixcol("#FDDBC7","black",0.100),
                                             
                                             mixcol("#FFFFFF","black",0.125),
                                             
                                             mixcol("#D1E5F0","black",0.100),
                                             mixcol("#92C5DE","black",0.075),
                                             mixcol("#4393C3","black",0.050),
                                             mixcol("#2166AC","black",0.025),
                                             
                                             "#053061"
))


# quite a bit brighter, neat for heatmap
col6 <- grDevices::colorRampPalette(c(       "#67001F", 
                                             mixcol("#B2182B","black",0.025 * 0.5),
                                             mixcol("#D6604D","black",0.050 * 0.5),
                                             mixcol("#F4A582","black",0.075 * 0.5),
                                             mixcol("#FDDBC7","black",0.100 * 0.5),
                                             
                                             mixcol("#FFFFFF","black",0.125 * 0.5),
                                             
                                             mixcol("#D1E5F0","black",0.100 * 0.5),
                                             mixcol("#92C5DE","black",0.075 * 0.5),
                                             mixcol("#4393C3","black",0.050 * 0.5),
                                             mixcol("#2166AC","black",0.025 * 0.5),
                                             
                                             "#053061"
))



# center with white - nicest if applicable
col2 <- grDevices::colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                      "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                      "#4393C3", "#2166AC", "#053061"))



palette_mnp_12.8_6 <- c(
  'A_IDH_LG' = mixcol('#f85eb6', 'gray80',0.15),
  'A_IDH_HG' = '#d34394ff',
  'O_IDH' = 'aquamarine3',
  'OLIGOSARC_IDH' = 'aquamarine4',
  'CTRL_CORPCAL' = '#7f8b45',
  'GBM_RTK_I' = '#51458b'
)



palette_mnp_12.8_7 <- c(
  'A_IDH_HG' = '#d34394ff',
  'O_IDH' = 'aquamarine3',
  'OLIGOSARC_IDH' = 'aquamarine4',
  'other' = 'gray80'
)



palette_mnp_12.8_AC_3 <- c(
  'A_IDH_LG' = mixcol('#f85eb6', 'gray80',0.15),
  'A_IDH_HG' = '#d34394ff',
  'other' = '#7f8b45'
  #'GBM_RTK_I' = '#51458b'
)


# palette_g2_g3 <- c(
#   'Grade 2' = mixcol(col2(10)[9], "white", 0.15),
#   'Grade 3' = mixcol(col2(10)[2], "white", 0.15)
# )
palette_g2_g3 <- c(
  'Grade 2' = mixcol( 'lightblue', 'lightgreen'),
  'Grade 3' = mixcol( 'darkblue', 'darkgreen')
)


# palette_p_r <- c(
#   'Primary' = mixcol(col2(10)[9], "white", 0.15),
#   'Recurrent' = mixcol(col2(10)[2], "white", 0.15)
# )

palette_p_r <- c(
  'Primary' = '#d9e39e', # #e3b29e
  'Recurrent' = '#3b4500'
)

palette_yes_no_1 <- c('Yes' = '#E49E27', 'No' = '#59B2E6')
palette_yes_no_2 <- c('Yes' = '#009E74', 'No' = '#CB75A4')


palette_infinium_signals <- c('methylated'='darkgreen', 'unmethylated'='red')


