#!/usr/bin/env R


# center with zero - nicest if applicable
col2 <- grDevices::colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                      "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                      "#4393C3", "#2166AC", "#053061"))


# center with gray
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

