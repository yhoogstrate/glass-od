install.packages('corrplot')
install.packages('stringr')
install.packages('dplyr')
install.packages('readxl')
install.packages('assertthat')
install.packages('ggplot2')
install.packages('tidyr')
install.packages('tibble')
install.packages('pathwork')
install.packages('extrafont')
install.packages("ggpubr")
install.packages('pbapply')
install.packages("rstudioapi")
install.packages('gghalves')









if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("minfi")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

