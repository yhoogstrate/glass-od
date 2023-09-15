#!/usr/bin/env


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


#source('scripts/load_chrom_sizes.R')
source('scripts/load_themes.R')



library(ggplot2)
library(patchwork)



# obtain data ----



