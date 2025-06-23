#!/usr/bin/env R

parse_mnp_CNVP_chr4_mean <- function(fn) {
    if(length(fn) > 1) {
      out <-  do.call(rbind, pbapply::pblapply(fn, parse_mnp_CNVP_chr4_mean))
    } else {
        
      dat <- read.table(fn, header=T) |> 
      dplyr::filter(chrom == "chr4") |> 
      dplyr::rename(array_sentrix_id = ID) |> 
      dplyr::summarise(array_sentrix_id = unique(array_sentrix_id), array_mnp_CNVP_v12.8_v5.2_CNVP_segment_chr4_weighted_mean = weighted.mean(seg.mean, num.mark))
  
      return(dat)
    }
}


test1 <- parse_mnp_CNVP_chr4_mean('data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/206467110176_R06C01__sample from 0099-R2__193135____brain_classifier_v12.8_sample_report__v1.1__131____run-298392/cnvp_v5.2/206467110176_R06C01.segments.seg')
test2 <- parse_mnp_CNVP_chr4_mean("data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/206137490041_R04C01__sample from P-38 [tmp-id]__174928____brain_classifier_v12.8_sample_report__v1.1__131____run-298742/cnvp_v5.2/206137490041_R04C01.segments.seg")

stopifnot(test1$array_mnp_CNVP_v12.8_v5.2_CNVP_segment_chr4_weighted_mean == 0.008)
stopifnot(round(test2$array_mnp_CNVP_v12.8_v5.2_CNVP_segment_chr4_weighted_mean, 8) == 0.04048089) # rounding is needed


test3 <- parse_mnp_CNVP_chr4_mean(c(
  'data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/206467110176_R06C01__sample from 0099-R2__193135____brain_classifier_v12.8_sample_report__v1.1__131____run-298392/cnvp_v5.2/206467110176_R06C01.segments.seg',
  "data/GLASS_OD/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/206137490041_R04C01__sample from P-38 [tmp-id]__174928____brain_classifier_v12.8_sample_report__v1.1__131____run-298742/cnvp_v5.2/206137490041_R04C01.segments.seg"
))

stopifnot(round(test3$array_mnp_CNVP_v12.8_v5.2_CNVP_segment_chr4_weighted_mean, 8) == c(0.008, 0.04048089))


rm(test1, test2, test3)


