#!/usr/bin/env R 


# all hq ----


data.mvalues.hq_samples <- readRDS("cache/mvalues.HQ_samples.Rds") 
data.mvalues.mask.hq_samples <- readRDS("cache/mvalues.HQ_samples.detP_mask.Rds")

stopifnot(colnames(data.mvalues.hq_samples) == colnames(data.mvalues.mask.hq_samples))
stopifnot(rownames(data.mvalues.hq_samples) == rownames(data.mvalues.mask.hq_samples))


