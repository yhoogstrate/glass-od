#!/usr/bin/env R 

# load ----


source('scripts/load_functions.R')


# all hq ----


data.mvalues.hq_samples <- readRDS("cache/mvalues.HQ_samples.Rds") 
data.mvalues.mask.hq_samples <- readRDS("cache/mvalues.HQ_samples.detP_mask.Rds")

stopifnot(colnames(data.mvalues.hq_samples) == colnames(data.mvalues.mask.hq_samples))
stopifnot(rownames(data.mvalues.hq_samples) == rownames(data.mvalues.mask.hq_samples))


## probe table - mask / detP counts ----

data.mvalues.probes <- data.mvalues.mask.hq_samples |> 
  is.na() |> 
  rowSums() |> 
  as.data.frame() |> 
  dplyr::rename(n_na = 1) |> 
  dplyr::mutate(good_probe = n_na == 0) |> 
  tibble::rownames_to_column('probe_id')

## probe table - probe annotation ----


tmp <- read.delim("~/mnt/neuro-genomic-1-ro/catnon/Methylation - EPIC arrays/EPIC.hg38.manifest.tsv.gz") |> 
  dplyr::rename(probe_id = probeID) |> 
  assertr::verify(data.mvalues.probes$probe_id %in% probe_id) |> 
  assertr::verify(!duplicated(probe_id))


data.mvalues.probes <- data.mvalues.probes |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))

rm(tmp)

## AC embryionic development genes ghisai ----



tmp <- c("MFAP2","OPRD1","CD101","MIR3681HG","AC064875.1","OSR1","LINC01121","AC007402.1","AC007402.2","TLX2","LINC01956","EN1","GRB14","HOXD11","HOXD10","HOXD-AS2","HOXD9","HOXD8","HOXD3","HOXD4",
"HAGLROS","IGFBP2","PAX3","TRIM71","AC112487.1","GPR156","SHOX2","IGF2BP2","LNX1","NMU","NKX6-1","HAND2","HAND2-AS1","ASB5","TENM3-AS1","NPR3","OTP","SPRY4","COL9A1","EYA4",
"AL139393.1","AL139393.3","TWIST1","IGF2BP3","HOXA2","HOXA3","HOXA5","HOXA6","HOXA7","HOXA9","HOXA10","HOXA11","HOXA11-AS","MNX1","NRG1","PLAT","GDF6","OSR2","PAX5","BX255923.2",
"AL353150.1","VAX1","EBF3","IGF2","IGF2-AS","AC090692.1","LINC02367","AC084816.1","AC009318.4","LRRK2","HOXC9","HOXC8","HOXC4","LINC01234","TBX5","GSX1","PAX9","LINC00648","C14orf39","SIX6",
"AC087473.1","ISL2","HAPLN3","LINC01579","GPR139","CRNDE","IRX5","AC126773.4","GJC1","TBX21","AC015909.3","CACNG5","GALR1","ZNF560","NLRP11","COL9A3","SIM2","SOX10","NR0B1","TMSB15A","IL13RA2")


data.mvalues.probes <- data.mvalues.probes |> 
  dplyr::mutate(catnon_embryionic_development = unlist(pbapply::pblapply(gene, split_match, gid=tmp)) |
                                                unlist(pbapply::pblapply(gene_HGNC, split_match, gid=tmp))  ) |> 
  (function(.) {
    print(sum(.$catnon_embryionic_development))
    assertthat::assert_that(sum(.$catnon_embryionic_development) == (3526))
    return(.)
  })()


rm(tmp)

## filter for good probes ----


data.mvalues.good_probes <- data.mvalues.probes |> 
  dplyr::filter(good_probe) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (694299))
    return(.)
  })() |> 
  dplyr::pull(probe_id)



