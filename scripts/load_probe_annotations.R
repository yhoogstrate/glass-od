#!/usr/bin/env R

# https://gdc.cancer.gov/content/improved-dna-methylation-array-probe-annotation

# EPIC ----

# https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/EPIC.hg38.manifest.tsv.gz
# + https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/archive/202209/EPIC.hg38.mask.tsv.gz
# + gencode 36 annotation https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/EPIC.hg38.manifest.gencode.v36.tsv.gz

# to be done
# + https://github.com/zhou-lab/KYCG_knowledgebase_EPIC


## manifest ----

metadata.cg_probes.epic <- read.table("data/Improved DNA Methylation Array Probe Annotation/EPIC/EPIC.hg38.manifest.tsv", header=T) |> 
  dplyr::rename(probe_id = Probe_ID)


## add masking ----


tmp <- read.table("data/Improved DNA Methylation Array Probe Annotation/EPIC/EPIC.hg38.mask.tsv", header=T) |> 
  dplyr::rename(probe_id = probeID)


stopifnot(intersect(colnames(tmp) , colnames(metadata.cg_probes.epic)) == c("probe_id"))


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','')) |> 
  dplyr::mutate(exclude = is.na(MASK_general) | MASK_general == T) # control probes (^ctl_.+$ instead of ^gc[0-9]+$) have no MASK_general

rm(tmp)


## add gene annotation ----


tmp <- read.table("data/Improved DNA Methylation Array Probe Annotation/EPIC/EPIC.hg38.manifest.gencode.v36.tsv", header=T) |> 
  dplyr::rename(probe_id = probeID) |> 
  dplyr::mutate(CpG_chrm = NULL, CpG_beg  = NULL, CpG_end  = NULL)


stopifnot(intersect(colnames(tmp) , colnames(metadata.cg_probes.epic)) == c("probe_id"))


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))

rm(tmp)



# 450K ----



# https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.tsv.gz
# + https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/archive/202209/HM450.hg38.mask.tsv.gz
# + gencode 36 annotation https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.gencode.v36.tsv.gz

# to be done
# + https://github.com/zhou-lab/KYCG_knowledgebase_EPIC


## manifest ----


metadata.cg_probes.450k <- read.table("data/Improved DNA Methylation Array Probe Annotation/450k/HM450.hg38.manifest.tsv", header=T) |> 
  dplyr::rename(probe_id = Probe_ID)


## add masking ----


tmp <- read.table("data/Improved DNA Methylation Array Probe Annotation/450k/HM450.hg38.mask.tsv", header=T) |> 
  dplyr::rename(probe_id = probeID)


stopifnot(intersect(colnames(tmp) , colnames(metadata.cg_probes.450k)) == c("probe_id"))


metadata.cg_probes.450k <- metadata.cg_probes.450k |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','')) |> 
  dplyr::mutate(exclude = is.na(MASK_general) | MASK_general == T) # control probes (^ctl_.+$ instead of ^gc[0-9]+$) have no MASK_general

rm(tmp)



## add gene annotation ----


tmp <- read.table("data/Improved DNA Methylation Array Probe Annotation/450k/HM450.hg38.manifest.gencode.v36.tsv", header=T) |> 
  dplyr::rename(probe_id = probeID) |> 
  dplyr::mutate(CpG_chrm = NULL, CpG_beg  = NULL, CpG_end  = NULL)


stopifnot(intersect(colnames(tmp) , colnames(metadata.cg_probes.450k)) == c("probe_id"))


metadata.cg_probes.450k <- metadata.cg_probes.450k |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))

rm(tmp)



# 27k ----


