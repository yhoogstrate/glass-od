#!/usr/bin/env R

# load ----

source('scripts/load_functions.R')

# info ----

# https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/infinium-methylationepic-manifest-column-headings.pdf

# https://gdc.cancer.gov/content/improved-dna-methylation-array-probe-annotation

# EPIC ----

# https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/EPIC.hg38.manifest.tsv.gz
# + https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/archive/202209/EPIC.hg38.mask.tsv.gz
# + gencode 36 annotation https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/EPIC.hg38.manifest.gencode.v36.tsv.gz

# to be done
# + https://github.com/zhou-lab/KYCG_knowledgebase_EPIC


## manifest ----

# 2022
metadata.cg_probes.epic <- read.table("data/Improved DNA Methylation Array Probe Annotation/EPIC/EPIC.hg38.manifest.tsv", header=T) |> 
  dplyr::rename(probe_id = Probe_ID) |> 
  assertr::verify(!duplicated(probe_id)) |> 
  dplyr::mutate(pos = round((CpG_beg + CpG_end )/2)) |> 
  dplyr::mutate(is_1P = CpG_chrm == 'chr1' & pos < 130 * 1000000) |> # rough margin
  dplyr::mutate(is_19Q = CpG_chrm == 'chr19' & pos > 23.5 * 1000000 ) # rough margin


# old manifest for probeCpGcnt & context35 ----
# these variables are only found in some strange files, presumably an old version of the manifest
# probeCpGcnt: the number of CpG in the probe.
# context35: the number of CpG in the [-35bp, +35bp] window.

tmp <- read.delim("~/mnt/neuro-genomic-1-ro/catnon/Methylation - EPIC arrays/EPIC.hg38.manifest.tsv.gz") |> 
  dplyr::rename(probe_id = probeID) |> 
  dplyr::select(probe_id, probeCpGcnt, context35)

metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','')) 




## add masking ----


tmp <- read.table("data/Improved DNA Methylation Array Probe Annotation/EPIC/EPIC.hg38.mask.tsv", header=T) |> 
  dplyr::rename(probe_id = probeID)


stopifnot(intersect(colnames(tmp) , colnames(metadata.cg_probes.epic)) == c("probe_id"))


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','')) |> 
  dplyr::mutate(exclude = is.na(MASK_general) | MASK_general == T) # control probes (^ctl_.+$ instead of ^gc[0-9]+$) have no MASK_general

rm(tmp)



## manifest from illumina ----

# https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html
# Infinium MethylationEPIC v1.0 B5 Manifest File (BPM Format)
# 172 MB
# Mar 13, 2020


tmp <- read.csv("data/Improved DNA Methylation Array Probe Annotation/EPIC/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip=7,header=T, sep=",") |> 
  dplyr::filter(grepl("^cg",IlmnID) & grepl("^cg",Name)) |> 
  dplyr::rename(probe_id = IlmnID) |> 
  dplyr::mutate(Name = NULL) |> 
  dplyr::rename(AlleleA_ProbeSeq_Illumina_manifest = AlleleA_ProbeSeq) |> 
  dplyr::rename(AlleleB_ProbeSeq_Illumina_manifest = AlleleB_ProbeSeq)


colnames(tmp)[colnames(tmp) %in% colnames(metadata.cg_probes.epic)]



metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','')) 


## G-CIMP ----

### IDH-mut ----


tmp <- read.delim('data/GLASS_NL/Methylation/Analysis/RF TCGA/IDHmtprobes.txt') |> 
  dplyr::rename(probe_id = 1)

metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::mutate(GCIMP_IDHmut_probe = probe_id %in% tmp$probe_id)

rm(tmp)



### IDH-wt ----


tmp <- read.delim('data/GLASS_NL/Methylation/Analysis/RF TCGA/IDHwtprobes.txt') |> 
  dplyr::rename(probe_id = 1)

metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::mutate(GCIMP_IDHwt_probe = probe_id %in% tmp$probe_id)

rm(tmp)


### PanGlioma ----


tmp <- read.delim('data/GLASS_NL/Methylation/Analysis/RF TCGA/pangliomaprobes.txt') |> 
  dplyr::rename(probe_id = 1)

metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::mutate(GCIMP_PanGlioma_probe = probe_id %in% tmp$probe_id)

rm(tmp)


### magic 90 ----


load("data/GLASS_NL/Methylation/Analysis/RF TCGA/GCIMPlow.90probes.Rda")
tmp <- GCIMPlow.probes |> 
  dplyr::rename(probe_id = probes)


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::mutate(GCIMP_low_90 = probe_id %in% tmp$probe_id)

rm(tmp)


## add gene annotation ----


tmp <- read.table("data/Improved DNA Methylation Array Probe Annotation/EPIC/EPIC.hg38.manifest.gencode.v36.tsv", header=T) |> 
  dplyr::rename(probe_id = probeID) |> 
  dplyr::mutate(CpG_chrm = NULL, CpG_beg  = NULL, CpG_end  = NULL)


stopifnot(intersect(colnames(tmp) , colnames(metadata.cg_probes.epic)) == c("probe_id"))


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))

rm(tmp)




# majority has no probe sequence?!
# plt |>
#   dplyr::filter(is.na(probe_edit_distance)) |>
#   head() |> 
#   dplyr::select(probe_id, AlleleA_ProbeSeq, AlleleB_ProbeSeq, probe_edit_distance)

## add sequence-context ----


if(!file.exists("cache/load_probe_annotation__sequence_contexts.Rds")) {
  
  tmp <- metadata.cg_probes.epic |> 
    dplyr::select(probe_id, Forward_Sequence) |> 
    dplyr::mutate(gc_sequence_context_1 = toupper(gsub("^.+(.[^A-Za-z][A-Za-z]{2}[^A-Za-z].).+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_2 = toupper(gsub("^.+(..[^A-Za-z][A-Za-z]{2}[^A-Za-z]..).+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_l1 = toupper(gsub("^.+(.)[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_l2 = toupper(gsub("^.+(.).[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_l3 = toupper(gsub("^.+(.)..[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_l4 = toupper(gsub("^.+(.)...[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_l5 = toupper(gsub("^.+(.)....[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_l6 = toupper(gsub("^.+(.).....[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_l7 = toupper(gsub("^.+(.)......[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_l8 = toupper(gsub("^.+(.).......[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_r1 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z](.).+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_r2 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z].(.).+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_r3 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z]..(.).+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_r4 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z]...(.).+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_r5 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z]....(.).+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_r6 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z].....(.).+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_r7 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z]......(.).+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(gc_sequence_context_r8 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z].......(.).+$","\\1", Forward_Sequence))) |> 
    dplyr::mutate(Forward_Sequence = NULL)
  
    saveRDS(tmp, file="cache/load_probe_annotation__sequence_contexts.Rds")
} else {
  tmp <- readRDS(file="cache/load_probe_annotation__sequence_contexts.Rds")
}


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))

rm(tmp)



## nearest CG ----



if(!file.exists("cache/load_probe_annotation__nearest_CG.Rds")) {
  
  closest_cg <- function(dna_string) {
    out <- gsub("[CG]","x",dna_string, fixed=T) |> 
      stringr::str_locate_all("CG") |> 
      as.data.frame() |> 
      dplyr::mutate(start = ifelse(start < 60, start - 60 ,start - 61), end = NULL) |> 
      dplyr::filter(start != 0)  |> 
      rbind(data.frame(start = 61)) |> 
      dplyr::arrange(abs(start)) |> 
      head(n=1) |> 
      dplyr::pull(start)
    
    return (as.numeric(out))
  }
  
  stopifnot(closest_cg("CTAAGTGCAGTCAGGATCTGTTAGTACAGTGGCTTTTGATGGAACAGCTGAGGCACACAT[CG]CGCGTGGCATGGACTCCGGGGCCGAACGCTCACGACCAAGACTTTTGCCCTTTTGAAATG") == 1)
  stopifnot(closest_cg("CTAAGTGCAGTCAGGATCTGTTAGTACAGTGGCTTTTGATGGAACAGCTGAGGCACACCG[CG]AACGTGGCATGGACTCCGGGGCCGAACGCTCACGACCAAGACTTTTGCCCTTTTGAAATG") == -1)
  stopifnot(closest_cg("CTAAGTGCAGTCAGGATCTGTTAGTACAGTGGCTTTTGATGGAACAGCTGAGGCACACAT[CG]ACGCGTGGCATGGACTCCGGGGCCGAACGCTCACGACCAAGACTTTTGCCCTTTTGAAAT") == 2)
  stopifnot(closest_cg("CTAAGTGCAGTCAGGATCTGTTAGTACAGTGGCTTTTGATGGAACAGCTGAGGCACACGA[CG]AACGTGGCATGGACTCCGGGGCCGAACGCTCACGACCAAGACTTTTGCCCTTTTGAAATG") == -2)
  
  tmp <- metadata.cg_probes.epic |> 
    #head(n=1000) |> 
    dplyr::select(probe_id, Forward_Sequence) |> 
    dplyr::rowwise() |> 
    dplyr::mutate(closest_CG = closest_cg(Forward_Sequence)) |> 
    dplyr::ungroup() |> 
    dplyr::mutate(Forward_Sequence = NULL)
  
  saveRDS(tmp, file="cache/load_probe_annotation__nearest_CG.Rds")
  
} else {
  tmp <- readRDS(file="cache/load_probe_annotation__nearest_CG.Rds")
}


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))


rm(tmp)


## is solo-WCGW ----
# >= 35 & == wCGw [https://www.nature.com/articles/s41467-022-34268-8]
# w: https://en.wikipedia.org/wiki/FASTA_format = [AT]


if(!file.exists("cache/load_probe_annotation__is_solo_WCGW.Rds")) {
  
  tmp <- metadata.cg_probes.epic |> 
    dplyr::select(probe_id, Forward_Sequence, closest_CG) |> 
    dplyr::mutate(tmp = gsub("^.+([ACTG])[^ACTG][CG]{2}[^ACTG]([ACTG]).+$","\\1:CG:\\2", Forward_Sequence)) |> 
    dplyr::mutate(is_solo_WCGW = (closest_CG >= 35) & (tmp %in% c("A:CG:A", "A:CG:T", "T:CG:A", "T:CG:T"))) |> 
    dplyr::select(probe_id, is_solo_WCGW)
  
  saveRDS(tmp, file="cache/load_probe_annotation__is_solo_WCGW.Rds")
  
} else {
  tmp <- readRDS(file="cache/load_probe_annotation__is_solo_WCGW.Rds")
}


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))


rm(tmp)




## repli-seq data ----


tmp <- read.table("output/tables/repli-seq_data.txt", header=T)
dplymetadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('', ''))



## AC embryionic development genes ghisai ----


tmp <- c("MFAP2","OPRD1","CD101","MIR3681HG","AC064875.1","OSR1","LINC01121","AC007402.1","AC007402.2","TLX2","LINC01956","EN1","GRB14","HOXD11","HOXD10","HOXD-AS2","HOXD9","HOXD8","HOXD3","HOXD4",
         "HAGLROS","IGFBP2","PAX3","TRIM71","AC112487.1","GPR156","SHOX2","IGF2BP2","LNX1","NMU","NKX6-1","HAND2","HAND2-AS1","ASB5","TENM3-AS1","NPR3","OTP","SPRY4","COL9A1","EYA4",
         "AL139393.1","AL139393.3","TWIST1","IGF2BP3","HOXA2","HOXA3","HOXA5","HOXA6","HOXA7","HOXA9","HOXA10","HOXA11","HOXA11-AS","MNX1","NRG1","PLAT","GDF6","OSR2","PAX5","BX255923.2",
         "AL353150.1","VAX1","EBF3","IGF2","IGF2-AS","AC090692.1","LINC02367","AC084816.1","AC009318.4","LRRK2","HOXC9","HOXC8","HOXC4","LINC01234","TBX5","GSX1","PAX9","LINC00648","C14orf39","SIX6",
         "AC087473.1","ISL2","HAPLN3","LINC01579","GPR139","CRNDE","IRX5","AC126773.4","GJC1","TBX21","AC015909.3","CACNG5","GALR1","ZNF560","NLRP11","COL9A3","SIM2","SOX10","NR0B1","TMSB15A","IL13RA2")


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::mutate(catnon_embryionic_development = unlist(pbapply::pblapply(genesUniq, split_match, gid=tmp))   ) |> 
  (function(.) {
    print(sum(.$catnon_embryionic_development))
    assertthat::assert_that(sum(.$catnon_embryionic_development) == (4310))
    return(.)
  })()

rm(tmp)


# Wies signature ----


tmp <- read.csv("data/GLASS_NL/Metadata/(Epi)genetic_data/ProbeSelection_IvR_FDR 1e-9 Delta 1_06072022.csv") |> 
  dplyr::rename(probe_id = Probe_ID) |> 
  dplyr::mutate(deep_significant = T) |> 
  dplyr::rename_with( ~ paste0("glass_nl_prim_rec__",.x), .cols=!matches("^probe_id$",perl = T))


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','')) |> 
  dplyr::mutate(glass_nl_prim_rec__deep_significant = ifelse(is.na(glass_nl_prim_rec__deep_significant), F, glass_nl_prim_rec__deep_significant))




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


