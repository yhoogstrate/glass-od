#!/usr/bin/env R

# load ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')
#source('scripts/load_palette.R')
#source('scripts/load_themes.R')


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
  dplyr::mutate(probe_type = dplyr::recode(type, "I"="I (red & green)", "II"="II (ligation Allele-A)")) |>  # II is 1 (probe) and I is 2 probes.. ugh!
  assertr::verify(!duplicated(probe_id)) |> 
  dplyr::mutate(pos = round((CpG_beg + CpG_end )/2)) |> 
  dplyr::mutate(is_1P = CpG_chrm == 'chr1' & pos < 130 * 1000000) |> # rough margin
  dplyr::mutate(is_19Q = CpG_chrm == 'chr19' & pos > 23.5 * 1000000 ) |>  # rough margin
  dplyr::rename(target_hg38 = target) # noticed an example where target == 'CA' - while in hg19 ref genome it should have been actually 'CG'







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
  #dplyr::filter(grepl("^cg",IlmnID) & grepl("^cg",Name)) |> 
  dplyr::rename(probe_id = IlmnID) |> 
  dplyr::mutate(Name = NULL) |> 
  dplyr::mutate(AlleleA_ProbeSeq = NULL) |> # identical to the other manifest
  dplyr::mutate(AlleleB_ProbeSeq = NULL) |> # identical to the other manifest
  dplyr::rename(Forward_Sequence_hg19 = Forward_Sequence)

colnames(tmp)[colnames(tmp) %in% colnames(metadata.cg_probes.epic)]


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','')) 


### fix orientation ---
##' both `mapFlag_A` and `mapYD_A` occasionally do not match
##' `Forward_Sequence_hg19`, possibly because reference genomes changed
##' while those are trivial for subsequence / context analysis


fn <- "cache/load_probe_annotations__mapped_orientation.Rds"
if(!file.exists(fn)) {
  
  tmp <- metadata.cg_probes.epic |> 
    dplyr::filter(grepl("^cg", probe_id)) |> 
    
    dplyr::mutate(needle = gsub("R","A", AlleleA_ProbeSeq)) |> 
    dplyr::filter(!is.na(Forward_Sequence_hg19)) |> 
    dplyr::mutate(haystack = gsub("[][]","",Forward_Sequence_hg19)) |> 
    dplyr::mutate(haystack_rc =  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(haystack)))) |> 
    
    dplyr::mutate(haystack = gsub("G","A", haystack)) |> 
    dplyr::mutate(haystack_rc = gsub("G","A", haystack_rc)) |> 
    
    dplyr::mutate(n_fwd    = !is.na(stringr::str_locate(haystack,    needle)[,2])) |> 
    dplyr::mutate(n_fwd_rc = !is.na(stringr::str_locate(haystack_rc, needle)[,2])) |> 
    
    dplyr::mutate(orientation_mapped = dplyr::case_when(
      n_fwd == T & n_fwd_rc == F ~ "reverse",
      n_fwd == F & n_fwd_rc == T ~ "forward",
      T ~ "error"
    )) |> 
    
    dplyr::select(probe_id, orientation_mapped) |> 
    
    assertr::verify(orientation_mapped != "error") # orientation, AlleleA_ProbeSeq as compared to Forward_Sequence_hg19

  saveRDS(tmp, file=fn)
  
} else {
  
  tmp <- readRDS(file = fn)
  
}


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','')) 


rm(tmp, fn)




## RC the Forward_Sequence_hg19 when direction is forward ----


fn <- "cache/load_probe_annotations__sequence_pre_target_post.Rds"
if(!file.exists(fn)) {

  tmp <- metadata.cg_probes.epic |> 
    
    dplyr::filter(!is.na(orientation_mapped)) |>  # some QC probes
    dplyr::filter(!is.na(Forward_Sequence_hg19)) |> 
    
    #dplyr::mutate(Forward_Sequence_hg19 = "AAAAAT[CG]GCCCC", orientation = "reverse") |> 
    dplyr::mutate(sequence_pre = ifelse(orientation_mapped == "reverse",
                                        gsub("\\x5B.+$","",Forward_Sequence_hg19),
                                        spgs::reverseComplement(gsub("^.+\\x5D","",Forward_Sequence_hg19), case="as is")
    )) |> 
    dplyr::mutate(sequence_target = ifelse(orientation_mapped == "reverse",
                                           gsub("^.+\\x5B(.+)\\x5D.+$","\\1",Forward_Sequence_hg19),
                                           spgs::reverseComplement(gsub("^.+\\x5B(.+)\\x5D.+$","\\1",Forward_Sequence_hg19), case="as is")
    )) |> 
    dplyr::mutate(sequence_post = ifelse(orientation_mapped == "reverse",
                                         gsub("^.+\\x5D","",Forward_Sequence_hg19),
                                         spgs::reverseComplement(gsub("\\x5B.+$","",Forward_Sequence_hg19), case="as is")
    )) |> 
    assertr::verify(nchar(sequence_pre) == 60) |> 
    assertr::verify(nchar(sequence_target) == 2) |> 
    assertr::verify(nchar(sequence_post) == 60) |> 
    dplyr::select(probe_id, sequence_pre, sequence_target, sequence_post) |> 
    assertr::verify(probe_id %in% metadata.cg_probes.epic$probe_id)
  
  saveRDS(tmp, file=fn)
  
} else {
  
  tmp <- readRDS(file = fn)
  
}


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','')) 


rm(tmp, fn)




## proceed with annotations ----


fn <- "cache/load_probe_annotations__probe_type_annotation.Rds"
if(!file.exists(fn)) {
  
  tmp <- metadata.cg_probes.epic |> 
    dplyr::filter(grepl("^cg", probe_id)) |> 
    dplyr::mutate(probe_type_orientation = dplyr::case_when(
      type == "I"  & orientation_mapped == "forward" ~ "I - forward (red & green)",
      type == "I"  & orientation_mapped == "reverse" ~ "I - reverse (red & green)",
      type == "II" & orientation_mapped == "forward" ~ "II - forward (Allele-A)",
      type == "II" & orientation_mapped == "reverse" ~ "II - reverse (Allele-A)",
      T ~ "error"
    )) |> 
    dplyr::select(probe_id, probe_type_orientation) |> 
    assertr::verify(probe_type_orientation != "error")
  
  
  #table(tmp$probe_type_orientation)
  
  
  saveRDS(tmp, file=fn)
  
} else {
  
  tmp <- readRDS(file = fn)
  
}


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','')) 


rm(tmp, fn)




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
    dplyr::filter(!is.na(sequence_pre)) |> 
    #head(n=10) |> 
    dplyr::mutate(gc_sequence_context_1_new = paste0(gsub("^.+(.)$","\\1",sequence_pre),  "[", sequence_target, "]", gsub("^(.).+$","\\1",sequence_post))) |> 
    dplyr::mutate(gc_sequence_context_1_old = toupper(gsub("^.+(.[^A-Za-z][A-Za-z]{2}[^A-Za-z].).+$","\\1", Forward_Sequence_hg19))) |> 
    
    dplyr::mutate(gc_sequence_context_2_new = paste0(gsub("^.+(..)$","\\1",sequence_pre),  "[", sequence_target, "]", gsub("^(..).+$","\\1",sequence_post))) |> 
    dplyr::mutate(gc_sequence_context_2_old = toupper(gsub("^.+(..[^A-Za-z][A-Za-z]{2}[^A-Za-z]..).+$","\\1", Forward_Sequence_hg19))) |> 
    
    # dplyr::mutate(gc_sequence_context_l1 = toupper(gsub("^.+(.)[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_l2 = toupper(gsub("^.+(.).[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_l3 = toupper(gsub("^.+(.)..[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_l4 = toupper(gsub("^.+(.)...[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_l5 = toupper(gsub("^.+(.)....[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_l6 = toupper(gsub("^.+(.).....[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_l7 = toupper(gsub("^.+(.)......[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_l8 = toupper(gsub("^.+(.).......[^A-Za-z][A-Za-z]{2}[^A-Za-z].+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_r1 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z](.).+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_r2 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z].(.).+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_r3 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z]..(.).+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_r4 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z]...(.).+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_r5 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z]....(.).+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_r6 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z].....(.).+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_r7 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z]......(.).+$","\\1", Forward_Sequence_hg19))) |> 
    # dplyr::mutate(gc_sequence_context_r8 = toupper(gsub("^.+[^A-Za-z][A-Za-z]{2}[^A-Za-z].......(.).+$","\\1", Forward_Sequence_hg19))) |> 
    
    dplyr::select(probe_id, gc_sequence_context_1_new, gc_sequence_context_2_new,
                            gc_sequence_context_1_old, gc_sequence_context_2_old) |> 
    
    assertr::verify(nchar(gc_sequence_context_1_new) == 6) |> 
    assertr::verify(nchar(gc_sequence_context_1_old) == 6) |> 
    assertr::verify(nchar(gc_sequence_context_2_new) == 8) |> 
    assertr::verify(nchar(gc_sequence_context_2_old) == 8) |> 
    assertr::verify(probe_id %in% metadata.cg_probes.epic$probe_id)
  
  
    saveRDS(tmp, file="cache/load_probe_annotation__sequence_contexts.Rds")
    
} else {
  
  tmp <- readRDS(file="cache/load_probe_annotation__sequence_contexts.Rds")
  
}


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))

rm(tmp)




## add overall G-C content ----
#' split per pre- and post sequence


fn <- "cache/load_probe_annotation__gc_content.Rds"

if(!file.exists(fn)) {
  
  tmp <- metadata.cg_probes.epic |> 
    dplyr::filter(!is.na(sequence_pre))  |> 
  
    dplyr::mutate(sequence_pre_48 = gsub("^.+(.{48})$","\\1",sequence_pre)) |> # probes are 50 long minus 1 or 2 CG bases (depending on probe type)
    assertr::verify(nchar(sequence_pre_48) == 48) |> 
    dplyr::mutate(sequence_post_48 = gsub("^(.{48}).+$","\\1",sequence_post)) |> # probes are 50 long minus 1 or 2 CG bases (depending on probe type)
    assertr::verify(nchar(sequence_post_48) == 48) |> 
    dplyr::mutate(sequence_48_flank = paste0(sequence_pre_48, sequence_target, sequence_post_48)) |> 
    dplyr::mutate(sequence_60_flank = paste0(sequence_pre, sequence_target, sequence_post)) |> 

    dplyr::mutate(sequence_pre_48_c_content   = as.numeric(Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_pre_48  ), letters="C") / Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_pre_48  ), letters="ACTG"))) |> 
    dplyr::mutate(sequence_post_48_c_content  = as.numeric(Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_post_48 ), letters="C") / Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_post_48 ), letters="ACTG"))) |> 
    dplyr::mutate(sequence_48_flank_c_content = as.numeric(Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_48_flank), letters="C") / Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_48_flank), letters="ACTG"))) |> 
    dplyr::mutate(sequence_60_flank_c_content = as.numeric(Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_60_flank), letters="C") / Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_60_flank), letters="ACTG"))) |> 
  
    dplyr::mutate(sequence_pre_48_g_content   = as.numeric(Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_pre_48  ), letters="G") / Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_pre_48  ), letters="ACTG"))) |> 
    dplyr::mutate(sequence_post_48_g_content  = as.numeric(Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_post_48 ), letters="G") / Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_post_48 ), letters="ACTG"))) |> 
    dplyr::mutate(sequence_48_flank_g_content = as.numeric(Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_48_flank), letters="G") / Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_48_flank), letters="ACTG"))) |> 
    dplyr::mutate(sequence_60_flank_g_content = as.numeric(Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_60_flank), letters="G") / Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence_60_flank), letters="ACTG"))) |> 
    
    dplyr::mutate(sequence_pre_48_gc_content   = sequence_pre_48_c_content   + sequence_pre_48_g_content) |> 
    dplyr::mutate(sequence_post_48_gc_content  = sequence_post_48_c_content  + sequence_post_48_g_content) |> 
    dplyr::mutate(sequence_48_flank_gc_content = sequence_48_flank_c_content + sequence_48_flank_g_content) |> 
    dplyr::mutate(sequence_60_flank_gc_content = sequence_60_flank_c_content + sequence_60_flank_g_content) |> 
    
    dplyr::select(probe_id, ends_with("_content"))
  
  saveRDS(tmp, file=fn)
  
} else {
  
  tmp <- readRDS(file=fn)
  
}


metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('',''))

rm(tmp)


## overall CG-C's vs independent C's ----



fn <- "cache/load_probe_annotation__CG_Cs_vs_independent_Cs.Rds"

if(!file.exists(fn)) {
  
  tmp <- metadata.cg_probes.epic |> 
    dplyr::filter(grepl("^cg", probe_id)) |> 
    dplyr::select(probe_id, sequence_pre, orientation_mapped, AlleleA_ProbeSeq, probe_type) |> 
    dplyr::mutate(tmp = gsub("^.+(.{50})$","\\1",sequence_pre)) |> 
    dplyr::mutate(n_CG = stringr::str_count(tmp, pattern = "CG")) |> 
    dplyr::mutate(n_CG_to_CR = stringr::str_count(AlleleA_ProbeSeq, pattern = "R")) |> 
    dplyr::mutate(n_CG_to_CA = n_CG - n_CG_to_CR) |> 
    dplyr::mutate(n_independent_A = stringr::str_count(gsub("CG","XX",tmp), pattern = "A")) |> 
    dplyr::mutate(n_independent_C = stringr::str_count(gsub("CG","XX",tmp), pattern = "C")) |> 
    dplyr::mutate(n_independent_T = stringr::str_count(gsub("CG","XX",tmp), pattern = "T")) |> 
    dplyr::mutate(n_independent_G = stringr::str_count(gsub("CG","XX",tmp), pattern = "G")) |> 
    dplyr::select(probe_id, n_CG, n_CG_to_CR, n_CG_to_CA, n_independent_A, n_independent_C, n_independent_T, n_independent_G)
  
  saveRDS(tmp, file=fn)
  
} else {
  
  tmp <- readRDS(file=fn)
  
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
    dplyr::select(probe_id, Forward_Sequence_hg19) |> 
    dplyr::rowwise() |> 
    dplyr::mutate(closest_CG = closest_cg(Forward_Sequence_hg19)) |> 
    dplyr::ungroup() |> 
    dplyr::mutate(Forward_Sequence_hg19 = NULL)
  
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
    dplyr::select(probe_id, Forward_Sequence_hg19, closest_CG) |> 
    dplyr::mutate(tmp = gsub("^.+([ACTG])[^ACTG][CG]{2}[^ACTG]([ACTG]).+$","\\1:CG:\\2", Forward_Sequence_hg19)) |> 
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


## Repli-Tali predictors NatCom 2022 ----
#' https://zenodo.org/records/7108429
#' 10.1038/s41467-022-34268-8

if(!file.exists('data/RepliTali_coefs.csv')) {
  download.file('https://raw.githubusercontent.com/jamieendicott/Nature_Comm_2022/main/RepliTali/RepliTali_coefs.csv','data/RepliTali_coefs.csv')
}

tmp <- read.csv("data/RepliTali_coefs.csv") |> 
  dplyr::filter(grepl("^cg", Coefficient)) |> 
  dplyr::rename(probe_id = Coefficient) |> 
  dplyr::mutate(RepliTali_coef = Value) |> 
  dplyr::mutate(RepliTali_coef_type = ifelse(Value < 0, "negative", "positive")) |> 
  dplyr::mutate(is_RepliTali_coef = T) |> 
  dplyr::select(probe_id, contains("RepliTali"))



metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
  dplyr::left_join(tmp, by=c('probe_id'='probe_id'), suffix=c('','')) |> 
  dplyr::mutate(is_RepliTali_coef = ifelse(is.na(is_RepliTali_coef),F, is_RepliTali_coef))
rm(tmp)


## AC embryionic development genes ghisai ----
#'@todo fix gewoon die lijst met probes zoals in de supplement

# tmp <- c("MFAP2","OPRD1","CD101","MIR3681HG","AC064875.1","OSR1","LINC01121","AC007402.1","AC007402.2","TLX2","LINC01956","EN1","GRB14","HOXD11","HOXD10","HOXD-AS2","HOXD9","HOXD8","HOXD3","HOXD4",
#          "HAGLROS","IGFBP2","PAX3","TRIM71","AC112487.1","GPR156","SHOX2","IGF2BP2","LNX1","NMU","NKX6-1","HAND2","HAND2-AS1","ASB5","TENM3-AS1","NPR3","OTP","SPRY4","COL9A1","EYA4",
#          "AL139393.1","AL139393.3","TWIST1","IGF2BP3","HOXA2","HOXA3","HOXA5","HOXA6","HOXA7","HOXA9","HOXA10","HOXA11","HOXA11-AS","MNX1","NRG1","PLAT","GDF6","OSR2","PAX5","BX255923.2",
#          "AL353150.1","VAX1","EBF3","IGF2","IGF2-AS","AC090692.1","LINC02367","AC084816.1","AC009318.4","LRRK2","HOXC9","HOXC8","HOXC4","LINC01234","TBX5","GSX1","PAX9","LINC00648","C14orf39","SIX6",
#          "AC087473.1","ISL2","HAPLN3","LINC01579","GPR139","CRNDE","IRX5","AC126773.4","GJC1","TBX21","AC015909.3","CACNG5","GALR1","ZNF560","NLRP11","COL9A3","SIM2","SOX10","NR0B1","TMSB15A","IL13RA2")
# 
# 
# metadata.cg_probes.epic <- metadata.cg_probes.epic |> 
#   dplyr::mutate(catnon_embryionic_development = unlist(pbapply::pblapply(genesUniq, split_match, gid=tmp))   ) |> 
#   (function(.) {
#     print(sum(.$catnon_embryionic_development))
#     assertthat::assert_that(sum(.$catnon_embryionic_development) == (4310))
#     return(.)
#   })()
# 
# rm(tmp)


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

# exp <- readRDS("cache/load_probe_annotation__sequence_contexts.Rds") |> 
#   dplyr::select(probe_id, gc_sequence_context_2_new) |> 
#   dplyr::mutate(gc_sequence_context_2 = gsub("[][]","",gc_sequence_context_2)) |> 
#   dplyr::filter(!is.na(gc_sequence_context_2)) |> 
#   dplyr::rename(sequence_context_stranded = gc_sequence_context_2) |> 
#   dplyr::mutate(sequence_rc = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(sequence_context_stranded)))) |> 
#   dplyr::mutate(sequence_context_unstranded = dplyr::case_when(
#     sequence_context_stranded < sequence_rc ~ paste0(sequence_context_stranded,"/",sequence_rc),
#     sequence_context_stranded > sequence_rc ~ paste0(sequence_rc,"/",sequence_context_stranded),
#     sequence_context_stranded == sequence_rc ~ paste0(sequence_context_stranded,"*",sequence_rc)
#   )) |> 
#   dplyr::mutate(sequence_rc = NULL) |> 
#   dplyr::tibble()
# 
# 
