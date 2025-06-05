#!/usr/bin/env R 


if(!exists('data.mvalues.probes')) { # cached M-values
  source('scripts/load_mvalues_hq_samples.R')
}


# load annotations necessary to calc
annot <- read.csv("data/Improved DNA Methylation Array Probe Annotation/EPIC/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip=7,header=T, sep=",") |> 
  dplyr::rename(probe_id = IlmnID) |> 
  dplyr::rename(Forward_Sequence_hg19 = Forward_Sequence) |> 
  dplyr::select(probe_id, AlleleA_ProbeSeq, Forward_Sequence_hg19)


# do the magic part 1: calc orientation
annot <- annot |> 
  dplyr::filter(grepl("^cg", probe_id)) |> 
  
  # calc orientation mapped
  dplyr::mutate(needle = gsub("R","A", AlleleA_ProbeSeq)) |> 
  dplyr::filter(!is.na(Forward_Sequence_hg19)) |> 
  dplyr::mutate(haystack = gsub("[][]","", Forward_Sequence_hg19)) |> 
  dplyr::mutate(haystack_rc =  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(haystack)))) |> 
  
  dplyr::mutate(haystack = gsub("G","A", haystack)) |> 
  dplyr::mutate(haystack_rc = gsub("G","A", haystack_rc)) |> 
  
  dplyr::mutate(n_fwd    = !is.na(stringr::str_locate(haystack,    needle)[,2])) |> 
  dplyr::mutate(n_fwd_rc = !is.na(stringr::str_locate(haystack_rc, needle)[,2])) |> 
  
  dplyr::mutate(orientation_mapped = dplyr::case_when(
    n_fwd == T & n_fwd_rc == F ~ "reverse",
    n_fwd == F & n_fwd_rc == T ~ "forward",
    T ~ "error"
  )) 




annot <- annot |> # magic part 2: using orientation determine sequence context 
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
  
  dplyr::mutate(gc_sequence_context_2_new = paste0(gsub("^.+(..)$","\\1",sequence_pre),  "[", sequence_target, "]", gsub("^(..).+$","\\1",sequence_post))) |> 
  
  dplyr::select(probe_id, gc_sequence_context_2_new)








