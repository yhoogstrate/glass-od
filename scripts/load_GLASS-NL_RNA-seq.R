#!/usr/bin/env R

# load GTF file (gene annot) ----


expression.glass.gtf <- read.delim('data/GLASS_NL/RNAseq/gencode.v34.primary_assembly.annotation.gtf',comment.char = "#",sep="\t",header=F) |> 
  dplyr::filter(V3 == "gene") |> 
  dplyr::mutate(gene_id = gsub("^.*gene_id[ ]+([^;]+);.+$","\\1", V9)) |> 
  dplyr::filter(grepl("_PAR_", gene_id) == F) |>   # these are odd equivalents of chrX positioned at chrY
  dplyr::mutate(ENSID = gsub("\\..+$","",gene_id)) |>  
  dplyr::mutate(gene_name = gsub("^.*gene_name[ ]+([^;]+);.+$","\\1", V9)) |> 
  dplyr::mutate(gene_type = gsub("^.*gene_type[ ]+([^;]+);.+$","\\1", V9)) |> 
  dplyr::mutate(gene_uid = paste0(ENSID , "_", gene_name))




# load count exonic features ----


expression.glass.exon <- read.delim('data/GLASS_NL/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt',skip=1,header=T) |>
  dplyr::rename_with(~ gsub("^X.+new.","", .x)) |>
  dplyr::rename_with(~ gsub(".Aligned.+bam$","", .x)) |>
  dplyr::rename_with(~ gsub(".","-", .x, fixed=T)) |>
  dplyr::filter(grepl("_PAR_", Geneid) == F) |>  # these are odd equivalents of chrX positioned at chrY
  dplyr::rename(gene_id = Geneid) |> 
  dplyr::left_join(expression.glass.gtf |> dplyr::select(gene_id, gene_uid),by=c('gene_id'='gene_id'))



## merge GTF and featureCounts per-gene stats ----


expression.glass.exon.metadata <- expression.glass.exon |> 
  dplyr::select(gene_id, Chr, Start, End, Strand, Length) |> 
  dplyr::left_join(expression.glass.gtf, by=c('gene_id' = 'gene_id')) |> 
  dplyr::rename(gene_chr = V1) |>  
  dplyr::rename(gene_strand = V7) |> 
  dplyr::mutate(gene_chr_center_loc = (V4 + V5) /  2) |> 
  dplyr::mutate(gene_loc = paste0(gene_chr, ":", round(gene_chr_center_loc / 1000000),"M")) 


# cleanup ---

expression.glass.exon <- expression.glass.exon |> 
  dplyr::select(-gene_id, -Chr, -Start, -End, - Strand, -Length) |>
  tibble::column_to_rownames('gene_uid')


