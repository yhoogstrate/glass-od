#!/usr/bin/env R


# load data ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')
#source('scripts/load_palette.R')
#source('scripts/load_themes.R')



if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}


if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}


if(!exists('data.protoemics.glass_nl')) {
  source('scripts/load_GLASS-NL_proteomics.R')
}



# 0. CpG to protein id ----

# 1stExon   3'UTR   5'UTR    Body ExonBnd TSS1500  TSS200 
# 207639   26759   20307   72027  308532    5483  110964   68775 


meta <- data.mvalues.probes |>
  dplyr::select(probe_id, UCSC_RefGene_Name, UCSC_RefGene_Group) |> 
  tidyr::separate_rows(
    UCSC_RefGene_Name,
    UCSC_RefGene_Group,
    sep = ";"
  ) |> 
  dplyr::filter(UCSC_RefGene_Name != "") |> 
  dplyr::distinct(probe_id, UCSC_RefGene_Name, UCSC_RefGene_Group, .keep_all = TRUE) |> 
  
  dplyr::filter(UCSC_RefGene_Group %in% c("TSS200","TSS1500","1stExon","5'UTR"))

tt_meth <- meta |> 
  dplyr::select(probe_id, UCSC_RefGene_Name) |> 
  dplyr::distinct(probe_id, UCSC_RefGene_Name, .keep_all = TRUE)






# 1a. load prot ----


metadata.prot <- glass_nl.metadata.array_samples |> 
  dplyr::filter(proteomics_sid %in% colnames(data.proteomics.glass_nl))


data.prot <- data.proteomics.glass_nl |> 
  dplyr::select(metadata.prot$proteomics_sid) |> 

  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('proteomics_sid') |> 
  dplyr::left_join(metadata.prot |> dplyr::select(proteomics_sid, "Sample_Name"), by=c('proteomics_sid'='proteomics_sid')) |> 
  tibble::column_to_rownames('Sample_Name') |> 
  dplyr::mutate(proteomics_sid = NULL) |> 
  t() |> 
  as.data.frame()


dim(data.prot)

data.prot <- data.prot |> 
  tibble::rownames_to_column('protein') |> 
  
  dplyr::rowwise() |>
  dplyr::filter(mean(dplyr::c_across(dplyr::where(is.numeric)), na.rm = TRUE) >= 6) |>
  dplyr::ungroup() |> 
  tibble::column_to_rownames('protein')


dim(data.prot)


# 1b. load meth data ----


metadata.meth <- glass_nl.metadata.array_samples |> 
  filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES)


data.meth <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  
  dplyr::select(metadata.meth$array_sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::left_join(metadata.meth |> dplyr::select(array_sentrix_id, "Sample_Name"), by=c('array_sentrix_id'='array_sentrix_id')) |> 
  tibble::column_to_rownames('Sample_Name') |> 
  dplyr::mutate(array_sentrix_id = NULL) |> 
  t() |> 
  as.data.frame()


dim(data.meth)



# 1c. load rna data ----


expression.glass.gtf <- read.delim('/home/youri/mnt/neuro-genomic-1-ro/glass/RNAseq/gencode.v34.primary_assembly.annotation.gtf',comment.char = "#",sep="\t",header=F) |> 
  dplyr::filter(V3 == "gene") |>
  dplyr::mutate(gene_id = gsub("^.*gene_id[ ]+([^;]+);.+$","\\1", V9)) |>
  dplyr::filter(grepl("_PAR_", gene_id) == F) |>  # these are odd equivalents of chrX positioned at chrY
  dplyr::mutate(ENSID = gsub("\\..+$","",gene_id)) |>
  dplyr::mutate(gene_name = gsub("^.*gene_name[ ]+([^;]+);.+$","\\1", V9)) |>
  dplyr::mutate(gene_type = gsub("^.*gene_type[ ]+([^;]+);.+$","\\1", V9)) |>
  dplyr::mutate(gene_uid = paste0(ENSID , "_", gene_name)) |> 
  dplyr::select(gene_id, gene_name)



expression.glass.exon <- read.delim('/home/youri/mnt/neuro-genomic-1-ro/glass/RNAseq/alignments/alignments-new/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt',skip=1,header=T) |>
  dplyr::rename_with(~ gsub( "^X.+new\\.", "", .x)) |>
  dplyr::rename_with(~ gsub( ".Aligned.+bam$", "", .x,)) |>
  dplyr::rename_with(~ gsub(".","-", .x,fixed=T)) |>
  dplyr::filter(grepl("_PAR_", Geneid) == F) |>
  dplyr::rename(gene_id = Geneid) |>
  dplyr::left_join(expression.glass.gtf |> dplyr::select(gene_id, gene_name),by=c('gene_id'='gene_id')) |>
  dplyr::select(-c(gene_id, Chr, Start, End, Strand, Length))

a <- expression.glass.exon |> 
  dplyr::group_by(gene_name) |>
  dplyr::summarise(
    # Sums the read counts across all columns that are numeric
    dplyr::across(where(is.numeric), sum),
    .groups = "drop" # Ungroup the data after summarization
  )

b <- a |> 
  dplyr::rowwise() |>
  dplyr::filter(mean(dplyr::c_across(dplyr::where(is.numeric)), na.rm = TRUE) >= 8) |>
  dplyr::ungroup() |> 
  tibble::column_to_rownames('gene_name')


expression.glass.exon = b



metadata.rna <- glass_nl.metadata.array_samples |> 
  dplyr::filter(!is.na(RNA_seq_sid)) |> 
  dplyr::filter(RNA_seq_sid %in% colnames(expression.glass.exon)) |> 
  dplyr::filter(!is.na(Sample_Name)) |> 
  dplyr::select(-starts_with("PC"))



data.rna <- expression.glass.exon |> 
  DESeq2::DESeqDataSetFromMatrix(data.frame(cond = as.factor(paste0('c',round(runif(ncol(expression.glass.exon)))+1) )), ~cond) |> 
  DESeq2::vst(blind=T) |> 
  SummarizedExperiment::assay() |> 
  as.data.frame(stringsAsFactors=F) |> 
  dplyr::select(metadata.rna$RNA_seq_sid) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('RNA_seq_sid') |> 
  dplyr::left_join(metadata.rna |> dplyr::select(RNA_seq_sid, "Sample_Name"), by=c('RNA_seq_sid'='RNA_seq_sid')) |> 
  tibble::column_to_rownames('Sample_Name') |> 
  dplyr::mutate(RNA_seq_sid = NULL) |> 
  t() |> 
  as.data.frame()


data.rna[1:8,1:8]



# 3a. correlation of per gene protein value / methylation value ----


metadata.combi <- dplyr::inner_join(
  metadata.prot |> dplyr::select(Sample_Name),
  metadata.meth |> dplyr::select(Sample_Name),
  by=c('Sample_Name'='Sample_Name')
) |> 
  dplyr::left_join(glass_nl.metadata.array_samples, by=c('Sample_Name'='Sample_Name'))


data.combi.prot <- data.prot |> dplyr::select(metadata.combi$Sample_Name)

data.combi.meth <- data.meth |> dplyr::select(metadata.combi$Sample_Name) |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::left_join(tt_meth |> dplyr::select(probe_id, UCSC_RefGene_Name), by=c('probe_id'= 'probe_id')) |> 
  dplyr::filter(!is.na(UCSC_RefGene_Name)) |> 
  dplyr::mutate(probe_id = NULL) |> 
  dplyr::group_by(UCSC_RefGene_Name) |>
  dplyr::summarise(
    across(where(is.numeric), mean, .names = "{.col}"),
    .groups = "drop"
  ) |> 
  tibble::column_to_rownames('UCSC_RefGene_Name')


stopifnot(colnames(data.combi.prot) == colnames(data.combi.meth))

isct <- intersect(rownames(data.combi.meth), rownames(data.combi.prot))

data.combi.meth <- data.combi.meth |>
  tibble::rownames_to_column("UCSC_RefGene_Name") |>
  dplyr::filter(UCSC_RefGene_Name %in% c(isct)) |> 
  dplyr::arrange(UCSC_RefGene_Name) |>
  tibble::column_to_rownames("UCSC_RefGene_Name")

# Get the order established by data.combi.meth
sorted_rownames <- rownames(data.combi.meth)

# Reorder data.combi.prot using the established order
data.combi.prot <- data.combi.prot[sorted_rownames, ]


stopifnot(colnames(data.combi.prot) == colnames(data.combi.meth))
stopifnot(rownames(data.combi.prot) == rownames(data.combi.meth))


c = cor(t(data.combi.prot) , t(data.combi.meth ), method="spearman")
d = diag(c)


title <- (meta |> 
            dplyr::pull(UCSC_RefGene_Group) |> 
            unique() |> 
            sort() |> 
            paste(collapse = ", "))
plot(sort(d),ylab="Per gene correlation mean methylation - proteomics", xlab = "Gene / protein",
     main=paste0(title, " x proteomics"), cex.main=0.7)
abline(h=0)
text(0.015*length(d), 0.1,paste0("median: ",round(median(d),3)), col="red", adj=0)
abline(h=median(d), col="gray",lty=2)




# 3b. correlation prot x RNA ----


metadata.combi <- dplyr::inner_join(
  metadata.prot |> dplyr::select(Sample_Name),
  metadata.rna |> dplyr::select(Sample_Name),
  by=c('Sample_Name'='Sample_Name')
) |> 
  dplyr::left_join(glass_nl.metadata.array_samples, by=c('Sample_Name'='Sample_Name'))



data.combi.prot <- data.prot |> dplyr::select(metadata.combi$Sample_Name)
data.combi.rna <- data.rna |> dplyr::select(metadata.combi$Sample_Name)

order <- data.frame(gene_id = intersect(rownames(data.combi.prot),rownames(data.combi.rna)))

data.combi.prot <- order |> 
  dplyr::left_join(data.combi.prot |> tibble::rownames_to_column('gene_id'), by=c('gene_id'='gene_id')) |> 
  tibble::column_to_rownames('gene_id')

data.combi.rna <- order |> 
  dplyr::left_join(data.combi.rna |> tibble::rownames_to_column('gene_id'), by=c('gene_id'='gene_id')) |> 
  tibble::column_to_rownames('gene_id')



stopifnot(colnames(data.combi.prot) == colnames(data.combi.rna))
stopifnot(rownames(data.combi.prot) == rownames(data.combi.rna))


c = cor(t(data.combi.prot), t(data.combi.rna), method="spearman")
print(dim(c))
d = diag(c)


plot(sort(d),ylab="Per gene correlation mean methylation - proteomics", xlab = "Gene / protein", 
     main="protein vs RNA", cex.main=0.7)
abline(h=0)
text(0.015*length(d), 0.2, paste0("median: ",round(median(d),3)), col="red", adj=0)
abline(h=median(d), col="gray",lty=2)




# 3c. DNA meth x RNA ----



metadata.combi <- dplyr::inner_join(
  metadata.rna |> dplyr::select(Sample_Name),
  metadata.meth |> dplyr::select(Sample_Name),
  by=c('Sample_Name'='Sample_Name')
) |> 
  dplyr::left_join(glass_nl.metadata.array_samples, by=c('Sample_Name'='Sample_Name'))


# 1. columns
data.combi.rna <- data.rna |> dplyr::select(metadata.combi$Sample_Name)

data.combi.meth <- data.meth |> dplyr::select(metadata.combi$Sample_Name) |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::left_join(tt_meth |> dplyr::select(probe_id, UCSC_RefGene_Name), by=c('probe_id'= 'probe_id')) |> 
  dplyr::filter(!is.na(UCSC_RefGene_Name)) |> 
  dplyr::mutate(probe_id = NULL) |> 
  dplyr::group_by(UCSC_RefGene_Name) |>
  dplyr::summarise(
    across(where(is.numeric), mean, .names = "{.col}"),
    .groups = "drop"
  ) |> 
  tibble::column_to_rownames('UCSC_RefGene_Name')


stopifnot(colnames(data.combi.meth) == colnames(data.combi.rna))




# 2. rows

isct <- intersect(rownames(data.combi.meth), rownames(data.combi.prot))

data.combi.rna <- data.combi.rna |>
  tibble::rownames_to_column("UCSC_RefGene_Name") |>
  dplyr::filter(UCSC_RefGene_Name %in% c(isct)) |> 
  dplyr::arrange(UCSC_RefGene_Name) |>
  tibble::column_to_rownames("UCSC_RefGene_Name")

# Get the order established by data.combi.meth
sorted_rownames <- rownames(data.combi.rna)

# Reorder data.combi.prot using the established order
data.combi.meth <- data.combi.meth[sorted_rownames, ]


stopifnot(colnames(data.combi.rna) == colnames(data.combi.meth))
stopifnot(rownames(data.combi.rna) == rownames(data.combi.meth))


c = cor(t(data.combi.rna), t(data.combi.meth), method="spearman")
d = diag(c)


title <- (meta |> 
            dplyr::pull(UCSC_RefGene_Group) |> 
            unique() |> 
            sort() |> 
            paste(collapse = ", "))
plot(sort(d),ylab="Per gene correlation mean methylation - rna", 
     xlab = "Gene / protein", 
     main=paste0(title," x RNA"), cex.main=0.7)
abline(h=0)
text(0.015*length(d), 0.1,paste0("median: ",round(median(d),3)), col="red", adj=0)
abline(h=median(d), col="gray",lty=2)







# 4. correlation of change (DMP) per protein / per methylation ----




