#!/usr/bin/env R

# load data ----

if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}


if(!exists('expression.glass.exon')) {
  source('scripts/load_GLASS-NL_RNA-seq.R')
}


# not even needed thus far
# source('scripts/load_mvalues_hq_samples.R')



# perform regression on CGC & RNA-seq of matching samples

meta <- glass_nl.metadata.array_samples |> 
  filter_GLASS_NL_idats(CONST_N_GLASS_NL_INCLUDED_SAMPLES) |> 
  dplyr::filter(!is.na(array_A_IDH_HG__A_IDH_LG_lr_v12.8 & !is.na(RNA_seq_sid))) |> 
  dplyr::select(Sample_Name, array_A_IDH_HG__A_IDH_LG_lr_v12.8, RNA_seq_sid) |> 
  dplyr::filter(!is.na(RNA_seq_sid)) |> 
  dplyr::rename(CGC = array_A_IDH_HG__A_IDH_LG_lr_v12.8) |> 
  dplyr::mutate(CGC = scale(CGC, center=T))


dat <- expression.glass.exon |> dplyr::select(meta$RNA_seq_sid)



# plot avg OR expression ----


# avg OR expression ----


plt <- dat |> 
  tibble::rownames_to_column('gene_uid') |> 
  dplyr::filter(grepl("_OR[0-9]", gene_uid)) |> 
  tibble::column_to_rownames('gene_uid') |> 
  #as.matrix() |> 
  rowMeans() |> 
  data.frame() |> 
  dplyr::rename(mean_read_count = 1) |> 
  dplyr::mutate(mean_read_count = round(mean_read_count,3)) |> 
  tibble::rownames_to_column('gene_uid')


ggplot(plt, aes(x = reorder(gene_uid, mean_read_count), y=mean_read_count)) +
  geom_point() +
  scale_y_continuous(trans = 'log1p')





# clean up low expressed  ----

rn = dat |>
  rowMeans()  |>
  data.frame () |>
  dplyr::rename(count = 1) |> 
  dplyr::filter(count >= 1.0) |> 
  tibble::rownames_to_column('gene_uid') |>
  dplyr::pull('gene_uid')

dat_trunc <- dat |> 
  tibble::rownames_to_column('gene_uid') |>
  dplyr::filter(gene_uid %in% rn) |> 
  tibble::column_to_rownames('gene_uid')





# deseq ----


dds <- DESeq2::DESeqDataSetFromMatrix(countData = dat_trunc, colData = meta, design= ~CGC)
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds) |> 
    as.data.frame() |> 
    dplyr::filter(!is.na(padj)) |> 
    dplyr::arrange(pvalue,padj) |> 
    tibble::rownames_to_column('gene_uid') |> 
    dplyr::left_join(expression.glass.exon.metadata , 
                       dplyr::select(gene_uid, gene_name, gene_type, gene_strand, gene_chr, gene_chr_center_loc, gene_loc),by=c('gene_uid'='gene_uid'))


# View(res)



ggplot(res, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point(data=subset(res, !grepl("^OR[0-9]", gene_name)), pch=16, cex=0.001, col="black") +
  geom_point(data=subset(res, grepl("^OR[0-9]", gene_name)), pch=16, cex=1, col="red") +
  theme_nature




plt <- res |> 
  dplyr::mutate(OR = ifelse(grepl("^OR[0-9]", gene_name), "Yes", "No"))


ggplot(plt, aes(x=OR, y=stat)) + 
  ggbeeswarm::geom_quasirandom(data=subset(plt, OR == "No"), pch=19,cex=0.01) +
  ggbeeswarm::geom_quasirandom(data=subset(plt, OR == "Yes"), pch=19,cex=0.7, col="red") +
  #ggpubr::stat_compare_means(method = "wilcox.test") +
  ggpubr::stat_compare_means(method = "t.test") +
  #geom_boxplot(width=0.75, outlier.shape=NA, outlier.color=NA, col="darkgray") +
  gghalves::geom_half_violin(fill=mixcol("#FFFFFF", "#EEEEEE"),side = "r", draw_quantiles = c(0.5), linewidth=theme_nature_lwd, col = palette_repel_text_col, width = 1.6) +
  theme_nature +
  labs(y="RNA-seq expression co-increase with scaled CGC")




