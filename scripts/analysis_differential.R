#!/usr/bin/env R


# load data ----


library(ggplot2)
library(minfi)
source('scripts/youri_gg_theme.R')


if(!exists('glass_od.metadata.patients')) {
  source('scripts/load_metadata.R')
}

if(!exists('glass_od.data.mvalues')) {
  source('scripts/load_mvalues.R')
}


# unpaired grading ----



metadata <- glass_od.metadata.idats |> 
  dplyr::filter(is.na(reason_excluded_patient)) |> 
  dplyr::filter(is.na(reason_excluded_resection)) |> 
  dplyr::filter(is.na(reason_excluded_resection_isolation)) |> 
  dplyr::filter(is.na(reason_excluded_array_sample)) |> 
  dplyr::filter(study_name == "GLASS-OD") |> 
  assertr::verify(!is.na(qc.pca.outlier))  |> 
  dplyr::filter(qc.pca.outlier == F) |> 
  dplyr::filter(!is.na(resection_tumor_grade)) |> 
  dplyr::mutate(grading = ifelse(resection_tumor_grade == 2, 0, 1))


data <- glass_od.data.mvalues |> 
  dplyr::select(metadata$sentrix_id)


stopifnot(metadata$sentrix_id == colnames(data))


design <- model.matrix(~0 + grading, data=metadata)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit)

top <- limma::topTable(fit, n=nrow(fit), sort="p") |> # adjust="BH"
  tibble::rownames_to_column('cg') |> 
  dplyr::mutate('A_IDH_LG_HG_predictor' = cg %in% probes)


ggplot(top, aes(x = logFC, y=-log10(adj.P.Val), col=A_IDH_LG_HG_predictor)) +
  geom_point(data = top |> dplyr::filter(A_IDH_LG_HG_predictor==F), pch=19,cex=0.1,alpha=0.5) +
  geom_point(data = top |> dplyr::filter(A_IDH_LG_HG_predictor==T), pch=19,cex=0.5)



# fully paired analysis ----
## obtain patients ----



# glass_od.metadata.idats$reason_exclusion_resection_isolation
array_samples_all <- glass_od.metadata.idats |> 
  dplyr::filter(is.na(reason_excluded_patient)) |> 
  dplyr::filter(is.na(reason_excluded_resection)) |> 
  dplyr::filter(is.na(reason_excluded_resection_isolation)) |> 
  dplyr::filter(is.na(reason_excluded_array_sample)) |> 
  dplyr::filter(!is.na(qc.pca.outlier ) & qc.pca.outlier == F)


array_samples_primary <- array_samples_all |> 
  dplyr::filter(resection_number == 1)


array_samples_recurrent <- array_samples_all |>  # take last recurrence, to max out the time effect
  dplyr::filter(resection_number > 1) |> 
  dplyr::arrange(desc(resection_id), desc(resection_isolation_id)) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::top_n(n=1)


metadata <- rbind(
  array_samples_primary |> dplyr::mutate(primary_recurrence_status = "primary")  ,
  array_samples_recurrent |> dplyr::mutate(primary_recurrence_status = "last_recurrence")
) |> 
  dplyr::filter(patient_id %in% intersect(array_samples_primary$patient_id, array_samples_recurrent$patient_id)) |> 
  dplyr::mutate(patient = as.factor(paste0("pat",patient_id)))


rm(array_samples_all, array_samples_primary, array_samples_recurrent)


## limma ----


data <- glass_od.data.mvalues |> 
  dplyr::select(metadata$sentrix_id)

#data <- data |> 
#  dplyr::slice_head(n=4000)

stopifnot(metadata$sentrix_id == colnames(data))


design <- model.matrix(~0 + patient + primary_recurrence_status, data=metadata)
fit <- limma::lmFit(data, design)

saveRDS(fit, file="cache/analysis_differential_fully-paired.Rds")




## dmpFinder ----

targets <- array_samples |> 
  dplyr::rename(Sample_Name = resection_isolation_id) |>
  dplyr::mutate(Array = gsub("^.+_","",sentrix_id)) |> 
  dplyr::rename(Slide = methylation_array_chip_id) |> 
  dplyr::mutate(Basename = gsub("_Grn.idat$","", channel_green)) |> 
  dplyr::select(Sample_Name, sentrix_id, status, Array,Slide, Basename) 




RGSet <- read.metharray.exp(targets = targets, force =T)


# manifest <- getManifest(RGSet)
# manifest
# 
# #MSet.illumina <- preprocessIllumina(RGSet, bg.correct = TRUE,
# #                                    normalize = "controls")
# 
# GRset.funnorm <- preprocessFunnorm(RGSet)
# 
# beta <- getBeta(GRset.funnorm)
# status  <- pData(GRset.funnorm)$status


# store m-values

proc <- preprocessNoob(RGSet, offset = 0, dyeCorr = TRUE, verbose = TRUE, dyeMethod="reference") 
proc.r.mvalue <- ratioConvert(proc , what = "M")

mvalues <- as.data.frame(assays(proc.r.mvalue)$M)

stats <- data.frame(mad = pbapply::pbapply(mvalues, 1, mad)) |> 
  dplyr::mutate(
    mean.a =
      pbapply::pbapply(mvalues |> 
                         dplyr::select(targets |> dplyr::filter(status == "primary") |> dplyr::mutate(id = gsub("^.+/","",Basename)) |>  dplyr::pull(id)),
                       1, mean, na.rm=T)
  ) |> 
  dplyr::mutate(
    mean.b =
      pbapply::pbapply(mvalues |> 
                         dplyr::select(targets |> dplyr::filter(status == "last_recurrence") |> dplyr::mutate(id = gsub("^.+/","",Basename)) |>  dplyr::pull(id)),
                       1, mean, na.rm=T)
  )




stats <- stats |> 
  dplyr::mutate(delta1 = mean.b - mean.a) |> 
  tibble::rownames_to_column('probe_id') 



beta <- getBeta(proc)
status  <- pData(proc)$status


dmp <- dmpFinder(beta, pheno = status  , type = "categorical") |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::left_join(stats,
                   by=c('probe_id'='probe_id'),
                   suffix = c('','')
                     )






ggplot(dmp |> 
         dplyr::filter(mad > 1) 
       , aes(x=delta1, y=-log10(pval))) +
  geom_point(pch=19,cex=0.1, alpha=0.1) +
  labs(x = "delta M-value primary recurrence [20 random matching pairs]") +
  youri_gg_theme



# odd pca?

pc <- prcomp(t(mvalues))
pc_plt <- pc$x |> 
  as.data.frame() |> 
  tibble::rownames_to_column('sentrix_id') |> 
  dplyr::left_join(targets, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))

ggplot(pc_plt, aes(x=PC1, y=PC2, label=Sample_Name) ) +
  #geom_point() +
  geom_text()

## 0006-R1
## 0019-R1




plt <- read.delim("~/mnt/neuro-genomic-1-ro/catnon/Methylation - EPIC arrays/EPIC.hg38.manifest.tsv.gz") |>
  dplyr::filter(MASK_general == F) |> 
  dplyr::mutate(pos = round((CpG_beg + CpG_end )/2)) |> 
  dplyr::select(CpG_chrm, pos, probeID, probe_strand  ) |> 
  dplyr::left_join(dmp, by=c('probeID' = 'probe_id'), suffix=c('','')) |> 
  dplyr::filter(!is.na(delta1) & !is.na(qval)) |> 
  dplyr::rename(chr = CpG_chrm) |> 
  dplyr::mutate(chr = factor(chr, levels=gtools::mixedsort(unique(as.character(chr))) )) |> 
  dplyr::mutate(significant = qval < 0.000001) |> 
  dplyr::filter(chr %in% c("chrM","chrX","chrY", "NA",NA) == F)  |> 
  dplyr::mutate(absdiff = abs(delta1))


ggplot(plt, aes(x=pos / 1000000,y=delta1, col=chr)) + 
  #geom_vline(xintercept=175, col="blue", alpha=0.5) + 
  facet_grid(cols = vars(chr), scales = "free", space="free") +
  geom_point(pch=19,cex=0.2, alpha=0.25) +
  #geom_point(data = subset(plt, significant==T), pch=21,cex=0.8,col='black',fill=NA, alpha=0.5) +
  geom_smooth(se=F,col="black", lwd=0.7, n=500) +
  youri_gg_theme +
  labs(x=NULL)


# maybe stranded adds sth?

ggplot(subset(plt, chr == "chr2"), aes(x=pos / 1000000,y=delta1, col=chr)) + 
  geom_vline(xintercept=176.125, col="blue", alpha=0.5) + 
  facet_grid(cols = vars(chr), scales = "free", space="free") +
  geom_point(pch=19,cex=0.2) +
  geom_point(data = subset(plt, significant==T & chr == "chr2"), pch=21,cex=0.8,col='black',fill=NA, alpha=0.5) +
  geom_smooth(se=F,col="black", lwd=0.7,n=1000, span=0.2) +
  youri_gg_theme +
  labs(x=NULL) + 
  xlim(170,180)



# fully paired analysis ----

