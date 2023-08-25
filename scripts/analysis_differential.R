#!/usr/bin/env R


# load data ----


library(ggplot2)
library(minfi)
source('scripts/youri_gg_theme.R')
source('scripts/load_functions.R')


if(!exists('data.mvalues.hq_samples')) {
  source('scripts/load_mvalues_hq_samples.R')
}



if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}



# analyses: GLASS-OD primary - last recurrence ----
## data: example ----

metadata.example <- data.frame(pid = (1:100) %% 50) |> 
  dplyr::arrange(pid) |> 
  dplyr::mutate(pid = factor(paste0("p", pid))) |> 
  dplyr::mutate(sid = factor(paste0("s",1:100))) |> 
  dplyr::mutate(condition = factor(c(rep("c1",75),rep("c2",25))))


replag <- function(c) {
  d <- c()
  
  x <- 1:length(c)
  y2 <- (x)*2
  y1 <- y2 - 1
  
  d[y1] <- c[x]
  d[y2] <- c[x]
  
  return(d)
}

reppos <- function(c, pl) {
  d = c()
  
  for(k in 1:length(c)) {
    p = levels(pl)[k]
    
    d[which(pl == p)] <- c[k]
  }
  
  return(d)
}


#reppos(runif(3), factor(c("p1","p4","p6","p6","p1","p4")))





metadata.example <- data.frame(
  condition = c("c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                "c1","c2",
                
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                "c1","c1",
                
                "c2","c2",
                "c2","c2",
                "c2","c2",
                "c2","c2",
                "c2","c2")
) |> 
  dplyr::mutate(pat = factor(replag(paste0("p",1:35)))) |> 
  dplyr::arrange(condition, pat)



data.example <- data.frame(
  `nodiff_A_001` = c(runif(45), runif(25)) + reppos(runif(35), metadata.example$pat),
  `nodiff_A_002` = c(runif(45), runif(25)) + reppos(runif(35), metadata.example$pat),
  `nodiff_A_003` = c(runif(45), runif(25)) + reppos(runif(35), metadata.example$pat),
  `nodiff_A_004` = c(runif(45), runif(25)) + reppos(runif(35), metadata.example$pat),
  `nodiff_A_005` = c(runif(45), runif(25)) + reppos(runif(35), metadata.example$pat),
  
  `nodiff_B_001` = c(runif(45), runif(25)) + reppos(rep(0, 35), metadata.example$pat),
  `nodiff_B_002` = c(runif(45), runif(25)) + reppos(rep(0, 35), metadata.example$pat),
  `nodiff_B_003` = c(runif(45), runif(25)) + reppos(rep(0, 35), metadata.example$pat),
  `nodiff_B_004` = c(runif(45), runif(25)) + reppos(rep(0, 35), metadata.example$pat),
  `nodiff_B_005` = c(runif(45), runif(25)) + reppos(rep(0, 35), metadata.example$pat),
  
  `nodiff_high_batch_A_001` = c(runif(45), runif(25)) + reppos(runif(35,10,20), metadata.example$pat),
  `nodiff_high_batch_A_002` = c(runif(45), runif(25)) + reppos(runif(35,10,20), metadata.example$pat),
  `nodiff_high_batch_A_003` = c(runif(45), runif(25)) + reppos(runif(35,10,20), metadata.example$pat),
  `nodiff_high_batch_A_004` = c(runif(45), runif(25)) + reppos(runif(35,10,20), metadata.example$pat),
  `nodiff_high_batch_A_005` = c(runif(45), runif(25)) + reppos(runif(35,10,20), metadata.example$pat),
  
  `C1_high_batch_A_001` = c(runif(45,1,2), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_002` = c(runif(45,1,2), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_003` = c(runif(45,1.5,3), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_004` = c(runif(45,1.5,3), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_005` = c(runif(45,10,12), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_006` = c(runif(45,7,8), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_007` = c(runif(45,0,5), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_008` = c(runif(45,0,5), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_009` = c(runif(45,0,5), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_high_batch_A_010` = c(runif(45,0,5), runif(25)) + reppos(runif(35), metadata.example$pat),
  
  
  `C1_low_batch_A_001` = c(runif(45,0,0.1), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_low_batch_A_002` = c(runif(45,0,0.1), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_low_batch_A_003` = c(runif(45,0,0.1), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_low_batch_A_004` = c(runif(45,0,0.1), runif(25)) + reppos(runif(35), metadata.example$pat),
  `C1_low_batch_A_005` = c(runif(45,0,0.1), runif(25)) + reppos(runif(35), metadata.example$pat)
  
) |> 
  t() |> 
  as.data.frame()


### test: indeed p.values and not f-test ----


design.example <- model.matrix(~pat + condition, data=metadata.example)
fit.example <- limma::lmFit(as.matrix(data.example), design.example)
fit.example <- limma::eBayes(fit.example,trend=T)
stats.example <- limma::topTable(fit.example, n=nrow(fit.example)) # adjust="BH", sort="p")

pval = as.data.frame(fit.example$p.value) |> 
  tibble::rownames_to_column('testid') 

ggplot(pval, aes(x=testid, y=conditionc2)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))


coef <- as.data.frame(fit.example$coefficients) |> 
  tibble::rownames_to_column('testid') 

ggplot(coef, aes(x=testid, y=conditionc2)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))


plot(coef$conditionc2, -log(pval$conditionc2))





## data: paired ----


metadata.p <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(163) |> 
  filter_primaries_and_last_recurrences(136) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::filter(dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p",patient_id))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 
  dplyr::mutate(batch = ifelse(isolation_person_name == "USA / Duke", "batch [US]", "batch [EU]")) |> 
  dplyr::select(sentrix_id, patient, pr.status, batch)


tmp <- rbind(
  metadata.p |> 
    dplyr::filter(pr.status == "primary") |> 
    dplyr::mutate(patient = as.character(patient)) |> 
    dplyr::mutate(patient = sample(patient, replace=F)) |> 
    assertr::verify(!duplicated(patient))
  ,
  metadata.p |> 
    dplyr::filter(pr.status == "recurrence") |> 
    dplyr::mutate(patient = as.character(patient)) |> 
    dplyr::mutate(patient = sample(patient, replace=F)) |> 
    assertr::verify(!duplicated(patient))
) |> 
  dplyr::mutate(patient_shuffled = as.factor(patient)) |> 
  dplyr::select(sentrix_id, patient_shuffled)

metadata.p <- metadata.p |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'),suffix=c('',''))
rm(tmp)


data.p <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.p$sentrix_id) |> 
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |> # this synthax, oh my
  dplyr::arrange(mad) |>
  dplyr::mutate(mad = NULL) 
  #dplyr::slice_head(n=125000)


### test: pat + condition  ----


design.p.p1 <- model.matrix(~factor(patient) + factor(pr.status), data=metadata.p)
fit.p.p1 <- limma::lmFit(data.p, design.p.p1)
fit.p.p1 <- limma::eBayes(fit.p.p1,trend=T)
stats.p.p1 <- limma::topTable(fit.p.p1,
                n=nrow(data.p),
                coef="factor(pr.status)recurrence",
                sort.by = "none",
                adjust.method="fdr") |> 
  tibble::column_to_rownames('probe_id') 


rm(design.p.p1, fit.p.p1)


sum(stats.p.p1$P.Value < 0.01)
sum(stats.p.p1$adj.P.Val < 0.01)



#plot(coef.p.p1, -log(pval.p.p1),pch=19,cex=0.2)




### test: pat_shuffled + condition ----


design.p.p1shuf <- model.matrix(~factor(patient_shuffled) + factor(pr.status), data=metadata.p)
fit.p.p1shuf <- limma::lmFit(data.p, design.p.p1shuf)
fit.p.p1shuf <- limma::eBayes(fit.p.p1shuf,trend=T)
stats.p.p1shuf <- limma::topTable(fit.p.p1shuf,
                                  n=nrow(data.p),
                                  coef="factor(pr.status)recurrence",
                                  sort.by = "none",
                                  adjust.method="fdr") |> 
  tibble::column_to_rownames('probe_id') 


rm(design.p.p1shuf, fit.p.p1shuf)


sum(stats.p.p1shuf$P.Value < 0.01)
sum(stats.p.p1shuf$adj.P.Val < 0.01)



### test: batch[US/EU] + condition ----

design.p.b1 <- model.matrix(~factor(batch) + factor(pr.status), data=metadata.p)
fit.p.b1 <- limma::lmFit(data.p, design.p.b1)
fit.p.b1 <- limma::eBayes(fit.p.b1,trend=T)
stats.p.b1 <- limma::topTable(fit.p.b1,
                              n=nrow(data.p),
                              coef="factor(pr.status)recurrence",
                              sort.by = "none",
                              adjust.method="fdr") |> 
  tibble::column_to_rownames('probe_id') 


rm(design.p.b1, fit.p.b1)


sum(stats.p.b1$P.Value < 0.01)
sum(stats.p.b1$adj.P.Val < 0.01)



### test:  removeBatchEffects() => unpaired ----
#' to see how big the effect of diverse primaries is

#t() |> #scale() |> #t() |>

dn.p <- data.p |>
  limma::removeBatchEffect(metadata.p$patient)


design.p.rmb.u <- model.matrix(~factor(pr.status), data=metadata.p)
fit.p.rmb.u <- limma::lmFit(dn.p, design.p.rmb.u)
fit.p.rmb.u <- limma::eBayes(fit.p.rmb.u,trend=T)
stats.p.rmb.u <- limma::topTable(fit.p.rmb.u,
                   n=nrow(data.p),
                   coef="factor(pr.status)recurrence",
                   sort.by = "none",
                   adjust.method="fdr") |> 
  tibble::column_to_rownames('probe_id') 


rm(dn.p, design.p.rmb.u, fit.p.rmb.u)



### test: condition (unpaired) ----
#' to see how big the effect of diverse primaries is


design.p.u <- model.matrix(~factor(pr.status), data=metadata.p)
fit.p.u <- limma::lmFit(data.p, design.p.u)
fit.p.u <- limma::eBayes(fit.p.u,trend=T)
stats.p.u <- limma::topTable(fit.p.u, 
               n=nrow(data.p),
               coef="factor(pr.status)recurrence",
               sort.by = "none",
               adjust.method="fdr") |> 
  tibble::column_to_rownames('probe_id') 

rm(design.p.u, fit.p.u)




## data: partially paired ----

metadata.pp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(163) |> 
  filter_primaries_and_last_recurrences(136) |> 
  
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"remainder")))) |> 
  dplyr::mutate(pr.status = factor(ifelse(resection_number == 1,"primary","recurrence"),levels=c("primary","recurrence"))) |> 
  dplyr::select(sentrix_id, patient, pr.status)


data.pp <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.pp$sentrix_id) |> 
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |> # this synthax, oh my
  dplyr::arrange(mad) |>
  dplyr::mutate(mad = NULL)
  #dplyr::slice_head(n=125000)




### test: pat + condition  ----


design.pp.p1 <- model.matrix(~factor(patient) + factor(pr.status), data=metadata.pp)
fit.pp.p1 <- limma::lmFit(data.pp, design.pp.p1)
fit.pp.p1 <- limma::eBayes(fit.pp.p1, trend=T)
stats.pp.p1 <- limma::topTable(fit.pp.p1,
                               n=nrow(data.pp),
                               coef="factor(pr.status)recurrence",
                               sort.by = "none",
                               adjust.method="fdr") |> 
  tibble::column_to_rownames('probe_id') 

rm(design.pp.p1, fit.pp.p1)


sum(stats.pp.p1$P.Value < 0.01)
sum(stats.pp.p1$adj.P.Val < 0.01)



### removeBatchEffects() => unpaired ----
#' to see how big the effect of diverse primaries is


dn.pp <- data.pp |> 
  limma::removeBatchEffect(metadata.pp$patient)


design.pp.rmb.u <- model.matrix(~factor(pr.status), data=metadata.pp)
fit.pp.rmb.u <- limma::lmFit(dn.pp, design.pp.rmb.u)
fit.pp.rmb.u <- limma::eBayes(fit.pp.rmb.u,trend=T)
stats.pp.rmb.u <- limma::topTable(fit.pp.rmb.u,
                                  n=nrow(data.pp),
                                  coef="factor(pr.status)recurrence",
                                  sort.by = "none",
                                  adjust.method="fdr") |> 
  tibble::column_to_rownames('probe_id') 

rm(dn.pp, design.pp.rmb.u, fit.pp.rmb.u)


### without pairing ----


design.pp.u <- model.matrix(~factor(pr.status), data=metadata.pp)
fit.pp.u <- limma::lmFit(data.pp, design.pp.u)
fit.pp.u <- limma::eBayes(fit.pp.u, trend=T)
stats.pp.u <- limma::topTable(fit.pp.u,
                              n=nrow(data.pp),
                              coef="factor(pr.status)recurrence",
                              sort.by = "none",
                              adjust.method="fdr") |> 
  tibble::column_to_rownames('probe_id') 

rm(design.pp.u, fit.pp.u)


sum(stats.pp.u$P.Value < 0.01)
sum(stats.pp.u$adj.P.Val < 0.01)



## plot: sort(pvals) ----



plt <- rbind(
  stats.p.b1 |> dplyr::mutate(type = "~batch + cond", paired = paste0("full [n=",nrow(metadata.p),"]")),
  stats.p.p1 |> dplyr::mutate(type = "~patient + cond", paired = paste0("full [n=",nrow(metadata.p),"]")),
  stats.p.p1shuf |> dplyr::mutate(type = "~patient_shuffled + cond", paired = paste0("full [n=",nrow(metadata.p),"]")),
  #stats.p.rmb.u |> dplyr::mutate(type = "removeBatchEffects() => cond", paired = paste0("full [n=",nrow(metadata.p),"]")),
  stats.p.u |> dplyr::mutate(type = "~cond", paired = paste0("full [n=",nrow(metadata.p),"]")),
  
  stats.pp.p1 |> dplyr::mutate(type = "~patient + cond", paired = paste0("partial [n=",nrow(metadata.pp),"]")),
  #stats.pp.rmb.u |> dplyr::mutate(type = "removeBatchEffects() => cond", paired = paste0("partial [n=",nrow(metadata.pp),"]")),
  stats.pp.u |> dplyr::mutate(type = "~cond", paired = paste0("partial [n=",nrow(metadata.pp),"]"))
) |> 
  dplyr::mutate(group = paste0(paired, ", ",type )) |> 
  dplyr::arrange(group, P.Value) |> 
  dplyr::group_by(group) |> 
  dplyr::mutate(x = 1:dplyr::n()) |> 
  dplyr::ungroup()


ggplot(plt,aes(x=x, y=P.Value, group = group, col=type, lty=paired)) +
  geom_line() +
  theme_bw() +
  labs(x= "probe, ordered by p-value") +
  geom_hline(yintercept = 0.01, col="red",lty=2,lwd=0.15)




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


# analyses: GLASS-OD g2 - g3 ----
#' @todo ASK what to do if first resection is G3 while last is G2 !!!

## data: partially paired ----

metadata.pp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(163) |> 
  filter_first_G2_and_last_G3(105) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"remainder")))) |> 
  assertr::verify(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(gr.status = ifelse(resection_tumor_grade == 2, "Grade2","Grade3"))



data.pp <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.pp$sentrix_id) |> 
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |> # this synthax, oh my
  dplyr::arrange(mad) |>
  dplyr::mutate(mad = NULL)




### test: pat + condition  ----


design.g23 <- model.matrix(~factor(patient) + factor(gr.status), data=metadata.pp)
fit.gr <- limma::lmFit(data.pp, design.g23)
fit.gr <- limma::eBayes(fit.gr, trend=T)
stats.gr <- limma::topTable(fit.gr,
                            n=nrow(data.pp),
                            coef="factor(gr.status)Grade3",
                            sort.by = "none",
                            adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.g23, fit.gr)


sum(stats.gr$P.Value < 0.01)
sum(stats.gr$adj.P.Val < 0.01)


plot(sort(stats.gr$P.Value),type="l")


# analyses: GLASS-OD AcCGAP ----
## data: partially paired ----

metadata.pp <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(163) |> 
  assertr::verify(!is.na(A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  assertr::verify(is.numeric(A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  dplyr::mutate(A_IDH_HG__A_IDH_LG_lr__lasso_fit = scale(A_IDH_HG__A_IDH_LG_lr__lasso_fit)) |> 
  dplyr::group_by(patient_id) |> 
  dplyr::mutate(is.paired = dplyr::n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(patient = as.factor(paste0("p",ifelse(is.paired,patient_id,"remainder"))))



data.pp <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.pp$sentrix_id) |> 
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |> # this synthax, oh my
  dplyr::arrange(mad) |>
  dplyr::mutate(mad = NULL)




### test: pat + condition  ----


design.lgc <- model.matrix(~factor(patient) + as.numeric(A_IDH_HG__A_IDH_LG_lr__lasso_fit), data=metadata.pp)
fit.lgc <- limma::lmFit(data.pp, design.lgc)
fit.lgc <- limma::eBayes(fit.lgc, trend=T)
stats.lgc <- limma::topTable(fit.lgc,
                            n=nrow(data.pp),
                            coef="as.numeric(A_IDH_HG__A_IDH_LG_lr__lasso_fit)",
                            sort.by = "none",
                            adjust.method="fdr") |> 
  tibble::rownames_to_column('probe_id') 

rm(design.lgc, fit.lgc)


sum(stats.lgc$P.Value < 0.01)
sum(stats.lgc$adj.P.Val < 0.01)

plot(sort(stats.lgc$P.Value),type="l")


# plot: g2-g3 x lgc ----



plt <- stats.gr |> 
  dplyr::left_join(stats.lgc,by=c('probe_id'='probe_id'),suffx=c('.gr','.lgc'))


plot(sort(stats.lgc$P.Value),type="l")


## powerplot ----
#' take intersect, then order


plt.pre <- stats.gr |> 
  dplyr::left_join(stats.lgc,by=c('probe_id'='probe_id'),suffix=c('.gr','.lgc')) |> 
  dplyr::mutate(col = ifelse(
    probe_id %in% (data.mvalues.probes |> dplyr::filter(catnon_embryionic_development) |>  dplyr::pull(probe_id)),
    "embrionic development", "other"
  ))

plt <- data.frame(
  p.gr = plt.pre |> dplyr::arrange(P.Value.gr) |> dplyr::pull(P.Value.gr),
  p.lgc = plt.pre |> dplyr::arrange(P.Value.lgc) |> dplyr::pull(P.Value.lgc)
) |> 
  dplyr::mutate(x = 1:dplyr::n()) |> 
  dplyr::mutate(delta = p.gr - p.lgc)


#plot(plt$x, plt$delta, type="l")
plot(plt$x, plt$p.gr, type="l")
lines(plt$x, plt$p.lgc)


plot(plt.pre$logFC.gr, plt.pre$logFC.lgc, pch=19,cex=0.1)
abline(h=0, col="red")
abline(v=0, col="red")


ggplot(plt.pre, aes(x=logFC.lgc, y=-log(P.Value.lgc), col=col)) + 
  geom_point(data = plt.pre |> dplyr::filter(col =="other"), pch=19, cex=0.1, alpha=0.1) + 
  geom_point(data = plt.pre |> dplyr::filter(col !="other"), pch=19, cex=0.5,alpha=0.8) +
  geom_vline(xintercept=0, col="red", lty=2, lwd=0.5) +
  geom_hline(yintercept=0, col="red", lty=2, lwd=0.5) +
  theme_bw()



# analyses: GLASS-OD + GLASS-NL LG/HG ----
#' multi-dataset mixed model to find changes specific to
#' OD
#' AC
#' general grading
## test pat + HG_AC + HG_OD + overall_HG ----


metadata.od <- glass_od.metadata.array_samples |> 
    filter_GLASS_OD_idats(163) |> 
    filter_first_G2_and_last_G3(105) |> 
    dplyr::select(sentrix_id, resection_id, patient_id, LG_HG_status, A_IDH_HG__A_IDH_LG_lr__lasso_fit) |>  # has to be WHO since OD has no methylation based split
    dplyr::rename(sample_id = resection_id) |> 
    dplyr::mutate(dataset = "GLASS-OD") |> 
    dplyr::mutate(LG_HG_status = factor(LG_HG_status, levels=c("LG","HG")))
metadata.ac <- glass_nl.metadata.array_samples |> 
      filter_GLASS_NL_idats(218) |> 
      filter_first_G2_and_last_G3(130) |> 
      dplyr::select(sentrix_id, Sample_Name, patient_id, LG_HG_status, A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV) |>  # has to be WHO since OD has no methylation based split
      dplyr::rename(sample_id = Sample_Name) |> 
      dplyr::rename(`A_IDH_HG__A_IDH_LG_lr__lasso_fit` = `A_IDH_HG__A_IDH_LG_lr__lasso_fit__10xCV`) |> 
      dplyr::mutate(dataset = "GLASS-NL") |> 
      dplyr::mutate(LG_HG_status = factor(LG_HG_status, levels=c("LG","HG")))

data.od <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.od$sentrix_id)

data.ac <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata.ac$sentrix_id)


design.od <- model.matrix(~factor(patient_id) + factor(LG_HG_status), data=metadata.od)
design.ac <- model.matrix(~factor(patient_id) + factor(LG_HG_status), data=metadata.ac)

fit.od <- limma::eBayes(limma::lmFit(data.od, design.od),trend=T)
fit.ac <- limma::eBayes(limma::lmFit(data.ac, design.ac),trend=T)

stats.od <- limma::topTable(fit.od,
                              n=nrow(data.od),
                              coef="factor(LG_HG_status)HG",
                              sort.by = "none", 
                              adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')
stats.ac <- limma::topTable(fit.ac,
                              n=nrow(data.ac),
                              coef="factor(LG_HG_status)HG",
                              sort.by = "none",
                              adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')


stopifnot(stats.od$probe_id == stats.ac$probe_id)

plt <- stats.od |> 
  dplyr::left_join(stats.ac, by=c('probe_id'='probe_id'),suffix=c('.od','.ac')) |> 
  dplyr::mutate(col1 = probe_id %in% (stats.b.hg_od |> dplyr::filter(adj.P.Val < 0.01) |> dplyr::pull(probe_id)))|> 
  dplyr::mutate(col2 = probe_id %in% (stats.b.hg_ac |> dplyr::filter(adj.P.Val < 0.01) |> dplyr::pull(probe_id)))


plt <- plt |> 
  dplyr::mutate(col3 = probe_id %in% (stats.b.hg_ac$probe_id[stats.b.hg_ac$adj.P.Val < 0.01]))

ggplot(plt, aes(x=`t.od`,y=`t.ac`, col=col1)) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  geom_point(pch=19, alpha=0.15,cex=0.02) +
  theme_bw()


ggplot(plt, aes(x=`t.od`,y=`t.ac`, col=col2)) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  geom_point(pch=19, alpha=0.15,cex=0.02) +
  theme_bw()


ggplot(plt, aes(x=`t.od`,y=`t.ac`, col=col3)) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  geom_point(data=subset(plt, col3==F),pch=19, alpha=0.15,cex=0.02) +
  geom_point(data=subset(plt, col3==T),pch=19, alpha=0.65,cex=0.05) +
  theme_bw()


ggplot(plt, aes(x=`t.od`,y=`t.ac`, col=col2)) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  geom_point(data=subset(plt, col2==F),pch=19, alpha=0.15,cex=0.02) +
  geom_point(data=subset(plt, col2==T),pch=19, alpha=0.65,cex=0.05) +
  theme_bw()


# x alzheimer quick test ----

plt <- stats.od |> 
  dplyr::left_join(data.probes.alzheimer, by=c('probe_id'='probe_id'),suffix=c('',''))

# ggplot(plt, aes(x=logFC, y=-adj.P.Val, col=is.na(FDR.p.value))) +
#   geom_point(data = subset(plt, is.na(FDR.p.value)),pch=19,cex=0.01) +
#   geom_point(data = subset(plt, !is.na(FDR.p.value)),pch=19,cex=0.15)

ggplot(plt, aes(x=Beta..difference, y=logFC, col=is.na(FDR.p.value), label=label)) +
  geom_point(data = subset(plt, is.na(FDR.p.value)),pch=19,cex=0.01) +
  geom_point(data = subset(plt, !is.na(FDR.p.value)),pch=19,cex=0.35)
#ggrepel::geom_text_repel(col="black",size=3)



plt <- stats.ac |> 
  dplyr::left_join(data.probes.alzheimer, by=c('probe_id'='probe_id'),suffix=c('',''))


# ggplot(plt, aes(x=logFC, y=-adj.P.Val, col=is.na(FDR.p.value))) +
#   geom_point(data = subset(plt, is.na(FDR.p.value)),pch=19,cex=0.01) +
#   geom_point(data = subset(plt, !is.na(FDR.p.value)),pch=19,cex=0.15)


# yes - concordance
ggplot(plt, aes(x=Beta..difference, y=logFC, col=is.na(FDR.p.value), label=label)) +
  geom_point(data = subset(plt, is.na(FDR.p.value)),pch=19,cex=0.01) +
  geom_point(data = subset(plt, !is.na(FDR.p.value)),pch=19,cex=0.35)
#ggrepel::geom_text_repel(col="black",size=3)






### subtests HG_AC specific ----


metadata <- rbind( metadata.od , metadata.ac) |> 
  dplyr::mutate(HG_OD = factor(ifelse(dataset == 'GLASS-OD' & LG_HG_status == "HG","HG_OD","non_HG_OD"),levels=c('non_HG_OD','HG_OD'))) |> 
  dplyr::mutate(HG_AC = factor(ifelse(dataset == 'GLASS-NL' & LG_HG_status == "HG","HG_AC","non_HG_AC"),levels=c('non_HG_AC','HG_AC'))) |> 
  dplyr::mutate(HG_stat = factor(dplyr::case_when(
    HG_OD == "HG_OD" ~ "HG_OD",
    HG_AC == "HG_AC" ~ "HG_AC",
    T ~ "LG"
  ), levels=c("LG","HG_AC","HG_OD")))

data <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata$sentrix_id)


#### A HG_AC -1 - HG_OD 1 ----

# HG_stat of HG_OD + HG_AC maakt niet veel uit
#model.matrix(~HG_OD + HG_AC, data=metadata) |> head()
#model.matrix(~HG_stat, data=metadata) |> head()


design <- model.matrix(~factor(patient_id) + factor(HG_OD) + factor(HG_AC), data=metadata)
fit <- limma::eBayes(limma::lmFit(data, design),trend=T)

#limma::makeContrasts(HG_OD - HG_AC - LG, levels=colnames(design))
contrast <- limma::makeContrasts(HG_statHG_OD - HG_statHG_AC - HG_statLG, levels=colnames(design))
fit2 <- limma::eBayes( limma::contrasts.fit( fit, contrast ) )
# 
# HG_statHG_OD - HG_statHG_AC
f = rownames(fit2$p.value)[fit2$p.value < 0.01]


#### B HG_AC 1 & HG_OD 1 ----

# HG_stat of HG_OD + HG_AC maakt niet veel uit
#model.matrix(~HG_OD + HG_AC, data=metadata) |> head()
#model.matrix(~HG_stat, data=metadata) |> head()


design <- model.matrix(~factor(patient_id) + factor(dataset) + factor(LG_HG_status), data=metadata)
fit <- limma::eBayes(limma::lmFit(data, design),trend=T)
stats.hg <- limma::topTable(fit,
                                 n=nrow(data),
                                 coef="factor(LG_HG_status)HG",
                                 sort.by = "none",
                                 adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')



design <- model.matrix(~factor(patient_id) + factor(HG_OD) + factor(HG_AC), data=metadata)
fit <- limma::eBayes(limma::lmFit(data, design),trend=T)
stats.a.hg_ac <- limma::topTable(fit,
                                 n=nrow(data),
                                 coef="factor(HG_AC)HG_AC",
                                 sort.by = "none",
                                 adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')
stats.a.hg_od <- limma::topTable(fit,
                                 n=nrow(data),
                                 coef="factor(HG_OD)HG_OD",
                                 sort.by = "none",
                                 adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')



design <- model.matrix(~factor(patient_id) + factor(LG_HG_status) + factor(HG_AC), data=metadata)
fit <- limma::eBayes(limma::lmFit(data, design),trend=T)
stats.b.hg_ac <- limma::topTable(fit,
                                 n=nrow(data),
                                 coef="factor(HG_AC)HG_AC",
                                 sort.by = "none",
                                 adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')
#plot(plt$t.ac, stats.b.hg_ac$t, pch=19,cex=0.01)


design <- model.matrix(~factor(patient_id) + factor(LG_HG_status) + factor(HG_OD), data=metadata)
fit <- limma::eBayes(limma::lmFit(data, design),trend=T)
stats.b.hg_od <- limma::topTable(fit,
                            n=nrow(data),
                            coef="factor(HG_OD)HG_OD",
                            sort.by = "none",
                            adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')
#plot(plt$t.ac, stats.b.hg_ac$t, pch=19,cex=0.01)
#plot(plt$t.od, stats.b.hg_od$t, pch=19,cex=0.01)


design <- model.matrix(~factor(patient_id) + factor(HG_AC), data=metadata)
fit <- limma::eBayes(limma::lmFit(data, design),trend=T)
stats.c.hg_ac <- limma::topTable(fit,
                                 n=nrow(data),
                                 coef="factor(HG_AC)HG_AC",
                                 sort.by = "none",
                                 adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')


design <- model.matrix(~factor(patient_id) + factor(HG_OD), data=metadata)
fit <- limma::eBayes(limma::lmFit(data, design),trend=T)
stats.c.hg_od <- limma::topTable(fit,
                                 n=nrow(data),
                                 coef="factor(HG_OD)HG_OD",
                                 sort.by = "none",
                                 adjust.method="fdr") |>
  tibble::rownames_to_column('probe_id')



plt <- stats.hg |>  dplyr::rename_with(~ paste0(.x, ".hg")) |>
  dplyr::left_join(
    stats.a.hg_ac |> dplyr::rename_with(~ paste0(.x, ".a.hg_ac")),
    by = c("probe_id.hg" = "probe_id.a.hg_ac")
  ) |>
  dplyr::left_join(
    stats.a.hg_od |> dplyr::rename_with(~ paste0(.x, ".a.hg_od")),
    by = c("probe_id.hg" = "probe_id.a.hg_od")
  ) |>
  dplyr::left_join(
    stats.b.hg_ac |> dplyr::rename_with(~ paste0(.x, ".b.hg_ac")),
    by = c("probe_id.hg" = "probe_id.b.hg_ac")
  ) |>
  dplyr::left_join(
    stats.b.hg_od |> dplyr::rename_with(~ paste0(.x, ".b.hg_od")),
    by = c("probe_id.hg" = "probe_id.b.hg_od")
  ) |>
  dplyr::left_join(
    stats.c.hg_ac |> dplyr::rename_with(~ paste0(.x, ".c.hg_ac")),
    by = c("probe_id.hg" = "probe_id.c.hg_ac")
  ) |>
  dplyr::left_join(
    stats.c.hg_od |> dplyr::rename_with(~ paste0(.x, ".c.hg_od")),
    by = c("probe_id.hg" = "probe_id.c.hg_od")
  )


c = plt |>
  dplyr::select(probe_id.hg, starts_with("t.")) |> 
  tibble::column_to_rownames('probe_id.hg') |> 
  as.matrix() |> 
  cor()

corrplot::corrplot(c, order="hclust")


sum(stats.hg$adj.P.Val < 0.01)

sum(stats.a.hg_od$adj.P.Val < 0.01)
sum(stats.a.hg_ac$adj.P.Val < 0.01)
plot(stats.a.hg_od$t, stats.a.hg_ac$t, pch=19,cex=0.01)

sum(stats.b.hg_od$adj.P.Val < 0.01)
sum(stats.b.hg_ac$adj.P.Val < 0.01)
plot(stats.b.hg_od$t, stats.b.hg_ac$t, pch=19,cex=0.01)

sum(stats.c.hg_od$adj.P.Val < 0.01)
sum(stats.c.hg_ac$adj.P.Val < 0.01)
plot(stats.c.hg_od$t, stats.c.hg_ac$t, pch=19,cex=0.01)


plt <- plt |> 
  dplyr::mutate(col1 = ifelse(adj.P.Val.a.hg_ac < 0.01, "signi AC", "AC-n")) |> 
  dplyr::mutate(col2 = ifelse(adj.P.Val.a.hg_od < 0.01, "signi OD", "OD-n")) |> 
  #dplyr::mutate(col1 = ifelse(P.Value.c.hg_ac < 0.01, "signi AC-only", "AC-n")) |> 
  #dplyr::mutate(col2 = ifelse(P.Value.c.hg_od < 0.01, "signi OD-only", "OD-n")) |> 
  dplyr::mutate(col3 = ifelse(adj.P.Val.b.hg_ac < 0.01, "signi AC-only", "AC-n")) |> 
  dplyr::mutate(col4 = ifelse(adj.P.Val.b.hg_od < 0.01, "signi OD-only", "OD-n")) |> 
  dplyr::mutate(col5 = paste0(col1, " - ", col2, " - ", col3, " - ", col4))


ggplot(plt, aes(x=t.a.hg_ac, y=t.a.hg_od, col=col5)) +
  geom_point(pch=19,cex=0.01,alpha=0.5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  xlim(-12,12) +
  ylim(-12,12) +
  theme_bw()

ggplot(plt, aes(x=t.b.hg_ac, y=t.b.hg_od, col=col5)) +
  geom_point(pch=19,cex=0.01,alpha=0.5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  xlim(-12,12) +
  ylim(-12,12) +
  theme_bw()

ggplot(plt, aes(x=t.c.hg_ac, y=t.c.hg_od, col=col5)) +
  geom_point(pch=19,cex=0.01,alpha=0.5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  xlim(-12,12) +
  ylim(-12,12) +
  theme_bw()


