#!/usr/bin/env R


source("scripts/load_functions.R")


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}

source("scripts/load_themes.R")
source("scripts/load_palette.R")




# EPIC ... 387-plate1 -----

# '/home/youri/mnt/neuro-genomic-1-ro/glass/GLASS_OD/DNA Methylation - EPIC arrays/EPIC2023-387-plate1'


tmp.glod <- list.files(path = "data/GLASS_OD/DNA Methylation - EPIC arrays/EPIC2023-387-plate1/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(array_filename = _) |> 
  dplyr::mutate(array_filename = paste0("data/GLASS_OD/DNA Methylation - EPIC arrays/EPIC2023-387-plate1/", array_filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (2 * 53))
    return(.)
  })()


tmp.other <- list.files(path = "/home/youri/mnt/neuro-genomic-1-ro/DNAm_arrays_overig/DNA Methylation - EPIC arrays/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(array_filename = _) |> 
  dplyr::mutate(array_filename = paste0("/home/youri/mnt/neuro-genomic-1-ro/DNAm_arrays_overig/DNA Methylation - EPIC arrays/", array_filename)) |> 
  dplyr::filter(!grepl("plate2", array_filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (2 * 4))
    return(.)
  })()


tmp.gsam <- list.files(path = "/home/youri/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays/EPIC2023-387-plate1/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(array_filename = _) |> 
  dplyr::mutate(array_filename = paste0("/home/youri/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays/EPIC2023-387-plate1/", array_filename)) |> 
  dplyr::filter(!grepl("plate2", array_filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (2 * 4))
    return(.)
  })()


tmp.mint <- list.files(path = "/home/youri/mnt/neuro-genomic-1-ro/MINT/Methylation_data/DNA Methylation - EPIC arrays/EPIC2023-387-plate1/", pattern = "_(Grn|Red).idat$", recursive = TRUE) |>
  data.frame(array_filename = _) |> 
  dplyr::mutate(array_filename = paste0("/home/youri/mnt/neuro-genomic-1-ro/MINT/Methylation_data/DNA Methylation - EPIC arrays/EPIC2023-387-plate1/", array_filename)) |> 
  dplyr::filter(!grepl("plate2", array_filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (2 * 3))
    return(.)
  })()




tmp <- rbind(tmp.glod, tmp.other, tmp.gsam, tmp.mint) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (2 * 64))
    return(.)
  })() |> 
  assertr::verify(file.exists(array_filename)) |>
  dplyr::mutate(array_sentrix_id = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", array_filename)) |>
  dplyr::mutate(array_sentrix_id = gsub("GSM[0-9]+_", "", array_sentrix_id)) |>
  dplyr::mutate(array_channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", array_filename)) |>
  tidyr::pivot_wider(id_cols = array_sentrix_id, names_from = array_channel, values_from = c(array_filename)) |>
  dplyr::rename(array_channel_green = Grn) |>
  dplyr::rename(array_channel_red = Red) |> 
  dplyr::mutate(
   resection_id = dplyr::case_when(
      array_sentrix_id == "207356840076_R02C01" ~ "110-R1",
      array_sentrix_id == "207356840171_R05C01" ~ "110-R2",
      array_sentrix_id == "207513910100_R04C01" ~ "108_02",
      array_sentrix_id == "207513910100_R01C01" ~ "108_96",
      array_sentrix_id == "207513910076_R08C01" ~ "108_03",
      
      T ~ array_sentrix_id
    )
  ) |> 
  dplyr::arrange(resection_id) |> 
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL, array_channel_red = NULL) |> 
  dplyr::mutate(Sample_Name = paste0(array_sentrix_id)) |>
  dplyr::mutate(Array = gsub("^.+_","", array_sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1", array_sentrix_id)) |> 
  dplyr::select(array_sentrix_id, resection_id,  Basename, Sample_Name, Array, Slide)


#RGSet <- minfi::read.metharray.exp(targets = tmp, force = T) #red/green channel together





s <- minfi::getSnpBeta(RGSet) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::left_join(tmp |> dplyr::select(array_sentrix_id, resection_id), by=c('array_sentrix_id'='array_sentrix_id')) |> 
  tibble::column_to_rownames('resection_id') |> 
  dplyr::mutate(array_sentrix_id = NULL)

d <- s |> 
  dist(diag=T, upper=T)


distplot(d) +
  theme(plot.background = element_rect(fill="white", colour=NA))




c <- s |> 
  t() |> 
  cor(method="spearman")

ggcorrplot(c, abs=F, reorder=T) +
  theme(plot.background = element_rect(fill="white", colour=NA))





# 0110 etc. ----

corrplot::corrplot(cor(t(d)))
func_ggcorrplot(cor(t(d)))


dat <- glass_od.metadata.array_samples |> 
  dplyr::filter(grepl("(0110|0064|0065|0066|0067|0108|0110)", resection_id)) |> 
  dplyr::arrange(resection_id) |> 
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL, array_channel_red = NULL) |> 
  dplyr::mutate(Sample_Name = paste0(array_sentrix_id)) |>
  dplyr::mutate(Array = gsub("^.+_","", array_sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1", array_sentrix_id)) |> 
  #filter_GLASS_OD_idats(14) |> 
  dplyr::select(array_sentrix_id, resection_id, Basename, Sample_Name, Array, Slide, patient_study_name, array_qc.pca.detP.outlier)




RGSet <- minfi::read.metharray.exp(targets = dat, force = T) #red/green channel together



s <- minfi::getSnpBeta(RGSet) |> 
  scale(center=F) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::left_join(dat |> dplyr::select(array_sentrix_id, resection_id, array_qc.pca.detP.outlier), by=c('array_sentrix_id'='array_sentrix_id')) |> 
  dplyr::mutate(resection_id = ifelse(array_qc.pca.detP.outlier, paste0(resection_id, "**"), resection_id)) |> 
  dplyr::mutate(array_qc.pca.detP.outlier = NULL) |> 
  tibble::column_to_rownames('resection_id') |> 
  dplyr::mutate(array_sentrix_id = NULL)

d <- s |> 
  dist(diag=T, upper=T)


fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS", type="p",pch=19,cex=1, col=grepl("0110|0108",row.names(s)) * 1 + 1)
#text(x, y, labels = row.names(s), cex=.7 , col=grepl("0110|0108",row.names(s)) * 1 + 1)


distplot(d)



plt <- d |> 
  tibble::rownames_to_column('resection_id') |> 
  tidyr::pivot_longer(cols=-c("resection_id"))

ggplot(plt, aes(x = resection_id, y=value)) +
  geom_point()


# 206467010147 [0065,0066] ----


dat <- glass_od.metadata.array_samples |> 
  dplyr::filter(grepl("206467010147", array_sentrix_id) ) |> 
  dplyr::arrange(resection_id) |> 
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL, array_channel_red = NULL) |> 
  dplyr::mutate(Sample_Name = paste0(array_sentrix_id)) |>
  dplyr::mutate(Array = gsub("^.+_","", array_sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1", array_sentrix_id)) |> 
  dplyr::select(array_sentrix_id, resection_id, Basename, Sample_Name, Array, Slide)

RGSet <- minfi::read.metharray.exp(targets = dat, force = T) #red/green channel together




d <- minfi::getSnpBeta(RGSet) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::left_join(dat |> dplyr::select(array_sentrix_id, resection_id), by=c('array_sentrix_id'='array_sentrix_id')) |> 
  tibble::column_to_rownames('resection_id') |> 
  dplyr::mutate(array_sentrix_id = NULL) |> 
  dist(diag=T, upper=T)


distplot(d)





# main dataset [all] ----


dat <- glass_od.metadata.array_samples |> 
  dplyr::filter(patient_study_name == "GLASS-OD" & arraychip_version == "EPICv1") |> 
  dplyr::arrange(resection_id) |> 
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL, array_channel_red = NULL) |> 
  dplyr::mutate(Sample_Name = paste0(array_sentrix_id)) |>
  dplyr::mutate(Array = gsub("^.+_","", array_sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1", array_sentrix_id)) |> 
  dplyr::select(array_sentrix_id, resection_id,  isolation_id, Basename, Sample_Name, Array, Slide)
#head(n=25)


RGSet <- minfi::read.metharray.exp(targets = dat, force = T) #red/green channel together





s <- minfi::getSnpBeta(RGSet) |> 
  scale(center=F) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::left_join(dat |> dplyr::select(array_sentrix_id, isolation_id), by=c('array_sentrix_id'='array_sentrix_id')) |> 
  tibble::column_to_rownames('isolation_id') |> 
  dplyr::mutate(array_sentrix_id = NULL)

d <- s |> 
  dist(diag=T, upper=T)


distplot(d) +
  theme(plot.background = element_rect(fill="white", colour=NA))


ggsave("output/figures/analysis_snp_fingerprinting__main_dataset.png", width=15, height=15)




c <- s |> 
  t() |> 
  cor(method="spearman")

ggcorrplot(c, abs=T, reorder=T) +
  theme(plot.background = element_rect(fill="white", colour=NA))

ggsave("output/figures/analysis_snp_fingerprinting__main_dataset_cor.png", width=15, height=15)



rm(dat)
rm(d)


# sd plot

sds = sapply(as.data.frame(t(s)), sd)
sds = sort(sds)

barplot(sds)



included <-  glass_od.metadata.array_samples |>
  filter_GLASS_OD_idats(212) |> 
  dplyr::pull(isolation_id)



plt <- data.frame(sds) |> 
  tibble::rownames_to_column('isolation_id') |> 
  dplyr::left_join(glass_od.metadata.array_samples, by=c('isolation_id'='isolation_id')) |> 
  dplyr::mutate(is_included = isolation_id %in% included) |> 
  dplyr::arrange(is_included, sds)


ggplot(plt, aes(x = reorder(isolation_id, sds), y=sds, fill=is_included)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5,hjust=1))




# main dataset [211 hq] ----


dat <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::arrange(resection_id) |> 
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL, array_channel_red = NULL) |> 
  dplyr::mutate(Sample_Name = paste0(array_sentrix_id)) |>
  dplyr::mutate(Array = gsub("^.+_","", array_sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1", array_sentrix_id)) |> 
  dplyr::select(array_sentrix_id, resection_id,  isolation_id, Basename, Sample_Name, Array, Slide)


RGSet <- minfi::read.metharray.exp(targets = dat, force = T) #red/green channel together



s <- minfi::getSnpBeta(RGSet) |> 
  scale(center=F) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::left_join(dat |> dplyr::select(array_sentrix_id, isolation_id), by=c('array_sentrix_id'='array_sentrix_id')) |> 
  tibble::column_to_rownames('isolation_id') |> 
  dplyr::mutate(array_sentrix_id = NULL)

d <- s |> 
  dist(diag=T, upper=T)


distplot(d) +
  theme(plot.background = element_rect(fill="white", colour=NA))


ggsave("output/figures/analysis_snp_fingerprinting__main_dataset_HQ.png", width=15, height=15,dpi=600)






c <- s |> 
  t() |> 
  cor(method="spearman")

ggcorrplot(c, abs=T, reorder=T) +
  theme(plot.background = element_rect(fill="white", colour=NA))


ggsave("output/figures/analysis_snp_fingerprinting__main_dataset_HQ_cor.png", width=15, height=15)



rm(dat)
rm(d)


# sd plot

sds = sapply(as.data.frame(t(s)), sd)
sds = sort(sds)

barplot(sds)



included <-  glass_od.metadata.array_samples |>
  filter_GLASS_OD_idats(212) |> 
  dplyr::pull(isolation_id)



plt <- data.frame(sds) |> 
  tibble::rownames_to_column('isolation_id') |> 
  dplyr::left_join(glass_od.metadata.array_samples, by=c('isolation_id'='isolation_id')) |> 
  dplyr::mutate(is_included = isolation_id %in% included) |> 
  dplyr::arrange(is_included, sds)


ggplot(plt, aes(x = reorder(isolation_id, sds), y=sds, fill=is_included)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5,hjust=1))




# baplot raw data example ----



dat <- glass_od.metadata.array_samples |> 
  dplyr::arrange(resection_id) |> 
  dplyr::filter(resection_id %in% c("0120-R1", "0108-R3","0065-R2")) |>
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL, array_channel_red = NULL) |> 
  dplyr::mutate(Sample_Name = paste0(array_sentrix_id)) |>
  dplyr::mutate(Array = gsub("^.+_","", array_sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1", array_sentrix_id))  |> 
  dplyr::mutate(resection_id = paste0(patient_id, "_", resection_id)) |> 
  dplyr::select(array_sentrix_id, resection_id, Basename, Sample_Name, Array, Slide)

RGSet <- minfi::read.metharray.exp(targets = dat, force = T) #red/green channel together


s <- minfi::getSnpBeta(RGSet) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::mutate(array_sentrix_id = gsub("GSM[0-9]+_","",array_sentrix_id))|> 
  dplyr::left_join(dat |> dplyr::select(array_sentrix_id, resection_id), by=c('array_sentrix_id'='array_sentrix_id')) |> 
  dplyr::arrange(resection_id) |> 
  tibble::column_to_rownames('resection_id') |> 
  dplyr::mutate(array_sentrix_id = NULL)




p1 <- data.frame(snp = s[1,]) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('snp')

p1 <- rbind(p1, p1 |> dplyr::mutate(`0065_0065-R2`=0.5))


p1 = ggplot(p1, aes(x=snp, y=`0065_0065-R2`, group=snp)) +
  geom_line() +
  ylim(0,1) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5,hjust=1))





p2 <- data.frame(snp = s[2,]) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('snp')

p2 <- rbind(p2, p2 |> dplyr::mutate(`0108_0108-R3`=0.5))


p2 = ggplot(p2, aes(x=snp, y=`0108_0108-R3`, group=snp)) +
  geom_line() +
  ylim(0,1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5,hjust=1))





p3 <- data.frame(snp = s[3,]) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('snp')

p3 <- rbind(p3, p3 |> dplyr::mutate(`0120_0120-R1`=0.5))


p3 = ggplot(p3, aes(x=snp, y=`0120_0120-R1`, group=snp)) +
  geom_line() +
  ylim(0,1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5,hjust=1))





library(patchwork)
p2 + p3





# main: X/Y plot ----
## complete data ----


dat |> 
  dplyr::filter(is.na(patient_sex)) |> 
  dplyr::pull(resection_id)



dat <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::mutate(col = ifelse(patient_id == "0124", resection_id, "other"))


ggplot(dat, aes(
  x=array_mnp_rsGender_12.8_dist,
  y=array_mnp_rsGender_12.8_chrYintensity,
  label=resection_id,
  #col=array_mnp_rsGender_11b4_predicted
  col=patient_sex
)) +
  geom_point() +
  geom_point(data= subset(dat, col != "other"), col="black") +
  ggrepel::geom_text_repel(data= subset(dat, col != "other"), col="black")




## analysis ----

dat <- glass_od.metadata.array_samples |> 
  dplyr::filter(patient_study_name == "GLASS-OD" & arraychip_version == "EPICv1") |> 
  dplyr::arrange(resection_id) |> 
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL, array_channel_red = NULL) |> 
  dplyr::mutate(Sample_Name = paste0(array_sentrix_id)) |>
  dplyr::mutate(Array = gsub("^.+_","", array_sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1", array_sentrix_id)) |> 
  dplyr::mutate(col = dplyr::case_when(
    grepl("0065|0066|0110", resection_id) ~ patient_id,
    T ~ "other"
  ))


dat |> 
  dplyr::filter(grepl("(0110|0064|0065|0066|0067|0108|0110)", resection_id)) |> 
  dplyr::select(array_sentrix_id, resection_id, Basename, Sample_Name, Array, Slide) |> 
  View()



ggplot(dat, aes(
  x=array_mnp_rsGender_12.8_dist,
  y=array_mnp_rsGender_12.8_chrYintensity,
  label=resection_id,
  #col=patient_sex,
  col=array_mnp_rsGender_11b4_predicted
)) +
  geom_point() +
  geom_point(data= subset(dat, col != "other"), col="black") +
  ggrepel::geom_text_repel(data= subset(dat, col != "other"), col="black")


# array_mnp_rsGender_12.8_AAwg
#  array_mnp_rsGender_12.8_ABwg"       
# "array_mnp_rsGender_12.8_BBwg"         
# "array_mnp_rsGender_12.8_dist
# "array_mnp_rsGender_12.8_AAchrX
# "array_mnp_rsGender_12.8_ABchrX
# "array_mnp_rsGender_12.8_BBchrX
# "array_mnp_rsGender_12.8_chrYintensity
# array_mnp_rsGender_12.8_predicted




# validation set ----



dat <- glass_od.metadata.array_samples |> 
  #filter_OD_validation_idats(CONST_N_OD_VALIDATION_INCLUDED_SAMPLES) |> 
  dplyr::filter(patient_study_name == "OD-validation" & arraychip_version == "EPICv1") |> 
  dplyr::arrange(resection_id) |> 
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL, array_channel_red = NULL) |> 
  dplyr::mutate(Sample_Name = paste0(array_sentrix_id)) |>
  dplyr::mutate(Array = gsub("^.+_","", array_sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1", array_sentrix_id))  |> 
  dplyr::mutate(resection_id = paste0(patient_id, "_", resection_id)) |> 
  dplyr::select(array_sentrix_id, resection_id, Basename, Sample_Name, Array, Slide)

RGSet <- minfi::read.metharray.exp(targets = dat, force = T) #red/green channel together


d <- minfi::getSnpBeta(RGSet) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::mutate(array_sentrix_id = gsub("GSM[0-9]+_","",array_sentrix_id))|> 
  dplyr::left_join(dat |> dplyr::select(array_sentrix_id, resection_id), by=c('array_sentrix_id'='array_sentrix_id')) |> 
  dplyr::arrange(resection_id) |> 
  tibble::column_to_rownames('resection_id') |> 
  dplyr::mutate(array_sentrix_id = NULL) |> 
  dist(diag=T, upper=T)

rm(dat)


distplot(d) +
  theme(plot.background = element_rect(fill="white", colour=NA))



ggsave("output/figures/analysis_snp_fingerprinting__validation_dataset.png", , width=9, height=9,dpi=600)


rm(d)



