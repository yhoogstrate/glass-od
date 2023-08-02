#!/usr/bin/env R

# load data ----


library(ggplot2)


if(!exists('glass_nl.metadata.resections') | !exists('glass_od.metadata.resections')) {
  source('scripts/load_metadata.R')
}


if(!exists('glass_nl.data.mvalues') | !exists('glass_od.data.mvalues')) {
  source('scripts/load_mvalues.R')
}




# merge ODG + AC with all HQ-probes ----


metadata.glass_od <- glass_od.metadata.idats |> 
  filter_GLASS_OD_idats(163) 

metadata.glass_nl <- glass_nl.metadata.idats |> 
  filter_GLASS_NL_idats(218)


metadata.combi <- data.frame(sentrix_id = c(metadata.glass_od$sentrix_id, metadata.glass_nl$sentrix_id)) |> 
  dplyr::left_join(metadata.glass_od |> 
                     dplyr::select(sentrix_id,
                                   resection_tumor_grade,
                                   
                                   mnp_predictBrain_v12.5_cal_class,
                                   mnp_predictBrain_v12.8_cal_class,
                                   mnp_predictBrain_v2.0.1_cal_class,
                                   
                                   mnp_rsGender_11b4_predicted,
                                   mnp_rsGender_12.8_predicted,
                                   
                                   resection_number,
                                   resection_id,
                                   patient_id,
                                   methylation_bins_1p19q_purity,
                                   
                                   isolation_person_name
                                   ),
                   by=c('sentrix_id' = 'sentrix_id')) |> 
  dplyr::left_join(metadata.glass_nl |> 
                     dplyr::select(sentrix_id, 
                                   Sample_Type, 
                                   Sample_Name, 
                                   Sample_Sex,
                                   
                                   #mnp_predictBrain_v12.5_cal_class,
                                   mnp_predictBrain_v12.8_cal_class
                                   #mnp_predictBrain_v2.0.1_cal_class
                                   
                                   ) , by=c('sentrix_id'='sentrix_id'), suffix=c('_od','_nl')) 

rm(metadata.glass_nl, metadata.glass_od)
gc()


tmp <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_GLASS-NL_combined.Rds") |> 
  assertr::verify(metadata$sentrix_id %in% sentrix_id)

metadata <- metadata |> 
  dplyr::left_join(tmp, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))

rm(tmp)




# 
# 
# plt.split <- rbind(plt |> dplyr::mutate(facet = "dataset", col=dataset)
#                    ,
#                    plt |> dplyr::mutate(facet = "grade", col=grade)
#                    ,
#                    plt |> dplyr::mutate(facet = "batch", col=batch)
#                    )
# 
# ggplot(plt.split, aes(x=PC4, y=PC3, col=col)) +
#   facet_grid(cols = vars(facet)) +
#   geom_point() +
#   theme_bw()
# 




### MDS over PCA ----

dim(data.combi)

d <- dist(pc$x[,1:50])
dim(d)

fit <- cmdscale(d, eig=T, k=2) # first 50 components?
fit


stopifnot(colnames(data.combi) == metadata.combi$sentrix_id)


plt <- metadata.combi |> 
  dplyr::left_join(
    fit$points |> 
      as.data.frame() |>
      tibble::rownames_to_column('sentrix_id') |> 
      dplyr::rename(MDS_PCA_coord1 = V1) |> 
      dplyr::rename(MDS_PCA_coord2 = V2)
    ,
    by=c('sentrix_id' = 'sentrix_id')
  ) |> 
  dplyr::left_join(
    pc |> purrr::pluck('x') |> as.data.frame() |> tibble::rownames_to_column('sentrix_id'),
    by=c('sentrix_id' = 'sentrix_id')
  ) |> 
  dplyr::mutate(dataset_x = ifelse(dataset == "GLASS-OD" & predictBrain_12.5_cal_class == "A_IDH_HG","GLASS-OD [misclass A_IDH_HG]", dataset))


# PC1 = quality [at least in glass od]
# PC2 = purity  [at least in glass od]

ggplot(plt, aes(x=PC1, y=PC2, col=puritygroups)) +
  geom_point() +
  theme_bw()


ggplot(plt, aes(x=PC1, y=PC2, col=puritygroups)) +
  geom_point() +
  theme_bw()


ggplot(plt, aes(x=PC10, y=PC8, col=sex)) +
  geom_point() +
  theme_bw()


relevant_components <-   (1:50) |> purrr::discard(function( xx ){ 
  return (xx %in% c(1,2,8,10) == T) 
}  )


ggplot(plt  |> dplyr::filter(grepl("147|171", Sample_Name)), aes(x=PC5, y=PC3, col=dataset_x, label=Sample_Name)) +
  geom_point() +
  theme_bw() + 
  geom_text(size = 2.5,col="black")


ggplot(plt, aes(x=PC10, y=PC17, col=dataset_x)) +
  geom_point() +
  theme_bw()




ggplot(plt, aes(x=MDS_PCA_coord1, y=MDS_PCA_coord2, col=grade)) +
  geom_point() +
  theme_bw()




ggplot(plt, aes(x=PC2, y=PC45, col=dataset)) + 
  geom_point() +
  theme_bw()


for(ppc in paste0("PC",1:50)) {
  w = wilcox.test(
    plt |> dplyr::filter(dataset == "GLASS-OD") |> dplyr::pull(paste0(ppc)),
    plt |> dplyr::filter(dataset != "GLASS-OD") |> dplyr::pull(paste0(ppc))
  )
  
  print(paste0(ppc, " -- ",w$p.value))
}


# find sex pc's: PC10 , PC8
for(ppc in paste0("PC",1:50)) {
  w = wilcox.test(
    plt |> dplyr::filter(sex == "F") |> dplyr::pull(paste0(ppc)),
    plt |> dplyr::filter(sex == "M") |> dplyr::pull(paste0(ppc))
  )
  
  print(paste0(ppc, " -- ",w$p.value))
}




plt |>
  dplyr::filter(dataset == "GLASS-OD") |>
  dplyr::pull(paste0("PC",ppc))


plt |>
  dplyr::filter(dataset == "GLASS-OD") |>
  dplyr::pull("PC1") |> 
  head()



### check probes in PC1 related to 1p ~ 19q

loading.pc1 <- pc$rotation |>
  as.data.frame() |> 
  dplyr::select('PC1') |> 
  dplyr::rename(PC1.loadings = PC1) |> 
  tibble::rownames_to_column("probeID")


plt <- read.delim("~/mnt/neuro-genomic-1-ro/catnon/Methylation - EPIC arrays/EPIC.hg38.manifest.tsv.gz") |>
  dplyr::filter(MASK_general == F) |> 
  dplyr::mutate(pos = round((CpG_beg + CpG_end )/2)) |> 
  dplyr::select(CpG_chrm, pos, probeID, probe_strand  ) |> 
  dplyr::left_join(loading.pc1, by=c('probeID' = 'probeID'), suffix=c('','')) |> 
  dplyr::rename(chr = CpG_chrm) |> 
  dplyr::mutate(chr = factor(chr, levels=gtools::mixedsort(unique(as.character(chr))) )) |> 
  dplyr::filter(chr %in% c("chrM","chrX","chrY", "NA",NA) == F)




ggplot(plt, aes(x = pos/1000000, y=PC1.loadings)) +
  facet_grid(cols = vars(chr), scales = "free", space="free") +
  geom_point(size=0.25,alpha=0.15) +
  geom_smooth(col="red") +
  theme_bw()



# ggplot(plt |> dplyr::filter(chr %in% c("chr3"), aes(x = pos/1000000, y=PC1.loadings)) +
#   facet_grid(cols = vars(chr), scales = "free", space="free") +
#   geom_point(size=0.25,alpha=0.15) +
#   geom_smooth(col="red") +
#   theme_bw()






### UMAP over PCA ----



plt$UMAP1 = NULL
plt$UMAP2 = NULL
rm(um, umd)



um <- M3C::umap(t(pc$x[,relevant_components]), seed=round(runif(1) * 10000))

umd <- um$data |> 
  as.data.frame() |> 
  dplyr::rename(UMAP1 = 1) |> 
  dplyr::rename(UMAP2 = 2) |> 
  tibble::rownames_to_column('sentrix_id')


plt <- plt |> 
  dplyr::left_join(umd, by=c('sentrix_id'='sentrix_id'),suffix=c('','')) |> 
  dplyr::mutate(puritygroups= cut(methylation_bins_1p19q_purity, breaks=4))


ggplot(plt, aes(x=UMAP1, y=UMAP2, col=dataset_x, label=resection_id)) +
  geom_point() +
  #geom_text(size=3, col="black") +
  theme_bw()



### tSNE over PCA ----


ts <- M3C::tsne(t(pc$x[,relevant_components]),perplex=10)#1 and 2 are purity and quality
rownames(ts$data) <- rownames(pc$x)

tsd <- ts$data |> 
  as.data.frame() |> 
  dplyr::rename(tSNE1 = 1) |> 
  dplyr::rename(tSNE2 = 2) |> 
  tibble::rownames_to_column('sentrix_id')


plt <- plt |> 
  dplyr::left_join(tsd, by=c('sentrix_id'='sentrix_id'),suffix=c('',''))


ggplot(plt |> dplyr::mutate(resection_number = as.factor(resection_number)), 
       aes(x=tSNE1, y=tSNE2, col=dataset_x, group=patient_id)) +
  geom_point() +
  #geom_line(data = plt |> dplyr::filter(dataset == "GLASS-OD")) +
  theme_bw()



## raw UMAP ----

# see if and how well this works given the feature size


## recursiveCorPlot ----





recursiveCorPlot::recursiveCorPlot(
  as.data.frame(pc$x[,relevant_components[1:25]]),
  metadata.combi |> tibble::column_to_rownames('sentrix_id') |>  dplyr::select(dataset, predictBrain_12.5_cal_class) |> 
    dplyr::mutate(`GLASS-OD` = dataset == "GLASS-OD") |> 
    dplyr::mutate(`GLASS-NL` = dataset == "GLASS-NL") |> 
    dplyr::mutate(`GLASS-OD [misclass A_IDH_HG]` = 
                    (dataset == "GLASS-OD" & 
                    `predictBrain_12.5_cal_class` == "A_IDH_HG")) |> 
    dplyr::mutate(dataset = NULL, predictBrain_12.5_cal_class = NULL)
    ,
  3,
  3
)


# merge ODG + AC with all HQ-probes not 1P / 19Q ----


probes.1P19Q <- read.delim("~/mnt/neuro-genomic-1-ro/catnon/Methylation - EPIC arrays/EPIC.hg38.manifest.tsv.gz") |>
  dplyr::filter(MASK_general == F) |> 
  dplyr::mutate(pos = round((CpG_beg + CpG_end )/2)) |> 
  dplyr::mutate(is_1P = CpG_chrm == 'chr1' & pos < 130 * 1000000) |> # rough margin
  dplyr::mutate(is_19Q = CpG_chrm == 'chr19' & pos > 23.5 * 1000000 ) |> # rough margin
  dplyr::filter(is_1P | is_19Q)


data.combi.non1P19Q <- data.combi |> 
  tibble::rownames_to_column('probeID') |> 
  dplyr::filter(probeID %in% probes.1P19Q$probeID == F) |> 
  tibble::column_to_rownames('probeID')



plt <- metadata.combi |> 
  dplyr::mutate(dataset_x = ifelse(dataset == "GLASS-OD" & predictBrain_12.5_cal_class == "A_IDH_HG","GLASS-OD [misclass A_IDH_HG]", dataset)) |> 
  dplyr::mutate(puritygroups= cut(methylation_bins_1p19q_purity, breaks=4))





## PCA ----

pc.non1P19Q <- prcomp(t(data.combi.non1P19Q))
dim(pc.non1P19Q$x)


relevant_components <-   (1:50) |> purrr::discard(function( xx ){ 
  return (xx %in% c(1,5,7,8) == T)  # grading is relevant , purity is not,
}  )



plt <- plt |> 
  dplyr::left_join(
    pc.non1P19Q |> purrr::pluck('x') |> as.data.frame() |> tibble::rownames_to_column('sentrix_id'),
    by=c('sentrix_id' = 'sentrix_id')
  )


# PC2 = purity
ggplot(plt, aes(x=PC1, y=PC2, col=puritygroups)) +
  geom_point() +
  theme_bw()

# PC4 = tumor grade 1p/19q?
ggplot(plt, aes(x=PC1, y=PC4, col=resection_tumor_grade)) +
  geom_point() +
  theme_bw()


# grade = 4
for(ppc in paste0("PC",1:50)) {
  w = wilcox.test(
    plt |> dplyr::filter(resection_tumor_grade == "2") |> dplyr::pull(paste0(ppc)),
    plt |> dplyr::filter(resection_tumor_grade == "3") |> dplyr::pull(paste0(ppc))
  )
  
  print(paste0(ppc, " -- ",w$p.value))
}



# find sex pc's: PC7 , PC8
for(ppc in paste0("PC",1:50)) {
  w = wilcox.test(
    plt |> dplyr::filter(sex == "F") |> dplyr::pull(paste0(ppc)),
    plt |> dplyr::filter(sex == "M") |> dplyr::pull(paste0(ppc))
  )
  
  print(paste0(ppc, " -- ",w$p.value))
}



### UMAP over PCA ----



plt$UMAP1 = NULL
plt$UMAP2 = NULL
rm(um, umd)



um <- M3C::umap(t(pc.non1P19Q$x[,relevant_components[1:25]]), seed=round(runif(1) * 10000))

umd <- um$data |> 
  as.data.frame() |> 
  dplyr::rename(UMAP1 = 1) |> 
  dplyr::rename(UMAP2 = 2) |> 
  tibble::rownames_to_column('sentrix_id')


plt <- plt |> 
  dplyr::left_join(umd, by=c('sentrix_id'='sentrix_id'),suffix=c('','')) 


ggplot(plt, aes(x=UMAP1, y=UMAP2, col=dataset_x, label=resection_id)) +
  geom_point() +
  #geom_text(size=3, col="black") +
  theme_bw()


### tSNE over PCA ----


ts <- M3C::tsne(t(pc.non1P19Q$x[,relevant_components[1:25]]),perplex=10)#1 and 2 are purity and quality
rownames(ts$data) <- rownames(pc.non1P19Q$x)

tsd <- ts$data |> 
  as.data.frame() |> 
  dplyr::rename(tSNE1 = 1) |> 
  dplyr::rename(tSNE2 = 2) |> 
  tibble::rownames_to_column('sentrix_id')


plt <- plt |> 
  dplyr::left_join(tsd, by=c('sentrix_id'='sentrix_id'),suffix=c('',''))


ggplot(plt |> dplyr::mutate(resection_number = as.factor(resection_number)), 
       aes(x=tSNE1, y=tSNE2, col=dataset_x, group=patient_id)) +
  geom_point() +
  #geom_vline(xintercept=15, col="red") +
  #geom_line(data = plt |> dplyr::filter(dataset == "GLASS-OD"), col="gray",lwd=0.6) +
  theme_bw()





recursiveCorPlot::recursiveCorPlot(
  as.data.frame(pc.non1P19Q$x[,relevant_components[1:25]]),
  metadata.combi |> tibble::column_to_rownames('sentrix_id') |>  dplyr::select(dataset, predictBrain_12.5_cal_class) |> 
    dplyr::mutate(`GLASS-OD` = dataset == "GLASS-OD" & grepl("USA",plt$isolation_person_name) == F) |> 
    dplyr::mutate(`GLASS-OD [USA]` = dataset == "GLASS-OD" & grepl("USA",plt$isolation_person_name) == T) |> 
    dplyr::mutate(`GLASS-NL` = dataset == "GLASS-NL") |>
    dplyr::mutate(`GLASS-OD [misclass A_IDH_HG]` = 
                    (dataset == "GLASS-OD" & 
                       `predictBrain_12.5_cal_class` == "A_IDH_HG")) |> 
    dplyr::mutate(dataset = NULL, predictBrain_12.5_cal_class = NULL)
  ,
  9,
  11
)



### different codel clusters ----
#' does not seem perfectly patient specific, but quite much

gA <- c(
  "201496850071_R02C01","203293640061_R08C01","205828590003_R02C01","203986510092_R01C01","203989100149_R04C01",
  "203989100149_R05C01","203989100149_R06C01","203989100149_R07C01","203989100149_R08C01","203990170017_R01C01",
  "204073520033_R02C01","204073520033_R03C01","204073520033_R04C01","204073520033_R05C01","204088040040_R02C01",
  "204088040075_R07C01","204798720018_R03C01","204787070003_R01C01","204787070003_R02C01","204787070003_R03C01",
  "204787070003_R04C01","204787070003_R05C01","204787070003_R06C01","204787070003_R08C01","204787070018_R01C01",
  "204787070018_R02C01","204787070018_R03C01","204787070018_R04C01","204787070018_R07C01","204808700073_R02C01",
  "204808700073_R03C01","204808700073_R04C01","204808700073_R05C01","204808700073_R06C01","204808700073_R08C01",
  "204808700074_R01C01","204808700074_R02C01","204808700074_R03C01","204808700074_R05C01","206467010147_R05C01",
  "206467010147_R07C01","206522890026_R02C01","206116800026_R01C01","206116800026_R02C01","206116800026_R03C01",
  "206116800026_R04C01","206116800026_R06C01","206116800026_R08C01","206116800035_R02C01","206116800035_R03C01",
  "206116800035_R04C01","206116800035_R05C01","206116800035_R06C01","206116800035_R07C01","206116800035_R08C01",
  "206116800056_R01C01","206116800056_R03C01","206116800056_R05C01","206116800056_R06C01","206116800056_R07C01",
  "206116800056_R08C01","206116800060_R01C01","206116800060_R02C01","206116800060_R03C01","206116800060_R04C01",
  "206116800060_R05C01","206116800060_R07C01","206116800060_R08C01","206116800140_R01C01","206116800140_R03C01",
  "206116800140_R04C01","206116800140_R05C01","206116800140_R06C01","206116800140_R07C01","206116800140_R08C01",
  "206119350032_R06C01","206119350032_R07C01","206119350032_R08C01","206119350033_R02C01","206119350033_R08C01",
  "206119350042_R01C01","206119350042_R02C01","206119350042_R05C01","206119350042_R07C01","206119350042_R08C01",
  "206137490041_R01C01","206137490041_R03C01","206137490041_R04C01","206137490041_R05C01","206137490041_R06C01",
  "206137490041_R07C01","206137490041_R08C01","206137490053_R01C01","206137490053_R06C01","206137490053_R07C01",
  "206137490053_R08C01","206137490057_R02C01","206137490057_R03C01","206137490057_R06C01","206137490057_R07C01",
  "206137490057_R08C01","206137490109_R01C01","206137490109_R02C01","206137490109_R04C01","206137490109_R05C01",
  "206137490109_R06C01","206137490109_R08C01","206467010068_R01C01","206467010068_R04C01","206467010068_R07C01",
  "206467010068_R08C01"
)

gB <- c(
  "205832320037_R07C01","206238130171_R03C01","204808700073_R01C01","206467010069_R01C01","206467010069_R04C01",
  "206467010069_R05C01","206467010069_R06C01","206467010069_R07C01","206467010069_R08C01","206467010089_R01C01",
  "206467010089_R04C01","206467010089_R05C01","206467010089_R06C01","206467010089_R07C01","206467010089_R08C01",
  "206467010147_R01C01","206467010147_R04C01","206467010180_R01C01","206467010180_R02C01","206467010180_R03C01",
  "206467010180_R04C01","206467010180_R05C01","206467010180_R06C01","206467010180_R07C01","206467010180_R08C01",
  "206522890026_R07C01","206522890026_R08C01","206522890076_R01C01","206522890076_R02C01","206522890076_R06C01",
  "206522890076_R07C01","206522890076_R08C01","206522890079_R01C01","206522890079_R02C01","206522890079_R03C01",
  "206522890079_R04C01","206522890079_R05C01","206522890079_R06C01","206522890079_R07C01","206522890079_R08C01",
  "206522890089_R01C01","206522890089_R02C01","206522890089_R03C01","206522890089_R04C01","206522890089_R05C01",
  "206522890089_R06C01","206522890089_R07C01","206522890089_R08C01","206119350042_R06C01","206467110176_R04C01",
  "206467110176_R05C01","206467110176_R06C01"
)


for(ppc in paste0("PC",1:50)) {
  w = wilcox.test(
    plt |> dplyr::filter(sentrix_id %in% gA) |> dplyr::pull(paste0(ppc)),
    plt |> dplyr::filter(sentrix_id %in% gB) |> dplyr::pull(paste0(ppc))
  )
  
  print(paste0(ppc, " -- ",w$p.value))
}


plt <- plt |> 
  dplyr::mutate(gGLOD = dplyr::case_when(
    sentrix_id %in% gA ~ "gA",
    sentrix_id %in% gB ~ "gB",
    T ~ as.character(NA),
  )) |> 
  dplyr::mutate(gGLOD = as.factor(gGLOD))


#' "PC1 -- 4.56532945081721e-15"
#' "PC2 -- 1.61787568692805e-17"
#' "PC3 -- 3.09760839505662e-16"
#' "PC5 -- 7.17639737347061e-24"
#' "PC7 -- 1.72055518902603e-08"
#' "PC10 -- 7.07185312968851e-06"
#' "PC11 -- 4.85890868287798e-08"
#' "PC25 -- 5.52534544987418e-05"

ggplot(plt |> dplyr::filter(dataset == "GLASS-OD"), aes(x=PC3, y=PC5, col=gsub("_.+","",isolation_person_name)))  +
  geom_point()





