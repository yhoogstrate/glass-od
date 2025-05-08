#!/usr/bin/env R

# load libs ----


library(ggplot2)
library(patchwork)


# load (meta)data ----


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}

if(!exists('glass_nl.metadata.array_samples')) {
  source('scripts/load_GLASS-NL_metadata.R')
}


metadata.glass_od <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(211) 

metadata.glass_nl <- glass_nl.metadata.array_samples |> 
  filter_GLASS_NL_idats(203)


metadata.combi <- data.frame(array_sentrix_id = c(metadata.glass_od$array_sentrix_id, 
                                                  metadata.glass_nl$array_sentrix_id)) |> 
  dplyr::left_join(metadata.glass_od |> 
                     dplyr::select(patient_id,
                                   
                                   resection_tumor_grade,
                                   resection_number,
                                   resection_id,
                                   
                                   isolation_person_name,
                                   
                                   array_sentrix_id,
                                   array_methylation_bins_1p19q_purity,
                                   array_mnp_predictBrain_v2.0.1_cal_class,
                                   array_mnp_predictBrain_v12.5_cal_class,
                                   array_mnp_predictBrain_v12.8_cal_class,


                                   array_methylation_bins_1p19q_purity
                     ) |> 
                     dplyr::rename(batch = isolation_person_name) |> 
                     dplyr::mutate(dataset = "GLASS-OD")
                   ,
                   by=c('array_sentrix_id' = 'array_sentrix_id')) |> 
  dplyr::left_join(metadata.glass_nl |> 
                     dplyr::select(array_sentrix_id, 
                                   Sample_Type, 
                                   Sample_Name, 
                                   Sample_Sex,
                                   
                                   #mnp_predictBrain_v12.5_cal_class,
                                   #mnp_predictBrain_v2.0.1_cal_class
                                   array_mnp_predictBrain_v12.8_cal_class,
                                   
                                   WHO_Classification2021
                     ) |> 
                     dplyr::mutate(dataset = "GLASS-NL")
                   , by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('_od','_nl'))  |> 
  dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = ifelse(!is.na(array_mnp_predictBrain_v12.8_cal_class_nl),array_mnp_predictBrain_v12.8_cal_class_nl,array_mnp_predictBrain_v12.8_cal_class_od),
                array_mnp_predictBrain_v12.8_cal_class_nl = NULL,
                array_mnp_predictBrain_v12.8_cal_class_od = NULL) |> 
  dplyr::mutate(dataset = ifelse(is.na(dataset_od), dataset_nl, dataset_od), dataset_nl = NULL, dataset_od = NULL) |> 
  dplyr::mutate(tumor_grade = dplyr::case_when(
    resection_tumor_grade == 3 ~ 3,
    resection_tumor_grade == 2 ~ 2,
    WHO_Classification2021 == "Astrocytoma, IDH-mutant, WHO grade 4" ~ 4,
    WHO_Classification2021 == "Astrocytoma, IDH-mutant, WHO grade 3" ~ 3,
    WHO_Classification2021 == "Astrocytoma, IDH-mutant, WHO grade 2" ~ 2,
    T ~ as.numeric(NA)
  ),WHO_Classification2021 = NULL, resection_tumor_grade = NULL) |> 
  dplyr::mutate(tumor_grade_h_l = dplyr::case_when(
    tumor_grade == 4 ~ "high grade",
    tumor_grade == 3 & dataset =="GLASS-OD" ~ "high grade",
    tumor_grade == 3 & dataset =="GLASS-NL" ~ "low grade",
    tumor_grade == 2 ~ "low grade",
    T ~ as.character(NA)
  )) |> 
  dplyr::mutate(batch_us = dplyr::case_when(
    is.na(batch) ~ "GLASS-NL [EU]",
    batch == "USA / Duke" ~ "GLASS-OD [batch US]",
    T ~ "GLASS-OD [batch EU]"
  )) |> 
  dplyr::left_join(
    rbind(
      metadata.glass_od |> 
        dplyr::select("array_sentrix_id", "array_mnp_predictBrain_v12.8_cal_class", starts_with("array_PC.GLASS_OD_NL_combined_excl_1P19Q."))
      ,
      metadata.glass_nl |> 
        dplyr::select("array_sentrix_id", "array_mnp_predictBrain_v12.8_cal_class", starts_with("array_PC.GLASS_OD_NL_combined_excl_1P19Q."))
    )
    ,
    by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')
  )

rm(metadata.glass_nl, metadata.glass_od)
gc()




# ODG + AC with all HQ-probes ----
### PCA itself ----

plt.split <- rbind(metadata.combi |> dplyr::mutate(facet = "dataset", col=dataset),
                   metadata.combi |> dplyr::mutate(facet = "grade", col=tumor_grade_h_l),
                   metadata.combi |> dplyr::mutate(facet = "batch", col=batch_us))

plt.split <- rbind(
                   metadata.combi |> dplyr::mutate(facet = "dataset", col=dataset),
                   #metadata.combi |> dplyr::mutate(facet = "grade", col=tumor_grade_h_l),
                   metadata.combi |> dplyr::mutate(facet = "mnp brain classifier", col=mnp_predictBrain_v12.8_cal_class)
                   #metadata.combi |> dplyr::mutate(facet = "batch", col=batch_us)
                   )


ggplot(plt.split, aes(x= PC.GLASS_OD_NL_combined.4, y= PC.GLASS_OD_NL_combined.3, col=col)) +
  facet_grid(cols = vars(facet)) +
  geom_point() +
  theme_bw()



plt <- metadata.combi |> 
  dplyr::mutate(col.grading = dplyr::case_when(
    dataset == "GLASS-OD" & mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","OLIGOSARC_IDH") ~ "OD [HG - A_IDH_HG & OLIGOSARC]",
    dataset == "GLASS-OD" & (mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","OLIGOSARC_IDH") == F) ~ "OD [LG ~ other]",
    
    dataset == "GLASS-NL" & mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") ~ "AC [HG - A_IDH_HG]",
    dataset == "GLASS-NL" & (mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == F) ~ "AC [LG ~ other]"
  )) |> 
  dplyr::mutate(col.grading = dplyr::case_when(
    dataset == "GLASS-OD" & mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","OLIGOSARC_IDH") ~ "HG",
    dataset == "GLASS-NL" & mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") ~ "HG",
    T ~ "LG"
  ))



p1 <- ggplot(plt, aes(x=PC.GLASS_OD_NL_combined.3, y=PC.GLASS_OD_NL_combined.4, col=batch_us, label=resection_id)) +
  geom_point(size=2.5) +
  geom_point(size=3.2, col="black",fill=NA, pch=21, alpha=0.75) +
  theme_bw()+ 
  theme(  legend.position = "bottom")

p2 <- ggplot(plt, aes(x=PC.GLASS_OD_NL_combined.3, y=PC.GLASS_OD_NL_combined.4, col=dataset, alpha=col.grading, label=resection_id)) +
  geom_point(size=2.5) +
  geom_point(data = subset(plt, col.grading == "HG"), size=3.2, col="black",fill=NA, pch=21) +
  scale_alpha_manual(values=c('HG'=1.0, 'LG'=0.45)) +
  #scale_color_manual(values=c('HG'='black', 'LG'='white')) +
  theme_bw() + 
  theme(  legend.position = "bottom")


p1 + p2



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


dat <- metadata.combi |>
  dplyr::select(`sentrix_id`, starts_with("PC.GLASS_OD_NL_combined.")) |> 
  tibble::column_to_rownames('sentrix_id') |> 
  t()

um <- M3C::umap(dat, seed=round(runif(1) * 10000))

umd <- um$data |> 
  as.data.frame() |> 
  dplyr::rename(UMAP1 = 1) |> 
  dplyr::rename(UMAP2 = 2) |> 
  tibble::rownames_to_column('sentrix_id')


plt <- metadata.combi |> 
  dplyr::left_join(umd, by=c('sentrix_id'='sentrix_id'),suffix=c('','')) |> 
  dplyr::mutate(puritygroups= cut(methylation_bins_1p19q_purity, breaks=4))


plt <- plt |> 
  dplyr::mutate(col.grading = dplyr::case_when(
    dataset == "GLASS-OD" & mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","OLIGOSARC_IDH") ~ "OD [HG - A_IDH_HG & OLIGOSARC]",
    dataset == "GLASS-OD" & (mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","OLIGOSARC_IDH") == F) ~ "OD [LG ~ other]",
    
    dataset == "GLASS-NL" & mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") ~ "AC [HG - A_IDH_HG]",
    dataset == "GLASS-NL" & (mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == F) ~ "AC [LG ~ other]"
  ))



ggplot(plt, aes(x=UMAP1, y=UMAP2, col=batch_us, label=resection_id)) +
  geom_point() +
  #geom_text(data = subset(plt, batch_us == "GLASS-OD [batch EU]") , size=3, col="black") +
  theme_bw()




p1 <- ggplot(plt, aes(x=UMAP1, y=UMAP2, col=batch_us, label=resection_id)) +
  geom_point() +
  #geom_text(data = subset(plt, batch_us == "GLASS-OD [batch EU]") , size=3, col="black") +
  theme_bw()

p2 <- ggplot(plt, aes(x=UMAP1, y=UMAP2, col=col.grading, label=resection_id)) +
  geom_point() +
  #geom_text(data = subset(plt, batch_us == "GLASS-OD [batch EU]") , size=3, col="black") +
  theme_bw()




### tSNE over PCA ----



dat <- metadata.combi |>
  dplyr::select(`sentrix_id`, starts_with("PC.GLASS_OD_NL_combined.")) |> 
  tibble::column_to_rownames('sentrix_id') |> 
  t()

ts <- M3C::tsne(dat,perplex=10)#1 and 2 are purity and quality
rownames(ts$data) <- colnames(dat)

tsd <- ts$data |> 
  as.data.frame() |> 
  dplyr::rename(tSNE1 = 1) |> 
  dplyr::rename(tSNE2 = 2) |> 
  tibble::rownames_to_column('sentrix_id')


plt <- plt |> 
  dplyr::left_join(tsd, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))

plt <- plt |> 
  dplyr::mutate(iso = tSNE1 < 0 & tSNE2 < 15 & dataset == "GLASS-OD")

ggplot(plt, 
       aes(x=tSNE1, y=tSNE2, col=batch_us, group=patient_id, label=resection_id)) +
  geom_point() +
  geom_text(data = subset(plt, batch_us == "GLASS-OD [batch EU]") , size=3, col="black") +
  theme_bw()


#sentrix_id resection_id
#1  201496850071_R02C01      0068-R2 - rare isoatie?

#2  203989100149_R04C01      0017-R2
#3  203989100149_R08C01      0018-R3

#4  204073520033_R02C01      0021-R1
#5  204073520033_R03C01      0021-R2

#6  204073520033_R04C01      0022-R1
#7  204073520033_R05C01      0022-R2

#8  204088040075_R07C01      0023-R2
#9  206116800035_R07C01      0029-R2
#10 206137490041_R04C01      0038-R2
#11 206137490041_R06C01      0037-R2


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


## ROC linear classifier ----
#' 10xCV predict OD of AC  ----
#' https://www.projectpro.io/recipes/plot-auc-roc-curve-r
#' https://machinelearningmastery.com/linear-classification-in-r/

# install.packages("caTools")    # For Logistic regression 
# install.packages('pROC')       # For ROC curve to evaluate model 
# install.packages('VGAM')
install.packages('cvAUC')
# https://stackoverflow.com/questions/41533811/roc-curve-in-linear-discriminant-analysis-with-r
library(caTools)
library(pROC)
library(VGAM)
library(MASS)
library(ROCR)
library(cvAUC)

data <- metadata.combi |> 
  dplyr::select(sentrix_id, dataset, paste0("PC.GLASS_OD_NL_combined.",1:50))  |>  # avoid curse of dimension
  dplyr::mutate(dataset = as.factor(dataset)) |> 
  dplyr::mutate(i = 1:dplyr::n()) |> 
  dplyr::mutate(slice = (i %% 10) + 1) |> 
  dplyr::mutate(i = NULL) |> 
  tibble::column_to_rownames('sentrix_id')




train <- subset(data, split == "TRUE") 
test <- subset(data, split == "FALSE") 

# fit.glm.gaus = glm(dataset ~ .,train , family="gaussian")     # we use the glm()-general linear model to create an instance of model
# fit.vglm <- vglm(dataset ~ ., family=multinomial, data=train)

fit.lda <- lda(dataset ~ ., data=train)
summary(fit.lda)
predictions.lda <- predict(fit.lda, test |> dplyr::select(paste0("PC.GLASS_OD_NL_combined.",1:50)))$class
table(predictions.lda, test$dataset)


pred <- predict(fit.lda, test |> dplyr::select(paste0("PC.GLASS_OD_NL_combined.",1:50)), type="response")
test_roc = roc(test$dataset ~ pred, plot = TRUE, print.auc = TRUE)


predd <- prediction(pred$posterior[,2], test$dataset) 
perf <- performance(predd,"tpr","fpr")
plot(perf,colorize=TRUE)
roc(test$dataset, pred$posterior[,2],plot=TRUE, legacy.axes = TRUE, 
    percent =TRUE, xlab="False Positive Percentage", ylab="True Positive Percentage")



# https://search.r-project.org/CRAN/refmans/cvAUC/html/cvAUC.html

fits <- c()
predicts <- c()
predictions <- c()
performances <- c()
tests <- c()
for(sslice in 1:10) {
  train <- data |> 
    dplyr::filter(slice != sslice) |> 
    dplyr::mutate(slice = NULL)
  
  test <- data |> 
    dplyr::filter(slice == sslice) |> 
    dplyr::mutate(slice = NULL)
  
  f.fit <- lda(dataset ~ ., data=train)
  p.predict <- predict(f.fit, test |> dplyr::select(paste0("PC.GLASS_OD_NL_combined.",1:50)))
  p.prediction <- prediction(p.predict$posterior[,2], test$dataset) 
  p.performance <- performance(p.prediction,"tpr","fpr")
  
  
  fits <- c(fits, f.fit)
  predicts <- c(predicts, list(p.predict))
  predictions <- c(predictions, p.prediction)
  performances <- c(performances, p.performance)
  tests <- c(tests, list(test))
}


q <- data.frame(
  predicted = c(
    predicts[[1]]$class,
    predicts[[2]]$class,
    predicts[[3]]$class,
    predicts[[4]]$class,
    predicts[[5]]$class,
    predicts[[6]]$class,
    predicts[[7]]$class,
    predicts[[8]]$class,
    predicts[[9]]$class,
    predicts[[10]]$class
  ),
  label = c(
    tests[[1]]$dataset,
    tests[[2]]$dataset,
    tests[[3]]$dataset,
    tests[[4]]$dataset,
    tests[[5]]$dataset,
    tests[[6]]$dataset,
    tests[[7]]$dataset,
    tests[[8]]$dataset,
    tests[[9]]$dataset,
    tests[[10]]$dataset
    
  ),
  name = c(
    rownames(tests[[1]]),
    rownames(tests[[2]]),
    rownames(tests[[3]]),
    rownames(tests[[4]]),
    rownames(tests[[5]]),
    rownames(tests[[6]]),
    rownames(tests[[7]]),
    rownames(tests[[8]]),
    rownames(tests[[9]]),
    rownames(tests[[10]])
  )
)
# 203989100149_R08C01 = 0018-R3
# q |> dplyr::filter(predicted!=label)

# p <- list(
#   predictions = list(
#     '1' = predicts[[1]]$posterior[,2],
#     '2' = predicts[[2]]$posterior[,2],
#     '3' = predicts[[3]]$posterior[,2],
#     '4' = predicts[[4]]$posterior[,2],
#     '5' = predicts[[5]]$posterior[,2],
#     '6' = predicts[[6]]$posterior[,2],
#     '7' = predicts[[7]]$posterior[,2],
#     '8' = predicts[[8]]$posterior[,2],
#     '9' = predicts[[9]]$posterior[,2],
#     '10' = predicts[[10]]$posterior[,2]
#   ),
#   labels = list(
#     '1' = tests[[1]]$dataset,
#     '2' = tests[[2]]$dataset,
#     '3' = tests[[3]]$dataset,
#     '4' = tests[[4]]$dataset,
#     '5' = tests[[5]]$dataset,
#     '6' = tests[[6]]$dataset,
#     '7' = tests[[7]]$dataset,
#     '8' = tests[[8]]$dataset,
#     '9' = tests[[9]]$dataset,
#     '10' = tests[[10]]$dataset
#   )
# )

p <- list(
  predictions = list(
    '1' = c(
          predicts[[1]]$posterior[,2],
          predicts[[2]]$posterior[,2],
          predicts[[3]]$posterior[,2],
          predicts[[4]]$posterior[,2],
          predicts[[5]]$posterior[,2],
          predicts[[6]]$posterior[,2],
          predicts[[7]]$posterior[,2],
          predicts[[8]]$posterior[,2],
          predicts[[9]]$posterior[,2],
          predicts[[10]]$posterior[,2]
    )
  ),
  labels = list(
    '1' = c(
          tests[[1]]$dataset,
          tests[[2]]$dataset,
          tests[[3]]$dataset,
          tests[[4]]$dataset,
          tests[[5]]$dataset,
          tests[[6]]$dataset,
          tests[[7]]$dataset,
          tests[[8]]$dataset,
          tests[[9]]$dataset,
          tests[[10]]$dataset
      )
  )
)



rocobj <- roc(p$labels$`1`, p$predictions$`1`)
auc <- round(auc(p$labels$`1`, p$predictions$`1`),5)

ggroc(rocobj, colour = 'red', size = 1) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ') - astrocytoma vs. oligodendroglioma')) +
  theme_bw()+ 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")





# ODG + AC with all HQ-probes not 1P / 19Q ----
## init plt ----


plt <- metadata.combi |> 
  dplyr::mutate(col.grading1 = dplyr::case_when(
    dataset == "GLASS-OD" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","OLIGOSARC_IDH") ~ "HG",
    dataset == "GLASS-NL" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") ~ "HG",
    T ~ "LG"
  )) |> 
  dplyr::mutate(col.grading2 = dplyr::case_when(
    dataset == "GLASS-OD" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","OLIGOSARC_IDH") ~ "OD [HG - A_IDH_HG & OLIGOSARC]",
    dataset == "GLASS-OD" & (array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","OLIGOSARC_IDH") == F) ~ "OD [LG ~ other]",
    
    dataset == "GLASS-NL" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") ~ "AC [HG - A_IDH_HG]",
    dataset == "GLASS-NL" & (array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG") == F) ~ "AC [LG ~ other]"
  ))



## PCA ----


plt <- glass_od.metadata.array_samples


for(ppc in paste0("PC",1:10)) {
  w = wilcox.test(
    plt |> dplyr::filter(isolation_person_name == "USA / Duke") |> dplyr::pull(paste0(ppc)),
    plt |> dplyr::filter(isolation_person_name != "USA / Duke") |> dplyr::pull(paste0(ppc))
  )
  
  print(paste0(ppc, " -- ",w$p.value))
}




p1 <- ggplot(plt, aes(x=PC1, y=PC2, col=isolation_person_name == "USA / Duke")) +
  geom_point() +
  theme_classic() +
  labs(col = "USA sample") +
  theme(legend.position = "bottom")

p2 <- ggplot(plt, aes(x=PC1, y=PC4, col=isolation_person_name == "USA / Duke")) +
  geom_point() +
  theme_classic() +
  labs(col = "USA sample") +
  theme(legend.position = "bottom")

p3 <- ggplot(plt, aes(x=PC2, y=PC4, col=isolation_person_name == "USA / Duke")) +
  geom_point() +
  theme_classic() +
  labs(col = "USA sample") +
  theme(legend.position = "bottom")

p1 + p2 + p3


### fig for export ----



plt <- metadata.combi |> 
  dplyr::mutate(col = 
                  dplyr::case_when(
                    dataset == "GLASS-OD" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","A_IDH","A_IDH_HG") ~ "v12.8 OD classified as AC",
                    dataset == "GLASS-OD" & array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","A_IDH","A_IDH_HG") == F ~ "v12.8 OD",
                    dataset == "GLASS-NL" & array_mnp_predictBrain_v12.8_cal_class %in% c("O_IDH","OLIGOSARC_IDH") ~ "v12.8 AC classified as OD",
                    dataset == "GLASS-NL" & array_mnp_predictBrain_v12.8_cal_class %in% c("O_IDH","OLIGOSARC_IDH") == F ~ "v12.8 AC",
                  )
  )



ggplot(plt, aes(x=array_PC.GLASS_OD_NL_combined_excl_1P19Q.3, y=array_PC.GLASS_OD_NL_combined_excl_1P19Q.4, col=col,
                label=resection_id)) +
  geom_point(data = subset(plt, grepl("misclass", col) == F), size=theme_nature_size/3) +
  geom_point(data = subset(plt, grepl("misclass", col) == T), size=theme_nature_size/3) +
  theme_nature + 
  scale_color_manual(values=c('v12.8 AC'= mixcol('lightblue','darkblue'),
                              'v12.8 OD'= mixcol('lightgreen','darkgreen'),
                              'v12.8 OD classified as AC'='red',
                              'v12.8 AC classified as OD'='orange')) +
  labs(col = "",
       x = "PC3", 
       y = "PC4",
       subtitle = format_subtitle("PCA on GLASS-NL + GLASS-OD combined, excl 1P / 19Q probes"))




ggsave("output/figures/vis_mixed_unsupervised_combined_no1p19_PCA.pdf", width = (8.5/3) * 0.95, height = 2.80)







# find dataset best segregating components (PC3 and PC4 usually)
for(ppc in paste0("PC.GLASS_OD_NL_combined_excl_1P19Q.",1:50)) {
  w = wilcox.test(
    plt |> dplyr::filter(dataset == "GLASS-NL") |> dplyr::pull(paste0(ppc)),
    plt |> dplyr::filter(dataset == "GLASS-OD") |> dplyr::pull(paste0(ppc))
  )
  
  print(paste0(ppc, " -- ",w$p.value))
}




## UMAP over 1st 50 PCA ----


dat <- metadata.combi |>
  dplyr::select(`sentrix_id`, starts_with("PC.GLASS_OD_NL_combined_excl_1P19Q.")) |> 
  tibble::column_to_rownames('sentrix_id') |> 
  t() |> 
  M3C::umap(seed=round(runif(1) * 10000)) |> 
  purrr::pluck('data') |> 
  as.data.frame() |> 
  dplyr::rename(UMAP1 = 1) |> 
  dplyr::rename(UMAP2 = 2) |> 
  tibble::rownames_to_column('sentrix_id')


plt <- plt |> 
  dplyr::left_join(dat, by=c('sentrix_id'='sentrix_id'),suffix=c('',''))




p1 <- ggplot(plt, aes(x=UMAP1, y=UMAP2, col=batch_us, label=resection_id)) +
  geom_point(size=2.5) +
  geom_point(size=3.2, col="black",fill=NA, pch=21, alpha=0.75) +
  theme_bw()+ 
  theme(  legend.position = "bottom")

p2 <- ggplot(plt, aes(x=UMAP1, y=UMAP2, col=dataset, alpha=col.grading1, label=resection_id)) +
  geom_point(size=2.5) +
  geom_point(data = subset(plt, col.grading1 == "HG"), size=3.2, col="black",fill=NA, pch=21) +
  scale_alpha_manual(values=c('HG'=1.0, 'LG'=0.45)) +
  #scale_color_manual(values=c('HG'='black', 'LG'='white')) +
  theme_bw() + 
  theme(  legend.position = "bottom")


p1 + p2



## LDA ----
#' 10xCV predict OD of AC  ----
#' https://www.projectpro.io/recipes/plot-auc-roc-curve-r
#' https://machinelearningmastery.com/linear-classification-in-r/

# install.packages("caTools")    # For Logistic regression 
# install.packages('pROC')       # For ROC curve to evaluate model 
# install.packages('VGAM')
# install.packages('cvAUC')
# https://stackoverflow.com/questions/41533811/roc-curve-in-linear-discriminant-analysis-with-r
library(caTools)
library(pROC)
library(VGAM)
library(MASS)
library(ROCR)
library(cvAUC)


lda.data <- metadata.combi |> 
  dplyr::select(array_sentrix_id, dataset, paste0("array_PC.GLASS_OD_NL_combined_excl_1P19Q.", 1:50))  |>  # avoid curse of dimension
  dplyr::mutate(dataset = as.factor(dataset)) |> 
  dplyr::mutate(i = 1:dplyr::n()) |> 
  dplyr::mutate(slice = (i %% 10) + 1) |> 
  dplyr::mutate(i = NULL) |> 
  tibble::column_to_rownames('array_sentrix_id')


#train <- subset(data, split == "TRUE") 
#test <- subset(data, split == "FALSE") 

# fit.glm.gaus = glm(dataset ~ .,train , family="gaussian")     # we use the glm()-general linear model to create an instance of model
# fit.vglm <- vglm(dataset ~ ., family=multinomial, data=train)

# fit.lda <- lda(dataset ~ ., data=train)
# summary(fit.lda)
# predictions.lda <- predict(fit.lda, test |> dplyr::select(paste0("PC.GLASS_OD_NL_combined.",1:50)))$class
# table(predictions.lda, test$dataset)
# 
# 
# pred <- predict(fit.lda, test |> dplyr::select(paste0("PC.GLASS_OD_NL_combined.",1:50)), type="response")
# test_roc = roc(test$dataset ~ pred, plot = TRUE, print.auc = TRUE)
# 
# 
# predd <- prediction(pred$posterior[,2], test$dataset) 
# perf <- performance(predd,"tpr","fpr")
# plot(perf,colorize=TRUE)
# roc(test$dataset, pred$posterior[,2],plot=TRUE, legacy.axes = TRUE, 
#     percent =TRUE, xlab="False Positive Percentage", ylab="True Positive Percentage")



# https://search.r-project.org/CRAN/refmans/cvAUC/html/cvAUC.html

lda.labels <- c()
lda.classes <- c()
lda.posterior <- data.frame()
for(sslice in 1:10) {
  lda.train <- lda.data |> 
    dplyr::filter(slice != sslice) |> 
    dplyr::mutate(slice = NULL)
  
  lda.test <- lda.data |> 
    dplyr::filter(slice == sslice) |> 
    dplyr::mutate(slice = NULL)
  
  lda.predict <- predict(MASS::lda(dataset ~ ., data=lda.train), lda.test |>  dplyr::mutate(slice = NULL) )
  ds <- lda.test$dataset
  names(ds) <- rownames(lda.test)
  
  lda.labels <- c(lda.labels, ds)
  lda.classes <- c(lda.classes, lda.predict$class)
  lda.posterior <- rbind(lda.posterior, lda.predict$posterior)

}

table(lda.classes , lda.labels)
metadata.combi |> dplyr::filter(array_sentrix_id %in% names(which(lda.classes != lda.labels))) |> 
  dplyr::select(!contains("_PC."))



plt <- plt |> 
  dplyr::left_join(
    lda.posterior |>
      dplyr::rename_with( ~ paste0("lda.posterior.", .x)) |> 
      tibble::rownames_to_column('array_sentrix_id')
    , by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))





### ROC + AUC ----


#rocobj <- pROC::roc(lda.labels, lda.posterior$`GLASS-OD`)
#auc <- round(pROC::auc(lda.labels, lda.posterior$`GLASS-NL`),5)


rocobj <- pROC::roc(lda.labels, lda.classes)
auc <- round(pROC::auc(lda.labels, lda.classes),4)

pROC::ggroc(rocobj, colour = 'red', size = 1) +
  theme_bw()+ 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +
  labs(subtitle = paste0('ROC Curve - LDA AC vs. OD ', '(AUC = ', auc, ')')) +
  theme_bw() +
  theme_cellpress


ggsave("output/figures/vis_mixed_unsupervised_combined_LDA_ROC.pdf", width = (8.5/3) * 0.95, height = 2.65)



### diff confidence/probabilities HG vs LG ----



plt <- metadata.combi |> 
  dplyr::left_join(
    lda.posterior |>
      dplyr::rename_with( ~ paste0("lda.posterior.", .x)) |> 
      tibble::rownames_to_column('array_sentrix_id')
    , by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('',''))





plt.p <- rbind(
  plt |>
    dplyr::filter(dataset == "GLASS-NL") |> 
    dplyr::mutate(lda.posterior = `lda.posterior.GLASS-OD`)
  ,
  plt |>
    dplyr::filter(dataset == "GLASS-OD") |> 
    dplyr::mutate(lda.posterior = `lda.posterior.GLASS-NL`)
) |> 
  dplyr::mutate(classification_status = ifelse(lda.posterior < 0.5, "match", "misclassification")) |> 
  dplyr::filter(!is.na(tumor_grade_h_l)) |> 
  dplyr::mutate(tumor_grade_h_l = factor(tumor_grade_h_l, levels=c("low grade", "high grade"))) |> 
  dplyr::filter(!is.na(tumor_grade_h_l))

table(plt.p$dataset)






ggplot(plt.p, aes(x=tumor_grade_h_l, y=-log(lda.posterior), col=classification_status)) +
  geom_hline(yintercept = 0, lwd=theme_cellpress_lwd) +
  #facet_grid(cols = vars(dataset)) +
  facet_wrap(~dataset, scales="free") +
  ggbeeswarm::geom_quasirandom(size=theme_nature_size/3) +
  ggpubr::stat_compare_means(aes(group=tumor_grade_h_l), label.x.npc=0.1, method = "wilcox.test", show_guide  = FALSE,  size=theme_cellpress_size, family=theme_nature_font_family) + 
  labs(y = "-log( LDA posterior probability )") +
  geom_hline(yintercept = -log(0.5), col="red", lwd=0.5,lty=2) +
  scale_color_manual(values=c('match'='black', 'misclassification' = 'red')) + 
  theme_nature +
  labs(x=NULL, col=NULL) +
  ylim(0, 72) 


ggsave("output/figures/vis_mixed_unsupervised__combined_posterior_wilcox.pdf", width=(8.5*0.95)/3, height=2.75)




plt.p <- plt.p |> 
  dplyr::mutate(col_mnp = array_mnp_predictBrain_v12.8_cal_class) |> 
  dplyr::mutate(col_mnp = gsub("^CTRL_.+$","CTRL_*", col_mnp)) |> 
  dplyr::mutate(col_mnp = gsub("^GBM_.+$","GBM_*", col_mnp))



ggplot(plt.p, aes(x=tumor_grade_h_l, y=-log(lda.posterior), col=col_mnp)) +
  geom_hline(yintercept = 0, lwd=theme_cellpress_lwd) +
  #facet_grid(cols = vars(dataset)) +
  facet_wrap(~dataset, scales="free") +
  ggbeeswarm::geom_quasirandom(size=theme_nature_size/3) +
  ggpubr::stat_compare_means(aes(group=tumor_grade_h_l), label.x.npc=0.1, method = "wilcox.test", show_guide  = FALSE,  size=theme_cellpress_size, family=theme_nature_font_family) + 
  labs(y = "-log( LDA posterior probability )") +
  geom_hline(yintercept = -log(0.5), col="red", lwd=0.5,lty=2) +
  scale_color_manual(values=c( "GBM_*"="blue", "CTRL_*"="#444444", palette_mnp_12.8_6)) + 
  theme_nature +
  labs(x=NULL, col=NULL) +
  ylim(0, 72) 



ggsave("output/figures/vis_mixed_unsupervised__combined_posterior_wilcox_MNP_col.pdf", width=(8.5*0.95)/3, height=2.75)




