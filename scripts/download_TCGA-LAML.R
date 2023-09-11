#!/usr/bin/env R


d <- read.table("data/TCGA-AML/TCGA_AML_450k_methylation_data.txt",header=T,stringsAsFactors = F)
m <- read.table("data/TCGA-AML/mutations.txt",header=T)


# download ----

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")


s <- m |> 
  dplyr::filter(!is.na(IDH1) | !is.na(IDH2))


# TARGET-AML
# TCGA-AB

library(TCGAbiolinks)
# https://rdrr.io/bioc/TCGAbiolinks/man/GDCquery.html
query1 <- GDCquery(
  project = "TARGET-AML",
  data.category = "DNA Methylation",
  data.format = "idat"
)




query <- GDCquery(
  project = "TARGET-AML",
  data.category = "Raw microarray data",
  data.type = "Raw intensities",
  experimental.strategy = "Methylation array",
  legacy = TRUE,
  file.type = ".idat",
  platform = "Illumina Human Methylation 450")

tmp <- getGDCprojects()
# range(d[6,10:135])


# https://www.genelibs.cn/ftp/raw_data_tcga/LAML/





query <- GDCquery(
  project = "TCGA-LAML",
  data.category = "DNA Methylation",
  file.type = "idat",
  legacy = F
)


GDCdownload(query = query)


# tmp <- getGDCprojects()
# 
# match.file.cases <- getResults(query,cols = c("cases","file_name"))
# match.file.cases$project <- "TCGA-LGG"


query <- TCGAbiolinks::GDCquery(
  project = "TCGA-LAML",
  legacy = F,
  data.category = "Clinical"
)




TCGAbiolinks::GDCdownload(query = query)











