#!/usr/bin/env R


# load data ----


library(ggplot2)


if(!exists('youri_gg_theme')) {
  source('scripts/youri_gg_theme.R')
}



if(!exists('glass_od.metadata.idats')) {
  source('scripts/load_metadata.R')
}



if(!exists('glass_od.data.mvalues')) {
  source('scripts/load_mvalues.R')
}




# load all samples, qc and corr w/ qc stats ----


metadata <- glass_od.metadata.idats |> 
  dplyr::filter(is.na(reason_excluded_patient)) |> 
  dplyr::filter(is.na(reason_excluded_resection)) |> 
  dplyr::filter(is.na(reason_excluded_resection_isolation)) |> 
  dplyr::filter(is.na(reason_excluded_array_sample)) |>  # those are typically the replicates
  assertr::verify(!duplicated(resection_id)) |> 
  assertr::verify(sentrix_id != "204808700074_R04C01") |>  # 0003-R3-repA
  dplyr::filter(study_name == "GLASS-OD")



# load m-values ----


data <- glass_od.data.mvalues |> 
    dplyr::select(all_of(metadata$sentrix_id))


#data.full = data
#data = data.full |> 
#  dplyr::slice_head(n=2000)

for(sdf in c(0.9,1.0,1.1,1.2,1.3,1.4,1.5)) {
  
    
  pca.raw <- data |>
    (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() # this synthax, oh my
  
  dim(pca.raw)
  
  pca.raw <- pca.raw |> 
    #dplyr::filter(mad > -9999) |> 
    #dim()
    dplyr::filter(mad > sdf) |>
    dplyr::mutate(mad = NULL)  |>  # remove the sd to obtain original vst matrix
    t() |>
    prcomp() |> 
    purrr::pluck('x')  |>   # take coordinates
    as.data.frame(stringsAsFactor=F) |>   # transform back from matrix to data.frame
    tibble::rownames_to_column('sentrix_id')
  
  
  ## find most variable features ---
  
  plt.pca <- metadata |>
    dplyr::left_join(pca.raw |>  dplyr::select(sentrix_id, paste0("PC",1:25)),
                     by=c('sentrix_id'='sentrix_id'),
                     suffix=c('','')
                     )

  #plt <- rbind(
    # plt |>
    #   dplyr::mutate(panel = "primary vs. non primary") |>
    #   dplyr::mutate(col = ifelse(resection_number == 1, "primary", "recurrence"))
    # ,
    #plt |>
     # dplyr::mutate(panel = "grading") |>
    #  dplyr::mutate(col = paste0("g",resection_tumor_grade)) |>
    #  dplyr::mutate(col = ifelse(col == "g4","g3", col))
  #)

  # 
  # ggplot(plt, aes(x=PC4, y=PC2, col=col)) + 
  #   facet_grid(cols = vars(panel), scales = "free", space="free")+ 
  #   geom_point() +
  #   theme_bw()
  # 
  
  
  for(pc in paste0("PC",1:6)) {
    wt <- wilcox.test(
      plt.pca |> 
        dplyr::filter(!is.na(resection_tumor_grade)) |> 
        dplyr::filter(resection_tumor_grade == 2) |> 
        dplyr::pull(pc),
      plt.pca |> 
        dplyr::filter(!is.na(resection_tumor_grade)) |> 
        dplyr::filter(resection_tumor_grade == 3) |> 
        dplyr::pull(pc)
    )
    
    print(paste0(pc, ": ", wt$p.value, "  sdf: ", sdf))
  }
  
  
  rm(plt, wt)
}




ggplot(plt, aes(x=PC4, y=PC3, col=col)) + 
  facet_grid(cols = vars(panel), scales = "free", space="free")+ 
  geom_point() +
  theme_bw()


dat <- plt |>
  dplyr::filter(!is.na(resection_tumor_grade)) |> 
  dplyr::filter(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(resection_tumor_grade = factor(paste0('G',resection_tumor_grade), levels=c('G2','G3')))


logistic_model <- glm(resection_tumor_grade ~ PC4, data=dat, family=binomial)

#Data frame with hp in ascending order
Predicted_data <- data.frame(PC4=seq(  min(dat$PC4), max(dat$PC4),len=500))

# Fill predicted values using regression model
Predicted_data$resection_tumor_grade = predict(
  logistic_model, Predicted_data, type="response")

# Plot Predicted data and original data points
plot((as.numeric(resection_tumor_grade)-1) ~ PC4, data=dat)
lines(as.numeric(resection_tumor_grade) ~ PC4, Predicted_data, lwd=2, col="green")




for(sdf in c(0.9,1.0,1.1,1.2,1.3,1.4,1.5)) {
  
  
  pca.raw <- data |>
    (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |>  # this synthax, oh my
    dplyr::filter(mad > sdf) |>
    dplyr::mutate(mad = NULL)  |>  # remove the sd to obtain original vst matrix
    t() |>
    prcomp() |> 
    purrr::pluck('x')  |>   # take coordinates
    as.data.frame(stringsAsFactor=F) |>   # transform back from matrix to data.frame
    tibble::rownames_to_column('sentrix_id')
  
  plt <- metadata |>
    dplyr::left_join(pca.raw |>  dplyr::select(sentrix_id, paste0("PC",1:25)),
                     by=c('sentrix_id'='sentrix_id'),
                     suffix=c('','')
    )
  
  
  
  
  for(pc in paste0("PC",1:8)) {
    wt <- wilcox.test(
      plt |> 
        dplyr::filter(!is.na(rs_gender_predicted)) |> 
        dplyr::filter(rs_gender_predicted == "F") |> 
        dplyr::pull(pc),
      plt |> 
        dplyr::filter(!is.na(rs_gender_predicted)) |> 
        dplyr::filter(rs_gender_predicted == "M") |> 
        dplyr::pull(pc)
    )
    
    
    print(paste0(pc, ": ", wt$p.value, "  sdf: ", sdf))
  }

}


# PCA from CATNON ----


pc <- readRDS('cache/prcomp_unsup_meth_catnon.Rds')
pc.features <- pc |>
  purrr::pluck('rotation') |>
  as.data.frame(stringsAsFactors=F) |>
  tibble::rownames_to_column('cg') |> 
  dplyr::pull('cg')
pc.dat <- data |> 
  tibble::rownames_to_column('cg') |>
  dplyr::filter(cg %in% pc.features) |> 
  tibble::column_to_rownames('cg')


dim(pc.dat)
dim(pc$rotation)

predict(na.pass(pc), na.pass(t(pc.dat))) |> 
  as.data.frame(stringsAsFactor=F) |>   # transform back from matrix to data.frame
  tibble::rownames_to_column('sentrix_id') |> 
  head(n=15)



plt <- metadata |>
  dplyr::left_join(pca.raw |>  dplyr::select(sentrix_id, paste0("PC",1:25)),
                   by=c('sentrix_id'='sentrix_id'),
                   suffix=c('','')
  ) |> 
  dplyr::mutate(col = paste0("g",resection_tumor_grade)) |>
  dplyr::mutate(col = ifelse(col == "g4","g3", col))



ggplot(plt, aes(x=PC1, y=PC2, col=col)) + 
  #facet_grid(cols = vars(panel), scales = "free", space="free")+ 
  geom_point() +
  theme_bw()





