#!/usr/bin/env R

# load ----


source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')


library(ggplot2)
library(patchwork)


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# GLASS-OD / OD ----


metadata <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(210)



## x-check batch effect isolation(s) ----


plt <- metadata |> 
  dplyr::mutate(col = ifelse(isolation_person_name == "USA / Duke", "Batch [US]", "Batch [EU]")) |> 
  dplyr::mutate(col2 = dplyr::case_when(
    isolation_id == "0100-R2-repA" ~ "US: DNA as is",
    isolation_id == "0100-R2-repB" ~ "US: DNA over column",
    
    isolation_id == "0101-R2-repA" ~ "US: DNA as is",
    isolation_id == "0101-R2-repB" ~ "US: DNA over column",
    
    isolation_id == "0102-R2-repA" ~ "US: DNA as is",
    isolation_id == "0102-R2-repB" ~ "US: DNA over column",
    
    isolation_id == "0104-R1-repA" ~ "US: DNA as is",
    isolation_id == "0104-R1-repB" ~ "US: DNA over column",
    
    resection_id == "101-R3" ~ "US: FFPE",
    resection_id == "102-R3" ~ "US: FFPE",
    
    resection_id == "0016-R2" ~ "EU: suspected frozen",
    
    T ~ col
  )) |> 
  dplyr::mutate(group = dplyr::case_when(
    isolation_id == "0100-R2-repA" ~ "0100-R2",
    isolation_id == "0100-R2-repB" ~ "0100-R2",
    
    isolation_id == "0101-R2-repA" ~ "0101-R2",
    isolation_id == "0101-R2-repB" ~ "0101-R2",
    
    isolation_id == "0102-R2-repA" ~ "0102-R2",
    isolation_id == "0102-R2-repB" ~ "0102-R2",
    
    isolation_id == "0104-R1-repA" ~ "0104-R1",
    isolation_id == "0104-R1-repB" ~ "0104-R1",
    
    T ~ resection_id
  ))


ggplot(plt, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=-array_PC2, col=col2, group=group)) +
  geom_line(lwd = 0.3, col="black") +
  geom_point(data = subset(plt, col2 %in% c("Batch [EU]", "Batch [US]")), alpha=0.26) +
  geom_point(data = subset(plt, col2 %in% c("Batch [EU]", "Batch [US]") == F)) +
  theme_cellpress





## Figure 2B: PCA ----

### new ----


plt.split <- rbind(
  metadata |> 
    dplyr::mutate(col = dplyr::case_when(
      is.na(isolation_material) ~ "NA",
      isolation_material == "ffpe" ~ "FFPE",
      isolation_material == "tissue" ~ "Tissue"
    )) |> 
    dplyr::mutate(facet = "A) Material")  |> 
    dplyr::mutate(stat = col)
  ,
  metadata |> 
    dplyr::mutate(col = ifelse(resection_number == 1, "primary", "recurrent")) |> 
    dplyr::mutate(facet = "B) Resection")  |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata |>
    dplyr::mutate(col = as.factor(paste0("Grade ",resection_tumor_grade))) |> 
    dplyr::mutate(facet = "C) Histological grade") |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata |> 
    dplyr::mutate(col = ifelse(array_mnp_predictBrain_v2.0.1_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", array_mnp_predictBrain_v2.0.1_cal_class)) |> 
    dplyr::mutate(facet = "D) MNP CNS Classifier 114b/2.0.1")  |> 
    dplyr::mutate(stat = 'NA')
   ,
  metadata |> 
    dplyr::mutate(col = ifelse(array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", array_mnp_predictBrain_v12.8_cal_class)) |> 
    dplyr::mutate(facet = "E) MNP CNS Classifier 12.8")  |> 
    dplyr::mutate(stat = 'NA')
)



# array_PC2 ~ 202 | array_GLASS_NL_g2_g3_sig ~ 218
ggplot(plt.split, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=array_GLASS_NL_g2_g3_sig, col=col)) + 
  facet_wrap(~facet, scales="free",ncol=5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label), col=stat),  size=theme_cellpress_size, cor.coef.name ="rho", show_guide = FALSE, label.y.npc = "bottom") +
  geom_point(size=theme_nature_size/3) +
  scale_color_manual(values=c(
    "Grade 2"= col3(11)[10], "O_IDH"=col3(11)[10], "primary"=col3(11)[10], 
    "Grade 3"="red", "A_IDH_HG"="red", "recurrent"="red",
    
    "FFPE"="darkgreen",
    "OLIGOSARC_IDH" = "orange",
    "NA" = "gray40",
    "other" = "purple",
    "Tissue"="brown"
  )) +
  labs(x="Astrocytoma CGC Lasso",y = "GLASS-NL signature (supervised primary - recurrent)", col="") +
  labs(subtitle=format_subtitle("Logit WHO grade")) +
  theme_nature


ggsave("output/figures/vis_LGC_x_PCA__scatter_A.pdf",width=8.5 * 0.975, height = 2.72)



### old ----

plt.split <- rbind(
  metadata |>
    dplyr::mutate(col = as.factor(paste0("Grade ",resection_tumor_grade))) |> 
    dplyr::mutate(facet = "Histological grade") |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata |> 
    dplyr::mutate(col = ifelse(array_mnp_predictBrain_v2.0.1_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", array_mnp_predictBrain_v2.0.1_cal_class)) |> 
    dplyr::mutate(facet = "MNP CNS Classifier 114b/2.0.1")  |> 
    dplyr::mutate(stat = 'NA')
  # ,
  # metadata |> 
  #   dplyr::mutate(col = ifelse(array_mnp_predictBrain_v12.5_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", array_mnp_predictBrain_v12.8_cal_class)) |> 
  #   dplyr::mutate(facet = "MNP CNS Classifier 12.5")  |> 
  #   dplyr::mutate(stat = 'NA')
  ,
  metadata |> 
    dplyr::mutate(col = ifelse(array_mnp_predictBrain_v12.8_cal_class %in% c("A_IDH_HG","O_IDH","OLIGOSARC_IDH") == F, "other", array_mnp_predictBrain_v12.8_cal_class)) |> 
    dplyr::mutate(facet = "MNP CNS Classifier 12.8")  |> 
    dplyr::mutate(stat = 'NA')
  ,
  metadata |> 
    dplyr::mutate(col = ifelse(isolation_person_name == "USA / Duke", "Batch [US]", "Batch [EU]")) |> 
    dplyr::mutate(facet = "Batch")  |> 
    dplyr::mutate(stat = col)
)




ggplot(plt.split, aes(x=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, y=array_PC2, col=col)) + 
  #facet_grid(cols = vars(facet), scales = "free", space="free") +
  facet_wrap(~facet, scales="free",ncol=5) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label), col=stat),  size=theme_cellpress_size, cor.coef.name ="rho", show_guide = FALSE) +
  geom_point(size=theme_cellpress_size/3) +
  scale_color_manual(values=c(
    "Grade 2"= col3(11)[10],"O_IDH"=col3(11)[10], 
    "Grade 3"="red","A_IDH_HG"="red",  
    
    "Batch [EU]"="darkgreen",
    "OLIGOSARC_IDH" = "orange",
    "NA" = "gray40",
    "other" = "purple",
    "Batch [US]"="brown"
  )) +
  labs(x="AcCGAP",y = "PC2", col="") +
  labs(subtitle=format_subtitle("Logit WHO grade")) +
  theme_cellpress





# pc loadings

pc <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_prcomp.Rds")


plot(sort(pc$rotation[,1]), type="l",col="gray")
lines(sort(pc$rotation[,2]), type="l", col="blue")
lines(sort(pc$rotation[,3]), type="l", col="gray")
lines(sort(pc$rotation[,4]), type="l", col="gray")
lines(sort(pc$rotation[,5]), type="l", col="gray")
lines(sort(pc$rotation[,6]), type="l", col="gray")
lines(sort(pc$rotation[,7]), type="l", col="gray")
lines(sort(pc$rotation[,8]), type="l", col="gray")
lines(sort(pc$rotation[,9]), type="l", col="gray")
lines(sort(pc$rotation[,10]), type="l", col="gray")
abline(h=0, col="red")


n = length(pc$rotation[,1])

s1 = sort(pc$rotation[,1])
fdiv1 = s1[2:n] - s1[1:(n-1)]

s2 = sort(pc$rotation[,2])
fdiv2 = s2[2:n] - s2[1:(n-1)]

s3 = sort(pc$rotation[,3])
fdiv3 = s3[2:n] - s3[1:(n-1)]

s4 = sort(pc$rotation[,4])
fdiv4 = s4[2:n] - s4[1:(n-1)]

s5 = sort(pc$rotation[,5])
fdiv5 = s5[2:n] - s5[1:(n-1)]


plot(sort(fdiv1, decreasing=T), type="l",col="gray",ylim=c(0,0.00001),xlim=c(0,2000))
lines(sort(fdiv2, decreasing=T), type="l",col="blue")
lines(sort(fdiv3, decreasing=T), type="l",col="gray")
lines(sort(fdiv4, decreasing=T), type="l",col="gray")
lines(sort(fdiv5, decreasing=T), type="l",col="gray")



#plot(sort(-pc$rotation[,2]), type="l")
#abline(h=0, col="red")

#plot(sort(pc$rotation[,3]), type="l")
#abline(h=0, col="red")

#plot(sort(pc$rotation[,4]), type="l")
#abline(h=0, col="red")

plot(sort(pc$rotation[,5]), type="l")
abline(h=0, col="red")





## Figure 2C: logit ----


stats <- metadata |> 
  dplyr::filter(resection_tumor_grade %in% c(2,3)) |> 
  dplyr::mutate(resection_tumor_grade__hg = ifelse(resection_tumor_grade == 3, 1 , 0)) |> 
  dplyr::mutate(resection = ifelse(resection_number == 1, "primary", "recurrent")) |> 
  dplyr::mutate(resection_recurrent = ifelse(resection == "primary", 0 , 1))



### logistic GLASS-NL x grade ----


stats <- stats |> 
  dplyr::mutate(covar = array_GLASS_NL_g2_g3_sig) |> 
  dplyr::mutate(covar_name = "GLASS-NL primary-recurrent signature")


model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit.restyled <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class, resection_recurrent) |> 
    dplyr::mutate(col = resection_tumor_grade__hg + 2) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_recurrent - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class, resection_recurrent) |> 
    dplyr::mutate(col = resection_tumor_grade__hg + 2) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_recurrent + 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(resection_tumor_grade__hg = 0) |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_recurrent) |> 
    dplyr::mutate(col = resection_recurrent)
)


p1 <- ggplot(plt.logit.restyled, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.restyled |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=2) +
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(0, 1), breaks=c(0, 1), labels=c("Primary", "Recurrent"), oob = scales::squish) +
  labs(col=NULL) +
  ggnewscale::new_scale_colour() +
  geom_line(data = plt.logit.restyled |> dplyr::filter(type == "data"),
            col="white",
            lwd=theme_nature_lwd * 2, alpha=0.65
  ) +
  geom_line(data = plt.logit.restyled |> dplyr::filter(type == "data"),
            aes(col=col),
            lwd=theme_nature_lwd) +
  scale_color_gradientn(colours = c("aquamarine3","aquamarine3","#d34394ff","#d34394ff"), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 2.50, 3), labels=c("Grade 2","","", "Grade 3"), oob = scales::squish) +
  theme_nature +
  annotate("text", y = modelr::seq_range(stats$covar, 8)[2], 
                   x = 0.1, label = paste0("p = ",format.pval(pval)), size=theme_nature_size) +
  labs(col=NULL, y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = "GLASS-NL MedMeth x WHO Grade") +
  scale_x_continuous(breaks = c(0,1), labels=c("Grade 2", "Grade 3")) + 
  theme(legend.box = "vertical") + # space dependent
  theme(legend.key.size = unit(0.6, 'lines'))
p1



### logistic GLASS-NL x resection ----


stats <- stats |> 
  dplyr::mutate(covar = array_GLASS_NL_g2_g3_sig) |> 
  dplyr::mutate(covar_name = "GLASS-NL primary-recurrent signature")


model <- glm(resection_recurrent ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$resection_recurrent = predict(model, Predicted_data, type="response")


plt.logit.restyled <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class, resection_recurrent) |> 
    dplyr::mutate(col = resection_tumor_grade__hg + 2) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_recurrent - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class, resection_recurrent) |> 
    dplyr::mutate(col = resection_tumor_grade__hg + 2) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = resection_recurrent + 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(resection_tumor_grade__hg = 0) |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_recurrent) |> 
    dplyr::mutate(col = resection_recurrent)
)


p2 <- ggplot(plt.logit.restyled, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.restyled |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=2) +
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(0, 1), breaks=c(0, 1), labels=c("Primary", "Recurrent"), oob = scales::squish) +
  labs(col=NULL) +
  ggnewscale::new_scale_colour() +
  geom_line(data = plt.logit.restyled |> dplyr::filter(type == "data"),
            col="white",
            lwd=theme_nature_lwd * 2, alpha=0.65
  ) +
  geom_line(data = plt.logit.restyled |> dplyr::filter(type == "data"),
            aes(col=col),
            lwd=theme_nature_lwd) +
  scale_color_gradientn(colours = c("aquamarine3","aquamarine3","#d34394ff","#d34394ff"), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 2.50, 3), labels=c("Grade 2","","", "Grade 3"), oob = scales::squish) +
  theme_nature +
  annotate("text", y = modelr::seq_range(stats$covar, 8)[2], x = 0.2, label = paste0("p = ",format.pval(pval)), size=theme_nature_size) +
  labs(col=NULL, y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = "GLASS-NL MedMeth x Resection") +
  scale_x_continuous(breaks = c(0,1),
                     labels=c("Primary", "Recurrent")) + 
  theme(legend.box = "vertical") + # space dependent
  theme(legend.key.size = unit(0.6, 'lines'))
p2 



### logistic AcCGAP x grade ----


stats <- stats |> 
  dplyr::mutate(covar = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(covar_name = "Astrocytoma CGC Lasso")


model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit.restyled <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class, resection_recurrent) |> 
    dplyr::mutate(col = resection_recurrent + 1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = ifelse(resection_tumor_grade__hg == 0, resection_tumor_grade__hg, resection_tumor_grade__hg - 0.175)) |> 
    dplyr::rename(y = covar),
  
  stats |>  # right point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class, resection_recurrent) |> 
    dplyr::mutate(col = resection_recurrent + 1) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = ifelse(resection_tumor_grade__hg == 1, resection_tumor_grade__hg, resection_tumor_grade__hg + 0.175)) |> 
    dplyr::rename(y = covar),
  
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(resection_recurrent = 0) |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_tumor_grade__hg) |> 
    dplyr::mutate(col = resection_tumor_grade__hg)
)


plt.logit.restyled.grade <- plt.logit.restyled


p3 <- ggplot(plt.logit.restyled, aes(x=x, y=y, group=group, col=col)) +
  geom_line(data = plt.logit.restyled |> dplyr::filter(type == "data"),
            aes(col=col),
            lwd=theme_nature_lwd) +
  scale_color_gradientn(colours = c("aquamarine3","aquamarine3","#d34394ff","#d34394ff"), na.value = "grey50", limits = c(1, 2), breaks=c(1, 1.50, 1.50, 2), labels=c("primary","","", "recurrent"), oob = scales::squish) +
  labs(col=NULL) +
  ggnewscale::new_scale_colour() +
  geom_line(data = plt.logit.restyled |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=1.0) +
  theme_nature +
  scale_fill_manual(values = palette_mnp_12.8_6) +
  annotate("text", y = modelr::seq_range(stats$covar, 16)[15], x = 0.375, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(0, 1), breaks=c(0, 1), labels=c("Grade 2", "Grade 3"), oob = scales::squish) +
  labs(col=NULL, y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = "CGC Lasso x WHO Grade") +
  scale_x_continuous(breaks = c(0,1), labels=c("Grade 2", "Grade 3")) + 
  theme(legend.box = "vertical") + # space dependent
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
p3




### logistic AcCGAP x resection ----


stats <- stats |> 
  dplyr::mutate(covar = array_A_IDH_HG__A_IDH_LG_lr__lasso_fit) |> 
  dplyr::mutate(covar_name = "Astrocytoma CGC Lasso")


model <- glm(resection_recurrent ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


expnd <- (max(stats$covar)-min(stats$covar)) * 0.075
Predicted_data <- data.frame(covar=modelr::seq_range(c(min(stats$covar) - expnd, max(stats$covar) + expnd), 500))
Predicted_data$resection_recurrent = predict(model, Predicted_data, type="response")


plt.logit.restyled <- rbind(
  stats |>  # left point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class, resection_recurrent) |> 
    dplyr::mutate(col = resection_tumor_grade__hg + 2) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    dplyr::mutate(x = ifelse(resection_recurrent == 0, resection_recurrent, resection_recurrent - 0.175)) |> 
    #dplyr::mutate(x = resection_recurrent - 0.06 ) |> 
    dplyr::rename(y = covar)
  ,
  stats |>  # right point line:
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class, resection_recurrent) |> 
    dplyr::mutate(col = resection_tumor_grade__hg + 2) |> 
    dplyr::mutate(group = paste0("id",1:dplyr::n())) |> 
    dplyr::mutate(type = "data") |> 
    #dplyr::mutate(x = resection_recurrent + 0.06 ) |> 
    dplyr::mutate(x = ifelse(resection_recurrent == 1, resection_recurrent, resection_recurrent + 0.175)) |>
    dplyr::rename(y = covar)
  ,
  Predicted_data |> # logit fit
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(resection_tumor_grade__hg = 0) |> 
    dplyr::mutate(group = "logit fit") |> 
    dplyr::mutate(type = "fit") |> 
    dplyr::rename(y = covar) |> 
    dplyr::mutate(x = resection_recurrent) |> 
    dplyr::mutate(col = resection_recurrent)
)

plt.logit.restyled.resection <- plt.logit.restyled

p4 <- ggplot(plt.logit.restyled, aes(x=x, y=y, group=group, col=col)) +
  # geom_line(data = plt.logit.restyled |> dplyr::filter(type == "data"),
  #           col="white",
  #           lwd=theme_nature_lwd * 2, alpha=0.65
  # ) +
  geom_line(data = plt.logit.restyled |> dplyr::filter(type == "data"),
            aes(col=col),
            lwd=theme_nature_lwd) +
  scale_color_gradientn(colours = c("aquamarine3","aquamarine3","#d34394ff","#d34394ff"), na.value = "grey50", limits = c(2, 3), breaks=c(2, 2.50, 2.50, 3), labels=c("Grade 2","","", "Grade 3"), oob = scales::squish) +
  labs(col=NULL) +
  ggnewscale::new_scale_colour() +
  geom_line(data = plt.logit.restyled |> dplyr::filter(type == "fit") ,
            aes(col=col),
            lwd=1.0) +
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(0, 1), breaks=c(0, 1), labels=c("Primary", "Recurrent"), oob = scales::squish) +
  theme_nature +
  annotate("text", y = modelr::seq_range(stats$covar, 16)[15], x = 0.375, label = paste0("p = ",format.pval(pval, digits=3)), size=theme_nature_size) +
  labs(col=NULL, y= stats |> dplyr::pull(covar_name) |> unique(), fill=NULL, x=NULL, subtitle = "CGC Lasso x Resection") +
  scale_x_continuous(breaks = c(0,1), labels=c("Primary", "Recurrent")) + 
  theme(legend.box = "vertical") + # space dependent
  theme(legend.key.size = unit(0.6, 'lines'))
p4


library(patchwork)
p3 + p4 + plot_layout(ncol=2)

ggsave("output/figures/vis_LGC_PCA_logistic.pdf", width=(8.5*0.95)*(1/3), height=3)




#### t-test ----

t.test(
  stats |>
    dplyr::filter(resection_tumor_grade == 2) |> 
    dplyr::pull(covar),
  stats |>
    dplyr::filter(resection_tumor_grade == 3) |> 
    dplyr::pull(covar)
)




### logistic PC2 ----



scaling_factor_x <- 0.0185 * 50 * (0.5/0.26)
stats <- stats |> 
  dplyr::mutate(covar = array_PC2) |> 
  dplyr::mutate(covar_name = "PC2")



model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


Predicted_data <- data.frame(covar=modelr::seq_range(stats$covar, 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  stats |> 
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data"),
  Predicted_data |> 
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2) |> 
  dplyr::mutate(xmin=covar - scaling_factor_x,
                xmax=covar + scaling_factor_x,
                xmin_wash = covar - (2.25 * scaling_factor_x), 
                xmax_wash = covar + (2.25 * scaling_factor_x), 
                ymin=`resection tumor grade` - 0.06,
                ymax=`resection tumor grade` + 0.06)



p2 <- ggplot(plt.logit, aes(x = covar, y=`resection tumor grade`)) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), aes(col=`resection tumor grade`),lwd=2) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin_wash,xmax=xmax_wash,ymin=ymin,ymax=ymax), fill=alpha("white",0.65)) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=`array_mnp_predictBrain_v12.8_cal_class`)) +
  annotate("text", x = modelr::seq_range(stats$covar, 8)[2], y = 2.5, label = paste0("p = ",format.pval(pval)), size=theme_cellpress_size) + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col=NULL, x = stats |> dplyr::pull(covar_name) |> unique(), fill=NULL) +
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 3), oob = scales::squish) +
  scale_fill_manual(values = palette_mnp_12.8_6) +
  theme_cellpress +
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
p2


#### t-test ----


t.test(
  stats |>
    dplyr::filter(resection_tumor_grade == 2) |> 
    dplyr::pull(covar),
  stats |>
    dplyr::filter(resection_tumor_grade == 3) |> 
    dplyr::pull(covar)
)




### logistic Zhang2017 Mortality risk clock ----




scaling_factor_x <- 0.002*5*(0.5/0.83)
stats <- stats |> 
  dplyr::mutate(covar = array_dnaMethyAge__ZhangY2017) |> 
  dplyr::mutate(covar_name = "Zhang 2017 Mortality risk clock")



model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


Predicted_data <- data.frame(covar=modelr::seq_range(stats$covar, 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  stats |> 
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data"),
  Predicted_data |> 
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2) |> 
  dplyr::mutate(xmin=covar - scaling_factor_x,
                xmax=covar + scaling_factor_x,
                xmin_wash = covar - (2.25 * scaling_factor_x), 
                xmax_wash = covar + (2.25 * scaling_factor_x), 
                ymin=`resection tumor grade` - 0.06,
                ymax=`resection tumor grade` + 0.06)



p3 <- ggplot(plt.logit, aes(x = covar, y=`resection tumor grade`)) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), aes(col=`resection tumor grade`),lwd=2) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin_wash,xmax=xmax_wash,ymin=ymin,ymax=ymax), fill=alpha("white",0.65)) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=`array_mnp_predictBrain_v12.8_cal_class`)) +
  annotate("text", x = modelr::seq_range(stats$covar, 8)[2], y = 2.5, label = paste0("p = ",format.pval(pval)), size=theme_cellpress_size) + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col=NULL, x = stats |> dplyr::pull(covar_name) |> unique(), fill=NULL) +
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 3), oob = scales::squish) +
  scale_fill_manual(values = palette_mnp_12.8_6) +
  theme_cellpress +
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
p3



#### t-test ----


t.test(
  stats |>
    dplyr::filter(resection_tumor_grade == 2) |> 
    dplyr::pull(covar),
  stats |>
    dplyr::filter(resection_tumor_grade == 3) |> 
    dplyr::pull(covar)
)




### logistic epiTOC2 ----



scaling_factor_x <- 16 * (0.5 / 0.58)
stats <- stats |> 
  dplyr::mutate(covar = array_epiTOC2_tnsc) |> 
  dplyr::mutate(covar_name = "epiTOC2")



model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


Predicted_data <- data.frame(covar=modelr::seq_range(stats$covar, 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  stats |> 
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data"),
  Predicted_data |> 
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2) |> 
  dplyr::mutate(xmin=covar - scaling_factor_x,
                xmax=covar + scaling_factor_x,
                xmin_wash = covar - (2.25 * scaling_factor_x), 
                xmax_wash = covar + (2.25 * scaling_factor_x), 
                ymin=`resection tumor grade` - 0.06,
                ymax=`resection tumor grade` + 0.06)



p4 <- ggplot(plt.logit, aes(x = covar, y=`resection tumor grade`)) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), aes(col=`resection tumor grade`),lwd=2) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin_wash,xmax=xmax_wash,ymin=ymin,ymax=ymax), fill=alpha("white",0.65)) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=`array_mnp_predictBrain_v12.8_cal_class`)) +
  annotate("text", x = modelr::seq_range(stats$covar, 8)[2], y = 2.5, label = paste0("p = ",format.pval(pval)), size=theme_cellpress_size) + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col=NULL, x = stats |> dplyr::pull(covar_name) |> unique(), fill=NULL) +
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 3), oob = scales::squish) +
  scale_fill_manual(values = palette_mnp_12.8_6) +
  theme_cellpress +
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
p4



#### t-test ----


t.test(
  stats |>
    dplyr::filter(resection_tumor_grade == 2) |> 
    dplyr::pull(covar),
  stats |>
    dplyr::filter(resection_tumor_grade == 3) |> 
    dplyr::pull(covar)
)




### logistic median.overall.methylation ----



scaling_factor_x <- 0.0185 * 0.16 * (0.5/0.48)
stats <- stats |> 
  dplyr::mutate(covar = array_median.overall.methylation) |> 
  dplyr::mutate(covar_name = "Median methylation")



model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


Predicted_data <- data.frame(covar=modelr::seq_range(stats$covar, 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  stats |> 
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data"),
  Predicted_data |> 
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2) |> 
  dplyr::mutate(xmin=covar - scaling_factor_x,
                xmax=covar + scaling_factor_x,
                xmin_wash = covar - (2.25 * scaling_factor_x), 
                xmax_wash = covar + (2.25 * scaling_factor_x), 
                ymin=`resection tumor grade` - 0.06,
                ymax=`resection tumor grade` + 0.06)



p3 <- ggplot(plt.logit, aes(x = covar, y=`resection tumor grade`)) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), aes(col=`resection tumor grade`),lwd=2) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin-(5 * scaling_factor_x),xmax=xmax+(5 * scaling_factor_x),ymin=ymin,ymax=ymax), fill=alpha("white",0.65)) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=`array_mnp_predictBrain_v12.8_cal_class`)) +
  annotate("text", x = modelr::seq_range(stats$covar, 8)[2], y = 2.5, label = paste0("p = ",format.pval(pval)), size=theme_cellpress_size) + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col=NULL, x = stats |> dplyr::pull(covar_name) |> unique(), fill=NULL) +
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 3), oob = scales::squish) +
  scale_fill_manual(values = palette_mnp_12.8_6) +
  theme_cellpress +
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
p3



### logistic median.glass_nl_supervised.methylation ----



scaling_factor_x <- 0.0185 * 0.22 * (0.5/0.39)
stats <- stats |> 
  dplyr::mutate(covar = array_median.glass_nl_supervised.methylation) |> 
  dplyr::mutate(covar_name = "Median methylation (AC primary - recurrence DMP's)")



model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


Predicted_data <- data.frame(covar=modelr::seq_range(stats$covar, 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  stats |> 
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data"),
  Predicted_data |> 
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2) |> 
  dplyr::mutate(xmin=covar - scaling_factor_x,
                xmax=covar + scaling_factor_x,
                xmin_wash = covar - (2.25 * scaling_factor_x), 
                xmax_wash = covar + (2.25 * scaling_factor_x), 
                ymin=`resection tumor grade` - 0.06,
                ymax=`resection tumor grade` + 0.06)



p4 <- ggplot(plt.logit, aes(x = covar, y=`resection tumor grade`)) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), aes(col=`resection tumor grade`),lwd=2) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin_wash,xmax=xmax_wash,ymin=ymin,ymax=ymax), fill=alpha("white",0.65)) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=`array_mnp_predictBrain_v12.8_cal_class`)) +
  annotate("text", x = modelr::seq_range(stats$covar, 8)[2], y = 2.5, label = paste0("p = ",format.pval(pval)), size=theme_cellpress_size) + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col=NULL, x = stats |> dplyr::pull(covar_name) |> unique(), fill=NULL) +
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 3), oob = scales::squish) +
  scale_fill_manual(values = palette_mnp_12.8_6) +
  theme_cellpress +
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
p4





### logistic hypoSC solo-WCGW ----



scaling_factor_x <- 0.00035 * (0.5 / 0.38)
stats <- stats |> 
  dplyr::mutate(covar = array_epiTOC2_hypoSC) |> 
  dplyr::mutate(covar_name = "solo-WCGW HypoClock score")



model <- glm(resection_tumor_grade__hg ~ covar, data = stats, family = binomial)
pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "covar") |> 
  dplyr::pull(`Pr(>|z|)`)


Predicted_data <- data.frame(covar=modelr::seq_range(stats$covar, 500))
Predicted_data$resection_tumor_grade__hg = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  stats |> 
    dplyr::select(resection_tumor_grade__hg, covar, array_mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data"),
  Predicted_data |> 
    dplyr::mutate(array_mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
) |> 
  dplyr::mutate(`resection tumor grade` = resection_tumor_grade__hg + 2) |> 
  dplyr::mutate(xmin=covar - scaling_factor_x,
                xmax=covar + scaling_factor_x,
                xmin_wash = covar - (2.25 * scaling_factor_x), 
                xmax_wash = covar + (2.25 * scaling_factor_x),
                ymin=`resection tumor grade` - 0.06,
                ymax=`resection tumor grade` + 0.06)




p6 <- ggplot(plt.logit, aes(x = covar, y=`resection tumor grade`)) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), aes(col=`resection tumor grade`),lwd=2) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin_wash,xmax=xmax_wash,ymin=ymin,ymax=ymax), fill=alpha("white",0.65)) +
  geom_rect(data = plt.logit |> dplyr::filter(type == "data"), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=`array_mnp_predictBrain_v12.8_cal_class`)) +
  annotate("text", x = modelr::seq_range(stats$covar, 8)[2], y = 2.5, label = paste0("p = ",format.pval(pval)), size=theme_cellpress_size) + 
  scale_y_continuous(breaks = c(2,3)) + 
  labs(col=NULL, x = stats |> dplyr::pull(covar_name) |> unique(), fill=NULL) +
  scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", limits = c(2, 3), breaks=c(2, 3), oob = scales::squish) +
  scale_fill_manual(values = palette_mnp_12.8_6) +
  theme_cellpress +
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
p6



### combined export ----

library(patchwork)

(p1 + p2) /
(p3 + p4) +
  plot_annotation(theme = theme_cellpress,
                  subtitle=format_subtitle("Logit WHO grade"))


ggsave("output/figures/vis_LGC_x_PCA__logits.pdf", width = (8.5 * 0.975), height = 4)



### multi LM ----




model <- lm(resection_tumor_grade__hg ~ PC2 + A_IDH_HG__A_IDH_LG_lr__lasso_fit, data = stats)
model.g <- glm(resection_tumor_grade__hg ~ PC2 + A_IDH_HG__A_IDH_LG_lr__lasso_fit, data = stats,family = binomial)
summary(model)
summary(model.g)
sjPlot::plot_model(model)


model2 <- lm(A_IDH_HG__A_IDH_LG_lr__lasso_fit ~ resection_tumor_grade__hg + batch_us, data = stats)
summary(model2)
sjPlot::plot_model(model2)

model3 <- lm(-PC2 ~ resection_tumor_grade__hg + batch_us, data = stats)
summary(model3)
sjPlot::plot_model(model3)


model3 <- lm(PC2 ~ resection_tumor_grade__hg, data = stats)
summary(model3)



# G-SAM ----



metadata <- gsam.metadata.array_samples |> 
  filter_GSAM_idats(73) |> 
  dplyr::mutate(resection__rec = ifelse(resection == "R1",0,1))


data <- data.mvalues.hq_samples |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  tibble::column_to_rownames('probe_id') |> 
  dplyr::select(metadata$sentrix_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (695840))
    return(.)
  })() |> 
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |> # this synthax, oh my
  dplyr::arrange(mad) |> 
  dplyr::mutate(mad = NULL)


## logistic LR & recurrence ----



ggplot(metadata, aes(x=resection, y=A_IDH_HG__A_IDH_LG_lr, group=patient)) +
  ggbeeswarm::geom_quasirandom() +
  geom_line(col="black", alpha=0.25, lwd=0.25) +
  theme_bw()


model <- glm(resection__rec ~  A_IDH_HG__A_IDH_LG_lr, 
             data = metadata,
             family = binomial)
summary(model)


pval <- model |>
  summary() |> 
  purrr::pluck('coefficients') |>
  as.data.frame() |>  
  tibble::rownames_to_column('coef') |>
  dplyr::filter(coef == "A_IDH_HG__A_IDH_LG_lr") |> 
  dplyr::pull(`Pr(>|z|)`)


Predicted_data <- data.frame(A_IDH_HG__A_IDH_LG_lr=modelr::seq_range(metadata$A_IDH_HG__A_IDH_LG_lr, 500))
Predicted_data$resection__rec = predict(model, Predicted_data, type="response")


plt.logit <- rbind(
  metadata |> 
    dplyr::select(resection__rec, A_IDH_HG__A_IDH_LG_lr, mnp_predictBrain_v12.8_cal_class) |> 
    dplyr::mutate(type = "data")
  ,
  Predicted_data |> 
    dplyr::mutate(mnp_predictBrain_v12.8_cal_class = "") |> 
    dplyr::mutate(type = "fit")
)


ggplot(plt.logit, aes(x = A_IDH_HG__A_IDH_LG_lr, y=`resection__rec`, col= mnp_predictBrain_v12.8_cal_class)) +
  geom_point(data = plt.logit |> dplyr::filter(type == "data"),pch="|",size=5) +
  geom_line(data = plt.logit |> dplyr::filter(type == "fit"), col="darkgreen") +
  theme_bw() + 
  annotate("text", x = 2.5, y = 0.5, label = paste0("p = ",format.pval(pval))) + 
  theme(legend.position = 'bottom') + 
  scale_y_continuous(breaks = c(0,1)) + 
  labs(col="mnp v12.8", y = "Resection (0 primary, 1 recurrent)")



## survival ----
### R1 ----


ggplot(
  metadata |> 
    dplyr::filter(resection == "R1"),
  aes(y= survivalDays , x=A_IDH_HG__A_IDH_LG_lr, shape=status)) + 
  geom_point()



### R2 ----

#plot(sort(stats$A_IDH_HG__A_IDH_LG_lr__lasso_fit)) + abline(h=4.65)


stats <-  metadata |> 
  dplyr::filter(resection == "R2") |> 
  dplyr::mutate(A_IDH_HG__A_IDH_LG__cut = cut(A_IDH_HG__A_IDH_LG_lr__lasso_fit, breaks=2)) |> 
  dplyr::mutate(LR_status_median = ifelse(A_IDH_HG__A_IDH_LG_lr__lasso_fit > 6.1, "LR high", "LR low")) |> 
  dplyr::mutate(LR_status_L1 = ifelse(A_IDH_HG__A_IDH_LG_lr__lasso_fit > 4.7, "LR high", "LR low")) |> 
  dplyr::mutate(LR_status_L2 = ifelse(A_IDH_HG__A_IDH_LG_lr__lasso_fit > 4.1, "LR high", "LR low"))


ggplot(stats, aes(y= survivalFromSecondSurgeryDays , x=A_IDH_HG__A_IDH_LG_lr, shape=status)) + 
  geom_point()


surv_object <- survival::Surv(time = stats$survivalFromSecondSurgeryDays, event=stats$os.event)
#surv_object <- survival::Surv(time = stats$survivalDays, event=stats$os.event)


fit1 <- survival::survfit(surv_object ~  A_IDH_HG__A_IDH_LG__cut , data = stats)
print(survminer::surv_pvalue(fit1)$pval)

survminer::ggsurvplot(fit1, data = stats, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                      xlab="Survival from recurrence")


fit2 <- survival::survfit(surv_object ~  LR_status_median , data = stats)
print(survminer::surv_pvalue(fit2)$pval)

survminer::ggsurvplot(fit2, data = stats, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                      xlab="Survival from recurrence")

fit3 <- survival::survfit(surv_object ~  LR_status_L1 , data = stats)
print(survminer::surv_pvalue(fit3)$pval)

survminer::ggsurvplot(fit3, data = stats, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                      xlab="Survival from recurrence")

fit4 <- survival::survfit(surv_object ~  LR_status_L2 , data = stats)
print(survminer::surv_pvalue(fit4)$pval)

survminer::ggsurvplot(fit4, data = stats, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                      xlab="Survival from recurrence")




fit_l <- survival::coxph(surv_object ~  A_IDH_HG__A_IDH_LG_lr__lasso_fit , data = stats)
summary(fit_l)




# GLASS-NL / AC ----


metadata <- glass_nl.metadata.array_samples |> 
  dplyr::filter(!qc.pca.detP.outlier) |>
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (218))
    return(.)
  })()


## PCA ----


masked <- readRDS('cache/mvalues.HQ_samples.detP_mask.Rds') |> 
  is.na() |> 
  rowSums() |> 
  as.data.frame() |> 
  dplyr::rename(count = 1) |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::filter(count > 0) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (66116))
    return(.)
  })() |> 
  dplyr::pull('probe_id')




data <- readRDS('cache/mvalues.HQ_samples.Rds') 

data <- data |> 
  dplyr::select(metadata$sentrix_id)

data <- data |> 
  tibble::rownames_to_column('probe_id')

data <- data |> 
  dplyr::filter(probe_id %in% (masked) == F) |> 
  tibble::column_to_rownames('probe_id')

dim(data)


data <- data |> 
  (function(.) dplyr::mutate(., mad =  apply( ., 1, stats::mad)) )() |> # this synthax, oh my
  dplyr::arrange(mad) |> 
  dplyr::mutate(mad = NULL)


data.all.pca.obj <- data |> 
  #dplyr::slice_head(n=100000) |> 
  t() |> 
  prcomp()

data.all.pca <- data.all.pca.obj |> 
  purrr::pluck('x') |> 
  as.data.frame(stringsAsFactors=F) |> 
  dplyr::select(paste0("PC",1:40)) |> 
  tibble::rownames_to_column('sentrix_id')


plt <- metadata |> 
  dplyr::left_join(data.all.pca, by=c('sentrix_id'='sentrix_id'), suffix=c('',''))


plt <- plt |> 
  dplyr::mutate(mnp_predictBrain_v12.8_cal_A_IDH_LG_odds = (mnp_predictBrain_v12.8_cal_A_IDH_LG/100) / (1 - (mnp_predictBrain_v12.8_cal_A_IDH_LG/100))) |> 
  dplyr::mutate(mnp_predictBrain_v12.8_cal_A_IDH_HG_odds = (mnp_predictBrain_v12.8_cal_A_IDH_HG/100) / (1 - (mnp_predictBrain_v12.8_cal_A_IDH_HG/100))) |> 
  dplyr::mutate(A_IDH_HG__A_IDH_LG_lr = log(mnp_predictBrain_v12.8_cal_A_IDH_HG_odds / mnp_predictBrain_v12.8_cal_A_IDH_LG_odds))


ggplot(plt, aes(x=A_IDH_HG__A_IDH_LG_lr, y=PC1)) +
  geom_point()




# detP mask filter
# order on MAD
# pick top 200.000
# prcomp


# plot(plt$PC2, plt$median.overall.methylation)
# plot(plt$PC2, plt$A_IDH_HG__A_IDH_LG_lr)
# plot(plt$PC2, plt$A_IDH_HG__A_IDH_LG_lr__lasso_fit)
# plot(plt$median.overall.methylation, plt$A_IDH_HG__A_IDH_LG_lr__lasso_fit)
# plt$A_IDH_HG__O_IDH_lr, plt$A_IDH_HG__A_IDH_LG_lr)


# c <- plt |>
#   dplyr::mutate(mnp_predictBrain_v12.8_cal_A_IDH_HG = -log(mnp_predictBrain_v12.8_cal_A_IDH_HG)) |> 
#   dplyr::select(
#   PC2,
#   A_IDH_HG__A_IDH_LG_lr,
#   A_IDH_HG__A_IDH_LG_lr__lasso_fit,
#   A_IDH_HGoligsarc__O_IDH_lr,
#   median.overall.methylation,
#   median.glass_nl_supervised.methylation,
#   A_IDH_HG__O_IDH_lr
#   #mnp_predictBrain_v12.8_cal_A_IDH_HG
# ) |> 
#   as.matrix()
# 
# 
corrplot::corrplot(abs(cor(c, method="spearman")), order="hclust")
corrplot::corrplot(abs(cor(c, method="kendall")), order="hclust")
# 
# 
# plot(plt$A_IDH_HGoligsarc__O_IDH_lr, plt$PC2)
# plot(plt$A_IDH_HG__O_IDH_lr, plt$PC2)
# 
# plot(plt$A_IDH_HGoligsarc__O_IDH_lr, plt$PC2)
# plot(plt$PC2)
# 
