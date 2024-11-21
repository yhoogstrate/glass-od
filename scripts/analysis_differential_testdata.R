#!/usr/bin/env R


# data: example ----


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








# data: interaction term ----

REPS <- 8

metadata <- data.frame(
  disease = c(rep("Astro", 2*REPS), rep("Oligo", 2*REPS)),
  status  = c(rep("primary", REPS), rep("recurrent", REPS), rep("primary", REPS), rep("recurrent", REPS))
) |>
  dplyr::mutate(recurrent_oligo = ifelse(disease == "Oligo" & status == "recurrent","Yes", "No")) |> 
  dplyr::mutate(recurrent_astro = ifelse(disease == "Astro" & status == "recurrent","Yes", "No")) |> 
  dplyr::mutate(disease = factor(disease, levels=c("Astro", "Oligo"))) |>
  dplyr::mutate(status = factor(status, levels=c("primary", "recurrent"))) |> 
  dplyr::mutate(recurrent_oligo = factor(recurrent_oligo, levels=c("No", "Yes"))) |> 
  dplyr::mutate(recurrent_astro = factor(recurrent_astro, levels=c("No", "Yes"))) |> 
  dplyr::mutate(recurrent_astro_oligo = dplyr::case_when(
    status == "recurrent" & disease == "Oligo" ~ -1,
    status == "recurrent" & disease == "Astro" ~ 1,
    T ~ 0
  ))




genes <- data.frame(
  flat1 = rnorm(4 * REPS),
  flat2 = rnorm(4 * REPS),
  flat3 = rnorm(4 * REPS),
  flat4 = rnorm(4 * REPS),
  flat5 = rnorm(4 * REPS),
  flat6 = rnorm(4 * REPS),
  flat7 = rnorm(4 * REPS),
  flat8 = rnorm(4 * REPS),
  flat9 = rnorm(4 * REPS),
  flat10 = rnorm(4 * REPS),
  flat11 = rnorm(4 * REPS),
  flat12 = rnorm(4 * REPS),
  
  prim_rec_1 = c(rnorm(REPS, 0), rnorm(REPS, 1.5), rnorm(REPS, 0), rnorm(REPS, 1.5)),
  prim_rec_2 = c(rnorm(REPS, 0), rnorm(REPS, 1.5), rnorm(REPS, 0), rnorm(REPS, 1.5)),
  
  
  prim_rec_oligo1 = c(rnorm(3 * REPS, 0), rnorm(REPS, 4.5)),
  prim_rec_oligo2 = c(rnorm(3 * REPS, 0), rnorm(REPS, 4.5)),
  
  prim_rec_astro1 = c(rnorm(REPS, 0), rnorm(REPS, 4.5), rnorm(2 * REPS, 0)),
  prim_rec_astro2 = c(rnorm(REPS, 0), rnorm(REPS, 4.5), rnorm(2 * REPS, 0))
  
  
  #prim_rec__prim_rec_astro1 = c(rnorm(REPS, 0), rnorm(REPS, 3),   rnorm(REPS, 0), rnorm(REPS, 1.5)),
  #prim_rec__prim_rec_astro2 = c(rnorm(REPS, 0), rnorm(REPS, 3), rnorm(REPS, 0), rnorm(REPS, 1.5)),
  
  
  #prim_rec__prim_rec_oligo1 = c(rnorm(REPS, 0), rnorm(REPS, 1.5), rnorm(REPS, 0), rnorm(REPS, 3)),
  #prim_rec__prim_rec_oligo2 = c(rnorm(REPS, 0), rnorm(REPS, 1.5), rnorm(REPS, 0), rnorm(REPS, 3))
) |>
  round(digits=2) |> 
  t() |>
  as.data.frame()



## prim - rec ----



design <- model.matrix(~status, data=metadata)
fit <- limma::lmFit(genes, design)
fit <- limma::eBayes(fit, trend=T)

limma::topTable(fit,
                n=nrow(genes),
                coef="statusrecurrent",
                sort.by = "none",
                adjust.method="fdr") |> 
  tibble::rownames_to_column('gene') |>
  dplyr::mutate(pval.stars = gtools::stars.pval(P.Value))





## interaction ----
# equal B's and P's, mirrored T's

design <- model.matrix(~status + status:disease, data=metadata)
design <- model.matrix(~status + disease +  status*disease, data=metadata)
fit <- limma::lmFit(genes, design)
fit <- limma::eBayes(fit, trend=T)
fit$coefficients |> colnames()


limma::topTable(fit,
                n=nrow(genes),
                coef="statusrecurrent",
                sort.by = "none",
                adjust.method="fdr") |> 
  tibble::rownames_to_column('gene') |>
  dplyr::mutate(pval.stars = gtools::stars.pval(P.Value))



limma::topTable(fit,
                n=nrow(genes),
                coef="statusrecurrent:diseaseOligo",
                sort.by = "none",
                adjust.method="fdr") |> 
  tibble::rownames_to_column('gene') |>
  dplyr::mutate(pval.stars = gtools::stars.pval(P.Value))




## type 2 [oligo] ----


design <- model.matrix(~status + recurrent_oligo , data=metadata)
fit <- limma::lmFit(genes, design)
fit <- limma::eBayes(fit, trend=T)
fit$coefficients |> colnames()


limma::topTable(fit,
                n=nrow(genes),
                coef="statusrecurrent",
                sort.by = "none",
                adjust.method="fdr") |> 
  tibble::rownames_to_column('gene') |>
  dplyr::mutate(pval.stars = gtools::stars.pval(P.Value))



limma::topTable(fit,
                n=nrow(genes),
                coef="recurrent_oligoYes",
                sort.by = "none",
                adjust.method="fdr") |> 
  tibble::rownames_to_column('gene') |>
  dplyr::mutate(pval.stars = gtools::stars.pval(P.Value))


## type 2 [astro] ----


design <- model.matrix(~status + recurrent_astro , data=metadata)
fit <- limma::lmFit(genes, design)
fit <- limma::eBayes(fit, trend=T)
fit$coefficients |> colnames()


limma::topTable(fit,
                n=nrow(genes),
                coef="statusrecurrent",
                sort.by = "none",
                adjust.method="fdr") |> 
  tibble::rownames_to_column('gene') |>
  dplyr::mutate(pval.stars = gtools::stars.pval(P.Value))



limma::topTable(fit,
                n=nrow(genes),
                coef="recurrent_astroYes",
                sort.by = "none",
                adjust.method="fdr") |> 
  tibble::rownames_to_column('gene') |>
  dplyr::mutate(pval.stars = gtools::stars.pval(P.Value))



## type 3 [1 & -1] ----


design <- model.matrix(~status + disease + recurrent_astro_oligo , data=metadata)
fit <- limma::lmFit(genes, design)
fit <- limma::eBayes(fit, trend=T)
fit$coefficients |> colnames()


limma::topTable(fit,
                n=nrow(genes),
                coef="statusrecurrent",
                sort.by = "none",
                adjust.method="fdr") |> 
  tibble::rownames_to_column('gene') |>
  dplyr::mutate(pval.stars = gtools::stars.pval(P.Value))



limma::topTable(fit,
                n=nrow(genes),
                coef="diseaseOligo",
                sort.by = "none",
                adjust.method="fdr") |> 
  tibble::rownames_to_column('gene') |>
  dplyr::mutate(pval.stars = gtools::stars.pval(P.Value))



limma::topTable(fit,
                n=nrow(genes),
                coef="recurrent_astro_oligo",
                sort.by = "none",
                adjust.method="fdr") |> 
  tibble::rownames_to_column('gene') |>
  dplyr::mutate(pval.stars = gtools::stars.pval(P.Value))




