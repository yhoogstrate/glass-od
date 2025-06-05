#!/usr/bin/env R

metadata <- data.frame(rn = paste0("s",1:30)) |> 
  dplyr::mutate(sid = rn) |> 
  tibble::column_to_rownames("rn") |> 
  dplyr::mutate(grade = as.factor(c(rep("G2",10), rep("G3", 20)))) |> 
  dplyr::mutate(TMZ1 = as.factor(c(rep("F",20), rep("T",10)))) |> 
  dplyr::mutate(rank = 1:dplyr::n())
metadata


data <- metadata |> 
  dplyr::mutate(cg001_int1 = c(runif(15,min=-5,max=-3),
                               runif(15,min=3,max=5)
  )) |> 
  dplyr::mutate(cg002_int2 = c(runif(10,min=-5,max=-3),
                               runif(10,min=-1,max=1),
                               runif(10,min=3,max=5)
  )) |> 
  dplyr::mutate(cg003_int3 = c(runif(10,min=-5,max=-3),
                               runif(10,min=-5,max=5),
                               runif(10,min=3,max=5)
  )) |> 
  dplyr::mutate(cg004_int1_mm = ifelse(1:dplyr::n()== 5, runif(15,min=3,max=5)  , cg001_int1) ) |> 
  dplyr::mutate(cg004_int2_mm = ifelse(1:dplyr::n()== 5, runif(15,min=3,max=5)  , cg002_int2) ) |> 
  dplyr::mutate(cg005_grad = c(runif(10,min=-5,max=-3),
                               runif(20,min=3,max=5)
  )) |> 
  dplyr::mutate(cg006_tmz1 = c(runif(20,min=-5,max=-3),
                               runif(10,min=3,max=5)
  )) |> 
  dplyr::select(starts_with("cg")) |> 
  t()



head(metadata)
head(data)



design <- model.matrix(~grade + TMZ1, data=metadata)
fit <- limma::lmFit(data, design)
fit <- limma::eBayes(fit, trend=T)
(fit$coefficients)


stats_grade <- limma::topTable(fit, n=nrow(data),
                               coef="gradeG3",
                               sort.by = "none", adjust.method="fdr") |> 
  tibble::rownames_to_column('cg_id') 


stats_TMZ1 <- limma::topTable(fit, n=nrow(data),
                               coef="TMZ1T",
                               sort.by = "none", adjust.method="fdr") |> 
  tibble::rownames_to_column('cg_id') 





plt1 <- metadata |>
  tidyr::pivot_longer(cols = c(grade, TMZ1), values_to = "val",names_to = "type")

p1 <- ggplot(plt1, aes(x=reorder(sid, rank), y=type, fill=val)) +
  geom_tile(col="black") +
  labs(x=NULL, y="metadata") +
  theme_minimal()



plt2 <- data |> 
  as.data.frame() |> 
  tibble::rownames_to_column('cg_id') |> 
  tidyr::pivot_longer(cols = -c("cg_id"), 
               names_to = "sid", 
               values_to = "value") |> 
  dplyr::mutate(rank = as.numeric(gsub("^s","", sid)))

  
p2 <- ggplot(plt2, aes(x=reorder(sid, rank), y=value)) +
  facet_grid(rows = vars(cg_id)) + 
  geom_line(group="1") +
  geom_point() +
  labs(x=NULL, y="example data (M-values)") +
  theme_minimal()



plt3 <- rbind(stats_grade |> 
  dplyr::select(cg_id, t, P.Value) |> 
  dplyr::mutate(fit = "Grade"),
stats_TMZ1 |> 
  dplyr::select(cg_id, t, P.Value) |> 
  dplyr::mutate(fit = "TMZ1")) 

plt3 <- rbind(
  plt3 |> dplyr::mutate(dot=T), 
  plt3 |> dplyr::mutate(t = 0, dot=F)
)


p3 <- ggplot(plt3, aes(x = fit, y=t, col=P.Value < 0.01)) +
  facet_grid(rows = vars(cg_id)) + 
  geom_line() +
  geom_point(data=subset(plt3, dot == T)) +
  labs(y = "fit's limma: ~grade + tmz (t-stat)", x=NULL) +
  theme_minimal()



library(patchwork)
(p1 + plot_spacer() + plot_layout(widths = c(4, 1))) / (p2 + p3 + plot_layout(widths = c(4, 1)))




