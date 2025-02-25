#!/usr/bin/env R

a <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_prcomp.Rds")
b <- readRDS("cache/analysis_unsupervised_PCA_GLASS-OD_x.Rds")


polycomb_probes <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & grepl("^cg", probe_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::rename(gene = GencodeCompV12_NAME) |> 
  dplyr::select(probe_id, gene) |>
  tibble::tibble() |> 
  tidyr::separate_longer_delim(gene, delim = ";") |> 
  dplyr::filter(gene %in% c(c(genes_polycomb_eed_homeobox, 
                              genes_polycomb_h3k27_homeobox,
                              genes_polycomb_prc2_homeobox,
                              genes_polycomb_suz12_homeobox
  ))) 


c = a$rotation[,2:3] |> 
  as.data.frame() |> 
  tibble::rownames_to_column("probe_id") |> 
  dplyr::mutate(polycomb_tf = probe_id %in% c(polycomb_probes$probe_id)) |> 
  dplyr::rename(`PC2 rotation` = PC2) |> 
  dplyr::rename(`PC3 rotation` = PC3)


ggplot(c, aes(x=PC2, y=polycomb_tf)) +
  geom_point(pch=19,cex=0.001, alpha=0.15)


ggplot(c, aes(x=PC3, y=polycomb_tf)) +
  geom_point(pch=19,cex=0.001, alpha=0.15)


ggplot(c, aes(PC2, fill=polycomb_tf)) +
  geom_density(alpha=0.25)
#ggpubr::stat_compare_means()

ggplot(c, aes(PC3, fill=polycomb_tf)) +
  geom_density(alpha=0.25)
#

ggplot(c, aes(x = polycomb_tf, y=`PC2 rotation`)) + 
  geom_violin(draw_quantiles = c(0.5), col="red") +
  ggpubr::stat_compare_means(method = "t.test", ref.group = ".all.")


ggplot(c, aes(x = polycomb_tf, y=`PC3 rotation`)) + 
  geom_violin(draw_quantiles = c(0.5), col="red") +
  ggpubr::stat_compare_means(method = "t.test", ref.group = ".all.")




factoextra::fviz_pca_biplot(a,
                            repel = F,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

factoextra::fviz_pca_biplot(a, geom="point",
                            label="none")

# plots ----



factoextra::fviz_pca_var(a, axes = c(3,2),
                         #geom="point",
                         label="none",
                         col.var = "red",
                         alpha.var  = 0.025,
                         select.var =
                           list(
                             name = polycomb_probes$probe_id #[1:1500]
                           )
) +
  xlim(-1,1) +
  ylim(-1,1) +
  coord_equal() +
  theme_nature

ggsave("output/figures/vis_PCA_loadings__fvis_pca_var__PC3_PC2.png", width=8.5*0.975/2, height=8.5*0.975/2 * 1.2, dpi=600)



factoextra::fviz_pca_var(a, axes = c(3,1),
                         #geom="point",
                         label="none",
                         col.var = "red",
                         alpha.var  = 0.025,
                         select.var =
                           list(
                             name = polycomb_probes$probe_id #[1:1500]
                           )
) +
  xlim(-1,1) +
  ylim(-1,1) +
  coord_equal() +
  theme_nature

ggsave("output/figures/vis_PCA_loadings__fvis_pca_var__PC3_PC1.png", width=8.5*0.975/2, height=8.5*0.975/2 * 1.2, dpi=600)



# facto_summarize ----


ct <- factoextra::facto_summarize(a, element ="var")
ct2 <- factoextra::facto_summarize(a, element ="var", axes=1:5)



plt <- ct2 |> 
  as.data.frame() |> 
  tibble::rownames_to_column('probe_id') |> 
  dplyr::mutate(polycomb_tf = probe_id %in% c(polycomb_probes$probe_id))



ggplot(plt, aes(x = polycomb_tf, y=`Dim.2`)) + 
  geom_violin(draw_quantiles = c(0.5), col="red") +
  ggpubr::stat_compare_means(method = "wilcox.test")


ggplot(plt, aes(x = polycomb_tf, y=`Dim.3`)) + 
  geom_violin(draw_quantiles = c(0.5), col="red") +
  ggpubr::stat_compare_means(method = "t.test")


ggplot(plt, aes(x = polycomb_tf, y=`Dim.3`)) + 
  geom_violin(draw_quantiles = c(0.5), col="red") +
  ggpubr::stat_compare_means(method = "t.test")




# test ----


data(iris)
d = prcomp(iris |>  dplyr::select( Sepal.Length, Sepal.Width, Petal.Length, Petal.Width ))





factoextra::fviz_pca_var(d, axes = c(3,2),
                         #geom="point",
                         label="none",
                         col.var = "red",
                         pointsize = 10,
                         alpha.var  = 0.5,
                         
                         
                         #arrowsize=12,
                         pointshape=3
                         
                         #select.var =
                        #   list(
                         #    name = polycomb_probes$probe_id[1:150]
                          # )
) +
  xlim(-1,1) +
  ylim(-1,1)




# contributions ----

# https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/