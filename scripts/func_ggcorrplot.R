#!/usr/bin/env R


if ("recursiveCorPlot" %in% rownames(installed.packages()) == F) {
  devtools::install_github("yhoogstrate/recursiveCorPlot")
}




ggcorrplot <- function(c,abs=F, reorder=T, cmethod="spearman") {
  if(abs) {
    d <- abs(cor(c, method=cmethod)) 
    h <- stats::hclust(stats::as.dist(1 - d), method = "ward.D2")
    d <- cor(c, method=cmethod)
  } else {
    d <- cor(c, method=cmethod)
    h <- stats::hclust(stats::as.dist(1 - d), method = "ward.D2")
  }
  
  print(d)
  
  if(reorder) {
    o <- h$labels[h$order] |>
      base::rev()
  } else {
    o <- h$labels
    #base::rev()
  }
  
  
  plt <- d |>
    base::as.data.frame() |>
    dplyr::select(dplyr::all_of(o)) |>
    base::t() |>
    base::as.data.frame() |>
    dplyr::select(dplyr::all_of(o)) |>
    base::t() |>
    base::as.matrix()
  
  o.join <- base::data.frame(name = o, i = 1:length(o))
  
  plt.gg <- reshape2::melt(plt) |>
    dplyr::rename(y = .data$`Var1`) |>
    dplyr::rename(x = .data$`Var2`) |>
    dplyr::mutate(x = as.factor(.data$`x`)) |>
    dplyr::mutate(y = as.factor(.data$`y`)) |>
    dplyr::left_join(o.join |> dplyr::rename(x.order = .data$`i`), by = c("x" = "name")) |>
    dplyr::left_join(o.join |> dplyr::mutate(i = dplyr::n() - .data$i + 1) |> dplyr::rename(y.order = .data$i), by = c("y" = "name"))
  
  
  
  ggplot2::ggplot(plt.gg, ggplot2::aes(
    x = .data$x.order,
    y = .data$y.order,
    radius = ((abs(.data$value) * 0.725) + 0.3) / 2 - 0.05,
    # [0.3 , 0.8] + 0.2 smoothened from lwd/border
    fill = .data$value,
    col = .data$value,
    label = .data$x
  )) +
    ggplot2::geom_tile(col = "gray", fill = "white", lwd = theme_cellpress_lwd / 2) +
    ggplot2::geom_rect(data=head(plt.gg, n=1), aes(xmin=0.5, xmax=base::length(o) + 0.5, ymin=0.5, ymax=base::length(o) + 0.5), col = "gray", fill = NA, lwd = theme_cellpress_lwd ) +
    ggplot2::scale_fill_gradientn(colours = col2(200), na.value = "grey50", limits = c(-1, 1), guide = "colourbar") + # guide = "colourbar",
    ggplot2::scale_color_gradientn(colours = col2(200), na.value = "grey50", limits = c(-1, 1), guide = "none") +
    recursiveCorPlot::geom_circle(radius.fixed = T) + # from THIS repo
    ggplot2::labs(y = NULL, x = NULL, main = NULL, fill="R / Rho") +
    ggplot2::coord_fixed() +
    ggplot2::scale_x_discrete(labels = NULL, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(name = NULL, breaks = base::length(o):1, labels = o, expand = c(0, 0)) +
    theme_cellpress +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.y = ggplot2::element_text(angle = 0, hjust = 1, vjust = 0.5), # used to be [3,6] reduce font size here, should become argument
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5, color = "gray80"),
      text = ggplot2::element_text(size = 13),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank()
    ) +
    theme(
      legend.position = 'right',
      legend.key.size = unit(0.3, 'lines')
    )
}



