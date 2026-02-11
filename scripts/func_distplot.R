#!/usr/bin/env R


# https://yunuuuu.github.io/ggalign/


distplot <- function(d) {

  o.join <- data.frame(name = attr(d, 'Labels')) |> 
    #dplyr::arrange(name) |> 
    dplyr::mutate(i = attr(d, 'Size'):1)
  
  if(!is.null(attr(d, "Facet"))) {
    o.join$Facet <- as.character(attr(d, "Facet"))
  }
  
  m <- as.matrix(d) |> 
    as.data.frame()
  
  
  plt.gg <- reshape2::melt(as.matrix(d), varnames = c("row", "col")) |> 
    dplyr::left_join(o.join |> dplyr::select(name, i) |> dplyr::rename(x.order = i), by=c('row' = 'name')) |> 
    dplyr::left_join(o.join |> dplyr::select(name, i) |> dplyr::rename(y.order = i), by=c('col' = 'name')) |>
    dplyr::mutate(size = (max(value) - value) / max(value)) |> 
    dplyr::mutate(size = 0.2 + (0.4*size) * 0.5)
  
  if(!is.null(attr(d, "Facet"))) {
    plt.gg <- plt.gg |> 
      dplyr::left_join(o.join |> dplyr::select(name, Facet), by=c('row' = 'name')) # or col?
  }
  
  
  ggplot2::ggplot(plt.gg, ggplot2::aes(
    x = .data$x.order,
    y = .data$y.order,
    radius = size,
    
    # [0.3 , 0.8] + 0.2 smoothened from lwd/border
    fill =value,
    col = value,
    label = .data$row
  )) +
    
    ggplot2::geom_tile(col = "darkgray", fill = "white", lwd = theme_nature_lwd/2) +
    ggplot2::geom_rect(data=head(plt.gg, n=1), aes(xmin=0.5,xmax=nrow(o.join)+0.5,ymin=0.5,ymax=nrow(o.join)+0.5), col = "gray", fill = NA, lwd = theme_nature_lwd) +
    ggplot2::scale_fill_gradientn(colours = col2(201)[1:101], na.value = "grey50", limits = c(0, max(plt.gg$value)), guide = "colourbar") + # guide = "colourbar",
    ggplot2::scale_color_gradientn(colours = col2(201)[1:101], na.value = "grey50", limits = c(0, max(plt.gg$value)), guide = "none") +
    recursiveCorPlot::geom_circle(radius.fixed = T) + # from THIS repo
    ggplot2::labs(y = NULL, x = NULL, main = NULL, fill="Eucl. dist B-value") +
    ggplot2::coord_fixed() +
    ggplot2::scale_x_discrete(labels = NULL, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(name = NULL, breaks = nrow(o.join):1, labels = o.join$name, expand = c(0, 0)) +
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
      legend.key.size = unit(0.3, 'lines'),
    )
    #theme(panel.spacing = unit(0.5, "lines")) # for facettes but this works poorly with coord_fixed
}



