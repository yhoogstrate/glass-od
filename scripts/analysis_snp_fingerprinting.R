#!/usr/bin/env R


#!/usr/bin/env R





func_snpdistplot <- function(d) {
  
  
  o.join <- data.frame(name = attr(d, 'Labels')) |> 
    dplyr::mutate(i = attr(d, 'Size'):1)
  

  m <- as.matrix(d) |> 
    as.data.frame()
  
  
  plt.gg <- melt(as.matrix(d), varnames = c("row", "col")) |> 
    dplyr::left_join(o.join |> dplyr::rename(x.order = i), by=c('row' = 'name')) |> 
    dplyr::left_join(o.join |> dplyr::rename(y.order = i), by=c('col' = 'name')) |>
    dplyr::mutate(size = (max(value) - value) / max(value)) |> 
    dplyr::mutate(size = 0.2 + (0.4*size) * 0.5)
  
  
  ggplot2::ggplot(plt.gg, ggplot2::aes(
    x = .data$x.order,
    y = .data$y.order,
    radius = size,
    
    # [0.3 , 0.8] + 0.2 smoothened from lwd/border
    fill =value,
    col = value,
    label = .data$row
  )) +
    ggplot2::geom_tile(col = "gray", fill = "white", lwd = theme_nature_lwd/2) +
    ggplot2::geom_rect(data=head(plt.gg, n=1), aes(xmin=0.5,xmax=nrow(o.join)+0.5,ymin=0.5,ymax=nrow(o.join)+0.5), col = "gray", fill = NA, lwd = theme_nature_lwd) +
    ggplot2::scale_fill_gradientn(colours = col2(201)[1:101], na.value = "grey50", limits = c(0, max(plt.gg$value)), guide = "colourbar") + # guide = "colourbar",
    ggplot2::scale_color_gradientn(colours = col2(201)[1:101], na.value = "grey50", limits = c(0, max(plt.gg$value)), guide = "none") +
    recursiveCorPlot::geom_circle(radius.fixed = T) + # from THIS repo
    ggplot2::labs(y = NULL, x = NULL, main = NULL, fill="Eucl. dist B-value") +
    ggplot2::coord_fixed() +
    ggplot2::scale_x_discrete(labels = NULL, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(name = NULL, breaks = nrow(o.join):1, labels = o.join$name, expand = c(0, 0)) +
    theme_nature +
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
}







dat <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) |> 
  dplyr::arrange(resection_id)


run <- dat |> 
  head(n=40) |> 
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL, array_channel_red = NULL) |> 
  dplyr::mutate(Sample_Name = paste0(array_sentrix_id)) |>
  dplyr::mutate(Array = gsub("^.+_","", array_sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1", array_sentrix_id)) |> 
  dplyr::select(Basename, Sample_Name, Array, Slide)

RGSet <- minfi::read.metharray.exp(targets = as.data.frame(run), force = T) #red/green channel together




d <- minfi::getSnpBeta(RGSet) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('array_sentrix_id') |> 
  dplyr::left_join(dat |> dplyr::select(array_sentrix_id, resection_id), by=c('array_sentrix_id'='array_sentrix_id')) |> 
  tibble::column_to_rownames('resection_id') |> 
  dplyr::mutate(array_sentrix_id = NULL) |> 
  dist(diag=T, upper=T)


func_snpdistplot(d)




