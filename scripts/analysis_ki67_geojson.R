#!/usr/bin/env R


mixcol_single <- function(c1, c2, ratio=0.5) {
  n <- names(c1)
  
  m1 <- col2rgb(c1)
  m2 <- col2rgb(c2)
  
  m <- (ratio * m2) + ((1 - ratio) * m1)
  r <- rgb(m[1,],m[2,],m[3,], maxColorValue = 255)
  names(r) <- n
  
  return(r)
}



mixcol <- function(c1, c2, ratio=c(0.5)) {
  out = c()
  
  for(rat in ratio) {
    out <- c(out, mixcol_single(c1, c2, rat))
  }
  
  return (out)
}




# center with gray - original also used in corrplot and some HI-publications
col3 <- grDevices::colorRampPalette(c(       "#67001F", 
                                             mixcol("#B2182B","black",0.025),
                                             mixcol("#D6604D","black",0.050),
                                             mixcol("#F4A582","black",0.075),
                                             mixcol("#FDDBC7","black",0.100),
                                             
                                             mixcol("#FFFFFF","black",0.125),
                                             
                                             mixcol("#D1E5F0","black",0.100),
                                             mixcol("#92C5DE","black",0.075),
                                             mixcol("#4393C3","black",0.050),
                                             mixcol("#2166AC","black",0.025),
                                             
                                             "#053061"
))




library(ggplot2)




parse_cells <- function(cell) {
  return (data.frame(
    id = cell$id,
    nuc_mean_x =  mean(t(data.frame(cell$nucleusGeometry$coordinates))[,1]),
    nuc_mean_y =  mean(t(data.frame(cell$nucleusGeometry$coordinates))[,2]),
    nuclear_dab_sd = cell$properties$measurements$`Nucleus: Hematoxylin OD std dev`,
    nuclear_dab_mean = cell$properties$measurements$`Nucleus: DAB OD mean`,
    nuclear_dab_max = cell$properties$measurements$`Nucleus: DAB OD max`,
    cellular_dab_sd = cell$properties$measurements$`Cell: Hematoxylin OD std dev`,
    cellular_dab_mean = cell$properties$measurements$`Cell: DAB OD mean`,
    cellular_dab_max = cell$properties$measurements$`Cell: DAB OD max`
    ))
}

dat <- rjson::fromJSON(file = "/home/youri/test2.json.geojson") # alt ctr C in qupath
dat2 <- dat$features
plt <- do.call(rbind, pbapply::pblapply(dat2, parse_cells))


plot(plt$nuclear_dab_mean, plt$nuclear_dab_sd)

plot(plt$nuclear_dab_mean, plt$cellular_dab_sd, pch=19,cex=0.01)
plot(plt$nuclear_dab_mean, plt$cellular_dab_mean, pch=19,cex=0.01)



ggplot(plt, aes(x=nuc_mean_x, y=-nuc_mean_y
                                                         #, col=nuclear_dab_mean
                                                         )) +
  geom_point(data=plt |> dplyr::filter(nuclear_dab_mean <= 0.5), col="darkgreen", alpha=0.25, cex=0.001) +
  geom_point(data=plt |> dplyr::filter(nuclear_dab_mean > 0.5), col="red", cex=0.45)
  #ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", oob = scales::squish)



plot(plt$nuclear_dab_mean, plt$nuclear_dab_sd, pch=19,cex=0.01)

