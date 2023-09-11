#!/usr/bin/env R


#!/usr/bin/env R

# requires to have Arial: - https://github.com/wch/extrafont/issues/32
#
# https://superuser.com/questions/1153990/anyone-know-how-to-install-arial-fonts-on-centos-7
# sudo yum install curl cabextract xorg-x11-font-utils fontconfig
# wget http://www.itzgeek.com/msttcore-fonts-2.0-3.noarch.rpm
# rpm -Uvh msttcore-fonts-2.0-3.noarch.rpm
#
# 1. update remotes / Rttf2pt1:
# install.packages('remotes')
# remotes::install_version("Rttf2pt1", version = "1.3.8")
#
# 2. install and update fonts: - https://fromthebottomoftheheap.net/2013/09/09/preparing-figures-for-plos-one-with-r/
# install.packages("extrafont")
# library("extrafont")
# font_import(')
# font_import(pattern='Arial')
#
# 3. check if Arial exists:
# extrafonts::loadfonts()
# sort(extrafont::fonts())

# pts to mm: 0.3527777778

theme_cellpress_lwd <- 0.5 / 2.14 # pt?
theme_cellpress_size <- 7 * (3.88 / 11) # ~2.47



theme_cellpress <- theme_bw() +
  theme(
  text =          element_text(size = 7, family = "Arial", face = "plain"),
  axis.text =     element_text(size = 7, family = "Arial", face = "plain"),
  axis.title.x =  element_text(size = 7, family = "Arial", face = "plain"),
  axis.title.y =  element_text(size = 7, family = "Arial", face = "plain"),
  axis.line =     element_line(linewidth = theme_cellpress_lwd),
  axis.ticks =    element_line(linewidth = theme_cellpress_lwd),
  
  strip.text =    element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1)),
  strip.text.x =  element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1)),
  strip.text.y =  element_text(size = 7, family = "Arial", face = "plain", margin=margin(1,1,1,1)),
  strip.background = element_blank(), # clean as possible

  legend.title =  element_text(size = 7, family = "Arial", face = "plain"),
  legend.text =   element_text(size = 7, family = "Arial", face = "plain"),
  
  legend.position = 'bottom',
  legend.margin   = margin(t=-2),
  legend.key.size = unit(0.2, 'lines'),
  legend.background = element_blank(),
  legend.box.background = element_blank(),
  
  plot.title =    element_text(size = 7, family = "Arial", face = "plain"), # `title` covers both title and subtitle
  plot.subtitle = element_text(size = 7, family = "Arial", face = "italic"),
  plot.caption =  element_text(size = 7, family = "Arial", face = "italic"),
  plot.background = element_blank(),

  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.spacing =      unit(0.1, "lines"), # facet_grid margin
  panel.border =       element_blank(), # no sqaure, but two lines instead (axis.line)
  panel.background = element_blank()
  
  

)


theme_cellpress_with_facet <- theme_cellpress +
  theme(
    #axis.ticks = element_blank(),
        panel.border = element_rect(linewidth=theme_cellpress_lwd, fill=NA)
        )

extrafont::loadfonts()


