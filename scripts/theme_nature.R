#!/usr/bin/env R

# nature:
# "We prefer the use of a ‘standard’ font, preferably 12-point Times New Roman"
# "serif" could be the alternative: quartzFonts()

# requires to have (Microsoft) Sans Serif: - https://github.com/wch/extrafont/issues/32
#
# https://superuser.com/questions/1153990/anyone-know-how-to-install-Times New Roman-fonts-on-centos-7
# sudo yum install curl cabextract xorg-x11-font-utils fontconfig
# wget http://www.itzgeek.com/msttcore-fonts-2.0-3.noarch.rpm
# sudo rpm -Uvh msttcore-fonts-2.0-3.noarch.rpm
#
# 1. update remotes / Rttf2pt1:
# install.packages('remotes')
# remotes::install_version("Rttf2pt1", version = "1.3.8")
# library(Rttf2pt1)
#
# 2. install and update fonts: - https://fromthebottomoftheheap.net/2013/09/09/preparing-figures-for-plos-one-with-r/
# probably this step has to be conducted under sudo? (sudo R)
# install.packages("extrafont")
# library("extrafont")
# font_import(pattern='Times New Roman')
# extrafont::font_import(pattern="Times New Roman.ttf", prompt=FALSE)
#
# sudo chmod -R 777 /usr/lib64/R/library/extrafontdb/metrics
# font_import(pattern='', prompt=FALSE)
#
# 3. check if Times New Roman exists:
# extrafont::loadfonts() 
# sort(extrafont::fonts())

# pts to mm: 0.3527777778



# For color palette's, have a look here:
# https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html



# For alternative theme(s), have a look here:
# https://github.com/koundy/ggplot_theme_Publication/blob/master/ggplot_theme_Publication-2.R

# font size range: 5 - 7

theme_nature_lwd <- 0.5 / 2.14 # pt?
theme_nature_size <- 6 * (3.88 / 11) # 
#theme_nature_size_annotate <- 6 * (3.88 / 11) # scaling_factor for ggplot::annotate(size=) to font size 6


theme_nature <- theme_bw() +
  theme(
    text =          element_text(size = 6, family = "Times New Roman", face = "plain"),
    axis.text =     element_text(size = 6, family = "Times New Roman", face = "plain", color="black"), # , angle=90, vjust =0.5
    axis.title.x =  element_text(size = 6, family = "Times New Roman", face = "plain", color="black"), # , vjust = -0.2
    axis.title.y =  element_text(size = 6, family = "Times New Roman", face = "plain", color="black"),
    axis.line =     element_line(linewidth = theme_nature_lwd),
    axis.ticks =    element_line(linewidth = theme_nature_lwd),
    
    strip.text =    element_text(size = 6, family = "Times New Roman", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.text.x =  element_text(size = 6, family = "Times New Roman", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.text.y =  element_text(size = 6, family = "Times New Roman", face = "plain", margin=margin(1,1,1,1), color="black"),
    strip.background = element_blank(), # clean as possible
    
    legend.title =  element_text(size = 6, family = "Times New Roman", face = "plain", color="black"),
    legend.text =   element_text(size = 6, family = "Times New Roman", face = "plain", color="black"),
    legend.position = 'bottom',
    legend.margin   = margin(t=-2),
    legend.key.size = unit(0.2, 'lines'),
    legend.key = element_blank(), # this should remove the white squares around the legend items, but seems to fail sometimes, probably due to the key.size above?
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    
    plot.title =      element_text(size = 6, family = "Times New Roman", face = "plain", color="black"), # `title` covers both title and subtitle
    plot.subtitle =   element_text(size = 6, family = "Times New Roman", face = "italic", color="#777777"),
    plot.caption =    element_text(size = 6, family = "Times New Roman", face = "italic", color="#777777"),
    plot.background = element_blank(),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing =      unit(0.1, "lines"), # facet_grid margin
    panel.border =       element_blank(), # no sqaure, but two lines instead (axis.line)
    panel.background =   element_blank()
  )


# Facette: use facet_wrap with scales="free" rather than facet_grid, for nice borders
# facet_wrap(~facet, scales="free",ncol=5)


# Exporting: - US letter size
# ggsave(..., width = (8.5 * 0.95), height = .. (<11.2))
# ggsave(..., width = (8.5 * 0.95) * (1/2), height = .. (<11.2)) # for two figs on one page width
# ggsave(..., width = (8.5 * 0.95) * (1/3), height = .. (<11.2)) # for three figs on one page width

library(extrafont)
loadfonts()
choose_font("Times New Roman")



#' legacy stuff:
#theme_nature_with_facet <- theme_nature +
#  theme(
#    #axis.ticks = element_blank(),
#        panel.border = element_rect(linewidth=theme_nature_lwd, fill=NA)
#        )

