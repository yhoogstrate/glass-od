#!/usr/bin/env



metadata <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) 


# 1. plt. resection nr ----

plt <- metadata

ggplot(plt, aes(x=resection_number, y=array_A_IDH_HG__A_IDH_LG_lr__lasso_fit , group=patient_id)) +
  geom_line(lwd=theme_cellpress_lwd, lty=1, col="darkgray",alpha=0.5) +
  geom_point(size = theme_cellpress_size/ 3) +
  stat_summary(aes(group=resection_number),fun = mean, geom = "crossbar", width = 0.5, color = "red", fatten = 1) +
  labs(x="Surgical intervention #", y="CGC''") +
  theme_cellpress


ggsave("output/figures/vis_CGC_aluvial.pdf", width= 1.15, height=1.75)






metadata <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES) 


# 2. plt. resection nr ----

plt <- metadata

ggplot(plt, aes(x=time_between_resection_and_array / 365, y=array_PC1 , col=isolation_material)) +
  #geom_line(lwd=theme_cellpress_lwd, lty=1, col="darkgray",alpha=0.5) +
  geom_point(size = theme_cellpress_size/ 3) +
  #stat_summary(aes(group=resection_number),fun = mean, geom = "crossbar", width = 0.5, color = "red", fatten = 1) +
  #labs(x="Surgical intervention #", y="CGC''") +
  theme_cellpress


ggsave("output/figures/vis_Quality_years.pdf", width= 1.15, height=2.7 * 0.75)




