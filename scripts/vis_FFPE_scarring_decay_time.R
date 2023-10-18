#!/usr/bin/env R

# load ----


source('scripts/load_palette.R')
source('scripts/load_functions.R')
source('scripts/load_themes.R')


library(ggplot2)


if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}



# plot x FFPE time ----


plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(210)


ggplot(plt, aes(x=-time_between_resection_and_array, array_epiTOC2_hypoSC, col =isolation_material)) +
  scale_y_reverse() +
  geom_point() + 
  theme_cellpress


# plot x patient age ----





ggplot(plt, aes(x=time_between_birth_and_resection, array_epiTOC2_hypoSC, col =isolation_material, label=isolation_id)) +
  scale_y_reverse() +
  geom_point() + 
  ggrepel::geom_text_repel(size=2.5, col="black") +
  theme_cellpress


plot(pp$age_at_resection, pp$time_between_birth_and_resection)
cor(pp$age_at_resection, as.numeric(pp$time_between_birth_and_resection))



