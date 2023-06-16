#!/usr/bin/env R


plt <- glass_od.metadata.idats |> 
  dplyr::filter(is.na(reason_excluded_array_sample)) |> 
  dplyr::filter(is.na(reason_excluded_patient)) |> 
  dplyr::filter(is.na(reason_excluded_resection_isolation)) |> 
  dplyr::filter(is.na(reason_excluded_resection))



plt.val <- plt |> 
  dplyr::mutate(y = mgmt_Estimated) |> 
  dplyr::mutate(type = "MGMT value")


plt.line <- rbind(
  plt |> 
    dplyr::mutate(y = mgmt_CI_Lower),
  plt |> 
    dplyr::mutate(y = mgmt_CI_Upper)
) |> 
  dplyr::mutate(type = "CI")


order <- plt |> 
  dplyr::select(resection_id, mgmt_Estimated) |> 
  dplyr::mutate(order = order(order(mgmt_Estimated, resection_id))) |> 
  dplyr::mutate(mgmt_Estimated = NULL)

plt.c <- rbind(plt.val, plt.line) |> 
  dplyr::left_join(order,by=c('resection_id'='resection_id')) |> 
  dplyr::arrange(order) |> 
  dplyr::mutate(mgmt_status = ifelse(is.na(mgmt_status), "unconfident", mgmt_status)) |> 
  dplyr::mutate(is.primary = ifelse(gsub("^.+\\-","",resection_id) == "R1", "primary", "recurrent"))


ggplot(plt.c, aes(x= reorder(resection_id, order), y=y, col=mgmt_status)) +
  
  facet_grid(cols = vars(is.primary), scales = "free", space="free") +
  
  geom_line(data = plt.c |> dplyr::filter(type == "CI"), lwd=1, col="gray80") +
  geom_point(data = plt.c |> dplyr::filter(type == "MGMT value")) +
  geom_hline(yintercept = plt.c$mgmt_Cutoff[1], lty=2, lwd=0.5) +
  labs(y = "predictMGMT estimated value", x="tumor resection", col = "predictMGMT status") +
  theme_light() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'bottom') +
  scale_color_manual(values=c('methylated' = hue_pal()(3)[2],
                              'unmethylated' = hue_pal()(3)[1],
                              'unconfident' = 'gray40'
                              ))

ggsave("output/figures/plt_mgmt_status.pdf", width=8.1, height=6)
