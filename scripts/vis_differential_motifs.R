#!/usr/bin/env R


# load ----


library(ggplot2)
library(patchwork)



source('scripts/load_constants.R')
source('scripts/load_functions.R')
source('scripts/load_palette.R')
source('scripts/load_themes.R')



if(!exists('data.mvalues.probes')) {
  source('scripts/load_mvalues_hq_samples.R')
}


#if(!exists('data.intensities.probes')) {
#  source('scripts/load_intensities_hq_samples.R')
#}




## plt.motifs.stranded ----


plt.motifs.stranded <- data.frame(# R way of thinking
  b1 = rep(c(rep("A",1),     rep("C",1),     rep("T",1),      rep("G",1)), 256/4),
  b2 = rep(c(rep("A",4),     rep("C",4),     rep("T",4),      rep("G",4)), 256/(4*4)),
  
  b3 = rep("C", 4*4),
  b4 = rep("G", 4*4),
  
  b5 = rep(c(rep("A",4*4),   rep("C",4*4),   rep("T", 4*4),   rep("G",4*4)), 256/(4*4*4)),
  b6 = rep(c(rep("A",4*4*4), rep("C",4*4*4), rep("T", 4*4*4), rep("G",4*4*4)), 256/(4*4*4*4))
) |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == (256))
    return(.)
  })() |> 
  dplyr::tibble() |> 
  dplyr::mutate(sequence_5p = paste0(b1, b2, b3, b4, b5, b6)) |> 
  dplyr::mutate(sequence_3p = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(sequence_5p)))) |> 
  assertr::verify(!duplicated(sequence_5p)) |> 
  assertr::verify(!duplicated(sequence_3p)) |> 
  dplyr::select(sequence_5p, sequence_3p) |> 
  dplyr::mutate(palindromic = sequence_5p == sequence_3p) |> 
  dplyr::mutate(
    oligo_sequence = tolower(dplyr::case_when(
      sequence_5p < sequence_3p ~ paste0(sequence_5p,"/",sequence_3p),
      sequence_5p > sequence_3p ~ paste0(sequence_3p,"/",sequence_5p),
      palindromic ~ paste0("*",sequence_5p,"/",sequence_3p)
    )))




### adds: DNMT1 - Adam NatCom 2020 ----
#' https://www.nature.com/articles/s41467-020-17531-8


tmp <- readxl::read_xlsx('data/DNMT1_Adam_et_al_2020_N2_data.xlsx') |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == 256)
    return(.)
  })() |> 
  dplyr::rename(sequence = D1) |> 
  dplyr::mutate(sequence = gsub(" ","", sequence)) |> 
  dplyr::rename(DNMT1_Adam_NatCom_2020_stranded = av) |>
  assertr::verify(!duplicated(sequence)) |> 
  assertr::verify(is.numeric(DNMT1_Adam_NatCom_2020_stranded))


plt.motifs.stranded <- plt.motifs.stranded |> 
  dplyr::left_join(tmp, by=c('sequence_5p'='sequence'), suffix=c('','')) |> 
  assertr::verify(!is.na(DNMT1_Adam_NatCom_2020_stranded))


rm(tmp)



### adds: DNMT1 - Adam NAR 2023 ----
#' https://academic.oup.com/nar/article-abstract/51/13/6622/7184160
#' "integrated with more data, but very similar"


tmp <- readxl::read_xlsx('data/DNMT1_Adam_NAR_2023_HM_combi_N2.xlsx') |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == 256)
    return(.)
  })() |> 
  dplyr::rename(sequence = Site) |> 
  dplyr::mutate(sequence = gsub(" ","", sequence)) |> 
  dplyr::rename(DNMT1_Adam_NAR_2023_stranded = HM) |>
  assertr::verify(!duplicated(sequence)) |> 
  assertr::verify(is.numeric(DNMT1_Adam_NAR_2023_stranded)) 


plt.motifs.stranded <- plt.motifs.stranded |> 
  dplyr::left_join(tmp, by=c('sequence_5p'='sequence'), suffix=c('','')) |> 
  assertr::verify(!is.na(DNMT1_Adam_NAR_2023_stranded))


rm(tmp)



### adds: DNMT3 - Gao NatCom 2020 ----
#' https://www.nature.com/articles/s41467-020-17109-4



tmp <- readxl::read_xlsx('data/41467_2020_17109_MOESM4_ESM.xlsx', skip=3) |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == (4096))
    return(.)
  })() |> 
  dplyr::rename(sequence_full = 1) |> 
  dplyr::rename(sequence_full_2 = 4) |> 
  dplyr::rename(sequence_full_3 = 7) |> 
  assertr::verify(sequence_full == sequence_full_2) |> 
  assertr::verify(sequence_full == sequence_full_3) |> 
  dplyr::rename(DNMT3A_Gao_NatCom_2020_stranded_3f = 2) |> 
  dplyr::rename(DNMT3B_Gao_NatCom_2020_stranded_3f = 5) |> 
  dplyr::mutate(sequence = gsub("^.(..) CG (..).$","\\1CG\\2", sequence_full)) |> 
  dplyr::select(sequence, DNMT3A_Gao_NatCom_2020_stranded_3f, DNMT3B_Gao_NatCom_2020_stranded_3f) |> 
  dplyr::group_by(sequence) |> 
  dplyr::summarise(DNMT3A_Gao_NatCom_2020_stranded = mean(DNMT3A_Gao_NatCom_2020_stranded_3f),
                   DNMT3B_Gao_NatCom_2020_stranded = mean(DNMT3B_Gao_NatCom_2020_stranded_3f)) |> 
  dplyr::ungroup() |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == (256))
    return(.)
  })() 


plt.motifs.stranded <- plt.motifs.stranded |> 
  dplyr::left_join(tmp, by=c('sequence_5p'='sequence'), suffix=c('','')) |> 
  assertr::verify(!is.na(DNMT3A_Gao_NatCom_2020_stranded)) |> 
  assertr::verify(!is.na(DNMT3B_Gao_NatCom_2020_stranded))


rm(tmp)



### adds: TET1/2 - Adam NatCom 2022 ----
#' https://www.nature.com/articles/s42003-022-03033-4



tmp <- readxl::read_xlsx('data/42003_2022_3033_MOESM4_ESM.xlsx', skip=3) |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == 256)
    return(.)
  })() |> 
  dplyr::rename(sequence = 1) |> 
  dplyr::mutate(sequence = gsub(" ", "" ,sequence)) |> 
  assertr::verify(!duplicated(sequence)) |> 
  dplyr::rename(TET1_5mC_Adam_NatCom_2022_stranded  = 2) |> 
  dplyr::rename(TET2_5mC_Adam_NatCom_2022_stranded  = 4) |> 
  dplyr::rename(TET1_5hmC_Adam_NatCom_2022_stranded = 3) |> 
  dplyr::rename(TET2_5hmC_Adam_NatCom_2022_stranded = 5) |> 
  assertr::verify(is.numeric(TET1_5mC_Adam_NatCom_2022_stranded)) |> 
  assertr::verify(is.numeric(TET2_5mC_Adam_NatCom_2022_stranded)) |> 
  assertr::verify(is.numeric(TET1_5hmC_Adam_NatCom_2022_stranded)) |> 
  assertr::verify(is.numeric(TET2_5hmC_Adam_NatCom_2022_stranded))


plt.motifs.stranded <- plt.motifs.stranded |> 
  dplyr::left_join(tmp, by=c('sequence_5p'='sequence'), suffix=c('','')) |> 
  assertr::verify(!is.na(TET1_5mC_Adam_NatCom_2022_stranded)) |> 
  assertr::verify(!is.na(TET2_5mC_Adam_NatCom_2022_stranded)) |> 
  assertr::verify(!is.na(TET1_5hmC_Adam_NatCom_2022_stranded)) |> 
  assertr::verify(!is.na(TET2_5hmC_Adam_NatCom_2022_stranded))


rm(tmp)




### adds: TET3 - Ravichandran SciAdv 2022 ----
#' https://www.science.org/doi/10.1126/sciadv.abm2427



tmp <- read.table("data/10.1126_sciadv.abm2427_Slope_summary.txt") |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == 256)
    return(.)
  })() |> 
  dplyr::rename(sequence = Motif) |> 
  tibble::rownames_to_column('sequence2') |> 
  assertr::verify(sequence == sequence) |> 
  dplyr::mutate(sequence2 = NULL) |> 
  dplyr::rename(TET3_Ravichandran_SciAdv_2022_stranded = ALL_84) |> 
  dplyr::tibble() |> 
  assertr::verify(!duplicated(sequence)) |> 
  assertr::verify(is.numeric(TET3_Ravichandran_SciAdv_2022_stranded)) |> 
  dplyr::select(sequence, TET3_Ravichandran_SciAdv_2022_stranded)


plt.motifs.stranded <- plt.motifs.stranded |> 
  dplyr::left_join(tmp, by=c('sequence_5p'='sequence'), suffix=c('','')) |> 
  assertr::verify(!is.na(TET3_Ravichandran_SciAdv_2022_stranded))


rm(tmp)




### adds: per probe outcomes (DMP etc.) ----


tmp <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |>
  dplyr::tibble() |> 
  dplyr::filter(!is.na(gc_sequence_context_2_new)) |> 
  dplyr::mutate(sequence_5p = toupper(gsub("[^a-zA-Z]","", gc_sequence_context_2_new))) |> 
  dplyr::select(sequence_5p, DMP__g2_g3__pp_nc__t, DMP__primary_recurrence__pp_nc__t) |> 
  dplyr::group_by(sequence_5p) |>
  dplyr::summarise(GLASS_OD__DMP__g2_g3__mean = mean(DMP__g2_g3__pp_nc__t),
                   GLASS_OD__DMP__g2_g3__median = median(DMP__g2_g3__pp_nc__t),
                   GLASS_OD__DMP__primary_recurrence__mean = mean(DMP__primary_recurrence__pp_nc__t),
                   GLASS_OD__DMP__primary_recurrence__median = median(DMP__primary_recurrence__pp_nc__t)
  ) |> 
  dplyr::ungroup()  |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == 256)
    return(.)
  })() |> 
  assertr::verify(sequence_5p %in% plt.motifs.stranded$sequence_5p)



plt.motifs.stranded <- plt.motifs.stranded |> 
  dplyr::left_join(tmp, by=c('sequence_5p'='sequence_5p')) |> 
  assertr::verify(!is.na(GLASS_OD__DMP__g2_g3__mean)) |> 
  assertr::verify(!is.na(GLASS_OD__DMP__g2_g3__median)) |> 
  assertr::verify(!is.na(GLASS_OD__DMP__primary_recurrence__mean)) |> 
  assertr::verify(!is.na(GLASS_OD__DMP__primary_recurrence__median))


rm(tmp)




## plt.motifs.unstranded [palindrom correction] ----



plt.motifs.unstranded <- plt.motifs.stranded |> 
  dplyr::select(
    `oligo_sequence`, `DNMT1_Adam_NatCom_2020_stranded`, `DNMT1_Adam_NAR_2023_stranded`,
    `DNMT3A_Gao_NatCom_2020_stranded`, `DNMT3B_Gao_NatCom_2020_stranded`, `TET1_5mC_Adam_NatCom_2022_stranded`,
    `TET1_5hmC_Adam_NatCom_2022_stranded`, `TET2_5mC_Adam_NatCom_2022_stranded`, `TET2_5hmC_Adam_NatCom_2022_stranded`,
    `TET3_Ravichandran_SciAdv_2022_stranded`
  ) |> 
  dplyr::group_by(oligo_sequence) |> 
  dplyr::summarise_all(mean) |> # all, except the one used to group_by
  dplyr::rename_with(~ gsub("_stranded$","_unstranded",.x))  |> 
  dplyr::filter(!duplicated(oligo_sequence)) |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == 136)
    return(.)
  })()




### adds: per probe outcomes (DMP etc.) ----


tmp <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |>
  dplyr::tibble() |> 
  assertr::verify(!is.na(gc_sequence_context_2_new)) |> 
  dplyr::select(gc_sequence_context_2_new, DMP__g2_g3__pp_nc__t, DMP__primary_recurrence__pp_nc__t) |> 
  dplyr::mutate(sequence_5p = gsub("[CG]","CG",gc_sequence_context_2_new, fixed=T), gc_sequence_context_2_new=NULL) |> 
  dplyr::mutate(sequence_3p = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(sequence_5p)))) |> 
  dplyr::mutate(palindromic = sequence_5p == sequence_3p) |> 
  dplyr::mutate(
    oligo_sequence = tolower(dplyr::case_when(
      sequence_5p < sequence_3p ~ paste0(sequence_5p,"/",sequence_3p),
      sequence_5p > sequence_3p ~ paste0(sequence_3p,"/",sequence_5p),
      palindromic ~ paste0("*",sequence_5p,"/",sequence_3p)
    ))) |> 
  dplyr::mutate(
    oligo_sequence = tolower(dplyr::case_when(
      sequence_5p < sequence_3p ~ paste0(sequence_5p,"/",sequence_3p),
      sequence_5p > sequence_3p ~ paste0(sequence_3p,"/",sequence_5p),
      palindromic ~ paste0("*",sequence_5p,"/",sequence_3p)
    ))) |> 
  
  dplyr::group_by(oligo_sequence) |>
  dplyr::summarise(GLASS_OD__DMP__g2_g3__mean = mean(DMP__g2_g3__pp_nc__t),
                   GLASS_OD__DMP__g2_g3__median = median(DMP__g2_g3__pp_nc__t),
                   GLASS_OD__DMP__primary_recurrence__mean = mean(DMP__primary_recurrence__pp_nc__t),
                   GLASS_OD__DMP__primary_recurrence__median = median(DMP__primary_recurrence__pp_nc__t)
                   ) |> 
  dplyr::ungroup()  |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == 136)
    return(.)
  })() |> 
  assertr::verify(oligo_sequence %in% plt.motifs.stranded$oligo_sequence)



plt.motifs.unstranded <- plt.motifs.unstranded |> 
  dplyr::left_join(tmp, by=c('oligo_sequence'='oligo_sequence')) |> 
  assertr::verify(!is.na(GLASS_OD__DMP__g2_g3__mean)) |> 
  assertr::verify(!is.na(GLASS_OD__DMP__g2_g3__median)) |> 
  assertr::verify(!is.na(GLASS_OD__DMP__primary_recurrence__mean)) |> 
  assertr::verify(!is.na(GLASS_OD__DMP__primary_recurrence__median))


rm(tmp)




# plots ----



pdf("output/figures/vis_differential__motif_tscores_x_enzyme_kinetics__corrplot.pdf", width=4,height=4)

plt <- plt.motifs |> 
  dplyr::select(-c(sequence_5p, sequence_3p, palindromic,
                   GLASS_OD__DMP__g2_g3__median, 
                   GLASS_OD__DMP__primary_recurrence__median, GLASS_OD__DMP__primary_recurrence__mean,
                   TET1_5hmC_Adam_NatCom_2022, TET2_5hmC_Adam_NatCom_2022
                   )) |> 
  dplyr::mutate(`-1 * TET3_Ravichandran_SciAdv_2022` = -1 * TET3_Ravichandran_SciAdv_2022, TET3_Ravichandran_SciAdv_2022=NULL) |> 
  dplyr::mutate(`-1 * GLASS_OD__DMP__g2_g3__mean` = -1 * GLASS_OD__DMP__g2_g3__mean, GLASS_OD__DMP__g2_g3__mean = NULL) |> 
  tibble::column_to_rownames('oligo_sequence') |> 
  cor(method="pearson")

colnames(plt) <- gsub("_"," ",colnames(plt))
rownames(plt) <- gsub("_"," ",rownames(plt))

corrplot::corrplot(plt, order="hclust", tl.cex=0.65 * (7/8),
                   cl.cex=0.65 * (7/8),
                   bg=NA,  tl.col = "black")


dev.off()



## motif: xx[CG]xx violin(s) stranded ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |>
  #dplyr::mutate(sequence_5p = toupper(paste0(gsub("^.+(..)$","\\1", sequence_pre), "CG", gsub("^(..).+$","\\1", sequence_post)))) |> 
  dplyr::mutate(sequence_5p = gsub("[^a-zA-Z]","",gc_sequence_context_2_new)) |> 
  dplyr::filter(!is.na(sequence_5p)) |> 
  dplyr::left_join(plt.motifs.stranded, by=c('sequence_5p'='sequence_5p'),suffix=c('',''))




plt <- plt |> 
  dplyr::group_by(sequence_5p) |> 
  dplyr::mutate(facet_rank = median(DMP__g2_g3__pp_nc__t)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(facet_name = paste0(
                    tolower(gsub("^(..)....$","\\1", sequence_5p)),
                    toupper(gsub("^..(..).+$","\\1", sequence_5p)),
                    tolower(gsub("^....(..)$","\\1", sequence_5p))
                  ))



ggplot(plt, aes(x=reorder(facet_name, -facet_rank), y=DMP__g2_g3__pp_nc__t)) +
  #facet_wrap(~reorder(facet_name, -facet_rank), scales="free_x", ncol=length(unique(plt$facet_name))) +
  #ggbeeswarm::geom_quasirandom(size=0.01, alpha=0.65) +
  #ggplot2::geom_violin(draw_quantiles = c(), linewidth=theme_cellpress_lwd, col = "white", fill="darkgray", adjust = 1.95) +
  geom_boxplot(width=0.75, outlier.shape=NA, outlier.color=NA, col="darkgray",
               coef=0.5, fill=NA,linewidth=theme_cellpress_lwd) +
  
  stat_summary(fun.y = median, fun.min = median, fun.max = median,
               geom = "crossbar", col="red", width=0.85, linewidth=theme_cellpress_lwd) +
  
  labs(x = NULL, y="Per probe t-score Grade 2 ~ Grade 3") +
  coord_cartesian(ylim = c(-6.75, 3.5)) + # soft clip
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, family = "mono")) +
  labs(subtitle=format_subtitle("Differential methylated motif overview"))


ggsave("output/figures/vis_differential__xxCGxx_violin.pdf", width = 11 * 0.975, height=3.75)





## motif: xx[CG]xx violin(s) unstranded ----



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |>
  #dplyr::mutate(sequence_context_stranded = toupper(paste0(gsub("^.+(..)$","\\1", sequence_pre), "CG", gsub("^(..).+$","\\1", sequence_post)))) |> 
  #dplyr::filter(!is.na(sequence_context_stranded)) |> 
  dplyr::mutate(sequence_5p = gsub("[^a-zA-Z]","",gc_sequence_context_2_new)) |> 
  dplyr::filter(!is.na(sequence_5p)) |> 
  dplyr::left_join(plt.motifs.stranded |> dplyr::select(sequence_5p, oligo_sequence), by=c('sequence_5p'='sequence_5p'), suffix=c('','')) |> 
  dplyr::left_join(plt.motifs.unstranded, by=c('oligo_sequence'='oligo_sequence'),suffix=c('',''))




plt <- plt |> 
  dplyr::group_by(oligo_sequence) |> 
  dplyr::mutate(facet_rank = median(DMP__g2_g3__pp_nc__t)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(facet_name = tolower(oligo_sequence))



ggplot(plt, aes(x=reorder(facet_name, -facet_rank), y=DMP__g2_g3__pp_nc__t)) +
  geom_boxplot(width=0.75, outlier.shape=NA, outlier.color=NA, col="darkgray",
               coef=0.5, fill=NA,linewidth=theme_nature_lwd) +
  
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", col="red", width=0.85, linewidth=theme_nature_lwd) +
  
  labs(x = NULL,
       y="Per probe t-score Grade 2 ~ Grade 3",
       caption=paste0("n=",length(unique(plt$facet_name))," sequence motifs, corrected palindromic motifs for strand (*)")) +
  coord_cartesian(ylim = c(-6.75, 3.5)) + # soft clip
  theme_nature + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, family = "mono")) +
  labs(subtitle=format_subtitle("Differential methylated motif overview"))



ggsave("output/figures/vis_differential_motifs__t_stat_grade_per_motif__unstranded.pdf", width=11 * 0.975, height=3.8)





## motif: xx[CG]xx violin(s) stranded x TET1----



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |>
  dplyr::mutate(sequence_context_stranded = toupper(paste0(gsub("^.+(..)$","\\1", sequence_pre), "CG", gsub("^(..).+$","\\1", sequence_post)))) |> 
  dplyr::filter(!is.na(sequence_context_stranded)) |> 
  dplyr::left_join(plt.motifs.stranded, by=c('sequence_context_stranded'='sequence_5p'),suffix=c('','')) |> 
  assertr::verify(is.numeric(TET1_5mC_Adam_NatCom_2022_stranded)) |> 
  dplyr::select(sequence_context_stranded, DMP__g2_g3__pp_nc__t, TET1_5mC_Adam_NatCom_2022_stranded) |> 
  dplyr::group_by(sequence_context_stranded) |> 
  dplyr::summarise(DMP__g2_g3__pp_nc__t = mean(DMP__g2_g3__pp_nc__t),
                   TET1_5mC_Adam_NatCom_2022_stranded = mean(TET1_5mC_Adam_NatCom_2022_stranded)) |> 
  dplyr::ungroup() |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 256) 
    return(.)
  })()




ggplot(plt, aes(y=TET1_5mC_Adam_NatCom_2022_stranded, x=DMP__g2_g3__pp_nc__t)) +
  geom_point() +
  #geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  labs(caption=paste0("n=",length(unique(plt$sequence_context_stranded))," motifs")) +
  theme_cellpress +
  theme(plot.background = element_rect(fill="white")) # png export




## motif: xx[CG]xx violin(s) stranded x TET2----



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |>
  dplyr::mutate(sequence_context_stranded = toupper(paste0(gsub("^.+(..)$","\\1", sequence_pre), "CG", gsub("^(..).+$","\\1", sequence_post)))) |> 
  dplyr::filter(!is.na(sequence_context_stranded)) |> 
  dplyr::left_join(plt.motifs.stranded, by=c('sequence_context_stranded'='sequence_5p'),suffix=c('','')) |> 
  assertr::verify(is.numeric(TET1_5mC_Adam_NatCom_2022_stranded)) |> 
  dplyr::select(sequence_context_stranded, DMP__g2_g3__pp_nc__t, TET2_5mC_Adam_NatCom_2022_stranded, palindromic) |> 
  dplyr::group_by(sequence_context_stranded) |> 
  dplyr::summarise(DMP__g2_g3__pp_nc__t = mean(DMP__g2_g3__pp_nc__t),
                   TET2_5mC_Adam_NatCom_2022_stranded = mean(TET2_5mC_Adam_NatCom_2022_stranded),
                   palindromic = all(palindromic)) |> 
  dplyr::ungroup() |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 256) 
    return(.)
  })()




ggplot(plt, aes(y=TET2_5mC_Adam_NatCom_2022_stranded, x=DMP__g2_g3__pp_nc__t, col=palindromic)) +
  geom_point() +
  #geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export




## motif: xx[CG]xx violin(s) stranded x TET3----



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |>
  dplyr::mutate(sequence_context_stranded = toupper(paste0(gsub("^.+(..)$","\\1", sequence_pre), "CG", gsub("^(..).+$","\\1", sequence_post)))) |> 
  dplyr::filter(!is.na(sequence_context_stranded)) |> 
  dplyr::left_join(plt.motifs.stranded, by=c('sequence_context_stranded'='sequence_5p'),suffix=c('','')) |> 
  assertr::verify(is.numeric(TET3_Ravichandran_SciAdv_2022_stranded)) |> 
  dplyr::group_by(sequence_context_stranded) |> 
  dplyr::summarise(DMP__g2_g3__pp_nc__t = mean(DMP__g2_g3__pp_nc__t),
                   TET3_Ravichandran_SciAdv_2022_stranded = mean(TET3_Ravichandran_SciAdv_2022_stranded)) |> 
  dplyr::ungroup() |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 256) 
    return(.)
  })()




ggplot(plt, aes(y=TET3_Ravichandran_SciAdv_2022_stranded, x=DMP__g2_g3__pp_nc__t)) +
  geom_point() +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress +
  theme(plot.background = element_rect(fill="white")) # png export



## motif: xx[CG]xx violin(s) unstranded x DNMT1 ----


plt <- plt.motifs.unstranded |> 
  dplyr::mutate(adjacent_cg = ifelse(grepl("cgcg", oligo_sequence),"Adjacent CG" ,"No adjacent CG"))


# DMP__g2_g3__pp_nc__t
# DMP__g2_g3__pp_nc_PC1__t
# DMP__PCs__pp_nc__PC1_t
ggplot(plt, aes(y=DNMT1_Adam_NAR_2023_unstranded, x=GLASS_OD__DMP__g2_g3__median, col=adjacent_cg, label=oligo_sequence)) +
  geom_point(size=theme_nature_size/3) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho", size=theme_nature_size,
                   label.x=-0.85) +
  scale_color_manual(
    values=c(
      'Adjacent CG' = mixcol( 'lightblue', 'lightgreen'),
      'No adjacent CG' = mixcol( 'darkblue', 'darkgreen')
    )
  ) +
  theme_nature +
  #ggrepel::geom_text_repel(data = subset(plt, grepl("cgcgcg|aacgtt|aacgtc|aacgtg", oligo_sequence)), col="black", size=theme_nature_size) +
  #coord_cartesian(ylim=c(0, NA))  +
  labs(x="T-score\nWHO grade 2 <> WHO grade 3",
       y="DNMT1 (Adam NAR 2023)",
       col="",
       caption=paste0("n=",length(unique(plt$oligo_sequence))," motifs (unstranded)"),
       subtitle=format_subtitle("TET1 x t-score WHO grade")
  )


ggsave("output/figures/vis_differential_motifs__DNMT1_x_grade__unstranded.pdf", width=8.5 * 0.975 * (1/5), height=2.18)



## motif: xx[CG]xx violin(s) unstranded x DNMT3A ----


plt <- plt.motifs.unstranded |> 
  dplyr::mutate(adjacent_cg = ifelse(grepl("cgcg", oligo_sequence),"Adjacent CG" ,"No adjacent CG"))


# DMP__g2_g3__pp_nc__t
# DMP__g2_g3__pp_nc_PC1__t
# DMP__PCs__pp_nc__PC1_t
ggplot(plt, aes(y=DNMT3A_Gao_NatCom_2020_unstranded, x=GLASS_OD__DMP__g2_g3__median, col=adjacent_cg, label=oligo_sequence)) +
  geom_point(size=theme_nature_size/3) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho", size=theme_nature_size,
                   label.x=-0.85) +
  scale_color_manual(
    values=c(
      'Adjacent CG' = mixcol( 'lightblue', 'lightgreen'),
      'No adjacent CG' = mixcol( 'darkblue', 'darkgreen')
    )
  ) +
  theme_nature +
  #ggrepel::geom_text_repel(data = subset(plt, grepl("cgcgcg|aacgtt|aacgtc|aacgtg", oligo_sequence)), col="black", size=theme_nature_size) +
  #coord_cartesian(ylim=c(0, NA))  +
  labs(x="T-score\nWHO grade 2 <> WHO grade 3",
       y="DNMT3A (Gao NatCom 2020)",
       col="",
       caption=paste0("n=",length(unique(plt$oligo_sequence))," motifs (unstranded)"),
       subtitle=format_subtitle("TET1 x t-score WHO grade")
  )


ggsave("output/figures/vis_differential_motifs__DNMT3A_x_grade__unstranded.pdf", width=8.5 * 0.975 * (1/5), height=2.18)



## motif: xx[CG]xx violin(s) unstranded x DNMT3B ----


plt <- plt.motifs.unstranded |> 
  dplyr::mutate(adjacent_cg = ifelse(grepl("cgcg", oligo_sequence),"Adjacent CG" ,"No adjacent CG"))


# DMP__g2_g3__pp_nc__t
# DMP__g2_g3__pp_nc_PC1__t
# DMP__PCs__pp_nc__PC1_t
ggplot(plt, aes(y=DNMT3B_Gao_NatCom_2020_unstranded, x=GLASS_OD__DMP__g2_g3__median, col=adjacent_cg, label=oligo_sequence)) +
  geom_point(size=theme_nature_size/3) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho", size=theme_nature_size,
                   label.x=-0.85) +
  scale_color_manual(
    values=c(
      'Adjacent CG' = mixcol( 'lightblue', 'lightgreen'),
      'No adjacent CG' = mixcol( 'darkblue', 'darkgreen')
    )
  ) +
  theme_nature +
  #ggrepel::geom_text_repel(data = subset(plt, grepl("cgcgcg|aacgtt|aacgtc|aacgtg", oligo_sequence)), col="black", size=theme_nature_size) +
  #coord_cartesian(ylim=c(0, NA))  +
  labs(x="T-score\nWHO grade 2 <> WHO grade 3",
       y="DNMT3B (Gao NatCom 2020)",
       col="",
       caption=paste0("n=",length(unique(plt$oligo_sequence))," motifs (unstranded)"),
       subtitle=format_subtitle("TET1 x t-score WHO grade")
  )


ggsave("output/figures/vis_differential_motifs__DNMT3B_x_grade__unstranded.pdf", width=8.5 * 0.975 * (1/5), height=2.18)





## motif: xx[CG]xx violin(s) unstranded x TET1----



plt <- plt.motifs.unstranded |> 
  dplyr::mutate(adjacent_cg = ifelse(grepl("cgcg", oligo_sequence),"Adjacent CG" ,"No adjacent CG"))


# DMP__g2_g3__pp_nc__t
# DMP__g2_g3__pp_nc_PC1__t
# DMP__PCs__pp_nc__PC1_t
ggplot(plt, aes(y=TET1_5mC_Adam_NatCom_2022_unstranded, x=GLASS_OD__DMP__g2_g3__median, col=adjacent_cg, label=oligo_sequence)) +
  geom_point(size=theme_nature_size/3) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho", size=theme_nature_size,
                   label.x=-0.85) +
  scale_color_manual(
    values=c(
      'Adjacent CG' = mixcol( 'lightblue', 'lightgreen'),
      'No adjacent CG' = mixcol( 'darkblue', 'darkgreen')
    )
  ) +
  theme_nature +
  #ggrepel::geom_text_repel(data = subset(plt, grepl("cgcgcg|aacgtt|aacgtc|aacgtg", oligo_sequence)), col="black", size=theme_nature_size) +
  coord_cartesian(ylim=c(0, NA))  +
  labs(x="T-score\nWHO grade 2 <> WHO grade 3",
       y="TET1 5mC (Adam NatCom 2022)",
       col="",
       caption=paste0("n=",length(unique(plt$oligo_sequence))," motifs (unstranded)"),
       subtitle=format_subtitle("TET1 x t-score WHO grade")
       )


ggsave("output/figures/vis_differential_motifs__TET1_x_grade__unstranded.pdf", width=8.5 * 0.975 * (1/5), height=2.18)




## motif: xx[CG]xx violin(s) unstranded x TET2----


plt <- plt.motifs.unstranded |> 
  dplyr::mutate(adjacent_cg = ifelse(grepl("cgcg", oligo_sequence),"Adjacent CG" ,"No adjacent CG"))




# DMP__g2_g3__pp_nc__t
# DMP__g2_g3__pp_nc_PC1__t
# DMP__PCs__pp_nc__PC2_t = mean(DMP__PCs__pp_nc__PC2_t),
ggplot(plt, aes(y=TET2_5mC_Adam_NatCom_2022_unstranded,
                x=GLASS_OD__DMP__g2_g3__median,
                col=adjacent_cg, label=oligo_sequence)) +
  geom_point(size=theme_nature_size/3) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho", size=theme_nature_size,
                   label.x=-0.85) +
  scale_color_manual(
    values=c(
      'Adjacent CG' = mixcol( 'lightblue', 'lightgreen'),
      'No adjacent CG' = mixcol( 'darkblue', 'darkgreen')
    )
  ) +
  theme_nature +
  #ggrepel::geom_text_repel(data = subset(plt, grepl("cgcgcg|aacgtt|aacgtc|aacgtg", oligo_sequence)), col="black", size=theme_nature_size) +
  coord_cartesian(ylim=c(0, NA))  +
  labs(x="T-score\nWHO grade 2 <> WHO grade 3",
       y="TET2 5mC (Adam NatCom 2022)",
       col="",
       caption=paste0("n=",length(unique(plt$oligo_sequence))," motifs (unstranded)"),
       subtitle=format_subtitle("TET2 x t-score WHO grade")
  )


ggsave("output/figures/vis_differential_motifs__TET2_x_grade__unstranded.pdf", width=8.5 * 0.975 * (1/5), height=2.18)





## motif: xx[CG]xx violin(s) unstranded x TET3----



plt <- plt.motifs.unstranded |> 
  dplyr::mutate(adjacent_cg = ifelse(grepl("cgcg", oligo_sequence),"Adjacent CG" ,"No adjacent CG"))




# DMP__g2_g3__pp_nc__t
# DMP__g2_g3__pp_nc_PC1__t
ggplot(plt, aes(y=-TET3_Ravichandran_SciAdv_2022_unstranded, x=GLASS_OD__DMP__g2_g3__median, col=adjacent_cg, label=oligo_sequence)) +
  geom_point(size=theme_nature_size/3) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho", size=theme_nature_size,
                   label.x=-0.85) +
  scale_color_manual(
    values=c(
      'Adjacent CG' = mixcol( 'lightblue', 'lightgreen'),
      'No adjacent CG' = mixcol( 'darkblue', 'darkgreen')
    )
  ) +
  theme_nature +
  #ggrepel::geom_text_repel(data = subset(plt, grepl("cgcgcg|aacgtt|aacgtc|aacgtg", oligo_sequence)), col="black", size=theme_nature_size) +
  labs(x="T-score\nWHO grade 2 <> WHO grade 3",
       y="-1 * TET3 (Ravichandran SciAdv 2022)",
       col="",
       caption=paste0("n=",length(unique(plt$oligo_sequence))," motifs (unstranded)"),
       subtitle=format_subtitle("TET3 x t-score WHO grade")
  )


ggsave("output/figures/vis_differential_motifs__TET3_x_grade__unstranded.pdf", width=8.5 * 0.975 * (1/5), height=2.18)







### motif: gc_sequence_context_1_ed ----


plt <- plt |> 
  dplyr::mutate(gc_sequence_context_1_edit = gsub("[CG]",":", gc_sequence_context_1, fixed=T)) |> 
  dplyr::mutate(gc_sequence_context_1_edit = dplyr::case_when(
    gc_sequence_context_1_edit == "A:A" ~ "A:A | T:T | A:T | T:A",
    gc_sequence_context_1_edit == "T:T" ~ "A:A | T:T | A:T | T:A",
    
    gc_sequence_context_1_edit == "A:T" ~ "A:A | T:T | A:T | T:A",
    gc_sequence_context_1_edit == "T:A" ~ "A:A | T:T | A:T | T:A",
    
    gc_sequence_context_1_edit == "G:G" ~ "G:G | C:C | C:G | G:C",
    gc_sequence_context_1_edit == "C:C" ~ "G:G | C:C | C:G | G:C",
    gc_sequence_context_1_edit == "C:G" ~ "G:G | C:C | C:G | G:C",
    gc_sequence_context_1_edit == "G:C" ~ "G:G | C:C | C:G | G:C",
    
    
    T ~ "G:* | C:*"
  ))
  #dplyr::filter(gc_sequence_context_1_edit %in% c("A:A | T:T", "A:T | T:A"))


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_1_edit)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  + 
  scale_colour_Publication() 



### motif: gc_sequence_context_2 ----


stat <- plt |> 
  dplyr::select(gc_sequence_context_2, DMP__g2_g3__pp_nc__t)

ggplot(stat, aes(x=gc_sequence_context_2, y=DMP__g2_g3__pp_nc__t)) +
  ggbeeswarm::geom_quasirandom(size=0.1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))




### m1 m2 m3 ----


pplt1 <- plt |> 
  dplyr::mutate(label = dplyr::case_when(
    gc_sequence_context_2 %in% c("AA[CG]TT") ~ "AA[CG]TT (*palin)",
    gc_sequence_context_2 %in% c("AA[CG]TG", "CA[CG]TT") ~ "AA[CG]TG/CA[CG]TT",
    gc_sequence_context_2 %in% c("AA[CG]TC", "GA[CG]TT") ~ "AA[CG]TC/GA[CG]TT",
    T ~ "other"
  )) 


p1 = ggplot(pplt1 , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=label)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(pplt1, label == "other"), pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(pplt1, label != "other"), pch=16, cex=0.15, alpha=0.15) + 
  
  geom_point(data=head(pplt1, n=1),pch=16, cex=1, alpha=0.8) + # for label
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y = "Per probe t-score Primary ~ Recurrence") +
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) +
  
  scale_color_manual(
    
    values = c(
      "AA[CG]TC/GA[CG]TT" = "darkblue",
               "AA[CG]TG/CA[CG]TT" = "darkblue",
               "AA[CG]TT (*palin)" = "darkblue",
               "other" = "#F8766D")
  )



pplt2 <- plt |> 
  dplyr::mutate(label = dplyr::case_when(
    grepl("CG$",gc_sequence_context_2) ~ "xx[CG]CG suffix",
    grepl("^CG",gc_sequence_context_2) ~ "CG[CG]xx prefix",
    T ~ "other"
  )) 

p2 = ggplot(pplt2 , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=label)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(pplt2, label == "other"), pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(pplt2, label != "other"), pch=16, cex=0.15, alpha=0.15) + 
  
  geom_point(data=head(pplt2, n=1),pch=16, cex=1, alpha=0.8) + # for label
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3",
       y = "Per probe t-score Primary ~ Recurrence") +
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) +
  
  scale_color_manual(
    
    values = c(
      "xx[CG]CG suffix" = "darkblue",
      "CG[CG]xx prefix" = "darkblue",
      "other" = "#F8766D")
  )


p1 + p2




### motif: gc_sequence_context_2_prefix ----


stat <- plt |> 
  dplyr::select(gc_sequence_context_2, DMP__g2_g3__pp_nc__t) |> 
  dplyr::mutate(gc_sequence_context_2_prefix = gsub("..$","", gc_sequence_context_2))

ggplot(stat, aes(x=gc_sequence_context_2_prefix, y=DMP__g2_g3__pp_nc__t)) +
  ggbeeswarm::geom_quasirandom(size=0.1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))


### motif: gc_sequence_context_2_suffix ----


stat <- plt |> 
  dplyr::select(gc_sequence_context_2, DMP__g2_g3__pp_nc__t) |> 
  dplyr::mutate(gc_sequence_context_2_suffix = gsub("^..","", gc_sequence_context_2))

ggplot(stat, aes(x=gc_sequence_context_2_suffix, y=DMP__g2_g3__pp_nc__t)) +
  ggbeeswarm::geom_quasirandom(size=0.1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))



### motif: gc_sequence_context_l1 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l1 %in% c('A','T'))) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_l2 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l2 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_l3 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l3 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 




### motif: gc_sequence_context_l4 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l4 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_l5 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l5 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 


### motif: gc_sequence_context_l6 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l6 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 


### motif: gc_sequence_context_l7 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l7 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 


### motif: gc_sequence_context_l8 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_l8 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_r1 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r1 %in% c('A','T'))) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) # resize colbox
  


### motif: gc_sequence_context_r2 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r2 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  theme(legend.key.size = unit(0.6, 'lines')) 



### motif: gc_sequence_context_r3 ----




ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r3 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox
  



### motif: gc_sequence_context_r4 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r4 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox




### motif: gc_sequence_context_r5 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r5 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox




### motif: gc_sequence_context_r6 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r6 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox



### motif: gc_sequence_context_r7 ----



ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r7 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox




### motif: gc_sequence_context_r8 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_r8 )) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox




## x: all epigenetic clocks ----


clocks <- plt |> 
  dplyr::select(contains("dnaMethyAge") | contains("epiTOC2")) |> 
  colnames()


for(clock in clocks) {
  
  clock_name <- gsub("__up_nc__t","",gsub("DMP__","", clock))
  
  plt.p <- plt |> 
    data.table::copy() |> # odd hack needed, because setnames also affects the former "glass_od.metadata.array_samples" object...
    data.table::setnames(old = c(clock), new = c("array_current_clock"))
    
  
  ggplot(plt.p, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=array_current_clock)) +
    geom_vline(xintercept=0, col="red") +
    geom_hline(yintercept=0, col="red") +
    
    geom_point(pch=16, cex=0.001, alpha=0.15) + 
    
    labs(col = gsub("_", " ",clock_name)) +
    
    theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
    
    theme(plot.background = element_rect(fill="white")) + # png export
    
    ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", oob = scales::squish)
  
  
  ggsave(paste0("output/figures/vis_differential__epiGenetic_clock__",clock_name,".png"), width = (8.5*0.975/2), height = 4.3)
  
  
  
  
}

## early ~ late repli ----

plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(probe_id %in% data.mvalues.good_probes) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(wgEncodeUwRepliSeqBjWaveSignalRep1)) |> 
  dplyr::filter(wgEncodeUwRepliSeqBjWaveSignalRep2 < 79) |>  # outliers that mess up linear statistics
  dplyr::filter(wgEncodeUwRepliSeqBg02esWaveSignalRep1 < 90) |> 
  dplyr::filter(wgEncodeUwRepliSeqBjWaveSignalRep1 < 95) |> 
  dplyr::filter(wgEncodeUwRepliSeqBjWaveSignalRep1 > 10) |> 
  dplyr::filter(wgEncodeUwRepliSeqNhekWaveSignalRep1 < 80)


ggplot(plt, aes(y= wgEncodeUwRepliSeqBjWaveSignalRep1, x=glass_nl_prim_rec__deep_significant)) +
  ggbeeswarm::geom_quasirandom(size=0.1) +
  ggpubr::stat_compare_means(aes(group=glass_nl_prim_rec__deep_significant), label.x.npc=0.1, method = "wilcox.test", show.lengend  = FALSE,  size=theme_cellpress_size)


wilcox.test(
  plt |>
    dplyr::filter(glass_nl_prim_rec__deep_significant == T) |> 
    dplyr::pull(wgEncodeUwRepliSeqBjWaveSignalRep1)
  ,
  plt |>
    dplyr::filter(glass_nl_prim_rec__deep_significant == F) |> 
    dplyr::pull(wgEncodeUwRepliSeqBjWaveSignalRep1)
)



ggplot(plt, aes(x=wgEncodeUwRepliSeqBg02esWaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqBg02esWaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqBjWaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqBjWaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqBjWaveSignalRep2, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white"))  # png export
  #geom_vline(xintercept=79, col="red", lwd=0.5)
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqBjWaveSignalRep2.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm06990WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm06990WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm12801WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm12801WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm12812WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm12812WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm12813WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm12813WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqGm12878WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqGm12878WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqHelas3WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqHelas3WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqHepg2WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqHepg2WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqHuvecWaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqHuvecWaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqImr90WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqImr90WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqK562WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqK562WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqMcf7WaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqMcf7WaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqNhekWaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqNhekWaveSignalRep1.png",width=8.5/2, height=8.5/2)


ggplot(plt, aes(x=wgEncodeUwRepliSeqSknshWaveSignalRep1, y=DMP__g2_g3__pp_nc__t)) +
  geom_point(pch=16,cex=0.0001,alpha=0.10) +
  geom_smooth(method='lm', formula= y~x, col="red", lwd=theme_cellpress_lwd) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)), col="1", cor.coef.name ="rho") +
  theme_cellpress + theme(plot.background = element_rect(fill="white")) # png export
ggsave("output/figures/vis_differential__g23__wgEncodeUwRepliSeqSknshWaveSignalRep1.png",width=8.5/2, height=8.5/2)

