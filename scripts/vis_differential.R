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


## plt.motifs ----




plt.motifs <- data.frame(# R way of thinking
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
    ))) |> 
  dplyr::filter(!duplicated(oligo_sequence)) |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == (136))
    return(.)
  })()


### adds: DNMT1 - Adam NatCom 2020 ----
#' https://www.nature.com/articles/s41467-020-17531-8


tmp <- readxl::read_xlsx('data/DNMT1_Adam_et_al_2020_N2_data.xlsx') |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == (256))
    return(.)
  })() |> 
  dplyr::rename(sequence_5p = D1) |> 
  dplyr::mutate(sequence_5p = gsub(" ","", sequence_5p)) |> 
  dplyr::rename(DNMT1_Adam_NatCom_2020_stranded = av) |>
  assertr::verify(!duplicated(sequence_5p)) |> 
  assertr::verify(is.numeric(DNMT1_Adam_NatCom_2020_stranded)) |> 
  dplyr::mutate(sequence_3p = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(sequence_5p)))) |> 
  assertr::verify(!duplicated(sequence_5p)) |> 
  assertr::verify(!duplicated(sequence_3p)) |> 
  dplyr::select(sequence_5p, sequence_3p, DNMT1_Adam_NatCom_2020_stranded) |> 
  dplyr::mutate(palindromic = sequence_5p == sequence_3p) |> 
  dplyr::mutate(
    oligo_sequence = tolower(dplyr::case_when(
      sequence_5p < sequence_3p ~ paste0(sequence_5p,"/",sequence_3p),
      sequence_5p > sequence_3p ~ paste0(sequence_3p,"/",sequence_5p),
      palindromic ~ paste0("*",sequence_5p,"/",sequence_3p)
    ))) |> 
  dplyr::group_by(oligo_sequence) |> 
  dplyr::mutate(DNMT1_Adam_NatCom_2020 = mean(DNMT1_Adam_NatCom_2020_stranded)) |> 
  dplyr::ungroup() |> 
  dplyr::select(oligo_sequence, DNMT1_Adam_NatCom_2020) |> 
  dplyr::distinct() |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == (136))
    return(.)
  })() |> 
  assertr::verify(oligo_sequence %in% plt.motifs$oligo_sequence)



plt.motifs <- plt.motifs |> 
  dplyr::left_join(tmp, by=c('oligo_sequence'='oligo_sequence'), suffix=c('','')) |> 
  assertr::verify(!is.na(DNMT1_Adam_NatCom_2020))


rm(tmp)


### adds: DNMT1 - Adam NAR 2023 ----
#' https://academic.oup.com/nar/article-abstract/51/13/6622/7184160
#' "integrated with more data, but very similar"



tmp <- readxl::read_xlsx('data/DNMT1_Adam_NAR_2023_HM_combi_N2.xlsx') |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == (256))
    return(.)
  })() |> 
  dplyr::rename(sequence_5p = Site) |> 
  dplyr::mutate(sequence_5p = gsub(" ","", sequence_5p)) |> 
  dplyr::rename(DNMT1_Adam_NAR_2023_stranded = HM) |>
  assertr::verify(!duplicated(sequence_5p)) |> 
  assertr::verify(is.numeric(DNMT1_Adam_NAR_2023_stranded)) |> 
  dplyr::mutate(sequence_3p = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(sequence_5p)))) |> 
  assertr::verify(!duplicated(sequence_5p)) |> 
  assertr::verify(!duplicated(sequence_3p)) |> 
  dplyr::select(sequence_5p, sequence_3p, DNMT1_Adam_NAR_2023_stranded) |> 
  dplyr::mutate(palindromic = sequence_5p == sequence_3p) |> 
  dplyr::mutate(
    oligo_sequence = tolower(dplyr::case_when(
      sequence_5p < sequence_3p ~ paste0(sequence_5p,"/",sequence_3p),
      sequence_5p > sequence_3p ~ paste0(sequence_3p,"/",sequence_5p),
      palindromic ~ paste0("*",sequence_5p,"/",sequence_3p)
    ))) |> 
  dplyr::group_by(oligo_sequence) |> 
  dplyr::mutate(DNMT1_Adam_NAR_2023 = mean(DNMT1_Adam_NAR_2023_stranded)) |> 
  dplyr::ungroup() |> 
  dplyr::select(oligo_sequence, DNMT1_Adam_NAR_2023) |> 
  dplyr::distinct() |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == (136))
    return(.)
  })() |> 
  assertr::verify(oligo_sequence %in% plt.motifs$oligo_sequence)


plt.motifs <- plt.motifs |> 
  dplyr::left_join(tmp, by=c('oligo_sequence'='oligo_sequence'), suffix=c('','')) |> 
  assertr::verify(!is.na(DNMT1_Adam_NAR_2023))




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
  dplyr::rename(DNMT3A_Gao_NatCom_2020_stranded = 2) |> 
  dplyr::rename(DNMT3B_Gao_NatCom_2020_stranded = 5) |> 
  dplyr::select(sequence_full, DNMT3A_Gao_NatCom_2020_stranded, DNMT3B_Gao_NatCom_2020_stranded) |> 
  dplyr::mutate(sequence_5p = gsub("^.(..) CG (..).$","\\1CG\\2", sequence_full)) |> 
  dplyr::mutate(sequence_3p = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(sequence_5p)))) |> 
  dplyr::mutate(palindromic = sequence_5p == sequence_3p) |> 
  dplyr::mutate(
    oligo_sequence = tolower(dplyr::case_when(
      sequence_5p < sequence_3p ~ paste0(sequence_5p,"/",sequence_3p),
      sequence_5p > sequence_3p ~ paste0(sequence_3p,"/",sequence_5p),
      palindromic ~ paste0("*",sequence_5p,"/",sequence_3p)
    ))) |> 
  dplyr::group_by(oligo_sequence) |> 
  dplyr::mutate(DNMT3A_Gao_NatCom_2020 = mean(DNMT3A_Gao_NatCom_2020_stranded)) |> 
  dplyr::mutate(DNMT3B_Gao_NatCom_2020 = mean(DNMT3B_Gao_NatCom_2020_stranded)) |> 
  dplyr::ungroup() |> 
  dplyr::select(oligo_sequence, DNMT3A_Gao_NatCom_2020, DNMT3B_Gao_NatCom_2020) |> 
  dplyr::distinct() |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == (136))
    return(.)
  })() |> 
  assertr::verify(oligo_sequence %in% plt.motifs$oligo_sequence)



plt.motifs <- plt.motifs |> 
  dplyr::left_join(tmp, by=c('oligo_sequence'='oligo_sequence'), suffix=c('','')) |> 
  assertr::verify(!is.na(DNMT3A_Gao_NatCom_2020)) |> 
  assertr::verify(!is.na(DNMT3B_Gao_NatCom_2020))


rm(tmp)

### adds: TET1/2 - Adam NatCom 2022 ----
#' https://www.nature.com/articles/s42003-022-03033-4



tmp <- readxl::read_xlsx('data/42003_2022_3033_MOESM4_ESM.xlsx', skip=3) |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == (256))
    return(.)
  })() |> 
  dplyr::rename(sequence_5p = 1) |> 
  dplyr::mutate(sequence_5p = gsub(" ", "" ,sequence_5p)) |> 
  assertr::verify(!duplicated(sequence_5p)) |> 
  dplyr::rename(TET1_5mC_Adam_NatCom_2022_stranded  = 2) |> 
  dplyr::rename(TET2_5mC_Adam_NatCom_2022_stranded  = 4) |> 
  dplyr::rename(TET1_5hmC_Adam_NatCom_2022_stranded = 3) |> 
  dplyr::rename(TET2_5hmC_Adam_NatCom_2022_stranded = 5) |> 
  assertr::verify(is.numeric(TET1_5mC_Adam_NatCom_2022_stranded)) |> 
  assertr::verify(is.numeric(TET2_5mC_Adam_NatCom_2022_stranded)) |> 
  assertr::verify(is.numeric(TET1_5hmC_Adam_NatCom_2022_stranded)) |> 
  assertr::verify(is.numeric(TET2_5hmC_Adam_NatCom_2022_stranded)) |> 
  dplyr::mutate(sequence_3p = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(sequence_5p)))) |> 
  assertr::verify(!duplicated(sequence_5p)) |> 
  assertr::verify(!duplicated(sequence_3p)) |> 
  dplyr::select(sequence_5p, sequence_3p, starts_with("TET")) |> 
  dplyr::mutate(palindromic = sequence_5p == sequence_3p) |> 
  dplyr::mutate(
    oligo_sequence = tolower(dplyr::case_when(
      sequence_5p < sequence_3p ~ paste0(sequence_5p,"/",sequence_3p),
      sequence_5p > sequence_3p ~ paste0(sequence_3p,"/",sequence_5p),
      palindromic ~ paste0("*",sequence_5p,"/",sequence_3p)
    ))) |> 
  dplyr::group_by(oligo_sequence) |> 
  dplyr::mutate(TET1_5mC_Adam_NatCom_2022 = mean(TET1_5mC_Adam_NatCom_2022_stranded)) |> 
  dplyr::mutate(TET2_5mC_Adam_NatCom_2022 = mean(TET2_5mC_Adam_NatCom_2022_stranded)) |> 
  dplyr::mutate(TET1_5hmC_Adam_NatCom_2022 = mean(TET1_5hmC_Adam_NatCom_2022_stranded)) |> 
  dplyr::mutate(TET2_5hmC_Adam_NatCom_2022 = mean(TET2_5hmC_Adam_NatCom_2022_stranded)) |> 
  dplyr::ungroup() |> 
  dplyr::select(oligo_sequence, TET1_5mC_Adam_NatCom_2022, TET2_5mC_Adam_NatCom_2022, TET1_5hmC_Adam_NatCom_2022, TET2_5hmC_Adam_NatCom_2022) |> 
  dplyr::distinct() |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == (136))
    return(.)
  })() |> 
  assertr::verify(oligo_sequence %in% plt.motifs$oligo_sequence)



plt.motifs <- plt.motifs |> 
  dplyr::left_join(tmp, by=c('oligo_sequence'='oligo_sequence'), suffix=c('','')) |> 
  assertr::verify(!is.na(TET1_5mC_Adam_NatCom_2022)) |> 
  assertr::verify(!is.na(TET2_5mC_Adam_NatCom_2022)) |> 
  assertr::verify(!is.na(TET1_5hmC_Adam_NatCom_2022)) |> 
  assertr::verify(!is.na(TET2_5hmC_Adam_NatCom_2022))


rm(tmp)




### adds: TET3 - Ravichandran SciAdv 2022 ----
#' https://www.science.org/doi/10.1126/sciadv.abm2427



tmp <- read.table("data/10.1126_sciadv.abm2427_Slope_summary.txt") |> 
  dplyr::rename(sequence_5p = Motif) |> 
  tibble::rownames_to_column('sequence') |> 
  assertr::verify(sequence_5p == sequence) |> 
  dplyr::mutate(sequence = NULL) |> 
  dplyr::rename(TET3_Ravichandran_SciAdv_2022_stranded = ALL_84) |> 
  dplyr::tibble() |> 
  assertr::verify(!duplicated(sequence_5p)) |> 
  assertr::verify(is.numeric(TET3_Ravichandran_SciAdv_2022_stranded)) |> 
  dplyr::mutate(sequence_3p = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(sequence_5p)))) |> 
  assertr::verify(!duplicated(sequence_5p)) |> 
  assertr::verify(!duplicated(sequence_3p)) |> 
  dplyr::select(sequence_5p, sequence_3p, starts_with("TET")) |> 
  dplyr::mutate(palindromic = sequence_5p == sequence_3p) |> 
  dplyr::mutate(
    oligo_sequence = tolower(dplyr::case_when(
      sequence_5p < sequence_3p ~ paste0(sequence_5p,"/",sequence_3p),
      sequence_5p > sequence_3p ~ paste0(sequence_3p,"/",sequence_5p),
      palindromic ~ paste0("*",sequence_5p,"/",sequence_3p)
    ))) |> 
  dplyr::group_by(oligo_sequence) |> 
  dplyr::mutate(TET3_Ravichandran_SciAdv_2022 = mean(TET3_Ravichandran_SciAdv_2022_stranded)) |> 
  dplyr::ungroup() |> 
  dplyr::select(oligo_sequence, TET3_Ravichandran_SciAdv_2022) |> 
  dplyr::distinct() |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == (136))
    return(.)
  })() |> 
  assertr::verify(oligo_sequence %in% plt.motifs$oligo_sequence)



plt.motifs <- plt.motifs |> 
  dplyr::left_join(tmp, by=c('oligo_sequence'='oligo_sequence'), suffix=c('','')) |> 
  assertr::verify(!is.na(TET3_Ravichandran_SciAdv_2022))


rm(tmp)



### adds: DMP outcomes GLASS-OD ----



tmp <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |>
  dplyr::tibble() |> 
  dplyr::filter(!is.na(gc_sequence_context_2)) |> 
  dplyr::select(gc_sequence_context_2, DMP__g2_g3__pp_nc__t, DMP__primary_recurrence__pp_nc__t) |> 
  dplyr::mutate(sequence_5p = gsub("[CG]","CG",gc_sequence_context_2,fixed=T), gc_sequence_context_2=NULL) |> 
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
  dplyr::mutate(GLASS_OD__DMP__g2_g3__mean = mean(DMP__g2_g3__pp_nc__t)) |> 
  dplyr::mutate(GLASS_OD__DMP__g2_g3__median = median(DMP__g2_g3__pp_nc__t)) |>
  dplyr::mutate(GLASS_OD__DMP__primary_recurrence__mean = mean(DMP__primary_recurrence__pp_nc__t)) |> 
  dplyr::mutate(GLASS_OD__DMP__primary_recurrence__median = median(DMP__primary_recurrence__pp_nc__t)) |> 
  dplyr::ungroup() |> 
  dplyr::select(oligo_sequence, GLASS_OD__DMP__g2_g3__mean, GLASS_OD__DMP__g2_g3__median, GLASS_OD__DMP__primary_recurrence__mean, GLASS_OD__DMP__primary_recurrence__median) |> 
  dplyr::distinct() |> 
  (function(.) {
    assertthat::assert_that(nrow(.) == (136))
    return(.)
  })() |> 
  assertr::verify(oligo_sequence %in% plt.motifs$oligo_sequence)


plt.motifs <- plt.motifs |> 
  dplyr::left_join(tmp, by=c('oligo_sequence'='oligo_sequence'), suffix=c('','')) |> 
  assertr::verify(!is.na(GLASS_OD__DMP__g2_g3__mean)) |> 
  assertr::verify(!is.na(GLASS_OD__DMP__g2_g3__median)) |> 
  assertr::verify(!is.na(GLASS_OD__DMP__primary_recurrence__mean)) |> 
  assertr::verify(!is.na(GLASS_OD__DMP__primary_recurrence__median))


rm(tmp)



# plots ----

## A: overall density ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


p_A <- ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.0085, col="black") + 
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3", y="Per probe t-score Primary ~ Recurrence") +
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)





## B: AcCGAP in OD ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


p_B <- ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t,
                       col=DMP__AcCGAP__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3", y="Per probe t-score Primary ~ Recurrence", col="Association A_IDH_LG ~ A_IDH_HG") +
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)




## C: GLASS-NL probes ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


p_C <- ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=glass_nl_prim_rec__deep_significant)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == F),pch=16, cex=0.001, alpha=0.15, show.legend=F) + 
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == T),pch=16, cex=0.1, alpha=0.35, show.legend=F) + 
  geom_point(data = head(plt, n=0),pch=16, show.legend=T) + 

  geom_vline(xintercept=0, col="red", alpha=0.1) +
  geom_hline(yintercept=0, col="red", alpha=0.1) +
  
  labs(x = "Per probe t-score Grade 2 ~ Grade 3", y="Per probe t-score Primary ~ Recurrence", col="Significant GLASS-NL\nPrimary ~ Recurrent") +
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  scale_color_manual(values=c('TRUE'= col3(3)[1], 'FALSE'='gray80'))


p_A + p_B + p_C


ggsave("/tmp/papbpc.png",width = 8.5*0.975, height=3.2, dpi=600)




## C: G-CIMP ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


### a: IDHmut ----


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=-10*GCIMP_IDHmut_probe)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, GCIMP_IDHmut_probe == F),pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, GCIMP_IDHmut_probe == T),pch=16, cex=0.1, alpha=0.35) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


### b: IDHwt ----



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=-10*GCIMP_IDHwt_probe)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, GCIMP_IDHwt_probe == F),pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, GCIMP_IDHwt_probe == T),pch=16, cex=0.1, alpha=0.35) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)



### c: PanGlioma ----



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=-10*GCIMP_PanGlioma_probe)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, GCIMP_PanGlioma_probe == F),pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, GCIMP_PanGlioma_probe == T),pch=16, cex=0.1, alpha=0.35) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


### d: G-CIMP low 90 ----



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=-10*GCIMP_low_90)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, GCIMP_low_90 == F),pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, GCIMP_low_90 == T),pch=16, cex=0.1, alpha=0.35) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


## X: RepliTali predictors ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=-10*RepliTali_coef)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, RepliTali_coef == F | is.na(RepliTali_coef)),pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, RepliTali_coef == T),pch=16, cex=1, alpha=0.35) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(RepliTali_coef) & RepliTali_coef != F)



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=RepliTali_coef)) +
  #geom_vline(xintercept=0, col="red") +
  #geom_hline(yintercept=0, col="red") +
  
  #geom_point(data = subset(plt, RepliTali_coef == F | is.na(RepliTali_coef)),pch=16, cex=0.001, alpha=0.15) + 
  #geom_point(data = subset(plt, RepliTali_coef == T),pch=16, cex=1, alpha=0.35) + 
  geom_point(pch=16, cex=1, alpha=0.35) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox


ggplot(plt, aes(x=DMP__primary_recurrence__pp_nc__t, y=RepliTali_coef)) +
  #geom_vline(xintercept=0, col="red") +
  #geom_hline(yintercept=0, col="red") +
  
  #geom_point(data = subset(plt, RepliTali_coef == F | is.na(RepliTali_coef)),pch=16, cex=0.001, alpha=0.15) + 
  #geom_point(data = subset(plt, RepliTali_coef == T),pch=16, cex=1, alpha=0.35) + 
  geom_point(pch=16, cex=1, alpha=0.35) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines'))  # resize colbox



## death clock ----



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::mutate(death_clock_probe = probe_id %in% c("cg03725309","cg25763716","cg13854219","cg25189904","cg15459165","cg19266329","cg24397007","cg23079012","cg27241845","cg06905155","cg16503724","cg19859270","cg02657160","cg14975410","cg05575921","cg14817490","cg21161138","cg12513616","cg20732076","cg06126421","cg15342087","cg01612140","cg25983901","cg12510708","cg00285394","cg01140244","cg23190089","cg07123182","cg26963277","cg18550212","cg10321156","cg25193885","cg07986378","cg23665802","cg04987734","cg19459791","cg00310412","cg23842572","cg19572487","cg01572694","cg08546016","cg18181703","cg03636183","cg24704287","cg11341610","cg14085840","cg26470501","cg05492306","cg25607249","cg01406381","cg07626482","cg03707168","cg08362785"))



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=-10*death_clock_probe, label=probe_id)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, death_clock_probe == F | is.na(death_clock_probe)),pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, death_clock_probe == T),pch=16, cex=1, alpha=0.35) + 
  
  #ggrepel::geom_text_repel(data = subset(plt, probe_id %in% 
  #                                         c("cg06905155","cg18181703","cg03636183","cg24704287","cg26470501")),nudge_y = 1) +
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)






## D: FFPE | Tissue ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__FFPE_or_FF__pp_nc__t)) + # 
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)



## E: FFPE decay time INTENSITY ----




plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == (693060)) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t,
                col=DPI__FFPE_decay_time__pp_nc__t)) +
  #facet_wrap(~probe_type, scales="free",ncol=2) +
  facet_wrap(~map, scales="free",ncol=2) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-8, 8), oob = scales::squish)


plt |> 
  dplyr::filter(type == "I") |> 
  dplyr::select(Forward_Sequence, AlleleB_ProbeSeq_Illumina_manifest,  AlleleA_ProbeSeq_Illumina_manifest) |> 
  View()


plt |> 
  dplyr::filter(type == "II") |> 
  dplyr::select(probe_id, Forward_Sequence, AlleleA_ProbeSeq_Illumina_manifest, mapFlag_A, mapQ_A, DMP__FFPE_decay_time__pp_nc__t) |> 
  dplyr::arrange(DMP__FFPE_decay_time__pp_nc__t) |> 
  View()




## E: FFPE decay time ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!grepl(" \\? ", probe_type_orientation))


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t,
                col=DMP__FFPE_decay_time__pp_nc__t)) +
  facet_wrap(~probe_type_orientation, ncol=4) + 
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)



a = plt |> 
  dplyr::select(CHR_hg38, gc_sequence_context_2, Forward_Sequence, DMP__FFPE_decay_time__pp_nc__t) |> 
  dplyr::arrange(DMP__FFPE_decay_time__pp_nc__t) 


b = plt |> 
  dplyr::select(CHR_hg38, gc_sequence_context_2, Forward_Sequence,DMP__FFPE_decay_time__pp_nc__t) |> 
  dplyr::arrange(-DMP__FFPE_decay_time__pp_nc__t)


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t,
                col=DMP__FFPE_decay_time__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)




### sandbox ----
#' uitzoeken wat 


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(abs(DMP__FFPE_decay_time__pp_nc__t) > 4.5) |> 
  dplyr::filter(grepl("?", probe_type_orientation))


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=c_content)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-8, 8), oob = scales::squish)



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=g_content)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(0, 1), oob = scales::squish)



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_content)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(0, 1), oob = scales::squish)



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=probeCpGcnt)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", trans = "log1p",
                                 breaks = 0:11,
                                 labels = c(0,"",2,"","","",6,"","","","",11)
  )  # , oob = scales::squish





ggplot(plt, aes(x=DMP__FFPE_decay_time__pp_nc__t, y=gc_content, col=probeCpGcnt)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", trans = "log1p",
                                 breaks = 0:11,
                                 labels = c(0,"",2,"","","",6,"","","","",11)
  )  # , oob = scales::squish





ggplot(plt, aes(x=DMP__FFPE_decay_time__pp_nc__t, y=gc_content, col=probeCpGcnt)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", trans = "log1p",
                                 breaks = 0:11,
                                 labels = c(0,"",2,"","","",6,"","","","",11)
  )  # , oob = scales::squish



ggplot(plt, aes(x=DMP__FFPE_decay_time__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=probeCpGcnt)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", trans = "log1p",
                                 breaks = 0:11,
                                 labels = c(0,"",2,"","","",6,"","","","",11)
  )  # , oob = scales::squish





## probe type




ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__FFPE_decay_time__pp_nc__t)) +
  facet_wrap(~type, scales="free",ncol=2) +

  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.5, alpha=0.75) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-7, 7), oob = scales::squish)




ggplot(plt, aes(x=type, y=probeCpGcnt)) +
  ggbeeswarm::geom_quasirandom(size=0.01)


ggplot(plt, aes(x=as.factor(probeCpGcnt), y=DMP__FFPE_decay_time__pp_nc__t)) +
  facet_wrap(~type, scales="free",ncol=2) +
  geom_violin(draw_quantiles = c(0.5), col="red")

#' wanneer een andere CpG wordt aangepast, wordt deze de-methylated?


## uberhaupt aantal C's



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(abs(DMP__epiTOC2_tnsc__up_nc__t) > 4.5 & !is.na(gc_content))


ggplot(plt, aes(x=probeCpGcnt, y=DMP__primary_recurrence__pp_nc__t, col=gc_content)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.1, alpha=0.75) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(0, 1), oob = scales::squish)




## uberhaupt aantal G's
## overall GC content
## overall intensity
## median b-value





## F: epiTOC2 / tsnc ----
#' normal epiTOC


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__epiTOC2_tnsc__up_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 

  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)




## G: epiTOC2 / [hypoSC!!] ----

plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()

ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=DMP__epiTOC2_hypoSC__pp_nc__t)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == F),pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt, glass_nl_prim_rec__deep_significant == T),pch=16, cex=0.1, alpha=0.35) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  ggplot2::scale_color_gradientn(colours = col3(200), na.value = "grey50", limits = c(-10, 10), oob = scales::squish)



## C content ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(c_content))


plt <- plt |> 
  dplyr::mutate(
    c_content_f = cut(c_content,
                      breaks=quantile(plt$c_content,
                      probs = seq(0,1, by=0.10),
                      na.rm=TRUE, 
                      names=TRUE,
                      include.lowest=TRUE, 
                      right = TRUE)))



p1 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=c_content)) +
  #facet_grid(cols = vars(Infinium_Design_Type))  +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, c_content_f == levels(plt$c_content_f)[1]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, c_content_f == levels(plt$c_content_f)[2]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, c_content_f == levels(plt$c_content_f)[3]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, c_content_f == levels(plt$c_content_f)[4]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, c_content_f == levels(plt$c_content_f)[5]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, c_content_f == levels(plt$c_content_f)[6]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, c_content_f == levels(plt$c_content_f)[7]), pch=16, cex=0.05, alpha=0.50) + 
  geom_point(data=subset(plt, c_content_f == levels(plt$c_content_f)[8]), pch=16, cex=0.1, alpha=0.55) + 
  geom_point(data=subset(plt, c_content_f == levels(plt$c_content_f)[9]), pch=16, cex=0.1, alpha=0.60) + 
  geom_point(data=subset(plt, c_content_f == levels(plt$c_content_f)[10]),  pch=16, cex=0.1, alpha=0.65) + 
  
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", 
                                 limits = c(0, 0.5), oob = scales::squish ) + 
  
  labs(col = "C content") +
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white")) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))



## G content ----



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(c_content))


plt <- plt |> 
  dplyr::mutate(
    g_content_f = cut(g_content,
                      breaks=quantile(plt$g_content,
                                      probs = seq(0,1, by=0.10),
                                      na.rm=TRUE, 
                                      names=TRUE,
                                      include.lowest=TRUE, 
                                      right = TRUE)))



p2 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=g_content)) +
  #facet_grid(cols = vars(Infinium_Design_Type))  +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, g_content_f == levels(plt$g_content_f)[1]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, g_content_f == levels(plt$g_content_f)[2]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, g_content_f == levels(plt$g_content_f)[3]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, g_content_f == levels(plt$g_content_f)[4]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, g_content_f == levels(plt$g_content_f)[5]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, g_content_f == levels(plt$g_content_f)[6]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, g_content_f == levels(plt$g_content_f)[7]), pch=16, cex=0.05, alpha=0.50) + 
  geom_point(data=subset(plt, g_content_f == levels(plt$g_content_f)[8]), pch=16, cex=0.1, alpha=0.55) + 
  geom_point(data=subset(plt, g_content_f == levels(plt$g_content_f)[9]), pch=16, cex=0.1, alpha=0.60) + 
  geom_point(data=subset(plt, g_content_f == levels(plt$g_content_f)[10]),  pch=16, cex=0.1, alpha=0.65) + 
  
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50",
                                 limits = c(0, 0.5), oob = scales::squish
                                 ) + 
  
  labs(col = "G content") +
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white")) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))



## GC content ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::filter(!is.na(c_content) & !is.na(g_content))


plt <- plt |> 
  dplyr::mutate(
    gc_content_f = cut(gc_content,
                      breaks=quantile(plt$gc_content,
                                      probs = seq(0,1, by=0.10),
                                      na.rm=TRUE, 
                                      names=TRUE,
                                      include.lowest=TRUE, 
                                      right = TRUE)))



p3 = ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_content)) +
  #facet_grid(cols = vars(Infinium_Design_Type))  +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, gc_content_f == levels(plt$gc_content_f)[1]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, gc_content_f == levels(plt$gc_content_f)[2]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, gc_content_f == levels(plt$gc_content_f)[3]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, gc_content_f == levels(plt$gc_content_f)[4]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, gc_content_f == levels(plt$gc_content_f)[5]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, gc_content_f == levels(plt$gc_content_f)[6]), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, gc_content_f == levels(plt$gc_content_f)[7]), pch=16, cex=0.05, alpha=0.50) + 
  geom_point(data=subset(plt, gc_content_f == levels(plt$gc_content_f)[8]), pch=16, cex=0.1, alpha=0.55) + 
  geom_point(data=subset(plt, gc_content_f == levels(plt$gc_content_f)[9]), pch=16, cex=0.1, alpha=0.60) + 
  geom_point(data=subset(plt, gc_content_f == levels(plt$gc_content_f)[10]),  pch=16, cex=0.1, alpha=0.65) + 
  
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", 
                                 limits = c(0, 1), oob = scales::squish
                                 ) + 
  
  labs(col = "GC content") +
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white")) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))




## H: probeCpGcnt ----


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=probeCpGcnt)) +
  #facet_grid(cols = vars(Infinium_Design_Type))  +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, probeCpGcnt == 0), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, probeCpGcnt == 1), pch=16, cex=0.0025, alpha=0.15) + 
  geom_point(data=subset(plt, probeCpGcnt == 2), pch=16, cex=0.005, alpha=0.20) + 
  geom_point(data=subset(plt, probeCpGcnt == 3), pch=16, cex=0.05, alpha=0.25) + 
  geom_point(data=subset(plt, probeCpGcnt == 4), pch=16, cex=0.01, alpha=0.30) + 
  geom_point(data=subset(plt, probeCpGcnt == 5), pch=16, cex=0.01, alpha=0.35) + 
  geom_point(data=subset(plt, probeCpGcnt == 6), pch=16, cex=0.01, alpha=0.40) + 
  geom_point(data=subset(plt, probeCpGcnt == 7), pch=16, cex=0.025, alpha=0.45) + 
  geom_point(data=subset(plt, probeCpGcnt == 8), pch=16, cex=0.05, alpha=0.50) + 
  geom_point(data=subset(plt, probeCpGcnt == 9), pch=16, cex=0.1, alpha=0.55) + 
  geom_point(data=subset(plt, probeCpGcnt == 10), pch=16, cex=0.1, alpha=0.60) + 
  geom_point(data=subset(plt, probeCpGcnt == 11), pch=16, cex=0.1, alpha=0.65) + 
  
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", trans = "log1p",
                                 breaks = 0:11,
                                 labels = c(0,"",2,"","","",6,"","","","",11)
  ) + # , oob = scales::squish
  
  labs(col = "CpG's per probe") +
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white")) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))


ggsave("output/figures/vis_differential__g23_prim-rec__probe_CpG_count.png", width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))




## I: closest_CG  ----

plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::mutate(closest_CG_25 = dplyr::case_when(
    closest_CG > 25 ~ 26,
    closest_CG < -25 ~ -26,
    T ~ closest_CG
  ))




ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=abs(closest_CG_25))) +
  
  facet_wrap(~abs(closest_CG_25) > 8, scales="free",ncol=5) +
  
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  ggplot2::scale_color_gradientn(colours = rev(col3(200))[1:175], na.value = "grey50", 
                                 #trans = "log1p",
                                 breaks = c(1,5,10,15,20,26),
                                 labels = c(1,"",10,"",20,"  25+")
  ) + # , oob = scales::squish
  
  labs(col = "distance to next CpG (bp)") +
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white")) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))



ggsave("output/figures/vis_differential__g23_prim-rec__closest_CG.png", width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))


### b: in quadrant ----


ggplot(plt, aes(y=DMP__g2_g3__pp_nc__t, x=as.factor(closest_CG))) +
  #geom_point(pch=16, cex=0.001, alpha=0.15, col="darkgray") +
  ggplot2::geom_violin(draw_quantiles = c(0.5), linewidth=theme_cellpress_lwd, col = "darkgray", adjust = 1.95) +
  theme_cellpress +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 




## J: solo-wCGw ----
# >= 35 & == wCGw [https://www.nature.com/articles/s41467-022-34268-8]


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=is_solo_WCGW)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, is_solo_WCGW==F) , pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=subset(plt, is_solo_WCGW==T) , pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=0) , pch=16) + # hack plotting empty data frame, but with right cex and alpha for guide
  

  labs(col = "is solo-WCGW", x = "Per probe t-score Grade 2 ~ Grade 3", y="Per probe t-score Primary ~ Recurrence") +
  scale_color_manual(values=c('TRUE'= col3(3)[1], 'FALSE'='gray80')) +

  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white")) +  # png export
  theme(legend.key.size = unit(0.6, 'lines'))

ggsave("output/figures/vis_differential__g23_prim-rec__solo-WCGW.png", width=(8.5 * 0.97 / 2), height=(8.5 * 0.97 / 2))



## K: discrete motifs ----
### a: motifs alone ----

plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |>
  dplyr::filter(!is.na(gc_sequence_context_2)) |> 
  
  dplyr::mutate(cg_prefix_suffix = grepl("^CG", gc_sequence_context_2) | grepl("CG$", gc_sequence_context_2)) |> 
  dplyr::mutate(top_affinity = gc_sequence_context_2 %in% c(
    "GA[CG]AG", "GA[CG]TC","CT[CG]TT","GA[CG]TG", "CT[CG]TT","CA[CG]TG","CA[CG]TT","GA[CG]TT","AA[CG]TT", 
    #            palin                                        palin                            palin
    "CT[CG]TC", "GA[CG]TC","CA[CG]AG","CA[CG]TC", "AA[CG]AG","CA[CG]TG","AA[CG]TG","AA[CG]TC","AA[CG]TT" )) |> 
  
  dplyr::mutate(col_facet1 = ifelse(cg_prefix_suffix, "CG prefix/suffix", "other")) |> 
  dplyr::mutate(col_facet2 = ifelse(top_affinity, "top affinity", "other")) 
  

plt <- rbind(
  plt |> 
    dplyr::mutate(col = col_facet1) |> 
    dplyr::mutate(facet = "CG prefix/suffix")
  ,
  plt |> 
    dplyr::mutate(col = col_facet2) |> 
    dplyr::mutate(facet = "top affinity")
) |> 
  dplyr::mutate(col = factor(col, levels=c(   "top affinity"  ,"CG prefix/suffix" ,"other"   )))


levels(plt$col)


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=col)) +
  facet_wrap(~facet, scales="free") +

  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, col == "other"), pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=subset(plt, col != "other"), pch=16, cex=0.001, alpha=0.15, show.legend = F) + 

  geom_point(data=head(plt, n=0),pch=16) + 
  
  labs(col = "sequence context", x = "Per probe t-score Grade 2 ~ Grade 3", y="Per probe t-score Primary ~ Recurrence") +
  scale_color_manual(values=c('top affinity'= col3(3)[1],
                              'other'='gray80',
                              'CG prefix/suffix'= col3(3)[3])) +
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white"))



### b: motifs x min closest_CG ----



## G: 1P ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()

ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=is_1P)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, !is_1P), pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=subset(plt,  is_1P), pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white"))




## G: 19Q ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()

ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=is_19Q)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, !is_19Q), pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=subset(plt,  is_19Q), pch=16, cex=0.001, alpha=0.15, show.legend = F) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white"))


## G: 1P | 19Q ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |> 
  dplyr::mutate(col = dplyr::case_when(
    CpG_chrm == "chr1" ~ CpG_chrm,
    CpG_chrm == "chr19" ~ CpG_chrm,
    T ~ "other"
  )) |> 
  dplyr::mutate(facet = ifelse(is_1P | is_19Q, "on 1p | 19q","other")) |> 
  dplyr::mutate( CpG_chrm = dplyr::case_when(
    CpG_chrm == "chr1" & is_1P ~ "chr1p",
    CpG_chrm == "chr1" & !is_1P ~ "chr1q",
  
    CpG_chrm == "chr19" & is_19Q ~ "chr19q",
    CpG_chrm == "chr19" & !is_19Q ~ "chr19p",
    T ~ CpG_chrm
  )) |> 
  dplyr::mutate(chr = factor(CpG_chrm, levels=gtools::mixedsort(unique(as.character(CpG_chrm))) ))



ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=probeCpGcnt)) +
  facet_grid(cols = vars(facet)) +

  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data=subset(plt, probeCpGcnt == 0), pch=16, cex=0.001, alpha=0.1) + 
  geom_point(data=subset(plt, probeCpGcnt == 1), pch=16, cex=0.0025, alpha=0.15) + 
  geom_point(data=subset(plt, probeCpGcnt == 2), pch=16, cex=0.005, alpha=0.20) + 
  geom_point(data=subset(plt, probeCpGcnt == 3), pch=16, cex=0.05, alpha=0.25) + 
  geom_point(data=subset(plt, probeCpGcnt == 4), pch=16, cex=0.01, alpha=0.30) + 
  geom_point(data=subset(plt, probeCpGcnt == 5), pch=16, cex=0.01, alpha=0.35) + 
  geom_point(data=subset(plt, probeCpGcnt == 6), pch=16, cex=0.01, alpha=0.40) + 
  geom_point(data=subset(plt, probeCpGcnt == 7), pch=16, cex=0.025, alpha=0.45) + 
  geom_point(data=subset(plt, probeCpGcnt == 8), pch=16, cex=0.05, alpha=0.50) + 
  geom_point(data=subset(plt, probeCpGcnt == 9), pch=16, cex=0.1, alpha=0.55) + 
  geom_point(data=subset(plt, probeCpGcnt == 10), pch=16, cex=0.1, alpha=0.60) + 
  geom_point(data=subset(plt, probeCpGcnt == 11), pch=16, cex=0.1, alpha=0.65) + 
  geom_point(data=head(plt, n=1),pch=16, cex=1, alpha=0.8) + 
  
  ggplot2::scale_color_gradientn(colours = rev(col3(200)), na.value = "grey50", trans = "log1p",
                                 breaks = 0:11,
                                 labels = c(0,"",2,"","","",6,"","","","",11)
  ) + # , oob = scales::squish
  
  labs(col = "CpG's per probe") +
  
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox
  theme(plot.background = element_rect(fill="white"))


ggplot(plt, aes(x = chr, y=DMP__g2_g3__pp_nc__t, fill=col)) +
  ggplot2::geom_violin(draw_quantiles = c(0.5), linewidth=theme_cellpress_lwd, col = "white", adjust = 1.95) +
  scale_fill_manual(values=c('chr1'='red','chr19'='blue','other'='darkgray')) +
  theme_cellpress +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


## L: PC1 ----

## M: PC2 ----

## N: PC3 ----



## X: Infinium_Design_Type ----
#' not need to be shown, effect is by CpG's per probe,
#' design type is defined by number of probes



plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })()


ggplot(plt, aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=Infinium_Design_Type)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(data = subset(plt, Infinium_Design_Type == "II"), pch=16, cex=0.001, alpha=0.15) + 
  geom_point(data = subset(plt,Infinium_Design_Type == "I"), pch=16, cex=0.001, alpha=0.15) +   
  
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) 


## Y: per oligo/motifs correlation ----


tmp.cor <- plt.motifs |> 
  dplyr::mutate(`TET3_Ravichandran_SciAdv_2022 * -1` = -1 * TET3_Ravichandran_SciAdv_2022) |> 
  dplyr::select(GLASS_OD__DMP__g2_g3__mean,
                TET1_5mC_Adam_NatCom_2022, TET2_5mC_Adam_NatCom_2022, `TET3_Ravichandran_SciAdv_2022 * -1`,
                DNMT1_Adam_NAR_2023, DNMT3A_Gao_NatCom_2020, DNMT3B_Gao_NatCom_2020) |> 
  cor(method = "spearman") |> 
  as.data.frame()  |> 
  tibble::rownames_to_column('enzyme') |> 
  dplyr::select(enzyme, GLASS_OD__DMP__g2_g3__mean) |> 
  dplyr::rename(cor_w_GLASS_OD__DMP__g2_g3__mean = GLASS_OD__DMP__g2_g3__mean) |> 
  dplyr::filter(enzyme != "GLASS_OD__DMP__g2_g3__mean")
  



plt <- plt.motifs |> 
  dplyr::mutate(`TET3_Ravichandran_SciAdv_2022 * -1` = -1 * TET3_Ravichandran_SciAdv_2022) |> 
  dplyr::select(oligo_sequence, GLASS_OD__DMP__g2_g3__mean,
                TET1_5mC_Adam_NatCom_2022, TET2_5mC_Adam_NatCom_2022, `TET3_Ravichandran_SciAdv_2022 * -1`,
                DNMT1_Adam_NAR_2023, DNMT3A_Gao_NatCom_2020, DNMT3B_Gao_NatCom_2020) |> 
  tidyr::pivot_longer(cols = -c(oligo_sequence, GLASS_OD__DMP__g2_g3__mean), names_to = "enzyme", values_to = "relative de-methylation affinity") |> 
  dplyr::left_join(tmp.cor, by=c('enzyme'='enzyme'), suffix=c('','') ) |> 
  dplyr::mutate(enzyme = gsub("_"," ", enzyme))




ggplot(plt, aes(x=GLASS_OD__DMP__g2_g3__mean, y=`relative de-methylation affinity`, label=oligo_sequence, col=cor_w_GLASS_OD__DMP__g2_g3__mean)) +
  facet_wrap(~enzyme, scales="free",ncol=3) +
  geom_point(size=theme_cellpress_size / 3) +
  ggpubr::stat_cor(method = "spearman", aes(label = after_stat(r.label)),  cor.coef.name ="rho", size=theme_cellpress_size, label.x.npc="right", hjust=1, col="black") +
  labs(x = "Mean t-score Grade 2 ~ Grade for probes with given sequence context",
       col = "Spearman correlation enzyme affinity with mean t-score") +
  theme_cellpress +
  ggplot2::scale_color_gradientn(colours = col4(200), na.value = "grey50", oob = scales::squish, 
                                 limits=c(-1,1),
                                 breaks = c(-1,-0.5,0,0.5,1),
                                 labels = c("-1","",0,"",1)
                                 )  +
  theme(legend.key.size = unit(0.6, 'lines'))


ggsave("output/figures/vis_differential__motif_tscores_x_enzyme_kinetics.pdf", width = 6.5, height=4.5)





## Z: motif corrplot ----


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



## motif: xx[CG]xx violin(s) ----


plt <- data.mvalues.probes |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED) 
    return(.)
  })() |> 
  dplyr::filter(detP_good_probe & !MASK_general) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == CONST_N_PROBES_UNMASKED_AND_DETP) 
    return(.)
  })() |>
  dplyr::filter(!is.na(gc_sequence_context_2))


plt.per.context <- plt |> 
  dplyr::select(gc_sequence_context_2) |> 
  dplyr::filter(!duplicated(gc_sequence_context_2)) |> 
  dplyr::mutate(gc_sequence_context_2_s = gsub("[CG]","CG",gc_sequence_context_2,fixed=T)) |> 
  dplyr::mutate(gc_sequence_context_2_rc = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(gc_sequence_context_2_s)))) |> 
  dplyr::mutate(facet_name = dplyr::case_when(
    gc_sequence_context_2_s < gc_sequence_context_2_rc ~ paste0(gc_sequence_context_2_s ,"/",gc_sequence_context_2_rc),
    gc_sequence_context_2_s > gc_sequence_context_2_rc ~ paste0(gc_sequence_context_2_rc,"/",gc_sequence_context_2_s ),
    gc_sequence_context_2_s == gc_sequence_context_2_rc ~ paste0("*",gc_sequence_context_2_s,"/",gc_sequence_context_2_s)
  ))
#dplyr::mutate(gc_sequence_context_2_s = NULL, gc_sequence_context_2_rc = NULL)


plt <- plt |> 
  dplyr::left_join(plt.per.context, by=c('gc_sequence_context_2'='gc_sequence_context_2'), suffix=c('','')) |> 
  dplyr::group_by(facet_name) |> 
  dplyr::mutate(facet_rank = median(DMP__g2_g3__pp_nc__t)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(facet_name = tolower(facet_name))


ggplot(plt, aes(x=reorder(facet_name, -facet_rank), y=DMP__g2_g3__pp_nc__t)) +
  #facet_wrap(~reorder(facet_name, -facet_rank), scales="free_x", ncol=length(unique(plt$facet_name))) +
  #ggbeeswarm::geom_quasirandom(size=0.01, alpha=0.65) +
  #ggplot2::geom_violin(draw_quantiles = c(), linewidth=theme_cellpress_lwd, col = "white", fill="darkgray", adjust = 1.95) +
  geom_boxplot(width=0.75, outlier.shape=NA, outlier.color=NA, col="darkgray",
              ,coef=0.5, fill=NA,linewidth=theme_cellpress_lwd) +
  
  stat_summary(fun.y = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5, col="red", width=0.85, linewidth=theme_cellpress_lwd) +
  
  labs(x = NULL, y="Per probe t-score Grade 2 ~ Grade 3") +
  coord_cartesian(ylim = c(-6.75, 3.5)) + # soft clip
  theme_cellpress + theme(legend.key.size = unit(0.6, 'lines')) + # resize colbox +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, family = "mono"))


ggsave("output/figures/vis_differential__xxCGxx_violin.pdf", width = 11 * 0.97, height=3.75)



### motif: gc_sequence_context_1 ----


ggplot(plt , aes(x=DMP__g2_g3__pp_nc__t, y=DMP__primary_recurrence__pp_nc__t, col=gc_sequence_context_1)) +
  geom_vline(xintercept=0, col="red") +
  geom_hline(yintercept=0, col="red") +
  
  geom_point(pch=16, cex=0.001, alpha=0.15) + 
  
  theme_cellpress + 
  
  theme(plot.background = element_rect(fill="white")) +  # png export
  
  theme(legend.key.size = unit(0.6, 'lines')) 


stat <- plt |> 
  dplyr::select(gc_sequence_context_1, DMP__g2_g3__pp_nc__t)

ggplot(stat, aes(x=gc_sequence_context_1, y=DMP__g2_g3__pp_nc__t)) +
  ggbeeswarm::geom_quasirandom(size=theme_cellpress_size/2) +
  theme_bw()


#### RC ----






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
  dplyr::filter(detP_good_probe & !MASK_general) |> 
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



# Figure S2A: corr covars per sample ----
## corr ----


plt <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(210) |> 
  dplyr::mutate(time_tissue_in_ffpe =  ifelse(isolation_material == "ffpe", time_between_resection_and_array, 0)) |> 
  dplyr::select(array_epiTOC2_tnsc, array_epiTOC2_hypoSC, contains("array_dnaMethyAge"), array_RepliTali,
                array_percentage.detP.signi, array_PC1, `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`,
                time_tissue_in_ffpe,
                array_GLASS_NL_g2_g3_sig, array_PC2, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit,
                
                array_PC3
  ) |> 
  dplyr::filter(!is.na(time_tissue_in_ffpe)) |> 
  
  dplyr::mutate(array_GLASS_OD_g2_g3_sig2 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig2 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_g2_g3_sig3 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig3 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_g2_g3_sig4 = NULL) |> 
  dplyr::mutate(array_GLASS_OD_prim_rec_sig4 = NULL) |> 
  
  dplyr::mutate(array_PC3 = -1 * array_PC3) |> 
  
  dplyr::mutate(AcCGAP = -1 * array_A_IDH_HG__A_IDH_LG_lr__lasso_fit, array_A_IDH_HG__A_IDH_LG_lr__lasso_fit = NULL) |> 
  dplyr::mutate(`-1 * array_dnaMethyAge__ZhangY2017` = -1 * array_dnaMethyAge__ZhangY2017 , array_dnaMethyAge__ZhangY2017 = NULL) |> 
  dplyr::mutate(`-1 * array_epiTOC2_hypoSC` = -1 * array_epiTOC2_hypoSC, array_epiTOC2_hypoSC = NULL) |>  # seems at inversed scale
  dplyr::mutate(`-1 * array_dnaMethyAge__LuA2019` = -1 * array_dnaMethyAge__LuA2019, array_dnaMethyAge__LuA2019 = NULL) |>  # seems at inversed scale
  dplyr::mutate(`-1 * QC: SPECIFICITY I GT MM 6` = -1 * `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`, `array_qc_SPECIFICITY_I_GT_Mismatch_6_(PM)_Red_smaller_NA_0.5`=NULL) |> 
  dplyr::select(-array_dnaMethyAge__YangZ2016) |>  # epiTOC1, near identical to epiTOC2
  dplyr::select(-array_dnaMethyAge__PCHorvathS2013) |>  # very similar to its 2018 equivalent
  dplyr::mutate(`array_qc_BISULFITE_CONVERSION_I_Beta_I-3_Beta_larger_0.12_0.18` = NULL) |> # contains N/A's
  dplyr::mutate(`array_qc_BISULFITE_CONVERSION_I_Beta_I-6_Beta_larger_0.2_0.3` = NULL)# contains N/A's


colnames(plt) <- gsub("array_","",colnames(plt))
colnames(plt) <- gsub("_"," ",colnames(plt))

corrplot::corrplot(cor(plt, method="spearman"), order="hclust", tl.cex=0.75, tl.pos="l")
#h = corrplot::corrplot(cor(plt, method="spearman"), order="hclust", tl.cex=0.75, tl.pos="l")

#p1 = ggcorrplot2::ggcorrplot(h$corr) +
#  theme(axis.title.x=element_blank(),
#        axis.text.x=element_blank(),
#        axis.ticks.x=element_blank())




