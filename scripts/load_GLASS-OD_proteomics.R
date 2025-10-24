#!/usr/bin/env R

# load data ----


source('scripts/load_constants.R')
source('scripts/load_functions.R')
#source('scripts/load_palette.R')


if(!exists('glass_od.metadata.proteomics')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# load proteomics raw per-peptide ----


#tmp <- readxl::read_xlsx('data/GLASS_OD/Protein - Tobias Weiss/GLODprot_raw_protein_matrix.xlsx')


# metadata.proteomics.glass_od <- tmp |> 
#     dplyr::select(Precursor.Id, Protein.Group, Protein.Ids,  Protein.Names, Genes,  First.Protein.Description, Proteotypic,  Stripped.Sequence, Modified.Sequence,  Precursor.Charge) |> 
#     dplyr::mutate(peptide_id = Precursor.Id) |> 
#     tibble::column_to_rownames("Precursor.Id")|> 
#     tibble::tibble()


# data.proteomics.glass_od <- tmp |> 
#   dplyr::select(Precursor.Id, contains("/scratch/")) |> 
#   tibble::column_to_rownames("Precursor.Id") |> 
#   dplyr::rename_with(~ gsub("/scratch/DIANN_A314/WU300344/","",.x)) |> 
#   dplyr::rename_with(~ gsub("\\.d$","",.x)) |> 
#   dplyr::mutate(
#     `20240215_003_S640489_GLODprot_pool` = NULL,
#     `20240215_036_S640489_GLODprot_pool` = NULL,
#     `20240215_069_S640489_GLODprot_pool` = NULL,
#     `20240215_102__S640489_GLODprot_pool` = NULL,
#     `20240215_135_S640489_GLODprot_pool` = NULL
#   ) |> 
#   dplyr::rename_with(~ gsub("^.+GLODprot_","GLODprot_",.x)) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(ncol(.) == CONST_N_GLASS_OD_ALL_PROTEOMICS)
#     assertthat::assert_that(all(colnames(.) %in% glass_od.metadata.proteomics$proteomics_id))
#     return(.)
#   })() |> 
#   dplyr::select(glass_od.metadata.proteomics$proteomics_id) |># same order
#   log2()
# 
# 
# rm(tmp)


# load proteomics per-protein ----

# 
# tmp <- read.csv('data/GLASS_OD/Protein - Tobias Weiss/WU300344_report.pg_matrix.tsv',
#                 check.names=F,
#                 sep="\t") |> 
#   dplyr::filter(!grepl("Y-FGCZCont", Protein.Ids)) |> 
#   dplyr::filter(!is.na(Genes)) |> 
#   dplyr::filter(Genes != "")
# 
# 
# metadata.proteomics.glass_od <- tmp |> 
#   dplyr::select(Protein.Group,
#                 Protein.Ids,   
#                 Protein.Names,  
#                 Genes, 
#                 First.Protein.Description ) |> 
#   dplyr::mutate(protein_id = Genes) |> 
#   tibble::column_to_rownames("Genes")|> 
#   tibble::tibble()
# 
# 
# 
# data.proteomics.glass_od <- tmp |> 
#   dplyr::select(Genes, contains("scratch")) |> 
#   tibble::column_to_rownames("Genes") |> 
#   dplyr::rename_with(~ gsub("/scratch/DIANN_A314/WU300344/","",.x)) |> 
#   dplyr::rename_with(~ gsub("\\.d$","",.x)) |> 
#   dplyr::mutate(
#     `20240215_003_S640489_GLODprot_pool` = NULL,
#     `20240215_036_S640489_GLODprot_pool` = NULL,
#     `20240215_069_S640489_GLODprot_pool` = NULL,
#     `20240215_102__S640489_GLODprot_pool` = NULL,
#     `20240215_135_S640489_GLODprot_pool` = NULL
#   ) |> 
#   dplyr::rename_with(~ gsub("^.+GLODprot_","GLODprot_",.x)) |> 
#   (function(.) {
#     print(dim(.))
#     assertthat::assert_that(ncol(.) == CONST_N_GLASS_OD_ALL_PROTEOMICS)
#     assertthat::assert_that(all(colnames(.) %in% glass_od.metadata.proteomics$proteomics_id))
#     return(.)
#   })() |> 
#   dplyr::select(glass_od.metadata.proteomics$proteomics_id) |>  # same order
#   log2()
# 
# 
# rm(tmp)
# 
# 

# load proteomics per-protein NORMALISED ----



tmp <- readxl::read_xlsx('data/GLASS_OD/Protein - Tobias Weiss/norm_prot_pgmatrix.xlsx') |> 
  dplyr::filter(grepl("_HUMAN$", Protein.Names))


tmp.2 <- read.csv('data/GLASS_OD/Protein - Tobias Weiss/WU300344_report.pg_matrix.tsv',
                   check.names=F,
                   sep="\t") |> 
  dplyr::filter(!is.na(Genes)) |> 
  dplyr::filter(Genes != "") |> 
  dplyr::filter(!grepl("Y-FGCZCont", Protein.Ids)) |> 
  dplyr::select(#Protein.Group, Protein.Ids,
                Protein.Names, Genes, First.Protein.Description)



metadata.proteomics.glass_od <- tmp |> 
  dplyr::select(Protein.Group,
                Protein.Ids,   
                Protein.Names,  
                isotopeLabel) |> 
  as.data.frame() |> 
  dplyr::left_join(tmp.2, by=c('Protein.Names'='Protein.Names')) |> 
  dplyr::filter(!is.na(Genes))




data.proteomics.glass_od <- tmp |> 
  dplyr::select(Protein.Names, contains("_GLODprot_")) |> 
  dplyr::filter(Protein.Names %in% metadata.proteomics.glass_od$Protein.Names) |> 
  dplyr::left_join(
    metadata.proteomics.glass_od |> dplyr::select(Protein.Names, Genes), by=c('Protein.Names'='Protein.Names')
  ) |> 
  dplyr::mutate(Protein.Names = NULL) |> 
  tibble::column_to_rownames("Genes") |> 
  dplyr::rename_with(~ gsub("\\.d$","",.x)) |> 
  dplyr::mutate(
    `20240215_003_S640489_GLODprot_pool` = NULL,
    `20240215_036_S640489_GLODprot_pool` = NULL,
    `20240215_069_S640489_GLODprot_pool` = NULL,
    `20240215_102__S640489_GLODprot_pool` = NULL,
    `20240215_135_S640489_GLODprot_pool` = NULL
  ) |> 
  dplyr::rename_with(~ gsub("^.+GLODprot_","GLODprot_",.x)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(ncol(.) == CONST_N_GLASS_OD_ALL_PROTEOMICS)
    assertthat::assert_that(all(colnames(.) %in% glass_od.metadata.proteomics$proteomics_id))
    return(.)
  })() |> 
  dplyr::select(glass_od.metadata.proteomics$proteomics_id) # same order


rm(tmp, tmp.2)




# add cell cycling annotation ----


clng <- c("AAAS","AATF","ABCB1","ABL1","ABRAXAS1","ABRAXAS2","ACTB","ACTL6A","ACTL6B","ACTR2","ACTR3","ACTR5",
          "ACTR8","ACVR1","ACVR1B","ADAM17","ADAMTS1","ADARB1","ADCYAP1","AFAP1L2","AGO4","AHCTF1","AHR","AICDA",
          "AIF1","AJUBA","AKAP8","AKAP8L","AKT1","AKT2","ALKBH4","ALMS1","ALOX15B","AMBRA1","ANAPC1","ANAPC10",
          "ANAPC11","ANAPC13","ANAPC15","ANAPC16","ANAPC2","ANAPC4","ANAPC5","ANAPC7","ANGEL2","ANK3","ANKFN1",
          "ANKK1","ANKLE1","ANKLE2","ANKRD17","ANKRD31","ANKRD53","ANLN","ANXA1","ANXA11","APBB1","APBB2","APC",
          "APEX2","APP","APPL1","APPL2","ARF1","ARF6","ARHGEF10","ARHGEF2","ARID1A","ARID1B","ARID2","ARL2","ARL3",
          "ARL8A","ARL8B","ARPP19","ASAH2","ASCL1","ASNS","ASPM","ASZ1","ATAD5","ATF2","ATF5","ATM","ATP2B4","ATR","ATRIP","ATRX","AUNIP","AURKA","AURKAIP1","AURKB","AURKC","AVEN","AVPI1","AXIN2","AZI2","BABAM1","BABAM2","BAG6","BAK1","BANF1","BANP","BAP1","BARD1","BAX","BAZ1B","BBS4","BCAT1","BCCIP","BCL2","BCL2L1","BCL2L11","BCL6","BCL7A","BCL7B","BCL7C","BCR","BECN1","BEX2","BEX4","BID","BIN1","BIN3","BIRC2","BIRC3","BIRC5","BIRC6","BIRC7","BLCAP","BLM","BMAL1","BMP2","BMP4","BMP7","BOD1","BOD1L1","BOD1L2","BOLL","BOP1","BORA","BRCA1","BRCA2","BRCC3","BRD4","BRD7","BRD8","BRDT","BRINP1","BRINP2","BRINP3","BRIP1","BRME1","BROX","BRSK1","BRSK2","BTBD18","BTC","BTG1","BTG2","BTG3","BTG4","BTN2A2","BTRC","BUB1","BUB1B","BUB3","C10orf90","C11orf80","C14orf39","C1orf146","C2CD3",
          "C6orf89","C9orf78","CABLES1","CABLES2","CACNB4","CACUL1","CALM1",
          "CALM2","CALM3","CALR","CAMK1","CAMK2A","CAMSAP3","CAPN3","CASP2","CASP3","CATSPERZ","CAV2",
          "CCAR1","CCAR2","CCDC102B","CCDC124","CCDC42","CCDC57","CCDC61","CCDC69","CCDC8","CCL2","CCNA1","CCNA2","CCNB1","CCNB1IP1","CCNB2","CCNB3","CCNC","CCND1","CCND2","CCND3","CCNDBP1","CCNE1","CCNE2","CCNF","CCNG1","CCNG2","CCNH","CCNI","CCNI2","CCNJ","CCNJL","CCNK","CCNL1","CCNL2","CCNO","CCNP","CCNQ","CCNT1","CCNT2","CCNY","CCNYL1","CCP110","CCPG1","CCSAP","CD28","CD2AP","CDC123",
          "CDC14A","CDC14B","CDC14C","CDC16","CDC20","CDC23","CDC25A","CDC25B","CDC25C","CDC26","CDC27","CDC34","CDC37","CDC42","CDC45","CDC5L","CDC6","CDC7","CDC73","CDCA2","CDCA3","CDCA5","CDCA8","CDK1","CDK10","CDK11A","CDK11B","CDK12","CDK13","CDK14","CDK15","CDK16","CDK17","CDK18","CDK19","CDK2","CDK20","CDK2AP1","CDK2AP2","CDK3","CDK4","CDK5","CDK5R1","CDK5R2","CDK5RAP1","CDK5RAP2","CDK5RAP3","CDK6","CDK7","CDK8","CDK9","CDKL1","CDKL2","CDKL3","CDKL4","CDKL5","CDKN1A","CDKN1B","CDKN1C","CDKN2A","CDKN2B","CDKN2C","CDKN2D","CDKN3","CDT1","CEBPA","CECR2","CENATAC","CENPA","CENPC","CENPE","CENPF","CENPH","CENPI","CENPJ","CENPK","CENPL","CENPM","CENPN","CENPO","CENPP","CENPQ","CENPS",
          "CENPT","CENPU","CENPV","CENPW","CENPX","CEP120","CEP126","CEP131","CEP135","CEP152","CEP164","CEP192","CEP250","CEP295","CEP295NL","CEP44","CEP55","CEP63","CEP68","CEP72","CEP76","CEP85","CEP97","CETN1","CETN2","CETN3","CFL1","CGREF1","CGRRF1","CHAF1A","CHAF1B","CHAMP1","CHD3","CHEK1","CHEK2","CHFR","CHMP1A","CHMP1B","CHMP2A","CHMP2B","CHMP3","CHMP4A","CHMP4B","CHMP4C","CHMP5","CHMP6","CHMP7","CHORDC1","CHTF18","CHTF8","CIAO1","CIAO2A","CIAO2B","CIB1","CINP","CIT","CITED2","CKAP2","CKAP5","CKS1B","CKS2","CLASP1","CLASP2","CLIP1","CLOCK","CLSPN","CLTA","CLTC","CLTCL1","CNPPD1","CNTD1","CNTLN","CNTRL","CNTROB","COPS5","CPSF3","CRADD","CREB3","CREBL2","CRLF3","CRNN","CROCC","CRY1","CSAG1","CSNK1A1","CSNK1D","CSNK2A1","CSNK2A2","CSNK2A3","CSPP1","CTBP1","CTC1","CTCF",
          "CTCFL","CTDNEP1","CTDP1","CTDSP1","CTDSP2",
          "CTDSPL","CTNNB1","CUL1","CUL2","CUL3","CUL4A","CUL4B","CUL5","CUL7","CUL9","CUZD1","CXCR5","CYLD","CYP1A1","CYP26B1","CYP27B1","DAB2IP","DACH1","DACT1","DAPK3",
          "DAZL","DBF4","DBF4B","DCAF13","DCDC1","DCLRE1A","DCTN1","DCTN2","DCTN3","DCTN6","DCUN1D3","DDB1","DDIAS","DDIT3","DDR2","DDRGK1","DDX11","DDX39B","DDX3X","DDX4","DEUP1","DGKZ","DIAPH3","DIRAS3","DIS3L2","DLG1","DLGAP5","DMAP1","DMC1","DMRT1","DMRTC2","DMTF1","DNA2","DNM2","DNMT3L","DONSON","DOT1L","DPF1","DPF2","DPF3","DR1","DRD2","DRD3","DRG1","DSCC1","DSN1","DTL","DTX3L","DUSP1","DUSP13B","DUSP3","DUX4","DYNC1H1","DYNC1LI1","DYNLT1","DYNLT3","DYRK3","E2F1","E2F2","E2F3","E2F4","E2F6","E2F7","E2F8","E4F1","ECD","ECRG4","ECT2","EDN1","EDN3","EDNRA","EFHC1","EGF","EGFR","EHMT2","EID1","EIF2AK4","EIF4E","EIF4EBP1","EIF4G1","EIF4G2","EME1","EME2","EML1","EML3","EML4","ENKD1","ENSA","ENTR1","EP300",
          "EP400","EPB41","EPB41L2","EPC1","EPC2","EPGN","EPM2A","EPS8","ERCC2","ERCC3","ERCC4","ERCC6","ERCC6L","EREG","ERH","ESCO1","ESCO2","ESPL1","ESX1","ETAA1","ETS1","EVI2B","EVI5","EXD1","EXO1","EXOC1","EXOC2","EXOC3","EXOC4","EXOC5","EXOC6","EXOC6B","EXOC7","EXOC8","EZH2","EZR","FAM107A","FAM110A","FAM32A","FAM83D","FAM9A","FAM9B",
          "FAM9C","FANCA","FANCD2","FANCI","FANCM","FAP","FBXL15","FBXL7","FBXO31","FBXO4","FBXO43","FBXO5","FBXO6","FBXO7","FBXW11","FBXW5","FBXW7","FEM1B","FEN1","FES","FGF10","FGF2","FGF8","FGFR1","FGFR2","FHL1","FIGN","FIGNL1","FIRRM","FKBP6","FLCN","FLNA","FMN2","FOSL1","FOXA1","FOXC1","FOXE3","FOXG1","FOXJ2","FOXJ3","FOXM1","FOXN3","FOXO4","FSD1","FZD3","FZD9","FZR1","GADD45A","GADD45B","GADD45G","GADD45GIP1","GAK","GAS1","GAS2","GATA3","GATA6","GBF1","GEM","GEN1","GFI1B","GIGYF2","GINS1","GINS3","GIPC1","GIT1","GJA1","GJC2","GLI1","GML","GMNC","GMNN","GNAI1","GNAI2","GNAI3","GOLGA2","GOLGA8B","GOLGA8F","GOLGA8S","GPER1","GPNMB","GPR132","GPR15LG","GPR3","GPSM1","GPSM2","GRK5","GSPT1","GSPT2","GTF2B","GTF2H1","GTPBP4","GTSE1","H1-8","H2AX","HACE1","HASPIN","HAUS1","HAUS2","HAUS3",
          "HAUS4","HAUS5","HAUS6","HAUS7","HAUS8","HBP1","HCFC1","HDAC3",
          "HDAC8","HECA","HECW2","HELLS","HEPACAM","HEPACAM2","HERC5","HES1","HEXIM1","HEXIM2","HFM1","HGF","HHEX","HINFP","HIPK2","HJURP",
          "HLA-G","HMCN1","HMG20B","HMGA2","HNRNPU","HORMAD1","HORMAD2","HOXA13","HOXC9","HOXD10","HPGD","HRAS","HSF1","HSF2BP","HSP90AB1","HSPA1A","HSPA1B","HSPA2","HTRA2","HTT","HUS1","HUS1B","HYAL1","ID2","ID4","IER3","IFFO1","IFNW1","IGF1","IGF1R","IGF2","IHO1","IK","IKZF1","IL10","IL1A","IL1B","INCA1",
          "INCENP","ING1","ING2","ING3","ING4","ING5","INHA","INHBA","INIP","INO80","INO80B","INO80C","INO80D","INO80E","INPPL1","INS","INSC","INSM1","INSM2","INSR","INTS13","INTS3","INTS7","IPO5","IPO7","IQGAP1","IQGAP2","IQGAP3","IRF1","IST1","ITGB1","ITGB3BP","JADE1","JADE2","JADE3","JTB","JUN","JUNB","JUND","KANK2","KASH5","KAT14","KAT2A","KAT2B","KAT5","KAT7","KATNA1","KATNB1","KCNA5","KCNH5","KCTD11","KCTD19","KDM8","KHDRBS1","KIAA0753","KIAA1614","KIF11","KIF13A","KIF14","KIF15","KIF18A","KIF18B","KIF20A","KIF20B","KIF22","KIF23","KIF25","KIF2A","KIF2B","KIF2C","KIF3A",
          "KIF3B","KIF4A","KIF4B","KIFC1","KIFC2","KIZ","KLF4","KLHDC3","KLHDC8B","KLHL13","KLHL18","KLHL21","KLHL22","KLHL42","KLHL9","KLK10","KLLN","KMT2E","KMT5A","KNL1","KNSTRN","KNTC1","KPNB1","KRT18","L3MBTL1","LATS1","LATS2","LCMT1","LEF1","LEP","LFNG","LGMN","LIF","LIG1","LIG3","LIG4","LIMK2","LIN54","LIN9","LIPA","LLGL1","LLGL2","LMLN","LMNA","LPIN1","LRP5","LRP6","LRRCC1","LSM10","LSM11","LSM14A","LZTS1","LZTS2","M1AP","MACROH2A1","MAD1L1","MAD2L1","MAD2L1BP","MAD2L2","MADD","MAEA","MAEL","MAGEA4","MAJIN","MAP10","MAP1S","MAP2K6","MAP3K11","MAP3K20","MAP3K7","MAP3K8","MAP4","MAP9","MAPK1","MAPK12","MAPK13","MAPK14","MAPK15","MAPK3","MAPK4","MAPK6","MAPK7",
          "MAPRE1","MAPRE2","MAPRE3","MARF1","MARK4","MARVELD1","MASTL","MAU2","MBIP","MBLAC1","MBTD1","MBTPS1","MBTPS2","MCIDAS","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MCM8","MCMBP","MCMDC2","MCPH1","MCRS1","MCTS1","MDC1","MDM1","MDM2","MDM4","MEAF6","MECOM","MECP2","MED1","MEI1","MEI4","MEIKIN","MEIOB","MEIOC","MEIOSIN","MEIS2","MELK","MEN1","MEPCE","METTL13","METTL3","MICAL3","MIIP","MIS12","MIS18A","MIS18BP1","MISP","MITD1","MKI67","MLF1","MLH1","MLH3","MMS19","MN1","MNAT1","MND1","MNS1","MNT","MOK","MORF4L1","MORF4L2","MOS","MOV10L1","MPLKIP","MRE11","MRFAP1L2","MRGBP","MRGPRX2","MRNIP","MRPL41","MSH2","MSH4","MSH5","MSH6","MSX1","MSX2","MTA3","MTBP","MTCL1","MUC1","MUS81","MX2","MYB","MYBBP1A","MYBL1","MYBL2","MYC","MYH10","MYH9","MYO16","MYO19","MYOCD","MYOG","MZT1","NAA10","NAA50","NAA60","NABP1","NABP2","NAE1","NANOGP8","NANOS2","NANOS3","NASP","NAT10","NBN","NCAPD2","NCAPD3","NCAPG","NCAPG2","NCAPH","NCAPH2","NCOR1",
          "NDC1","NDC80","NDE1","NDEL1","NDP","NEDD1","NEDD9","NEK1","NEK10","NEK11","NEK2","NEK3","NEK4","NEK6","NEK7","NEK9","NES","NEUROG1","NF2","NFIA","NFIB","NFRKB","NHERF1","NIN","NIPBL","NKX3-1","NLE1","NLRP2B","NLRP5","NME6","NOLC1","None","None","None","None","None","None","None","None","None","None","NOP53","NOX5","NPAT","NPM1","NPM2",
          "NPPC","NPR2","NR2E1","NR2F2","NR3C1","NR4A1","NRDE2","NSFL1C","NSL1","NSMCE2","NSUN2","NTMT1",
          "NUBP1","NUDC","NUDT15","NUDT16","NUDT6","NUF2","NUGGC","NUMA1","NUP214","NUP37","NUP43","NUP62","NUPR1","NUPR2","NUSAP1","OBSL1","ODF2","OIP5","OOEP","OPN1LW","OPN1MW","OR1A2","OR2A4","ORC1","ORC4","OSGIN2","OVOL1","OVOL2","P3H4","PABIR1","PAF1","PAFAH1B1","PAGR1","PAK4","PARD3","PARD3B","PARD6A","PARD6B","PARD6G","PARP3","PARP9","PAX6","PAXIP1","PBK","PBRM1","PBX1","PCID2","PCLAF","PCM1","PCNA",
          "PCNP","PCNT","PDCD2L","PDCD6IP","PDE3A","PDE4DIP","PDGFB","PDGFRB","PDS5A","PDS5B","PDXP","PELO","PER2","PES1","PHACTR4","PHB2","PHF10","PHF13","PHF8","PHGDH","PHIP","PHOX2B","PIAS1","PIBF1","PIDD1","PIK3C3","PIK3R4","PIM1","PIM2","PIM3","PIMREG","PIN1","PINX1","PIWIL1","PIWIL2","PIWIL3","PIWIL4","PKD1","PKD2","PKHD1","PKIA","PKMYT1","PKN2","PKP4","PLAGL1","PLCB1","PLCG2","PLD6","PLEC","PLK1","PLK2","PLK3","PLK4","PLK5","PLRG1","PLSCR1","PMF1-BGLAP","PMF1","PML","PNPT1","POC1B","POC5","POGZ","POLA1","POLDIP2","POLE","POU4F1","PPM1A","PPM1D",
          "PPM1G","PPME1","PPP1CA","PPP1CB","PPP1CC","PPP1R10","PPP1R12A","PPP1R13B","PPP1R15A","PPP1R1C","PPP1R35","PPP1R9B","PPP2CA","PPP2CB","PPP2R1A","PPP2R2D","PPP2R3B","PPP2R5B","PPP2R5C","PPP3CA","PPP5C","PPP6C","PRAP1","PRC1","PRCC","PRDM11","PRDM5","PRDM7","PRDM9","PRICKLE1","PRKACA","PRKCA","PRKCB","PRKCD","PRKCE","PRKDC","PRMT2","PRMT5","PRNP","PROX1","PRPF19","PRPF40A","PRR11","PRR19","PRR5","PSMA8","PSMC3IP","PSMD10","PSMD13","PSME1","PSME2","PSME3","PSMG2","PSRC1","PTCH1","PTEN","PTGS2","PTK6","PTP4A1","PTPA","PTPN11","PTPN3","PTPN6","PTPRC","PTPRK","PTTG1","PTTG2","PUM1","PUM2","PYHIN1","RAB11A","RAB11FIP3","RAB11FIP4","RAB35","RAB6C","RABGAP1","RACGAP1","RACK1","RAD1","RAD17","RAD18","RAD21","RAD21L1","RAD23A","RAD50","RAD51","RAD51AP1","RAD51B","RAD51C","RAD51D",
          "RAD54B","RAD54L","RAD9A","RAD9B","RAE1","RALA","RALB","RAN","RANBP1","RARA","RASA1","RASSF1","RASSF2","RASSF4","RB1","RB1CC1","RBBP4","RBBP8","RBL1","RBL2","RBM14","RBM38","RBM46","RBM7","RCBTB1","RCC1","RCC2","RDX","REC114","REC8","RECQL5","REEP3","REEP4","RFPL1","RFWD3","RGCC","RGS14","RGS2","RHEB","RHNO1","RHOA","RHOB","RHOC","RHOU","RIF1","RINT1","RIOK2","RIOK3","RIPOR2","RMDN1","RMI1","RMI2","RNASEH2B",
          "RNF103-CHMP3","RNF112","RNF167","RNF2","RNF212","RNF212B","RNF4","RNF8","ROCK1","ROCK2","ROPN1B","RPA1","RPA2","RPA3","RPA4","RPL10L","RPL23","RPL24","RPL26","RPRD1B","RPRM","RPS15A","RPS27L","RPS3","RPS6KA1","RPS6KA2","RPS6KA3","RPS6KB1","RPTOR","RRM1","RRM2","RRM2B","RRP8","RRS1","RSPH1","RTEL1","RTF2","RTKN","RTTN","RUNX3","RUVBL1","RUVBL2","RXFP3","SAPCD2","SASS6","SBDS","SCAND3","SCRIB","SDCBP","SDCCAG8","SDE2","SEH1L","SENP5","SENP6","SEPTIN1","SEPTIN10","SEPTIN11","SEPTIN12","SEPTIN14",
          "SEPTIN2","SEPTIN3","SEPTIN4","SEPTIN5","SEPTIN6","SEPTIN7","SEPTIN8","SEPTIN9","SERTAD1","SETD2","SETDB2","SETMAR","SFN","SFRP1","SGF29","SGO1","SGO2","SGSM3","SH2B1","SH3GLB1","SHB","SHCBP1L","SHOC1","SIAH1","SIAH2","SIK1","SIN3A","SIPA1","SIRT1","SIRT2","SIRT7","SIX3","SKA1","SKA2","SKA3","SKIL","SKP2","SLC16A1","SLC25A31","SLC25A5","SLC26A8","SLC39A5","SLC6A4","SLF1","SLF2","SLFN11","SLX4",
          "SMARCA2","SMARCA4","SMARCA5","SMARCAD1","SMARCB1","SMARCC1","SMARCC2","SMARCD1","SMARCD2","SMARCD3","SMARCE1","SMC1A","SMC1B","SMC2","SMC3","SMC4","SMC5","SMC6","SMIM22","SMOC2","SMPD3","SND1","SNRK","SNX18","SNX33","SNX9","SOGA1","SON","SOX15","SOX2","SOX9","SPAG5","SPAG8","SPAST","SPATA17","SPATA22","SPC24","SPC25","SPDL1","SPDYA","SPDYC","SPDYE1","SPDYE10","SPDYE11","SPDYE12","SPDYE13","SPDYE14","SPDYE15","SPDYE16","SPDYE17","SPDYE18","SPDYE2","SPDYE21","SPDYE2B","SPDYE3","SPDYE4","SPDYE5","SPDYE6","SPDYE8","SPDYE9","SPECC1L-ADORA2A","SPECC1L","SPHK1","SPICE1","SPIN1","SPIN2A","SPIN2B","SPIRE1","SPIRE2","SPO11","SPOUT1","SPRY1","SPRY2","SPTBN1","SRA1","SRC","SRPK1","SRPK2","SSNA1","SSTR5","SSX2IP","STAG1","STAG2","STAG3","STAMBP","STARD9","STAT3","STAT5B","STEAP3","STIL","STK10","STK11","STK33","STMN1","STOX1","STRA8","STRADA","STRADB","STXBP4","SUGT1","SUN1","SUN2","SUSD2","SUV39H1","SUV39H2","SVIL","SYCE1","SYCE1L","SYCE2","SYCE3","SYCP1","SYCP2","SYCP2L","SYCP3","SYF2","TACC1","TACC2","TACC3","TADA2A","TADA3","TAF1","TAF10","TAF1L","TAF2","TAF6","TAL1","TAOK1","TAOK2","TAOK3","TARDBP","TAS1R2","TAS2R13","TBCD","TBCE","TBRG1","TBRG4","TBX1","TBX2","TBX20","TBX3","TCF3","TCIM","TDRD1","TDRD12","TDRD9","TDRKH","TELO2","TENT4A","TENT4B","TENT5B","TERB1","TERB2","TERF1","TERF2","TERT","TESMIN","TET2","TEX11","TEX12","TEX14","TEX15","TEX19","TFAP4","TFDP1","TFDP2",
          "TFDP3","TFPT","TGFA","TGFB1","TGFB2","TGFBR1","TGM1","THAP1","THAP5","THAP9","THOC1","THOC5","TICRR","TIMELESS","TIPIN","TIPRL","TK1","TLE6","TLK1","TLK2","TM4SF5","TMEM14B","TMEM250","TMEM67","TMEM8B","TMOD3","TMPRSS11A","TNF","TNKS","TOM1L1","TOM1L2","TOP1","TOP2A","TOP2B","TOP3A","TOP3B","TOPBP1","TP53","TP53BP1","TP53BP2","TP53I13","TP53INP1","TP63","TP73","TPD52L1","TPPP","TPR","TPRA1","TPX2","TRAPPC12","TREX1","TRIAP1","TRIM21","TRIM32","TRIM35","TRIM36","TRIM37","TRIM39","TRIM71","TRIM75","TRIOBP","TRIP13","TRNP1","TRRAP","TSC1","TSC2","TSC22D2","TSG101","TSPYL2","TSSK4","TTBK1","TTC19","TTC28","TTI1","TTI2","TTK","TTL","TTLL12","TTN","TTYH1","TUBA1A","TUBA1B","TUBA1C","TUBA3C","TUBA3D","TUBA3E","TUBA4A","TUBA4B","TUBA8","TUBAL3","TUBB","TUBB1","TUBB2A","TUBB2B","TUBB3","TUBB4A","TUBB4B","TUBB6","TUBB8","TUBB8B","TUBD1","TUBE1",
          "TUBG1","TUBG2","TUBGCP2","TUBGCP3","TUBGCP4","TUBGCP5","TUBGCP6","TUSC2","TXLNG","TXNIP","TXNL4A","TXNL4B","UBA3","UBB","UBD","UBE2A","UBE2B","UBE2C","UBE2E2","UBE2I","UBE2L3","UBE2S","UBR2","UBXN2B","UCHL5","UFL1","UHMK1","UHRF1","UHRF2","UIMC1","UNC119","UPF1","URGCP","USH1C","USP16","USP17L2","USP19","USP2","USP22","USP26","USP28","USP29","USP3","USP33","USP37","USP39","USP44","USP47","USP51","USP8","USP9X","UTP14C","UVRAG","UXT","VASH1","VCP","VPS4A","VPS4B","VPS72","VRK1","WAC","WAPL","WASHC5","WASL","WBP2NL","WDHD1","WDR12","WDR5","WDR6","WDR62","WDR76","WEE1","WEE2","WIZ","WNK1","WNT10B","WNT4","WNT5A","WRAP73","WRN","WTAP","XIAP","XPC","XPO1","XRCC2","XRCC3","YEATS2","YEATS4","YTHDC2","YTHDF2","YWHAE","YY1","YY1AP1","ZBED3","ZBTB17","ZBTB49","ZC3H12D","ZC3HC1","ZCWPW1","ZFP36L1","ZFP36L2","ZFP42","ZFYVE19","ZFYVE26","ZMPSTE24","ZMYND11","ZNF16","ZNF207","ZNF268","ZNF318","ZNF324","ZNF503","ZNF541","ZNF655","ZNF703","ZNF830","ZNRD2","ZNRF4","ZPR1","ZSCAN21","ZW10","ZWILCH","ZWINT","ZZZ3",
          "KI67","MKI67")



metadata.proteomics.glass_od <- metadata.proteomics.glass_od |> 
  dplyr::mutate(cellcycling_go_ = Genes %in% clng)


rm(clng)


# GO:996023 Collagen containing ECM ----




# add fibronectin / FN1 ----


metadata.proteomics.glass_od <- metadata.proteomics.glass_od |> 
  dplyr::mutate(is_fibronectin_fn1 = Genes %in% c("FN1", "FINC", "FIBRONETIN", "FIBRONETIN1", "FN", "CIG", "GFND2", "LETS", "MSF"))


#metadata.proteomics.glass_od |>
#  dplyr::filter(is_fibronectin_fn1) 



# DPA: CGC ----

fn <- "cache/analysis_differential_proteomics__GLASS-OD__stats.cgc.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(protein_id, logFC, t, P.Value, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DPA__GLASS_OD__CGC__", .x), .cols=!matches("^protein_id$", perl = T))
  
  
  metadata.proteomics.glass_od <- metadata.proteomics.glass_od |> 
    dplyr::left_join(tmp, by=c('Genes'='protein_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DPA GLASS-OD x CGC is missing")
}

rm(fn)



# DPA: prim - rec (naive) ----

fn <- "cache/analysis_differential_proteomics__GLASS-OD__stats.prim-rec.naive.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(protein_id, logFC, t, P.Value, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DPA__GLASS_OD__prim-rec__naive__", .x), .cols=!matches("^protein_id$", perl = T))
  
  
  metadata.proteomics.glass_od <- metadata.proteomics.glass_od |> 
    dplyr::left_join(tmp, by=c('Genes'='protein_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DPA GLASS-OD x prim-rec (naive) is missing")
}

rm(fn)



# DPA: prim - rec (pat corrected) ----

fn <- "cache/analysis_differential_proteomics__GLASS-OD__stats.prim-rec.pat-corrected.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(protein_id, logFC, t, P.Value, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DPA__GLASS_OD__prim-rec__pat_corrected__", .x), .cols=!matches("^protein_id$", perl = T))
  
  
  metadata.proteomics.glass_od <- metadata.proteomics.glass_od |> 
    dplyr::left_join(tmp, by=c('Genes'='protein_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DPA GLASS-OD x prim-rec (naive) is missing")
}

rm(fn)




# DPA: grade (naive) ----

fn <- "cache/analysis_differential_proteomics__GLASS-OD__stats.grade.naive.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(protein_id, logFC, t, P.Value, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DPA__GLASS_OD__grade__naive__", .x), .cols=!matches("^protein_id$", perl = T))
  
  
  metadata.proteomics.glass_od <- metadata.proteomics.glass_od |> 
    dplyr::left_join(tmp, by=c('Genes'='protein_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DPA GLASS-OD x prim-rec (naive) is missing")
}

rm(fn)



# DPA: grade (pat corrected) ----

fn <- "cache/analysis_differential_proteomics__GLASS-OD__stats.grade.pat-corrected.Rds"
if(file.exists(fn)) {
  
  tmp <- readRDS(fn) |> 
    dplyr::select(protein_id, logFC, t, P.Value, adj.P.Val) |> 
    dplyr::rename_with(~paste0("DPA__GLASS_OD__grade__pat_corrected__", .x), .cols=!matches("^protein_id$", perl = T))
  
  
  metadata.proteomics.glass_od <- metadata.proteomics.glass_od |> 
    dplyr::left_join(tmp, by=c('Genes'='protein_id'), suffix=c('','') )
  
  rm(tmp)
  
} else {
  warning("DPA GLASS-OD x prim-rec (naive) is missing")
}

rm(fn)





