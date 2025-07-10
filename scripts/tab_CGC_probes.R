#!/usr/bin/env R

fit = readRDS("cache/LGC_predictor_probe_based_lm.Rds")

coef_best <- coef(b, s = "lambda.min")
coef(b, labeldat = b$lambda)

model <- rbind(
fit$beta |> 
  as.matrix() |> 
  as.data.frame() |> 
  dplyr::rename(value = s0) |> 
  tibble::rownames_to_column('coefficient') |> 
  dplyr::filter(value != 0) |> 
  dplyr::arrange(value) |> 
  dplyr::mutate( type = "weight")
,
data.frame(
  coefficient = c(''),
  value = c(fit$a0),
  type = "intercept"
)
)




model <- model |> 
  dplyr::left_join( 
    data.mvalues.probes |> 
      dplyr::select(probe_id,
                    
                    GCIMP_IDHmut_probe, 
                    
                    UCSC_RefGene_Name,  GencodeCompV12_NAME, CHR_hg38, Start_hg38,  End_hg38, Strand_hg38, gc_sequence_context_2_new
                    
      ) 

    ,   by=c('coefficient'='probe_id'), suffix = c('','') )



model$UCSC_RefGene_Name <- sapply(model$UCSC_RefGene_Name, function(x) {
  x |> 
    strsplit(split = ";") |>          # base strsplit returns a list
    unlist() |> 
    unique() |> 
    sort() |>                         # optional, remove if you want original order
    paste(collapse = ";")
})


model$GencodeCompV12_NAME <- sapply(model$GencodeCompV12_NAME, function(x) {
  x |> 
    strsplit(split = ";") |>          # base strsplit returns a list
    unlist() |> 
    unique() |> 
    sort() |>                         # optional, remove if you want original order
    paste(collapse = ";")
})


write.table(model, file=paste0("output/tables/tab_CGC_probes.txt"), sep="\t", quote = F, row.names = F)




