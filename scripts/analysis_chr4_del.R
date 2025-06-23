#!/usr/bin/env R

# 1. load ----


library(ggplot2)

if(!exists('glass_od.metadata.array_samples')) {
  source('scripts/load_GLASS-OD_metadata.R')
}


# 2. select samples ----


metadata <- glass_od.metadata.array_samples |> 
  filter_GLASS_OD_idats(CONST_N_GLASS_OD_INCLUDED_SAMPLES)


# 3. calc weighted CNV ----

out <- parse_mnp_CNVP_chr4_mean(metadata$array_mnp_CNVP_v12.8_v5.2_CNVP_segments)


write.table(out, file="cache/analysis_chr4_del__weighted_cnv_offsets.Rds")



# 4. define cut-off boundaries ----


marks = data.frame() |> 
  rbind(
    data.frame(resection_id = "0123-R1", chr4hd = "Yes"),
    data.frame(resection_id = "0055-R1", chr4hd = "Yes"),
    data.frame(resection_id = "0112-R1", chr4hd = "Yes"), # obviously present and obviously subclonal
    data.frame(resection_id = "0018-R2", chr4hd = "No"),  # partial arm
    data.frame(resection_id = "0042-R3", chr4hd = "Yes"),
    data.frame(resection_id = "0024-R2", chr4hd = "No"),
    data.frame(resection_id = "0109-R1", chr4hd = "Inconclusive"), # borderline
    data.frame(resection_id = "0010-R2", chr4hd = "Yes"),
    data.frame(resection_id = "0053-R1", chr4hd = "Inconclusive"), # poor quality sample
    data.frame(resection_id = "0065-R2", chr4hd = "Yes"),
    data.frame(resection_id = "0043-R2", chr4hd = "No"), # borderline to unvisible
    data.frame(resection_id = "0022-R1", chr4hd = "Yes"), # obviously present and obviously subclonal
    data.frame(resection_id = "0075-R2", chr4hd = "Yes"),
    data.frame(resection_id = "0043-R1", chr4hd = "Yes"),
    data.frame(resection_id = "0052-R2", chr4hd = "No"), 
    data.frame(resection_id = "0077-R2", chr4hd = "Yes"), # marginally present and obviously subclonal
    data.frame(resection_id = "0037-R2", chr4hd = "Inconclusive"),   # Not clear enough, more no than yes
    data.frame(resection_id = "0090-R3", chr4hd = "No"), # clearly no
    data.frame(resection_id = "0040-R1", chr4hd = "Yes"), # seems present and obviously subclonal
    data.frame(resection_id = "0002-R1", chr4hd = "Yes"), # obviously present and obviously subclonal
    data.frame(resection_id = "0006-R2", chr4hd = "Inconclusive"), # borderline
    data.frame(resection_id = "0018-R3", chr4hd = "No"), # low purity
    data.frame(resection_id = "0058-R3", chr4hd = "No"),
    data.frame(resection_id = "0108-R3", chr4hd = "No"), # low purity
    data.frame(resection_id = "0013-R1", chr4hd = "Yes"), # obviously present and obviously subclonal
    data.frame(resection_id = "0068-R3", chr4hd = "No"),
    data.frame(resection_id = "0017-R3", chr4hd = "No"),
    data.frame(resection_id = "0030-R1", chr4hd = "No"),
    data.frame(resection_id = "0094-R1", chr4hd = "Inconclusive"), # has an arm loss
    data.frame(resection_id = "0052-R1", chr4hd = "No"),
    data.frame(resection_id = "0124-R2", chr4hd = "Yes")
  )

plt <- metadata |> 
  dplyr::left_join(out, by=c('array_sentrix_id'='array_sentrix_id')) |>  
  dplyr::select(array_sentrix_id, array_methylation_bins_1p19q_purity, array_mnp_CNVP_v12.8_v5.2_CNVP_segment_chr4_weighted_mean, resection_id) |> 
  dplyr::left_join(marks, by=c('resection_id'='resection_id'))


# impute
plt <- plt |> 
  dplyr::mutate(chr4hd = ifelse(is.na(chr4hd) & array_mnp_CNVP_v12.8_v5.2_CNVP_segment_chr4_weighted_mean > -0.05, "No", chr4hd)) |> 
  dplyr::mutate(chr4hd = ifelse(is.na(chr4hd) & array_mnp_CNVP_v12.8_v5.2_CNVP_segment_chr4_weighted_mean < -0.1, "Yes", chr4hd)) |> 
  dplyr::mutate(chr4hd = ifelse(is.na(chr4hd) & resection_id %in% c("0082-R1", "0094-R2"), "Yes", chr4hd))



ggplot(plt, aes(x=array_methylation_bins_1p19q_purity, y=array_mnp_CNVP_v12.8_v5.2_CNVP_segment_chr4_weighted_mean, label=resection_id, col=chr4hd)) +
  geom_point() +
  ggrepel::geom_text_repel(size=2.5) +
  scale_y_continuous(limits = c(-0.3, NA)) +
  scale_color_manual(values = c(`Yes`= 'darkgreen', `No` = 'red', `Inconclusive` = mixcol('yellow', 'orange'), `NA` = 'darkgray'))


metadata |> 
  dplyr::filter(resection_id == "0124-R2") |> 
  dplyr::select(array_mnp_CNVP_v12.8_v5.2_CNVP_segments) |> 
  as.data.frame()


df = plt |> dplyr::filter(chr4hd != "Inconclusive") |> 
  dplyr::mutate(label = ifelse(chr4hd == "Yes", 1,0)) |> 
  dplyr::rename(x1 = array_methylation_bins_1p19q_purity) |> 
  dplyr::rename(x2 = array_mnp_CNVP_v12.8_v5.2_CNVP_segment_chr4_weighted_mean)

model <- glm(label ~ x1 + x2, data = df, family = "binomial")





# Extract coefficients
coef <- coef(model)
b0 <- coef[1]  # Intercept
b1 <- coef[2]  # Coefficient for x1
b2 <- coef[3]  # Coefficient for x2

# Create a data.frame for the boundary line
x1_vals <- seq(0, max(df$x1), length.out = 100)
x2_vals <- -(coef[1] + coef[2] * x1_vals) / coef[3]
boundary <- data.frame(x1 = x1_vals, x2 = x2_vals)




sprintf("x2 = -(%.3f + %.3f * x1) / %.3f", b0, b1, b2)
sprintf("x1 = -(%.3f + %.3f * x2) / %.3f", b0, b2, b1)



# Plot data and decision boundary
ggplot(df, aes(x = x1, y = x2, color = as.factor(label))) +
  geom_point() +
  geom_line(data = boundary, aes(x = x1, y = x2), color = "black") +
  labs(color = "Class") +
  theme_minimal() +
  xlim(0, 0.7) +
  
  geom_hline(yintercept=-b0/b2, lty=2,lwd=0.5,col="red") +
  geom_abline(intercept=-b0/b2, slope = -b1/b2, lty=2,lwd=0.5,col="red")

sprintf("Decision boundary: weighted CNV offset = %.4f + %.4f * purity", -b0/b2, -b1/b2)


# 5. apply ----


out <- plt |> 
  dplyr::mutate(logit_prob_chr4hd = 

predict(model, plt |>
          dplyr::rename(x1 = array_methylation_bins_1p19q_purity) |> 
          dplyr::rename(x2 = array_mnp_CNVP_v12.8_v5.2_CNVP_segment_chr4_weighted_mean),
        type = "response")
) |> 
  dplyr::mutate(array_mnp_CNVP_v12.8_v5.2_CNVP_logit_chr4hd = ifelse(logit_prob_chr4hd > 0.5, "Yes", "No")) |> 
  dplyr::select(array_sentrix_id, array_mnp_CNVP_v12.8_v5.2_CNVP_logit_chr4hd)

saveRDS(out, file = "cache/analysis_chr4_del__labels.Rds")

