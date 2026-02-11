

library(survival)
library(ggplot2)
library(dplyr)
library(patchwork)



svvl_forest <- function(model_or_data,
                        estimate = NULL,
                        lower = NULL,
                        upper = NULL,
                        label = NULL,
                        group = NULL,
                        null_line = NULL,
                        point_size = 1.5,
                        color = "black",
                        shape = 15,
                        xlim = NULL,
                        font_size = theme_cellpress_size,
                        exponentiate = TRUE) {
  
  require(ggplot2)
  require(broom)
  
  # Check if input is a coxph model
  if (inherits(model_or_data, "coxph")) {
    # Extract tidy results from coxph model
    data <- broom::tidy(model_or_data, conf.int = TRUE, exponentiate = exponentiate) |> 
      dplyr::mutate(term = gsub('`','',term))
    
    # Create formatted labels with HR and CI (without p-value)
    data$label <- paste0(
      data$term, " ",
      "HR: ", sprintf("%.2f", data$estimate),
      " (", sprintf("%.2f", data$conf.low), 
      "-", sprintf("%.2f", data$conf.high), ")"
    )
    
    # Create p-value labels
    data$p_label <- ifelse(data$p.value < 0.00001, "p<0.00001", 
                           sprintf("p=%.4f", data$p.value))
    
    # Extract concordance index from model
    concordance <- summary(model_or_data)$concordance
    c_index <- concordance[1]
    c_se <- concordance[2]
    caption_text <- sprintf("Concordance Index: %.3f (SE: %.3f)", c_index, c_se)
    
    # Set defaults for model input
    if (is.null(null_line)) null_line <- if(exponentiate) 1 else 0
    if (is.null(xlab)) xlab <- if(exponentiate) "Hazard Ratio" else "Coefficient"
    
    estimate_col <- rlang::sym("estimate")
    lower_col <- rlang::sym("conf.low")
    upper_col <- rlang::sym("conf.high")
    label_col <- rlang::sym("label")
    
  } else {
    # Input is a data frame
    data <- model_or_data
    data$p_label <- NULL  # No p-values for data frame input
    caption_text <- NULL
    
    # Set defaults for data frame input
    if (is.null(null_line)) null_line <- 0
    if (is.null(xlab)) xlab <- "Effect Size"
    
    estimate_col <- rlang::enquo(estimate)
    lower_col <- rlang::enquo(lower)
    upper_col <- rlang::enquo(upper)
    label_col <- rlang::enquo(label)
  }
  
  group_col <- rlang::enquo(group)
  
  # Add row number for alternating background
  #data$row_num <- seq_len(nrow(data))
  data <- data |> dplyr::mutate(row_num = order(order(term)))
  
  # Create base plot with alternating backgrounds
  p <- ggplot(data, aes(y = reorder(label, -row_num))) +
    geom_rect(aes(xmin = 0.0, xmax = Inf, 
                  ymin = as.numeric(factor(!!label_col)) - 0.5,
                  ymax = as.numeric(factor(!!label_col)) + 0.5,
                  fill = row_num %% 2 == 0),
              alpha = 0.3, show.legend = FALSE) +
    scale_fill_manual(values = c("TRUE" = "gray90", "FALSE" = "white")) +
    geom_vline(xintercept = null_line, linetype = "dashed", 
               color = "gray60", linewidth = theme_cellpress_lwd) +
    geom_errorbarh(aes(xmin = !!lower_col, xmax = !!upper_col),
                   height = 0.2, linewidth = theme_cellpress_lwd * 2, color = color) +
    geom_point(aes(x = !!estimate_col), 
               size = point_size, shape = shape, color = color) +
    {if(!is.null(data$p_label)) geom_text(data=subset(data, !is.na(p_label)),
      aes(x = !!upper_col, label = p_label), 
                                          hjust = -0.1,
                                          size = font_size,
                                          color = "gray30")} +
    labs(x = "Hazard ratio", y = NULL, caption = caption_text) +
    #scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    #theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(hjust = 0),
      plot.title = element_text(hjust = 0.5)
    )

  
  # Add x-axis limits if specified
  #if (!is.null(xlim)) {
  #  p <- p + coord_cartesian(xlim = xlim)
  #}
  
  # Add faceting if group is specified
  #if (!rlang::quo_is_null(group_col)) {
  #  p <- p + facet_grid(rows = vars(!!group_col), scales = "free_y", space = "free_y")
  #}
  
  return(p)
}







