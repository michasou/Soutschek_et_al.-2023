#new plot volcano function

plotVolcano_new <- function(
    df,
    sig = 0.05,
    text.size = 20,
    anno = NULL,
    xlim_min = -5, xlim_max = 5,
    ylim_max = NA,
    colors = c("blue", "red", "grey"),
    text_color = c("black", "white"),
    ybreaks = c(0, 2, 5, 10, 25, 269),
    yaxis_mode = c("FDR", "neglog10"),
    ...
){
  library(ggplot2)
  library(ggrepel)
  library(ggnewscale)
  library(scales)
  
  yaxis_mode <- match.arg(yaxis_mode)
  
  # Convert DFrame if needed
  if (inherits(df, "DFrame")) df <- as.data.frame(df)
  
  # Significance flag
  df$`Sign.` <- ifelse(df$FDR >= sig, "not",
                       ifelse(df$logFC > 0, "up", "down"))
  df$`Sign.` <- factor(df$`Sign.`, levels = c("up", "down", "not"))
  
  if (is.null(df$names)) df$names <- row.names(df)
  df$dir <- ifelse(df$logFC < 0, "downr", "upr")
  
  # Avoid log issues
  df$FDR[df$FDR < 1e-300] <- 1e-300
  df$neglog10FDR <- -log10(df$FDR)
  
  # Annotation subset
  if (!is.null(anno)) {
    anno <- df[df$names %in% anno, ]
  }
  
  # Optional reverse-log transform
  reverselog_trans <- function(base = 10) {
    trans <- function(x) -log(x, base)
    inv   <- function(x) base^(-x)
    trans_new(
      paste0("reverselog-", base),
      trans, inv,
      domain = c(1e-300, 1)
    )
  }
  
  # Core plot
  p <- ggplot(
    df,
    aes(
      x = logFC,
      y = if (yaxis_mode == "neglog10") neglog10FDR else FDR
    )
  ) +
    geom_point(aes(alpha = 0.7, color = `Sign.`)) +
    scale_color_manual(values = c(
      "up"   = colors[1],
      "down" = colors[2],
      "not"  = colors[3]
    )) +
    labs(alpha = NULL)
  
  # ----- Y-AXIS OPTIONS -----
  if (yaxis_mode == "neglog10") {
    p <- p + scale_y_continuous(
      name   = expression(-log[10](FDR)),
      breaks = ybreaks,
      labels = function(x) format(10^(-x), scientific = TRUE),
      expand = expansion(mult = c(0.02, 0.05)),
      limits = if (!is.na(ylim_max)) c(0, ylim_max) else NULL
    )
  } else if (yaxis_mode == "FDR") {
    p <- p + scale_y_continuous(
      name   = "FDR",
      trans  = reverselog_trans(base = 10),
      breaks = 10^(-ybreaks),
      labels = scales::scientific,
      expand = expansion(mult = c(0.02, 0.05)),
      limits = if (!is.na(ylim_max)) c(10^(-ylim_max), 1) else NULL
    )
  }
  # Annotation layer
  if (!is.null(anno) && nrow(anno) > 0) {
    p <- p +
      ggnewscale::new_scale_color() +
      ggrepel::geom_label_repel(
        data = anno,
        aes(label = names, color = dir),
        segment.colour = "black",
        show.legend = FALSE,
        max.overlaps = 100,
        force_pull = 1,
        ...
      ) +
      scale_color_manual(values = c("upr" = text_color[1], "downr" = text_color[2]))
  }
  
  
  #more infos
  p <- p +  xlab("log2-Fold Change") +
    theme_classic(base_size = text.size) +
    theme(
      legend.title = element_text(size = text.size - 1),
      legend.text  = element_text(size = text.size - 3),
      axis.title   = element_text(size = text.size)
    ) +
    guides(alpha = "none")
  
  if (!all(is.na(c(xlim_min, xlim_max)))) {
    p <- p + coord_cartesian(xlim = c(xlim_min, xlim_max))
  }
  
  p
}