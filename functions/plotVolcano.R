# Plot volcano function

plotVolcano <- function(df,sig = 0.05, text.size = 20, log.base = exp(1), anno, 
                        colors = c("blue","red","grey"),
                        text_color = c("black","white"),
                        ...){
  
  # trans function
  reverselog_trans <- function(base = log.base) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(-Inf, Inf))
  }
  
  if(class(df) == "DFrame") df <- as.data.frame(df)
  df$`Sign.` <- ifelse(df$FDR >= sig, "not",ifelse(df$logFC > 0,"up","down"))
  df$`Sign.` <- factor(df$`Sign.`, levels = c("up", "down", "not"))
  if(is.null(df$names)) df$names <- row.names(df)
  df$dir <- ifelse(df$logFC < 0, "downr","upr")
  
  # add little FDR value
  df$FDR <- ifelse(df$FDR < 1e-300,df$FDR + 1e-300,df$FDR)
  
  p <- ggplot(df, aes(logFC, FDR)) +
    geom_point(aes(alpha = 0.7, color=`Sign.`)) + 
    scale_color_manual(values = c("up" = colors[1],"down" = colors[2],
                                  "not" = colors[3])) +
    scale_y_continuous(trans = reverselog_trans()) +
    new_scale_color() +
    geom_label_repel(segment.colour = "black", 
                    aes(label=ifelse(names %in% anno, names, "")),
                    show.legend = FALSE,
                    max.overlaps = 100,
                    force_pull = 1,
                    ...) +
    scale_color_manual(values = c("upr" = text_color[1],"downr"  = text_color[2])) +
    scale_fill_manual(values = c("up" = colors[1],"down" = colors[2],
                                 "not" = colors[3])) +
    xlab("log2-Fold Change") + ylab("FDR") + 
    theme_classic(base_size = text.size) +
    labs(alpha = NULL) +
    theme(legend.title=element_text(size=text.size-1),
          legend.text=element_text(size=text.size-3), 
          axis.title =element_text(size=text.size)) + 
    guides(alpha = "none")
  
  
  
  
  
  p
}