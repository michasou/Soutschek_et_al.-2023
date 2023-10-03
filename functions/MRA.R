dea2viper <- function(dea, regulon, pleiotropy = TRUE, verbose=FALSE, weights=c(2, 1)){
  sig <- sign(dea$logFC)*-log10(dea$FDR)
  names(sig) <- row.names(dea)
  regulon <- regulon[intersect(names(regulon), row.names(dea))]
  res <- msviper(sig, regulon, pleiotropy=pleiotropy, verbose=verbose)
  res <- as.data.frame(res$es[c("nes","size","p.value")])
  res$FDR <- p.adjust(res$p.value)
  res <- res[order(res$p.value),]
  colnames(res) <- paste0("activity.",colnames(res))
  colnames(dea) <- paste0("expression.", colnames(dea))
  res <- cbind(factor=row.names(res), res, dea[row.names(res),])
  res$agg.FDR <- apply(res[,c("activity.FDR","expression.FDR")], 1, 
                   weights=weights, FUN=aggregation::lancaster)
  res$agg.FDR[which(sign(res$activity.nes) != sign(res$expression.logFC))] <- 1
  res
}