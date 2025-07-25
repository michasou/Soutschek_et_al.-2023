setwd("/mnt/schratt/enrichMir_datasets/STAR/")
lf <- list.files(pattern="ReadsPerGene\\.out\\.tab$")
names(lf) <- gsub("ReadsPerGene.out.tab","",lf,fixed=TRUE)
m <- sapply(lf, FUN=function(x) read.delim(x, header=T)[,4])
row.names(m) <- read.delim(lf[[1]], header=TRUE)[,1]

cd <- read.delim("../SraRunTable.txt", header=TRUE, sep=",", row.names = 1)
cd <- cd[colnames(m),c("Cell_Line","Batch","miRNA")]
cd <- cbind(cd, t(m[1:3,]))
m <- m[-1:-3,]

library(SummarizedExperiment)
se <- SummarizedExperiment( list(counts=m), colData=cd)

g <- rtracklayer::import.gff("/reference/Homo_sapiens/Ensembl/GRCh38.p13/Annotation/Release_98-2019-09-3/genes.gtf")
g2 <- g@elementMetadata[,c("gene_id","gene_name")]
g2 <- g2[!duplicated(g2),]
row.names(g2) <- g2$gene_id

rowData(se)$symbol <- g2[row.names(se),"gene_name"]
se <- se[rowSums(assay(se)>=10)>1,]

se2 <- SEtools::aggSE(se, by="symbol")

saveRDS(se, "geneId.SE.rds")
saveRDS(se2, "symbols.SE.rds")
