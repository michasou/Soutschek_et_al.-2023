#' clusterEnrichment
#'
#' @param clusters A named list of sets of genes of interest
#' @param sets A named list of reference genesets (must use the same gene
#' identifiers as `clusters`)
#' @param universe The optional universe (by default, the union of clusters is
#' used)
#' @param minSetSize The minimum set size (in universe) to consider
#' @param threshold The p-value threshold
#' @param family The model family for deviance calculation.
#'
#' @return A data.frame
clusterEnrichment <- function(clusters, sets, universe=NULL, minSetSize=10, threshold=0.05, family=c("binomial","poisson")){
    family <- match.arg(family)
    if(is.null(universe)) universe <- unlist(clusters)
    sets <- lapply(sets, y=universe, FUN=intersect)
    clusters <- lapply(clusters, y=universe, FUN=intersect)
    sets <- sets[sapply(sets,length)>=minSetSize]
    sl <- sapply(sets, length)
    cs <- sapply(clusters, length)
    universe <- length(universe)

    b <- sapply(clusters, FUN=function(x){
        sapply(sets, FUN=function(y){
            length(intersect(x,y))
        })
    })

    # # OLD chi^2 version
    # p <- as.data.frame(t(apply(b, 1, p=cs/sum(cs), FUN=function(x,p){
    #     test <- chisq.test(x, p=p)
    #     w <- which.min(test$expected)
    #     max.dev <- test$expected
    #     max.dev[w] <- sum(test$observed)-max.dev[w]
    #     V <- sqrt( test$statistic / sum( max.dev^2 / test$expected ) )
    #     en <- test$observed/test$expected
    #     c( p.value=test$p.value, stat=as.numeric(test$statistic), FDR=NA_real_,
    #        cramer.V=as.numeric(V), en )
    # })))
    # p$FDR <- p.adjust(p$p.value, method="fdr")


    dev <- sapply(1:nrow(b),sz=cs,mod=family,FUN=function(g, sz, mod){
        x <- b[g,]
        expected <- sz*sl[g]/universe
        enr <- log1p(x)-log1p(expected)
        p=sum(x)/sum(sz)
        # calculation of
        if(mod=="binomial"){
            term1<-sum(x*log(x/(sz*p)), na.rm=TRUE)
            nx<-sz-x
            term2<-sum(nx*log(nx/(sz*(1-p))), na.rm=TRUE)
            dev <- 2*(term1+term2)
        }else{
            dev <- 2*sum(x*log(x/(sz*p)),na.rm=TRUE)-2*sum(x-sz*p)
        }
        pval <- pchisq(dev, length(x)-1, lower.tail=FALSE)
        c(deviance=dev, p.value=pval, FDR=NA_real_, enr)
    })
    dev <- as.data.frame(t(dev))
    dev$FDR <- p.adjust(dev$p.value, method="fdr")
    colnames(dev)[4:ncol(dev)] <- paste("enrichment",colnames(dev)[4:ncol(dev)],sep=".")
    row.names(dev) <- row.names(b)
    dev <- dev[dev$p.value<threshold,]

    dev[order(dev$p.value),]
}

plotClusterEnrichment <- function(dev, k=5, top=NULL, trim.names = 50L,
                                  color=colorRampPalette(c("blue","black", "yellow"))(50),
                                  ...){
    dev <- dev[order(dev$p.value),]
    if(is.null(top)){
        d <- 1-cor(t(dev[,grep("enrichment",colnames(dev))]))
        bg <- cluster_louvain(knn.graph(d,k))
        dev$cluster <- NA_integer_
        dev[bg$names,"cluster"] <- bg$membership

        s1 <- row.names(dev)[unique(apply(dev[,grep("enrichment",colnames(dev))],2,which.max))]
        s2 <- row.names(dev)[!duplicated(dev$cluster)]
        s <- union(s1,s2)
    }else{
        s <- row.names(dev)[seq_len(min(top,nrow(dev)))]
    }
    deve <- dev[s,grep("enrichment",colnames(dev))]
    colnames(deve) <- gsub("enrichment\\.","",colnames(deve))
    label_names <- str_trunc(row.names(deve), trim.names)
    #pheatmap(deve, color=color, border_color = NA, labels_row = label_names, ...)
    ComplexHeatmap::Heatmap(deve, col=color, row_labels = label_names,name = "enrichment", ...)
}


#' knn.graph
#'
#' Builds a KNN graph from a distance matrix
#'
#' @param d A distance matrix
#' @param k The number of neighbors
#'
#' @return A igraph object
knn.graph <- function(d, k=5){
    d <- as.matrix(d)
    diag(d) <- 0
    for(i in seq_len(nrow(d))){
        d[i,-order(d[i,])[seq_len(k)]] <- 0
    }
    igraph::graph_from_adjacency_matrix(d, "undirected", weighted=TRUE)
}
