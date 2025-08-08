######### SE equality tester #########
# tgermade
# 04.02.20


#' SEequality
#'
#' @param se1 SE object
#' @param se2 SE object
#'
#' @return text with detailed information about the equality of the SE objects
#'
SEequality <- function(se1, se2){

  # libraries
  
  suppressPackageStartupMessages({
    library(edgeR)
    library(SummarizedExperiment)
    library(SEtools)
  })

  # tests
  
  ## assay dimensions
  dimensions <- all(dim(se1) == dim(se2))
  ## rownames
  tx_names <- all(rownames(se1) %in% rownames(se2)) & all(rownames(se2) %in% rownames(se1))
  ## row order
  tx_order <- all(rownames(se1) == rownames(se2))
  ## assay number
  nr_assays <- length(assays(se1)) == length(assays(se2))
  ## assay content
  if(nr_assays){
    tx_content <- all(
      sapply(1:length(assays(se1)), FUN=function(i) all(assays(se1)[[i]]==assays(se2)[[i]]))
    )
  } else tx_content <- FALSE
  ## rowData column number
  nr_rowdata <- length(rowData(se1)) == length(rowData(se2))
  ## rowData column names
  cols_rowdata <- all(colnames(rowData(se1)) == colnames(rowData(se2)))
  ## rowData dimensions
  if(nr_rowdata & length(rowData(se1))>0){
    dim_rowdata <- all(
      sapply(1:length(rowData(se1)), FUN=function(i)
        all(dim(rowData(se1)[[i]]) == dim(rowData(se2)[[i]])))
    )
  }
  ## colData dimensions
  dims_coldata <- all(dim(colData(se1)) == dim(colData(se2)))
  ## colData colnames
  colnames_coldata <- all(colnames(colData(se1)) %in% colnames(colData(se2))) &
    all(colnames(colData(se2)) %in% colnames(colData(se1)))
  ## colData columns
  cols_coldata <- all(
    sapply(1:length(colData(se1)), FUN=function(i)
      all(se1[[i]]==se2[[i]]))
    )
  
  
  # output
  
  cat(" same assay dimensions:", dimensions, "\n",
      "same assay rownames:", tx_names, "\n",
      "same assay rowname order:", tx_order, "\n",
      "same assay content:", tx_content, "\n")
  
  if(length(colnames(rowData(se1))) >= 1){
    cat(" same number rowData columns:", nr_rowdata, "\n",
        "same rowData column names:", cols_rowdata, "\n")
  } else cat(" no rowData columns\n")
  
  if(exists("dim_rowdata")){
    cat(" same rowData column dimensions:", dim_rowdata, "\n")
  } else cat(" no equality info about rowData columns\n")
  
  cat(" same colData dimensions:", dims_coldata, "\n",
      "same colData names:", colnames_coldata, "\n",
      "same colData content:", cols_coldata, "\n")
}
