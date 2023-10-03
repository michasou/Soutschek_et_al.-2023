suppressPackageStartupMessages({
  library(Biostrings)
  library(openxlsx)
})

get_mismatches <- function(seq1, seq2) {
  pa <- pairwiseAlignment(seq1, seq2)
  n_mismatches <- length(pa@pattern@mismatch@unlistData)
  n_indels <- length(pa@pattern@indel@unlistData)
  return(n_mismatches + n_indels)
}

orgs <- c(
  "human" = "hg38",
  "chimpanzee" = "panTro6",
  "gorilla" = "gorGor6",
  "orangutan" = "ponAbe3",
  "macaque" = "macFas5",
  "marmoset" = "calJac3",
  "mouse" = "mm39",
  "rat" = "rn6", 
  "dog" = "canFam5",
  "cow" = "bosTau9",
  #  "opossum" = "monDom5", ## different format... 
  "chicken" = "galGal6")

## load data
WB <- loadWorkbook("../res/miRNA_orthologs.xlsx")
df_orthologs <- read.xlsx(WB, sheet = "human miR with orthologs")

## create dataframe for mismatches
df_mismatches <- df_orthologs[c("miRNA_id", "Chromosome", "Strand", "Start", "End")]

## count mismatches for orthologs
hsa_seqs <- readDNAStringSet("../data/mirnas/hsa.fa")


for(org in orgs){
  if(org == "hg38") next
  if(org == "panTro6") next
  message(paste("Counting hits for", org, "..."))
  hits_seqs_file <- file.path("../data/hits", paste("hsa_", org, ".fa", sep=""))
  hits_seqs <- readDNAStringSet(hits_seqs_file)
  df_mismatches[[org]] <- sapply(df_mismatches$miRNA_id, function(miRNA_id){
    if(miRNA_id %in% names(hits_seqs)) get_mismatches(hits_seqs[miRNA_id], hsa_seqs[miRNA_id])
    else NA
  })
  message("...done!")
}
sheetName <- "mismatches in orthologs"
if(!sheetName %in% WB$sheet_names) addWorksheet(WB, sheetName)
writeData(WB, sheet = sheetName, df_mismatches)
saveWorkbook(WB, "../res/miRNA_orthologs.xlsx", overwrite = TRUE)

