### Script that performs (reciprocal) blast in order to identify human specific miRNAs

reverse_blast <- FALSE
genome_path <- "/mnt/reference/reference/other_genomes"

suppressPackageStartupMessages({
  library(BSgenome)
  library(GenomicFeatures)
  library(reticulate)
  library(tidyverse)
  library(openxlsx)
})


get_mirna_seqs <- function(genome, organism, onlyPrimary = TRUE) {
  gff <- paste(organism, "gff3", sep = ".")
  grl <- rtracklayer::import(file.path("../data/mirnas/", gff)) # load annotation
  grl <- grl[seqnames(grl) %in% seqlevels(genome)]
  if(onlyPrimary) {
    grl <- grl[grl$type=="miRNA_primary_transcript"] # filter for primary transcripts
    outfile <- paste("../data/mirnas/", organism, ".fa", sep = "")
  } else {
    outfile <- paste("../data/mirnas/", organism, "_all.fa", sep = "")
  }
  seqs <- suppressWarnings(BSgenome::getSeq(genome, grl)) # extract sequences
  names(seqs) <- mcols(grl)$Name # name sequences
  seqs <- seqs[!duplicated(names(seqs))] # remove duplicates
  writeXStringSet(seqs, outfile)
}


## Downloading from mirbase
organisms = c("hsa", "ptr", "ggo", "ppy", "mmu", "rno", "cfa", "bta", "gga")
for(organism in organisms){
  org_file <- paste(organism, "gff3", sep = ".")
  dest_file <- file.path("../data/mirnas", org_file)
  if(!file.exists(dest_file))
    download.file(file.path("ftp://mirbase.org/pub/mirbase/CURRENT/genomes", org_file), dest_file)
}

## get genomes and create blast databases
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

for(org in orgs){
  system2("./get_genome.sh", c(org, genome_path))
}

## get miRNA sequences for human
## only primary sequences:
if(!file.exists("../data/mirnas/hsa.fa")) {
  get_mirna_seqs(genome_hsa, "hsa") # generate fasta file with miRNAs
  rm(genome_hsa)
}

## all sequences:
if(!file.exists("../data/mirnas/hsa_all.fa")) {
  if(!exists("genome_hsa")) genome_hsa <- readDNAStringSet("../data/genomes/hg38/hg38.fa")  
  genome_hsa <- readDNAStringSet("../data/genomes/hg38/hg38.fa")
  get_mirna_seqs(genome_hsa, "hsa", onlyPrimary = FALSE) # generate fasta file with miRNAs
  rm(genome_hsa)
}

## create data frames with orthologs and mismatches
gr_hsa <- rtracklayer::import("../data/mirnas/hsa.gff3")
gr_hsa <- gr_hsa[gr_hsa$type=="miRNA_primary_transcript"]
df_orthologs <- as.data.frame(gr_hsa) %>%
  transmute(miRNA_id = Name, Chromosome = seqnames, Strand = strand, Start = start, End = end) %>%
  arrange(miRNA_id)
df_mismatches <- as.data.frame(gr_hsa) %>%
  transmute(miRNA_id = Name, Chromosome = seqnames, Strand = strand, Start = start, End = end) %>%
  arrange(miRNA_id)

## reciprocal blast for all organisms
blast_res_names <- c("qseqid", "seqid", "sstart", "send", "strand", "sseq", "qlen", "qstart", "qend", "length", "evalue")
mirna_file <- "../data/mirnas/hsa.fa"
message("Performing reciprocal blast for all organisms...")
for(org in orgs) {
  if(org == "hg38") next
  
  ## blast hsa miRNAs against org genome
  message(paste("blast hsa miRNAs against", org, "genome..."))
  db_path <- file.path(genome_path, org, org)
  system2("./rec_blast.sh", c(mirna_file, db_path, 1, 1))
  
  ## load blast hits from results file
  hits_file <- file.path("../data/hits", paste("hsa", org, sep="_"))
  hits <- read.table(hits_file)
  names(hits) <- blast_res_names
  
  ## filter and convert to GRanges
  hits_filt <- hits %>%
    filter(length/qlen > 0.6 & length/qlen < 1.3) 
  
  hits_gr <- hits_filt %>%
    mutate(start = pmin(sstart, send), end = pmax(sstart, send)) %>%
    mutate(strand = ifelse(strand == "plus", "+", "-")) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    sortSeqlevels() %>%
    sort()
  names(hits_gr) <- hits_gr$qseqid
  
  message("...done!")
  
  ## extract hit sequences from org genome
  message("extracting hit sequences...")
  genome_path_org <- file.path(genome_path, org, paste(org, ".fa", sep = ""))
  genome <- readDNAStringSet(genome_path_org)
  hits_seqs <- getSeq(genome, hits_gr)
  rm(genome)
  hits_seqs_file <- file.path("../data/hits", paste("hsa_", org, ".fa", sep=""))
  writeXStringSet(hits_seqs, hits_seqs_file)
  
  has_match <- df_orthologs$miRNA_id %in% hits_filt$qseqid
  
  message("...done!")
  if(reverse_blast){
    ##  blast hits against hsa genome 
    message(paste("blast", org, "sequences against hsa genome..."))
    system2("./rec_blast.sh", c(hits_seqs_file, "../data/genomes/hg38/hg38", 1, 1))
    reverse_hits_file <- file.path("../data/hits", paste("hsa", org, "hg38", sep="_"))
    reverse_hits <- read.table(reverse_hits_file)
    names(reverse_hits) <- blast_res_names
    
    reverse_hits_filt <- reverse_hits %>%
      filter(length/qlen > 0.6 & length/qlen < 1.3) 
    
    reverse_hits_gr <- reverse_hits_filt %>%
      mutate(start = pmin(sstart, send), end = pmax(sstart, send)) %>%
      mutate(strand = ifelse(strand == "plus", "+", "-")) %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
      sortSeqlevels() %>%
      sort()
    names(reverse_hits_gr) <- reverse_hits_gr$qseqid
    message("...done!")
    
    ol <- findOverlaps(gr_hsa, reverse_hits_gr)
    df_ol <- data.frame(n1 = gr_hsa$Name[ol@from], n2= names(reverse_hits_gr)[ol@to])
    df_ol_match <- df_ol[df_ol$n1 == df_ol$n2,]
    
    has_reverse_match <- df_orthologs$miRNA_id %in% df_ol_match$n1
    
    message(paste("Found", sum(has_match), "matches and", sum(has_reverse_match), "reverse matches in", length(gr_hsa), "sequences."))
    
    has_match <- has_match & mas_reverse_match
  } else {
    message(paste("Found", sum(has_match), "matches in", length(gr_hsa), "sequences."))
  }
  df_orthologs[[org]] <- has_match
}

df <- read.delim("../res/miRNA_conservation.csv")

write.table(df_orthologs, "../res/miRNA_conservation.csv", quote=FALSE, row.names = FALSE, sep = "\t")
non_humans <- orgs[orgs!=orgs["human"]]
with_ortholog <- apply(df_orthologs[non_humans], 1, any)
write.xlsx(list(df_orthologs[with_ortholog,], df_orthologs[!with_ortholog,]),
           file = "../res/miRNA_conservation.xlsx",
           sheetName = c("human miR with orthologs", "human miR without orthologs"))





