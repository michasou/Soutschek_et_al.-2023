#!/bin/bash

source /conda/etc/profile.d/conda.sh
conda activate conservation

S1=`basename $1 .fa`

# Add name to outfile
echo ">>$2" >> $3

# convert alignment from fasta to clustal format
readseq $1 -all -format=clustal -output=tmp.aln &>/dev/null

# run RNAz
RNAz tmp.aln >> $3
echo >> $3

rm tmp.aln