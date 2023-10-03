#!/bin/bash

source /conda/etc/profile.d/conda.sh 
conda activate conservation

S1=`basename $1 .fa`
S2=`basename $2`

blastn -db $2 \
  -query $1 \
  -dust no \
  -soft_masking false \
  -max_target_seqs $3 \
  -max_hsps $4 \
  -outfmt '6 qseqid sseqid sstart send sstrand sseq qlen qstart qend length evalue' \
  -out ../data/hits/"$S1"_"$S2" \
  -num_threads 4