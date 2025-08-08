#!/bin/bash

source /conda/etc/profile.d/conda.sh 
conda activate conservation

echo $1
echo $2

if [ -f $2/$1/$1.fa ]; then
  echo "$1 genome already exists!"
else
  echo "downloading $1 genome..."
  rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/$1/bigZips/$1.fa.gz ../data/genomes/$1/
  echo "...done!"
  echo "uncompressing $1 genome..."
  gunzip $2/$1/$1.fa.gz
  echo "...done!"
fi

if [ -f ../data/genomes/$1/$1.nsq ]; then
  echo "$1 database already exists!"
else
  echo "creating $1 database..."
  makeblastdb -in $2/$1/$1.fa -parse_seqids -dbtype nucl -out ../data/genomes/$1/$1
  echo "...done!"
fi