#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=5120


ref=/nfs/nas22/fs2202/biol_bc_kleele_2/Michael/Bioinformatics/Reference/Homo_Sapiens/Ensembl/GRCh38.p14/Release_112_2024
salmon=~/Salmon/salmon-latest_linux_x86_64/bin/salmon

cd /nfs/nas22/fs2202/biol_bc_kleele_2/Michael/Bioinformatics/Internal_Data/2019_soutschek_Ngn2Timecourse_Human_RiboZero/

for f in ./NovaSeq_20190930_NOV203_o5949_DataDelivery/*.fastq.gz; do
  f2=`basename $f _R1.fastq.gz`;
  f3=${f2/20190930.A-/};
  echo $f3;
  $salmon quant -p 24 -i $ref/salmon_index \
    -l A -r $f --validateMappings --dumpEq \
    -g $ref/tx2gene -o ./Salmon_Ens112/salmon/$f3
done





