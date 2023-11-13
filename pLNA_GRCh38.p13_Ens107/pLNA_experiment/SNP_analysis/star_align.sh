#! /bin/bash

ref=/mnt/reference/reference/Homo_sapiens/Ensembl/GRCh38.p13/Annotation/Release_107-2022/

cd /mnt/schratt/internData/2022_soutschek_Ngn2pLNAs_polyARNA/11_2022/X204SC22021444-Z01-F010/01.RawData/

for f in ./*_*/; do
  f2=`basename $f`
  echo $f2;
  
  STAR --genomeDir $ref/star_index \
  --genomeLoad NoSharedMemory \
  --readFilesIn "$f"*_1.fq.gz "$f"*_2.fq.gz \
  --readFilesCommand zcat \
  --runThreadN 12 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \
  --alignSJDBoverhangMin 1 \
  --alignSJoverhangMin 8 \
  --outFilterMatchNmin 30 \
  --outFilterMismatchNmax 10 \
  --outFilterMismatchNoverLmax 0.05 \
  --outFilterMultimapNmax 50 \
  --twopassMode Basic \
  --chimSegmentMin 15 \
  --chimJunctionOverhangMin 15 \
  --chimScoreMin 15 \
  --chimScoreSeparation 10 \
  --chimOutType Junctions SeparateSAMold\
  --outFileNamePrefix ../../star/$f2. \
  --outSAMtype BAM Unsorted \
  --alignEndsProtrude 3 ConcordantPair \
  --outSAMattributes All \
  --outMultimapperOrder Random > ../../star/$f2.log

samtools sort -l 9 -m 3500M -@ 12 ../../star/$f2.Aligned.out.bam -o ../../star/$f2.sorted.bam &&
samtools index ../../star/$f2.sorted.bam &&
rm ../../star/$f2.Aligned.out.bam

done





