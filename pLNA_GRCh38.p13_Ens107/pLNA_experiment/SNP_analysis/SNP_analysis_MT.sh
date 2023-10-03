#https://www.ebi.ac.uk/sites/ebi.ac.uk/files/content.ebi.ac.uk/materials/2014/140217_AgriOmics/dan_bolser_snp_calling.pdf

cd /mnt/schratt/internData/2022_soutschek_Ngn2pLNAs_polyARNA/11_2022/star_indexed
ref=/mnt/schratt/internData/2022_soutschek_Ngn2pLNAs_polyARNA/11_2022/assembly/Homo_sapiens.GRCh38.dna.primary_assembly.fa

for f in ./*sorted.bam ; do
  f2=`basename $f`;
  samtools index ./$f
  samtools view -b $f "MT" > ./MT_$f2; samtools index ./MT_$f2
done

#index genome
samtools faidx $ref


#SNPs
samtools mpileup -uD -f $ref \
   ./MT_*sorted.bam > Results_samtools.raw.bcf
   
   
bcftools view Results_samtools.raw.bcf > Results_SNPs.bcf


#filter
#bcftools view Results_SNPs.bcf | vcfutils.pl varFilter -D100 > Results_SNPs_filtered.bcf
