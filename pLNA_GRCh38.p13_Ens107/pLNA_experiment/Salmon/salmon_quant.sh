ref=/mnt/reference/reference/Homo_sapiens/Ensembl/GRCh38.p13/Annotation/Release_107-2022/
salmon=/common/salmon-1.8.0/bin/salmon

cd /mnt/schratt/internData/2022_soutschek_Ngn2pLNAs_polyARNA/11_2022/X204SC22021444-Z01-F010/01.RawData/

for f in ./*_*/; do
  f2=`basename $f`;
  echo $f2;
  $salmon quant -p 12 -i $ref/salmon1.8.0 \
    -l A -1 "$f"*_1.fq.gz -2 "$f"*_2.fq.gz --validateMappings --gcBias \
    -g $ref/tx2gene -o ../../salmon/$f2
done

