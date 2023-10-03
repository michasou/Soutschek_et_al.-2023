ref=/mnt/reference/reference/Homo_sapiens/Ensembl/GRCh38.p13/Annotation/Release_107-2022/

wget http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz -P $ref
wget http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P $ref

gzip -dk $ref/Homo_sapiens.GRCh38.107.gtf.gz
gzip -dk $ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir $ref/star_index \
--genomeFastaFiles $ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile $ref/Homo_sapiens.GRCh38.107.gtf \
--sjdbOverhang 149