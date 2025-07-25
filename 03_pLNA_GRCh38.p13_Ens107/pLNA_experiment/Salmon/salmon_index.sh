folder=/mnt/reference/reference/Homo_sapiens/Ensembl/GRCh38.p13/Annotation/Release_107-2022/
salmon=/common/salmon-1.8.0/bin/salmon

$salmon index -t $folder/transcripts_plus_decoys.fa.gz -i $folder/salmon1.8.0 --decoys $folder/decoys.txt -p 6 -k 31 
