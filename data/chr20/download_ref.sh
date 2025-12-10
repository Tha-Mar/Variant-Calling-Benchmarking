#these are the commands i used to download and index the reference

#this downloads chromosome 20 reference genome from ensemble
wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz

#decompress the file
gunzip Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz

#used mv to rename the file 
mv Homo_sapiens.GRCh38.dna.chromosome.20.fa chr20.fa

#index with samtools and bwa
samtools faidx chr20.fa
bwa-mem2 index chr20.fa

