#create sequence dict for gatk
gatk CreateSequenceDictionary \
  -R data/chr20/chr20.fa \
  -O data/chr20/chr20.dict

#create simulated reads using wgsim
wgsim -N 6551823 -1 150 -2 150 -r 0.001 -R 0.15 -S 42 data/chr20/chr20.fa \
    data/raw-sim/R1.fq data/raw-sim/R2.fq > data/truth/truth_raw.txt

# convert wgsim output to vcf
echo -e "##fileformat=VCFv4.2\n##contig=<ID=20,length=64444167>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > data/truth/truth.vcf
awk 'NR > 3 {
    chrom=$1; pos=$2; ref=$3; alt=$4;
    if (alt == "M") alt="C";
    if (alt == "W") alt="T";
    if (alt == "Y") alt="T";
    if (alt == "R") alt="G";
    if (alt == "K") alt="T";
    if (alt == "S") alt="C";
    if (alt == "-") next;
    print chrom "\t" pos "\t.\t" ref "\t" alt "\t.\tPASS\t.";
}' data/truth/truth_raw.txt >> data/truth/truth.vcf
