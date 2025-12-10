# Variant Calling Benchmarking

Benchmarking BCFtools, FreeBayes, and GATK HaplotypeCaller on simulated chr20 data.

## Setup

Need conda. Install everything with:

```bash
conda create -n variantcall -c bioconda -c conda-forge \
    bcftools freebayes samtools bwa-mem2 gatk4 matplotlib numpy python=3.10

conda activate variantcall
```

## Running

```bash
bash run_pipeline.sh
```

Takes about 2-3 hours. FreeBayes is slow.

## What it does

1. Downloads chr20 from Ensembl (release 112)
2. Indexes it with samtools, BWA-MEM2, and GATK
3. Simulates 6.5M read pairs with wgsim (30x coverage, 0.1% mutation rate, 15% indels)
4. Converts wgsim truth output to VCF
5. Aligns reads with BWA-MEM2
6. Calls variants with all three tools
7. Compares to truth set (TP/FP/FN)
8. Plots results

## Results

Check `results/plots/` for figures and `results/metrics/accuracy_summary.txt` for numbers.

BCFtools is fast and has good recall but lots of false positives. FreeBayes is extremely slow and calls way too many variants on simulated data. GATK is somewhere in between.

## Directory structure

```
data/chr20/          reference genome and indices
data/raw-sim/        simulated fastq files
data/truth/          truth VCF from wgsim
output/bwa-mem2/     aligned BAM
results/vcf/         called variants
results/metrics/     runtime and accuracy stats
results/plots/       comparison plots
```

## Important notes

- GATK needs read groups in the BAM - already handled in alignment step
- GATK needs .dict file - created during indexing
- wgsim outputs IUPAC codes, conversion script turns them into single bases
- FreeBayes doesn't multithread well, expect ~2 hours
- Simulated data is easier than real data, accuracy numbers are optimistic

## Manual steps if you don't want to use the master script

Download reference:
```bash
cd data/chr20
wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
mv Homo_sapiens.GRCh38.dna.chromosome.20.fa chr20.fa
samtools faidx chr20.fa
gatk CreateSequenceDictionary -R chr20.fa
bwa-mem2 index chr20.fa
cd ../..
```

Simulate reads:
```bash
wgsim -N 6551823 -1 150 -2 150 -r 0.001 -R 0.15 -S 42 \
    data/chr20/chr20.fa data/raw-sim/R1.fq data/raw-sim/R2.fq > data/truth/truth_raw.txt
```

Convert truth to VCF (wgsim uses IUPAC codes):
```bash
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
```

Align:
```bash
bwa-mem2 mem -t 8 data/chr20/chr20.fa \
    -R '@RG\tID:sim\tSM:sim_sample\tPL:ILLUMINA\tLB:lib1\tPU:unit1' \
    data/raw-sim/R1.fq data/raw-sim/R2.fq > output/bwa-mem2/sim-alignment.sam
samtools sort -@ 8 -o output/bwa-mem2/sim-alignment.bam output/bwa-mem2/sim-alignment.sam
samtools index output/bwa-mem2/sim-alignment.bam
```

Run callers, evaluate, and plot:
```bash
bash scripts/run_callers.sh
bash scripts/evaluate_calls.sh
python scripts/plot_results.py
```

## Troubleshooting

`GATK fails with read group error` - make sure you used -R flag during alignment

`GATK fails with dict error` - run `gatk CreateSequenceDictionary -R data/chr20/chr20.fa`

`command not found` - activate the conda environment

## Citations

BCFtools: Li H. (2011) Bioinformatics. 27(21):2987-93

FreeBayes: Garrison E, Marth G. (2012) arXiv:1207.3907

GATK: McKenna A, et al. (2010) Genome Res. 20:1297-303

BWA-MEM2: Vasimuddin Md, et al. (2019) IEEE IPDPS
