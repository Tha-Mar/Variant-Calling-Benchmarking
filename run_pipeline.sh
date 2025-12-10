#!/usr/bin/env bash
set -e

mkdir -p data/chr20 data/raw-sim data/truth output/bwa-mem2 results/vcf results/metrics results/eval results/plots

# download reference if needed
if [[ ! -f data/chr20/chr20.fa ]]; then
    cd data/chr20
    wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
    gunzip Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
    mv Homo_sapiens.GRCh38.dna.chromosome.20.fa chr20.fa
    cd ../..
fi

# index reference
if [[ ! -f data/chr20/chr20.fa.fai ]]; then
    samtools faidx data/chr20/chr20.fa
fi

if [[ ! -f data/chr20/chr20.fa.0123 ]]; then
    bwa-mem2 index data/chr20/chr20.fa
fi

# simulate reads
if [[ ! -f data/raw-sim/R1.fq ]]; then
    bash scripts/simulate_reads.sh
fi

# align
if [[ ! -f output/bwa-mem2/sim-alignment.bam ]]; then
    bash scripts/alignment-index.sh
fi

# call variants
bash scripts/run_callers.sh

# evaluate
bash scripts/evaluate_calls.sh > results/metrics/accuracy_summary.txt

# plot
python scripts/plot_results.py
