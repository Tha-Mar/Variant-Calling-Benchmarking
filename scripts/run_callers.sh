#!/usr/bin/env bash
set -euo pipefail

#run bcftools, this also captures the runtime and additional metrics using the -v flag
#this creates a pileup of the alignment and then calls the variants
#this uses 4 threads for the mpileup and call steps
/usr/bin/time -v bash -c '
  bcftools mpileup --threads 4 \
      -f data/chr20/chr20.fa \
      -Ou output/bwa-mem2/sim-alignment.bam |
  bcftools call --threads 4 -mv -Oz -o results/vcf/bcftools.vcf.gz
' 2> results/metrics/bcftools.time.log

bcftools index -t results/vcf/bcftools.vcf.gz

#run freebayes, this also captures the runtime and additional metrics using the -v flag
#this uses the freebayes algorithm to call the variants
#freebayes is single threaded by default and cant run with multiple threads unless you use freebayes-parallel
/usr/bin/time -v \
  freebayes -f data/chr20/chr20.fa output/bwa-mem2/sim-alignment.bam \
  > results/vcf/freebayes.vcf \
  2> results/metrics/freebayes.time.log

#compress the vcf file and index it using bgzip and tabix
bgzip -f results/vcf/freebayes.vcf
tabix -p vcf results/vcf/freebayes.vcf.gz

#run gatk haplotype caller, this also captures the runtime and additional metrics using the -v flag
#this uses the gatk haplotype caller algorithm to call the variants
#ran using 4 threads
/usr/bin/time -v \
  gatk HaplotypeCaller \
    -R data/chr20/chr20.fa \
    -I output/bwa-mem2/sim-alignment.bam \
    -O results/vcf/gatk.vcf.gz \
    --native-pair-hmm-threads 4 \
  2> results/metrics/gatk.time.log
  #gatk is outputted as a compressed file already 
