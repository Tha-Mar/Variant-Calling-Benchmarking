#use bwa-mem2 to align reads and create sam file 
bwa-mem2 mem -t 8 data/chr20/chr20.fa \
-R '@RG\tID:sim\tSM:sim_sample\tPL:ILLUMINA\tLB:lib1\tPU:unit1' \
data/raw-sim/R1.fq data/raw-sim/R2.fq > output/bwa-mem2/sim-alignment.sam

#sort the alignment with samtools 
samtools sort -@ 8 -o output/bwa-mem2/sim-alignment.bam output/bwa-mem2/sim-alignment.sam

#index alignment with samtools 
samtools index output/bwa-mem2/sim-alignment.bam
