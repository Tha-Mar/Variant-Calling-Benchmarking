#create simulated reads using wgsim 

#this command generates 6,551,823 forward and reverse reads at 150 length each to achieve approximately 30x coverage.
#it sets a mutation rate at .001 with 15% being indels and 85% substitutions
wgsim -N 6551823 -1 150 -2 150 -r 0.001 -R 0.15 -S 42 data/chr20/chr20.fa \
    data/raw-sim/R1.fq data/raw-sim/R2.fq > data/truth/truth.vcf
