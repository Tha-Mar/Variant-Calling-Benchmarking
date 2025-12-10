
# Evaluate bcftools

mkdir -p results/eval

awk 'NR>3 {print $1":"$2":"$4":"$5}' data/truth/truth.vcf \
  | sort > results/eval/truth.keys

#  I need to make key files for each caller. need the keys because isec wasnt working with gatk vcf

# bcftools 
bcftools view -H results/vcf/bcftools.vcf.gz \
  | awk '{print $1":"$2":"$4":"$5}' \
  | sort > results/eval/bcftools.keys

# freebayes
bcftools view -H results/vcf/freebayes.vcf.gz \
  | awk '{print $1":"$2":"$4":"$5}' \
  | sort > results/eval/freebayes.keys

# gatk
bcftools view -H results/vcf/gatk.vcf.gz \
  | awk '{print $1":"$2":"$4":"$5}' \
  | sort > results/eval/gatk.keys


# Compute TP / FN / FP for each tool

echo "BCFtools"
echo "TP (truth ∩ bcftools):"
comm -12 results/eval/truth.keys results/eval/bcftools.keys | wc -l

echo "FN (truth - bcftools):"
comm -23 results/eval/truth.keys results/eval/bcftools.keys | wc -l

echo "FP (bcftools - truth):"
comm -13 results/eval/truth.keys results/eval/bcftools.keys | wc -l


echo "Freebayes"
echo "TP (truth ∩ freebayes):"
comm -12 results/eval/truth.keys results/eval/freebayes.keys | wc -l

echo "FN (truth - freebayes):"
comm -23 results/eval/truth.keys results/eval/freebayes.keys | wc -l

echo "FP (freebayes - truth):"
comm -13 results/eval/truth.keys results/eval/freebayes.keys | wc -l


echo "GATK"
echo "TP (truth ∩ gatk):"
comm -12 results/eval/truth.keys results/eval/gatk.keys | wc -l

echo "FN (truth - gatk):"
comm -23 results/eval/truth.keys results/eval/gatk.keys | wc -l

echo "FP (gatk - truth):"
comm -13 results/eval/truth.keys results/eval/gatk.keys | wc -l