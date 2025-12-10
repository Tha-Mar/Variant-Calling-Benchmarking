
# Evaluate bcftools

mkdir -p results/eval

awk 'NR>3 {print $1":"$2":"$4":"$5}' data/truth/truth.vcf \
  | sort > results/eval/truth.keys
# i used bcftools isec before which worked find with BCFtools and freebayes but had issues with GATK 
#  I need to make key files for each caller. need the keys because isec wasnt working with gatk vcf
#this allows me to extract the calls to compare them to the truth set 
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

comm -12 results/eval/truth.keys results/eval/bcftools.keys | wc -l

comm -23 results/eval/truth.keys results/eval/bcftools.keys | wc -l

comm -13 results/eval/truth.keys results/eval/bcftools.keys | wc -l


comm -12 results/eval/truth.keys results/eval/freebayes.keys | wc -l

comm -23 results/eval/truth.keys results/eval/freebayes.keys | wc -l

comm -13 results/eval/truth.keys results/eval/freebayes.keys | wc -l


comm -12 results/eval/truth.keys results/eval/gatk.keys | wc -l

comm -23 results/eval/truth.keys results/eval/gatk.keys | wc -l

comm -13 results/eval/truth.keys results/eval/gatk.keys | wc -l