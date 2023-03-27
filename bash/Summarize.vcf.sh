VCF=$1
RES=$2

# depth indv/locus
vcftools --vcf $VCF --out $RES --depth
vcftools --vcf $VCF --out $RES --site-mean-depth
vcftools --vcf $VCF --out $RES --geno-depth

# missing data indv/locus
vcftools --vcf $VCF --out $RES --missing-indv
vcftools --vcf $VCF --out $RES --missing-site

# allele freq/indv freq buden
vcftools --vcf $VCF --out $RES --indv-freq-burden
vcftools --vcf $VCF --out $RES --freq2
vcftools --vcf $VCF --out $RES --singletons
vcftools --vcf $VCF --out $RES --012

# heterozygosity per individual
vcftools --vcf $VCF --out $RES --het

# SNP call quality
vcftools --vcf $VCF --out $RES --site-quality
