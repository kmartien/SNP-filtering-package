FNAME=$1
RES=$2

#2.1 Remove genotypes with minDP < 5*
#  Includes only genotypes greater than or equal to the "--minDP" value
vcftools --vcf $FNAME.vcf --out $FNAME.minDP5 --minDP 5 --recode --recode-INFO-all

#2.2 Remove loci with minQ < 20*
#  Includes only loci with QUAL > this threshold (calculated across individuals).
vcftools --vcf "$FNAME.minDP5.recode.vcf" --out "$FNAME.minDP5minQ20" --minQ 20 --recode --recode-INFO-all

#2.3 Remove loci with meanDP > 15
# Includes only loci with mean read depth across all individuals greater than or equal to the "min-meanDP" value
vcftools --vcf $FNAME.minDP5minQ20.recode.vcf --out $FNAME.minDP5minQ20meanDP15 --min-meanDP 15 --recode --recode-INFO-all

#2.4 Remove loci with minor allele count < 3*
#  Include only sites with Minor Allele Count greater than or equal to the "--mac" value. Allele count is simply the number of times that allele appears over all individuals at that site.
vcftools --vcf $FNAME.minDP5minQ20meanDP15.recode.vcf --out $FNAME.minDP5minQ20meanDP15mac3 --mac 3 --recode --recode-INFO-all

#2.5 Remove monomorphic sites*
vcftools --vcf $FNAME.minDP5minQ20meanDP15mac3.recode.vcf --out $RES --min-alleles 2 --recode --recode-INFO-all
