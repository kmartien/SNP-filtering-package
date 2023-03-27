FNAME=$1
RES=$2

conda activate trimming

vcfallelicprimitives vcf/$FNAME.recode.vcf --keep-info --keep-geno > vcf/$RES.vcf

conda deactivate