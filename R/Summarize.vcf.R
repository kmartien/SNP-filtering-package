summarize.vcf <- function(vcf.dir, res.dir, fname, res.name){

  vcf <- paste0(vcf.dir,"/",fname)
  res <- paste0(res.dir,"/",res.name)
  print(res)
  
# depth indv/locus
  vcftools.idepth(vcf, res.fn = res)
  vcftools.siteDepth(vcf, res)
  vcftools.genoDepth(vcf, res)

#  system2(command = "vcftools", args = c("--vcf ", vcf, "--out ", res, "--depth"),
#          stdout = FALSE)
#  system2(command = "vcftools", args = c("--vcf ", vcf, "--out ", res, "--site-mean-depth"),
#          stdout = FALSE)
#  system2(command = "vcftools", args = c("--vcf ", vcf, "--out ", res, "--geno-depth"),
#          stdout = FALSE)
  
# missing data indv/locus
  vcftools.imiss(vcf, res)
  vcftools.siteMiss(vcf, res)
  
#  system2(command = "vcftools", args = c("--vcf ", vcf, "--out ", res, "--missing-indv"),
#          stdout = FALSE)
#  system2(command = "vcftools", args = c("--vcf ", vcf, "--out ", res, "--missing-site"),
#          stdout = FALSE)

# allele freq/indv freq buden
  vcftools.freqBurden(vcf, res)
  vcftools.freq2(vcf, res)
  vcftools.singletons(vcf, res)
  vcftools.012(vcf, res)
#  system2(command = "vcftools", args = c("--vcf ", vcf, "--out ", res, "--indv-freq-burden"),
#          stdout = FALSE)
#  system2(command = "vcftools", args = c("--vcf ", vcf, "--out ", res, "--freq2"),
#          stdout = FALSE)
#  system2(command = "vcftools", args = c("--vcf ", vcf, "--out ", res, "--singletons"),
#          stdout = FALSE)
#  system2(command = "vcftools", args = c("--vcf ", vcf, "--out ", res, "--012"),
#          stdout = FALSE)
  
# heterozygosity per individual
  vcftools.het(vcf, res)
#  system2(command = "vcftools", args = c("--vcf ", vcf, "--out ", res, "--het"),
#          stdout = FALSE)

# SNP call quality
  vcftools.siteQual(vcf, res)
#  system2(command = "vcftools", args = c("--vcf ", vcf, "--out ", res, "--site-quality"),
#          stdout = FALSE)

visualize.vcf.ind(res.dir, res.name)
visualize.vcf.loci(res.dir, res.name)

}