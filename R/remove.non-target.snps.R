library(vcfR)
library(vcftoolsR)
library(tidyverse)
library(dplyr)
source("R/0 Load vcf2genos functions.R")

description = "GTseq.prod.gt.20000.corrected.targets"
date = format(Sys.time(), "%Y%b%d")

vcf <- read.vcfR(paste0("vcf/", description, ".vcf"))

tgt <- create.tgt(vcf)
loci <- unique(tgt$locus)

non.target.loci <- filter(tgt, POS != 130) %>% select(c(locus, CHROM, POS)) %>%
  distinct()

write.csv(non.target.loci[,c(2,3)], file="non.target.locs.to.remove.csv", row.names=FALSE)
vcftools.rmPOS(paste0("vcf/", description), 
                 paste0("vcf/", description, ".targets.only"),
                 "non.target.locs.to.remove.txt")

fname <- paste0(description, ".targets.only.recode")
vcf <- read.vcfR(paste0("vcf/", fname, ".vcf"))
