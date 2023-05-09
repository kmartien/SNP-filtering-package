library(vcfR)
library(pinfsc50)
library(reshape2)
library(ggplot2)
library(strataG)
library(tidyverse)
library(dplyr)
#library(easyGgplot2)
library(grid)
library(gridExtra)
library(gdata)
library(Hmisc)
library(Mnov.GTseq.data)
source("R/0 Load vcf2genos functions.R")
source("R/geno.smry.by.locus.R")
source("R/remove.rejected.loci.R")

description = "GTseq.val.all.final.filtered.recode"
date = format(Sys.time(), "%Y%b%d")

vcf <- read.vcfR(paste0("vcf/", description, ".vcf"))

###### (2)VIOLIN AND (3)HEATMAP PLOTS, (4)COUNT BIALLELIC  ##############################
violin.plot(vcf, paste0("results/", description, "_", date))
heat.map(vcf, paste0("results/", description, "_", date))
count.of.biallelic.loci <- sum(is.biallelic(vcf))

###### (5)CONVERT VCF TO TIDY ###########################################################
tgt <- create.tgt(vcf)

# Next two lines replaces | with / in phased loci. NOTE: TRY MICROHAPLOT!!!
tgt$gt_GT_alleles<-sapply(tgt$gt_GT_alleles, gsub, pattern="|", replace="/", fixed = TRUE)
tgt$gt<-sapply(tgt$gt, gsub, pattern="|", replace="/", fixed = TRUE)

# Replace underscores in locus names with hyphens so that I can separate CHROM and POS later
tgt$CHROM<-sapply(tgt$CHROM, gsub, pattern="_", replace="-", fixed = TRUE)
tgt$locus <- paste0(tgt$CHROM,"_",tgt$POS)

###### (6) CREATE DATAFRAME, GET DEPTH BY LOCUS AND #####################################
###### (7) CREATE DATAFRAME WITH COUNTS AND RATIOS, SUMMARIZE GENO TYPES BY LOCUS #######
geno.df <- create.gtype.df(tgt,paste0("results/", description, "_", date))
loc.dpth <- depth.by.locus(tgt,paste0("results/", description, "_", date))
create.gtype.cts.ratios.df(tgt,paste0("results/", description, "_", date))
geno.sum <- geno.smry.by.locus(tgt, paste0("results/", description, "_", date))
#geno.sum$pct.R.homo <- geno.sum$R.HOMO / sum(geno.sum$R.HOMO, geno.sum$HET, geno.sum$A.HOMO)
names(geno.sum)[1] <- "locus"


###### (8) GENERATE MULTIPLOT_LOCI_EDIT.CSV AND PLOT ALLELE COUNTS BY LOCUS #############
#mplot_loci_edit <- readRDS("microhaplot/GTseq.prod.20000.FLASh.locus_annotation.csv")
mplot_loci_edit <- geno.sum %>% mutate(pct.R.homo = R.HOMO/(R.HOMO+HET+A.HOMO))
mplot_loci_edit$status <- "Accept"
mplot_loci_edit$status[which(mplot_loci_edit$pct.R.homo > 0.95 | mplot_loci_edit$pct.R.homo < 0.05)] <- "Reject"
mplot_loci_edit$status[which(mplot_loci_edit$A.HOMO == 0)] <- "Reject"
mplot_loci_edit$status[which(mplot_loci_edit$HET == 0)] <- "Reject"
mplot_loci_edit$minDP <- 10
#mplot_loci_edit$max.ar.hm <- 0.25
#mplot_loci_edit$min.ar.hz <- 0.43
mplot_loci_edit$minAR <-  .15
mplot_loci_edit$maxAR <- .35
# Nate Campbell uses maxAR = 0.33, minAR = 0.1, but different definitions

tgt <- remove.rejected.loci(tgt, mplot_loci_edit)
mplot_loci_edit <- filter(mplot_loci_edit, status == "Accept") %>% select(locus, minDP, minAR, maxAR, status)

###### (8) PLOT ALLELE COUNTS BY LOCUS ###################################################
plot.allele.ct.by.locus(tgt=tgt, locus.ARs=mplot_loci_edit, fn.prefix = paste0("results/", description,"_", date))

save(tgt, mplot_loci_edit, file = paste0("data/", description, "_", date, "_vcf2tgt.rda"))

###### NEXT: USE ALLELE COUNT PLOTS TO ADJUST AR AND minDP VALUES IN #####################
###### "multiplot_loci_edit.csv", THEN MOVE ON TO XXXX.R. ################################

##################### Adjusting min and maxAR ###################################
#################################################################################
#increasing the maxAR brings the blue lines closer to the green line, slope = 1
#increasing the minAR brings the red lines closer to the blue lines
#decreasing them widens the angle

#slope1(upper blue line) is (1-maxAR)/maxAR
#slope2(lower blue line) is maxAR/(1-maxAR)
#slope3(upper red line) is (1-minAR)/minAR
#slope4(lower red line) is minAR/(1-minAR)

# outside of red lines are homozygotes, inside blue lines are heterozygotes,
#in between red and blue are NA's
