library(tidyverse)
library(dplyr)

project <- "GTseq.prod.all_val.loci.new.20readsMin"
load(paste0("data/", project, ".tgt.rda"))
#labels <- read.table(paste0("microhaplot/non.humpback.samples.label.txt"))
#names(labels) <- c("fname", "Indiv", "stratum")

tgt$questionable.hap <- FALSE

# Identify haplotypes with Ns in them
haps.with.ns <- sapply(1:nrow(tgt), function(i){
  ifelse(length(grep("N", tgt$haplo.1[i],)) > 0 || length(grep("N", tgt$haplo.2[i],)) > 0, TRUE, FALSE) 
})
tgt$questionable.hap[which(haps.with.ns)] <- TRUE

# Identify haplotypes with Xs in them
haps.with.Xs <- sapply(1:nrow(tgt), function(i){
  ifelse(length(grep("X", tgt$haplo.2[i],)) > 0 || length(grep("X", tgt$haplo.2[i],)) > 0, TRUE, FALSE) 
})
tgt$questionable.hap[which(haps.with.Xs)] <- TRUE

# identify individuals with more than 2 haplotypes after filtering
tgt$questionable.hap[which(!is.na(tgt$haplo.3))] <- TRUE

genos.to.check <- filter(tgt, questionable.hap == TRUE)

write.csv(genos.to.check, file = paste0("data-raw/", project, ".genos.to.check.csv"))

#tgt <- left_join(tgt, select(labels, c(Indiv, stratum)))
#geno.table <- left_join(select(labels, c(Indiv, stratum)), geno.table)

geno.table$num.genos <- do.call(rbind, lapply(1:nrow(geno.table), function(i){
  length(which(!is.na(geno.table[i,2:ncol(geno.table)])))/2
}))

genos.per.locus <- do.call(rbind, lapply(2:(ncol(geno.table)-1), function(l){
  length(which(!is.na(geno.table[1:nrow(geno.table),l])))
}))
genos.per.locus.names <- do.call(rbind,lapply(strsplit(colnames(geno.table)[2:(ncol(geno.table)-1)], split = "_"), function(loc){paste0("Mnov_gtseq_", loc[4])}))
genos.per.locus <- distinct(data.frame(genos.per.locus.names, genos.per.locus))
names(genos.per.locus) <- c("locus", "genos")

loci <- unique(genos.per.locus$locus)
haps.per.locus <- do.call(rbind, lapply(loci, function(loc){
  haps <- filter(genos.filtered, locus == loc)
  length(unique(haps$haplo))
}))
haps.per.locus <- data.frame(loci, as.integer(haps.per.locus))
names(haps.per.locus) <- c("locus", "num.haps")
loc.sum <- left_join(haps.per.locus, genos.per.locus)
write.csv(loc.sum, file = paste0("results/", project, "locus.summary.csv"))
save(loc.sum, file = paste0("data/", project, "locus.summary.rda"))
