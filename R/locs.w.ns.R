library(dplyr)

geno.file <- read.csv("results/GTseq_val_all-reported_diploid_haplotype.csv")


loci <- unique(geno.file$locus)

geno.file$contains.ns <- sapply(1:nrow(geno.file), function(i){
  ifelse(length(grep("N", geno.file$haplotype.1[i],)) > 0, TRUE, FALSE) 
})

length(which(geno.file$contains.ns == TRUE))
locs.w.ns <- filter(geno.file, contains.ns == TRUE)
locs.w.ns <- unique(locs.w.ns$locus)

rejected.loci <- locs.w.ns

geno.file.filtered <- geno.file[-which(geno.file$locus %in% rejected.loci),]
  
missing.data.ind <- data.frame(table(geno.file.filtered$indiv.ID)) %>%
  mutate(missing = length(loci)-Freq)
names(missing.data.ind) <- c("labID", "genos", "missing")
length(which(missing.data.ind$missing > round(length(loci)*.1)))

rejected.inds <- missing.data.ind$labID[which(missing.data.ind$missing > round(length(loci)*.1))]

geno.file.filtered <- geno.file.filtered[-which(geno.file.filtered$indiv.ID %in% rejected.inds),]
num.inds <- length(unique(geno.file.filtered$indiv.ID))

missing.data.loc <- data.frame(table(geno.file.filtered$locus)) %>%
  mutate(missing = num.inds-Freq)
length(which(missing.data.loc$missing > round(num.inds * 0.2)))

