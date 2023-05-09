library(dplyr)

project <- "GTseq.val.all"

min.allele.balance.het <- 3/7
min.AR.het <- 3/10
max.AR.homo <- 2/10
min.genos.per.ind <- 307

locus_annotation <- readRDS(paste0("data/", project, "-microhaplot-locus_annotation.rda"))
genos <- readRDS(paste0("/Users/Shared/KKMDocuments/Documents/Github.Repos/Shiny/microhaplot/", project, ".rds"))[,-1]
genos <- mutate(genos, id.loc = paste0(id, "-", locus))
genos$to.remove <- FALSE

# remove haplotypes that are NA
genos$to.remove[which(is.na(genos$haplo))] <- TRUE

# remove haplotypes with rank > 2 and AR < min.allele.balance.het 
genos$to.remove[which(genos$rank > 2 & genos$allele.balance < min.allele.balance.het)] <- TRUE


# check loci for which some individuals have more than one haplotype, decide which to drop
excess.haps <- filter(genos, rank > 2) %>% filter(to.remove == FALSE)
inds.to.check <- filter(genos, id.loc %in% excess.haps$id.loc)
locs.to.check <- data.frame(table(inds.to.check$locus)) %>% mutate(Status = "Accept")
names(locs.to.check)[1] <- "locus"
write.csv(locs.to.check, file = paste0("data-raw/", project, ".locs.to.check.csv"), row.names = FALSE)

### SEE IF I CAN FIGURE OUT WHAT'S GOING ON 155, 219, AND 302
locs.to.drop <- read.csv(file = paste0("data-raw/", project, ".locs.to.check.csv")) %>%
  filter(Status %in% c("Reject", "??")) %>% select(locus)

genos$to.remove[which(genos$locus %in% locs.to.drop$locus)] <- TRUE

# remove genotypes for individuals with fewer than 10 reads at a locus
inds <- unique(genos$id)
loci <- unique(genos$locus)
tot.depth <- data.frame(do.call(rbind, lapply(inds, function(i){
  df <- filter(genos, id == i) %>% filter(to.remove == FALSE)
  do.call(rbind, lapply(loci, function(l){
    haps <- filter(df, locus == l)
    data.frame(i, l, sum(haps$depth))
  }))
})))
names(tot.depth) <- c("id", "locus", "tot.depth")

genos.to.drop <- filter(tot.depth, tot.depth < 10) %>% mutate(id.loc = paste0(id, "-", locus))
genos$to.remove[which(genos$id.loc %in% genos.to.drop$id.loc)] <- TRUE

# remove genotypes for individuals that have fewer than 5 reads for its major haplotype at a locus
genos.to.drop <- filter(genos, depth < 5) %>% filter(rank == 1) %>% filter(to.remove == FALSE)
genos$to.remove[which(genos$id.loc %in% genos.to.drop$id.loc)] <- TRUE

# remove filtered genotypes
genos.filtered <- filter(genos, to.remove == FALSE) 
loci <- unique(genos.filtered$locus)

# calculate haplotype freqs at remaining loci
hapfreqs <- lapply(loci, function(l){
  df <- filter(genos.filtered, locus == l)
  hapfreq <- table(df$haplo)
  return(sort(hapfreq, decreasing = TRUE))
})
names(hapfreqs) <- loci

# create tgt-like data structure
tgt <- do.call(rbind, lapply(unique(genos.filtered$id.loc), function(x){
  df <- filter(genos.filtered, id.loc == x)
  loc <- df$locus[1]
  ind <- df$id[1]
  haplo.x <- df$haplo[1]
  depth.x <- df$depth[1]
  if(nrow(df) == 1){
    haplo.y <- haplo.x
    depth.y <- 0
  } else {
    haplo.y <- df$haplo[2]
    depth.y <- df$depth[2]
  }
  data.frame(locus = loc, Indiv = ind, haplo.x = haplo.x, haplo.y = haplo.y, depth.x = depth.x, depth.y = depth.y)
}))
tgt$allele.balance <- tgt$depth.y/tgt$depth.x

# summarize missing data
missing.data.ind <- data.frame(table(tgt$Indiv)) %>%
  mutate(missing = length(loci)-Freq)
names(missing.data.ind) <- c("labID", "genos", "missing")
length(which(missing.data.ind$genos > min.genos.per.ind))
rejected.inds <- missing.data.ind$labID[which(missing.data.ind$genos < min.genos.per.ind)]

tgt <- tgt[-which(tgt$Indiv %in% rejected.inds),]
num.inds <- length(unique(tgt$Indiv))

missing.data.loc <- data.frame(table(tgt$locus)) %>%
  mutate(missing = num.inds-Freq)
length(which(missing.data.loc$missing > round(num.inds * 0.2)))

tgt$CHROM <- tgt$locus
tgt$POS <- NA
tgt$REF <- NA
tgt$ALT <- NA
tgt$QUAL <- NA
tgt$FILTER <- NA
tgt$DP <- tgt$ID <- tgt$gt_GT <- tgt$gt_GL <- NA
tgt$gt_AD <- paste0(tgt$depth.x,",",tgt$depth.y)
tgt$gt_RO <- tgt$depth.x
tgt$A2 <- tgt$depth.y
tgt$gt_GT_alleles <- paste0(tgt$haplo.x,"/",tgt$haplo.y)
tgt$ratio.x <- tgt$depth.x/(tgt$depth.x + tgt$depth.y)
tgt$ratio.y <- tgt$depth.y/(tgt$depth.x + tgt$depth.y)
tgt$gt <- tgt$gt_GT_alleles
tgt$minAR <- min.AR.het
tgt$maxAR <- max.AR.homo
