library(vcfR)
library(microhaplot)

run.label <- "non.humpback.samples"


# for your dataset: customize the following paths
# sam.path <- "data-raw/sam.files/GTseq.val.all"
# label.path <- "microhaplot/RunMS43.label.txt"
# vcf.path <- "vcf/GTseq.val.all.final.filtered.recode.vcf"
# app.path <- "~/Shiny/microhaplot"

sam.path <- paste0("data-raw/sam.files/", run.label)
label.path <- file.path("microhaplot", paste0(run.label, ".label.txt"))
vcf.path <- paste0("vcf/", run.label, ".POS130.vcf")
out.path <- "results/microhaplot"
app.path <- "/Users/Shared/KKMDocuments/Documents/Github.Repos/Shiny/microhaplot"

vcf <- read.vcfR(vcf.path)
#locus.ignore.file <- read.csv(paste0("microhaplot/",run.label, ".locus_annotation.csv"))

# I've prepped the data, so can just jump straight to running the Shiny app
haplo.read.tbl <- prepHaplotFiles(run.label = "non.humpback.samples.POS130",
                                  sam.path = sam.path,
                                  out.path = out.path,
                                  label.path = label.path,
                                  vcf.path = vcf.path,
                                  app.path = app.path,
                                  n.jobs = 2) # assume running on dual core






runShinyHaplot(app.path)
