library(ctwas)

args = commandArgs(trailingOnly=TRUE)
if (length(args) <5) {
  stop(" 6 arguments:
       * zscore file name
       * ld genotype file (multiple file names in a .txt file)
       * weight
       * config file name (.R)
       * out file name
       * outputdir", call.=FALSE)
}

zdf.SNP <- data.table::fread(args[1], header = T)
data.table::setnames(zdf.SNP, old =  't.value', new =  'z')
zdf.SNP <- zdf.SNP[, c("id", "z")]
ld_pgenfs <- read.table(args[2], header = F, stringsAsFactors = F)[,1]
weight <- args[3]
outname <- args[5]
outputdir <- args[6]


ld_regions_custom <- "/home/simingz/ctwas/inst/extdata/example_regions.bed"
thin <- 1
ncore <- 1
prob_single <- 0.8

source(args[4]) # config

# get gene z score
ld_exprfs <- paste0(outputdir, "/", outname, "_chr", 1:22, ".expr.gz")
load(file = paste0(outputdir, "/", outname, "_zdf.Rd"))
zdf <- rbind(zdf.SNP, zdf.gene)
rm(zdf.SNP, zdf.gene); gc()
save.image("temp.Rd")
print("done")
# run ctwas_rss
ctwas_rss(zdf, ld_pgenfs, ld_exprfs, ld_regions = "EUR", ld_regions_custom = ld_regions_custom, thin = thin, outputdir = outputdir, outname = outname, ncore = ncore, prob_single = prob_single)
