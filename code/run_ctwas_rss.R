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
ld_exprfs <- vector()
zdf.gene <- NULL
for (i in 1:22){
  ld_pgenf <- ld_pgenfs[i]
  res <- impute_expr_z(zdf.SNP, ld_pgenf = ld_pgenf, weight = weight,
                       method = "lasso", outputdir = outputdir,
                       outname = outname)
  ld_exprfs[i] <- res$ld_exprf
  zdf.gene <- rbind(zdf.gene, res$zdf)
}
zdf <- rbind(zdf.SNP, zdf.gene)
save(zdf.gene, file = paste0(outputdir, "/", outname, "_zdf.Rd"))

# run ctwas_rss
ctwas_rss(zdf, ld_pgenfs, ld_exprfs, ld_regions = "EUR", ld_regions_custom = ld_regions_custom, thin = thin, outputdir = outputdir, outname = outname, ncore = ncore, prob_single = prob_single)
