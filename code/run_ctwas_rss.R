library(ctwas)

args = commandArgs(trailingOnly=TRUE)
if (length(args) <7) {
  stop(" 7 arguments:
       * zscore file name
       * ld genotype file (multiple file names in a .txt file)
       * weight
       * config file name (.R)
       * out expr z file name
       * out file name
       * outputdir", call.=FALSE)
}

zdf.SNP <- data.table::fread(args[1], header = T)
data.table::setnames(zdf.SNP, old =  't.value', new =  'z')
zdf.SNP <- zdf.SNP[, c("id", "z")]
ld_pgenfs <- read.table(args[2], header = F, stringsAsFactors = F)[,1]
weight <- args[3]
outname.e <- args[5]
outname <- args[6]
outputdir <- args[7]


ld_regions_custom <- "/home/simingz/ctwas/inst/extdata/example_regions.bed"
thin <- 1
ncore <- 1
prob_single <- 0.8

source(args[4]) # config

# get gene z score
if (file.exists(paste0(outputdir, "/", outname.e, "_zdf.Rd"))){
  ld_exprfs <- paste0(outputdir, "/", outname.e, "_chr", 1:22, ".expr.gz")
  load(file = paste0(outputdir, "/", outname.e, "_zdf.Rd"))
} else {
  ld_exprfs <- vector()
  zdf.gene <- NULL
  for (i in 1:22){
    ld_pgenf <- ld_pgenfs[i]
    res <- impute_expr_z(zdf.SNP, ld_pgenf = ld_pgenf, weight = weight,
                         method = "lasso", outputdir = outputdir,
                         outname = outname.e)
    ld_exprfs[i] <- res$ld_exprf
    zdf.gene <- rbind(zdf.gene, res$zdf)
  }
  save(zdf.gene, file = paste0(outputdir, "/", outname.e, "_zdf.Rd"))
}

zdf <- rbind(zdf.SNP, zdf.gene)
rm(zdf.SNP, zdf.gene); gc()

# run ctwas_rss
ctwas_rss(zdf, ld_pgenfs, ld_exprfs, ld_regions = "EUR", ld_regions_custom = ld_regions_custom, thin = thin, outputdir = outputdir, outname = outname, ncore = ncore, prob_single = prob_single)
