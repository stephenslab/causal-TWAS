# will generate locuszoom type of view for the region.

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3 ) {
  stop(" at least 3 arguments must be supplied:
       * expr file (Rd file)
       * phenotype file (Rd file)
       * out file name
       * region file (optional, chrnumber<tab>start<tab>end)", call.= FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "gwas.R"))

load(args[2])
pheno <- phenores$Y

outname <- paste0(args[3], ".exprgwas.txt")

load(args[1]) # exprres
anno <- cbind(exprres$chrom, exprres$p0, exprres$p1)
colnames(anno) <- c("chr", "p0", "p1")
exprres$expr <- exprres$expr[, order(anno[,"p0"], anno[, "chr"])]
anno <- anno[order(anno[,"p0"], anno[, "chr"]),]

if (length(args) == 3){

  geno <- exprres$expr

  GWAA(geno, pheno, snpname = colnames(geno), anno = anno, outname, family = gaussian, ncore = 1, nSplits = 1, compress = T)

} else {
  regions <- read.table(args[4], stringsAsFactors = F)

  for (i in 1:nrow(regions)){
    chr <- regions[i, 1]
    p0 <- regions[i, 2]
    p1 <- regions[i, 3]
    geno <- exprres$expr[ , (exprres$chrom == chr & (exprres$p0 > p0 | exprres$p1 < p1))]

    GWAA(geno, pheno, snpname = colnames(geno), anno = anno, outname, family = gaussian, ncore = 1, nSplits = 1, outname = outname, compress = T)
  }
}





