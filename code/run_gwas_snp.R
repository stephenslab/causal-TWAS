# will generate locuszoom type of view for the region.

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3 ) {
  stop(" at least 3 arguments must be supplied:
       * geno file (Rd file)
       * phenotype file (Rd file)
       * out file name
       * region file (optional, chrnumber<tab>start<tab>end)", call.= FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "gwas.R"))

load(args[2])
pheno <- phenores$Y

outname <- paste0(args[3], ".snpgwas.txt")

load(args[1]) # dat
snpname = dat$snp[,1]
anno <- cbind(dat$chr, dat$pos, dat$pos) # EPACT format

if (length(args) == 3){

  geno <- dat$G

  GWAA(geno, pheno, snpname = snpname, anno = anno, outname, family = gaussian, ncore = 5, nSplits = 10, compress = T)

} else {
  regions <- read.table(args[4], stringsAsFactors = F)

  for (i in 1:nrow(regions)){
    chr <- regions[i, 1]
    p0 <- regions[i, 2]
    p1 <- regions[i, 3]
    geno <- dat$G[ , (dat$chr == chr & dat$pos > p0 & dat$pos < p1)]

    GWAA(geno, pheno, snpname = snpname, anno = anno, outname, family = gaussian, ncore = 1, nSplits = 1, outname = outname, compress = T)
  }
}





