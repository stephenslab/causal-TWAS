# will generate locuszoom type of view for the region.
library(tools)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3 ) {
  stop(" at least 3 arguments must be supplied:
       * geno file (Rd file or txt file)
       * phenotype file (Rd file or txt file)
       * out file name
       * region file (optional, chrnumber<tab>start<tab>end)", call.= FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "gwas.R"))
source(paste0(codedir,"input_reformat.R"))

pfile <- args[1]
if (file_ext(pfile) == "txt"){
  pfiles <- read.table(pfile, header = F, stringsAsFactors = F)[,1]
  phenofiles <- read.table(args[2], header = F, stringsAsFactors = F)[,1]
  outnames <- paste0(args[3], "-B", 1:length(pfiles), ".snpgwas.txt")
  combine <- T
} else {
  pfiles <- pfile
  phenofiles <- args[2]
  outnames <- paste0(args[3], ".snpgwas.txt")
  combine <- F
}
pfileRds <- paste0(drop_ext(pfiles), ".FBM.Rd")


for (b in 1:length(pfiles)){
  load(phenofiles[b])
  pheno <- phenores$Y

  outname <- outnames[b]

  # load genotype data
  load(pfileRds[b])
  dat$G <- dat$G[]

  snpname = dat$snp[,1]
  anno <- cbind(dat$chr, dat$pos, dat$pos) # EPACT format

  if (length(args) == 3){

    geno <- dat$G

    GWAA(geno, pheno, snpname = snpname, anno = anno, outname, family = gaussian, ncore = 3, nSplits = 100, compress = T)

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
}

# combine all results into 1 file for batch input
outdflist <- list()
if (isTRUE(combine)) {
  for (b in 1:length(pfiles)){
    outdflist[[b]] <- read.table(paste0(outnames[b], ".gz"), header = T, stringsAsFactors = F, comment.char = "")
  }
  outdf <- do.call(rbind, outdflist)

  outname <-  paste0(args[3], ".snpgwas.txt")
  write.table(outdf, file = outname , row.names=F, col.names=T, sep="\t", quote = F)

  system(paste0("bgzip ", outname))
  system(paste0("tabix -p bed ", outname, ".gz"))

}



