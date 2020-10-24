# will generate locuszoom type of view for the region.
library(tools)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3 ) {
  stop(" at least 3 arguments must be supplied:
       * expr file (Rd file or txt file)
       * phenotype file (Rd file or txt file)
       * out file name
       * region file (optional, chrnumber<tab>start<tab>end)", call.= FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "gwas.R"))
source(paste0(codedir,"input_reformat.R"))

pfile <- args[1]
if (file_ext(pfile) == "txt"){
  pfiles <- read.table(pfile, header =F, stringsAsFactors = F)[,1]
  outnames <- paste0(args[3], "-B", 1:length(pfiles), ".exprgwas.txt")
  combine <- T
} else {
  pfiles <- pfile
  outnames <- paste0(args[3], ".exprgwas.txt")
  combine <- F
}

phenofile <- args[2]
load(phenofile)
pheno <- phenores$Y

for (b in 1:length(pfiles)){

  load(pfiles[b]) # exprres

  outname <- outnames[b]

  anno <- cbind(exprres$chrom, exprres$p0, exprres$p1)
  colnames(anno) <- c("chr", "p0", "p1")

  if (length(args) == 3){

    geno <- exprres$expr

    GWAA(geno, pheno, snpname = exprres$gnames, anno = anno, outname, family = gaussian, ncore = 1, nSplits = 1, compress = T)

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
}

# combine all results into 1 file for batch input
outdflist <- list()
if (isTRUE(combine)) {
  for (b in 1:length(pfiles)){
      outdflist[[b]] <- read.table(paste0(outnames[b], ".gz"), header = T, stringsAsFactors = F, comment.char = "")
  }
  outdf <- do.call(rbind, outdflist)

  outname <-  paste0(args[3], ".exprgwas.txt")
  write.table(outdf, file = outname , row.names=F, col.names=T, sep="\t", quote = F)
  system(paste0("bgzip ", outname))
}



