library(susieR)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
  stop("at least 6 arguments must be supplied:
       * genotype file name
       * expr Rd
       * pheno Rd
       * param txt
       * region file (chr, p0, p1, name, optional: pip, ifcausal, ... )
       * out file name
       * L (default 1, L in susie, optional)", call.=FALSE)
}

# pay attention to filter <- T/F

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir,"input_reformat.R"))


# genotype plink to R data file
pfile <- args[1]
pfileRd <- paste0(pfile, ".Rd")

if (!file.exists(pfileRd)) {
  print("format converting ...")
  plink2Rd(pfile, pfileRd); gc()
} # pfile.traw.gz & .raw.gz should exist

# load genotype data and scale
load(pfileRd)
print("scaling ...")
dat$G <- scaleRcpp(dat$G)

# load expression Rd file, variable: exprres
load(args[2])
dat$expr <-  scaleRcpp(exprres$expr)

# load phenotype Rd file, variable: phenores
load(args[3])

# run susie using given priors
param <- read.table(args[4], header = T, row.names = 1)

prior.gene <- param["gene.pi1", "estimated"]
prior.SNP <- param["snp.pi1", "estimated"]

regions <- read.table(args[5], stringsAsFactors = F, header =T)

filter <- T
if (filter) {
  regions <- regions[regions$ifcausal == 1 | regions$PIP > 0.3, ]
  flank = 500000
  regions$p0 <- regions$p0 - flank
  regions$p1 <- regions$p1 + flank
}

outname <- args[6]
L = 1
if (length(args) == 7){
  L = as.numeric(args[7])
  outname <- paste0(outname, ".L", L)
}

for (i in 1:nrow(regions)){
  chr <- regions[i, "chr"]
  p0 <- regions[i, "p0"]
  p1 <- regions[i, "p1"]
  name <- regions[i, "name"]

  susieres <- list()
  print(name)

  idx.gene <- exprres$chrom == chr & exprres$p0 > p0 & exprres$p1 < p1
  idx.SNP <- dat$chr == chr & dat$pos > p0 & dat$pos < p1

  X.gene <- exprres$expr[ , idx.gene, drop = F]
  X.SNP <-  dat$G[ , idx.SNP, drop = F]

  prior <- c(rep(prior.gene, dim(X.gene)[2]), rep(prior.SNP, dim(X.SNP)[2]))

  susieres[["susie"]] <- susie(cbind(X.gene, X.SNP), phenores$Y, L=L, prior_weights = prior)

  anno.gene <- cbind(colnames(X.gene), exprres$chrom[idx.gene],  exprres$p0[idx.gene])
  anno.SNP <- cbind(dat$snp[idx.SNP,], dat$chr[idx.SNP,], dat$pos[idx.SNP,])
  susieres[["anno"]] <- rbind(anno.gene, anno.SNP)

  save(susieres, file = paste(outname, name, "susieres.Rd", sep = "-"))

  outdf <- cbind(susieres[["anno"]], susieres[["susie"]]$pip)
  colnames(outdf) <- c("name", "chr", "pos", "pip")

  write.table( outdf , file= paste(outname, name, "susieres.txt", sep = "-")  , row.names=F, col.names= T, sep="\t", quote = F)
}




