library(susieR)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
  stop("at least 6 arguments must be supplied:
       * genotype file name
       * expr Rd
       * pheno Rd
       * param txt
       * region file (chr, p0, p1, name, optional: pip, ifcausal, ... )
       * out file name
       * L (default 1, L in susie, optional)
       ", call.=FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir,"input_reformat.R"))

nonoverlapping <- T
PIPfilter <- 0.1
chunksize = 500000

# load genotype data
pfile <- args[1]
pfileRd <- paste0(drop_ext(pfile), ".FBM.Rd")

if (!file.exists(pfileRd)) {
  print("format converting from pgen to FBM...")
  pgen2fbm(pfile, select = NULL, scale = T, type = "double")
}

load(pfileRd)

# load expression Rd file, variable: exprres
load(args[2])
dat$expr <-  scaleRcpp(exprres$expr)

# load phenotype Rd file, variable: phenores
load(args[3])
phenores$Y <-  phenores$Y - mean(phenores$Y)

# run susie using given priors
param <- read.table(args[4], header = T, row.names = 1)

prior.gene <- param["gene.pi1", "estimated"]
prior.SNP <- param["snp.pi1", "estimated"]

regions <- read.table(args[5], stringsAsFactors = F, header =T)

if (nonoverlapping){
  stpos <- min(regions$p0) - chunksize/2
  edpos <- max(regions$p1) + chunksize/2
  p0 <- stpos + 0: floor((edpos - stpos)/chunksize) * chunksize
  p1 <- stpos + 1: ceiling((edpos - stpos)/chunksize) * chunksize

} else {
  regions <- regions[regions$PIP > PIPfilter, ]
  regions$p0 <- regions$p0 - chunksize/2
  regions$p1 <- regions$p1 + chunksize/2
}

outname <- args[6]
L = 3
if (length(args) == 7){
  L = as.numeric(args[7])
  outname <- paste0(outname, ".L", L)
}

outlist <- list()

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
  X.SNP <-  dat$G[ , idx.SNP]

  prior <- c(rep(prior.gene, dim(X.gene)[2]), rep(prior.SNP, dim(X.SNP)[2]))

  susieres[["susie"]] <- susie(cbind(X.gene, X.SNP), phenores$Y, L=L, prior_weights = prior)

  anno.gene <- cbind(colnames(X.gene), exprres$chrom[idx.gene],  exprres$p0[idx.gene])
  anno.SNP <- cbind(dat$snp[idx.SNP,], dat$chr[idx.SNP,], dat$pos[idx.SNP,])
  susieres[["anno"]] <- rbind(anno.gene, anno.SNP)

  save(susieres, file = paste(outname, name, "susieres.Rd", sep = "-"))

  outdf <- cbind(susieres[["anno"]], susieres[["susie"]]$pip)
  colnames(outdf) <- c("name", "chr", "pos", "pip")

  outlist[[name]] <- outdf[outdf[, "name"] == name, ]
}

if (nonoverlapping) {
  outgdf <- do.call(rbind, outlist)
  write.table( outgdf , file= paste(outname, "susieres.txt", sep = "-")  , row.names=F, col.names= T, sep="\t", quote = F)
}

saveRDS(outlist, file = paste(outname, "susieres.rds", sep = "-"))


