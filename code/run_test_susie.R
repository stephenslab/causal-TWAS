library(susieR)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
  stop("4 arguments must be supplied:
       * genotype file name
       * expr Rd
       * pheno Rd
       * prior.R (will be sourced)
       * region file (optional, chrnumber<tab>start<tab>end)
       * out file name", call.=FALSE)
}

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

source(args[4]) # priors

outname <- args[6]

regions <- read.table(args[5], stringsAsFactors = F)

for (i in 1:nrow(regions)){
  chr <- regions[i, 1]
  p0 <- regions[i, 2]
  p1 <- regions[i, 3]
  X.gene <- exprres$expr[ , (exprres$chrom == chr & (exprres$p0 > p0 | exprres$p1 < p1))]
  X.SNP <-  dat$G[ , (dat$chr == chr & dat$pos > p0 & dat$pos < p1)]

  prior <- c(rep(prior.gene, dim(X.gene)[2]), rep(prior.SNP, dim(X.SNP)[2]))

  res <- susie(cbind(X.gene, X.SNP), phenores$Y, L=1, prior_weights = prior)
  save(res, file = paste(outname,chr,p0,p1,"susieres.Rd", sep = "-"))
}





