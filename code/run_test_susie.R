library(mr.ash.alpha)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("4 arguments must be supplied:
       * genotype file name
       * expr Rd
       * pheno Rd
       * prior.R (will be sourced)
       * out file name", call.=FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir,"input_reformat.R"))
source(paste0(codedir,"mr.ash2.R"))

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

outname <- args[5]

prior.SNP <- 1E-3
prior.gene <- 0.1
prior <- c(rep(prior.SNP, p-1), prior.gene)
res <- susie(X,y,L=1, prior_weights=prior)






