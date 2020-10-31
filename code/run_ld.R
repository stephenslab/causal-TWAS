library(logging)
library(tools)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("3 arguments must be supplied:
       * genotype file name
       * expr Rd
       * pheno Rd
       * out file name", call.=FALSE)
}

addHandler(writeToFile, file="run_ld.R.log", level='DEBUG')
loginfo('script started ... ')

loginfo("input arg 1 (genotype): %s ", args[1])
loginfo("input arg 2 (expr): %s ", args[2])
loginfo("input arg 3 (expr): %s ", args[3])
loginfo("input arg 4 (outname): %s ", args[4])


codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir,"ld.R"))

pfile <- args[1]
if (file_ext(pfile) == "txt"){
  pfiles <- read.table(pfile, header =F, stringsAsFactors = F)[,1]
  efiles <- read.table(args[2], header =F, stringsAsFactors = F)[,1]
} else {
  pfiles <- pfile
  efiles <- args[2]
}
pfileRds <- paste0(pfiles, ".unscaled.FBM.Rd")


# load pheno Rd file, variable: pheno
load(args[3])

outname <- args[4]

load(pfileRds[1])
load(efiles[1])

geno.t <- list()

geno.t[[1]] <- .eqtl_geno("MRPS21")

name <- "rs4839015"
temp <- dat$G[, dat$snp == name, drop = F]
colnames(temp) <- name
geno.t[[name]] <- temp

rsqlist <- list()

for (b in 1:length(pfiles)){
  # load genotype data
  load(pfileRds[b])

  # load expression Rd file, variable: exprres
  load(efiles[b])

  causalg.eqtl <-lapply(exprres$gnames[phenores$batch[[b]]$param$idx.cgene], .eqtl_geno)
  cg <- do.call(cbind, causalg.eqtl)

  csnp <- dat$G[ ,phenores$batch[[b]]$param$idx.cSNP]
  colnames(csnp) <- dat$snp[phenores$batch[[b]]$param$idx.cSNP, 1]

  rsqlist[[b]] <- lapply(geno.t, ld_max, y = cbind(cg, csnp), stats = "R.squared")
}

loginfo("ld for genes and causal genes/SNPs done.")



