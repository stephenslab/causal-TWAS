library(data.table)
library(logging)
library(tools)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("4 arguments must be supplied:
       * genotype file name (or multiple file names in a .txt file)
       * impute expr Rd file (or multiple file names in a .txt file)
       * param R file
       * out file name", call.=FALSE)
}

addHandler(writeToFile, file="run_simulate_phenotype.R.log", level='DEBUG')
loginfo('script started ... ')

loginfo("input arg 1 (genotype): %s ", args[1])
loginfo("input arg 2 (imputed expression): %s ", args[2])
loginfo("input arg 3 (param): %s ", args[3])
loginfo("input arg 4 (outname): %s ", args[4])

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir,"input_reformat.R"))
source(paste0(codedir,"simulate_phenotype.R"))


outname <- args[4]

# prepare input files. Genotype and expression are already scaled.
pfile <- args[1]

if (file_ext(pfile) == "txt"){
  pfiles <- read.table(pfile, header = F, stringsAsFactors = F)[,1]
  efiles <- read.table(args[2], header =F, stringsAsFactors = F)[,1]
} else {
  pfiles <- pfile
  efiles <- args[2]
}
pfileRds <- paste0(pfiles, ".FBM.Rd")


# set sigma based on pve
source(args[3])
set.seed(SED)

nlist <- list()
for ( b in 1:length(pfileRds)){
  load(pfileRds[b])
  load(efiles[b])
  nlist[[b]] <- c(ncol(dat$G), ncol(exprres$expr))
}
ndf <- do.call(rbind, nlist)
M <- sum(ndf[,1])
J <- sum(ndf[,2])
N <- nrow(dat$G)
expr.meanvar <- 1 # expression has been scaled

J.c <- round(J * pi_beta)
M.c <- round(M * pi_theta)

sigma_beta <- sqrt(pve.expr / (J.c * expr.meanvar * (1 - pve.snp - pve.expr)))
sigma_theta <- sqrt(pve.snp / (M.c * (1 - pve.snp - pve.expr)))

# Simulate phenotype
loginfo("start to simulate phenotype ...")
phenores <- list("batch" = list())
for ( b in 1:length(pfileRds)){
  # load genotype data
  pfileRd <- pfileRds[b]
  load(pfileRd)

  # load expression data
  load(efiles[b])
  dat$expr <-  exprres$expr

  loginfo("generating phenotype file: %s", paste0(outname, "-pheno.Rd"))
  phenores[["batch"]][[b]] <- simulate_phenotype(dat, mode = "snp-expr",
                                 pi_beta =  pi_beta,
                                 pi_theta = pi_theta,
                                 sigma_beta = sigma_beta,
                                 sigma_theta =sigma_theta)
}

Y <- matrix(rowSums(do.call(cbind, lapply(phenores[["batch"]], '[[', "Y.g"))) + rnorm(N), ncol = 1)
phenores$Y <- Y
var.y <- var(Y)
phenores$param$pve.snp.truth <-  sum(unlist(lapply(lapply(phenores[["batch"]], '[[', "param"), '[[', "var.snp")))/var.y
phenores$param$pve.gene.truth <- sum(unlist(lapply(lapply(phenores[["batch"]], '[[', "param"), '[[', "var.gene")))/var.y

save(phenores, file = paste0(outname, "-pheno.Rd"))

logwarn(str(summary(warnings())))
loginfo("simulation done.")
