
library(logging)
library(tools)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("4 arguments must be supplied:
       * genotype file name (or multiple file names in a .txt file)
       * weight dir
       * param R file
       * out file name", call.=FALSE)
}

addHandler(writeToFile, file="run_simulate_data.R.log", level='DEBUG')
loginfo('script started ... ')

loginfo("input arg 1 (genotype): %s ", args[1])
loginfo("input arg 2 (weight): %s ", args[2])
loginfo("input arg 3 (param): %s ", args[3])
loginfo("input arg 4 (outname): %s ", args[4])


codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir,"input_reformat.R"))
source(paste0(codedir,"simulate_phenotype.R"))

# load genotype data
pfile <- args[1]
outname <- args[4]
if (file_ext(pfile) == "txt"){
  pfiles <- read.table(pfile, header =F, stringsAsFactors = F)[,1]
  outnames <- paste0(outname, "-B", 1:length(pfiles))
} else {
  outnames <- outname
  pfiles <- pfile
}
pfileRds <- paste0(drop_ext(pfiles), ".FBM.Rd")

loginfo("start to impute expression ...")
nlist <- list()
for ( b in 1:length(pfileRds)){
  pfileRd <- pfileRds[b]
  outname <- outnames[b]

  load(pfileRd)

  dat$G <- dat$G[] # TODO: revise other functions in this script to take FBM as input.

  # impute expression
  loginfo("generating imputed expression file: %s",  paste0(outname, "-cis-expr.Rd"))

  weight <- args[2]
  exprres <- cis_expr(dat, weight, method = "lasso", checksnps = F)
  save(exprres, file = paste0(outname, "-cis-expr.Rd"))
  nlist[[b]] <- c(ncol(dat$G), ncol(exprres$expr), sum(apply(exprres$expr, 2, var))) # genotype is scaled
}
rm(dat);gc()

# set sigma based on pve
source(args[3])

ndf <- do.call(rbind, nlist)
M <- sum(ndf[,1])
J <- sum(ndf[,2])
expr.meanvar <- 1 # expression will be scaled

J.c <- round(J * pi_beta)
M.c <- round(M * pi_theta)

sigma_beta <- sqrt(pve.expr / (J.c * expr.meanvar * (1 - pve.snp - pve.expr)))
sigma_theta <- sqrt(pve.snp / (M.c * (1 - pve.snp - pve.expr)))

loginfo("start to simulate phenotype ...")
phenoreslist <- list()
for ( b in 1:length(pfileRds)){
  pfileRd <- pfileRds[b]
  outname <- outnames[b]

  load(pfileRd)
  dat$G <- dat$G[] # TODO: revise other functions in this script to take FBM as input.
  N <- nrow(dat$G)

  load(paste0(outname, "-cis-expr.Rd"))
  dat$expr <-  scaleRcpp(exprres$expr)

  loginfo("generating phenotype file: %s", paste0(outname, "-pheno.Rd"))
  phenores <- simulate_phenotype(dat, mode = "snp-expr",
                                 pi_beta =  pi_beta,
                                 pi_theta = pi_theta,
                                 sigma_beta = sigma_beta,
                                 sigma_theta =sigma_theta)

  phenoreslist[[b]] <- phenores
}

Y <- matrix(rowSums(do.call(cbind, lapply(phenoreslist, '[[', "Y.g"))) + rnorm(N), ncol = 1)
var.y <- var(Y)

for ( b in 1:length(outnames)){
  outname <- outnames[b]
  phenores <- phenoreslist[[b]]
  phenores$Y <- Y
  phenores$param$pve.snp.truth <-  phenores$param$var.snp/var.y
  phenores$param$pve.gene.truth <- phenores$param$var.gene/var.y
  save(phenores, file = paste0(outname, "-pheno.Rd"))
}

logwarn(str(summary(warnings())))
loginfo("simulation done.")
