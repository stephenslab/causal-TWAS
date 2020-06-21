
library(logging)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("4 arguments must be supplied:
       * genotype file name
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

# genotype plink to R data file
pfile <- args[1]
pfileRd <- paste0(pfile, ".Rd")

if (!file.exists(pfileRd)) {
  loginfo("genotype format converting ...")
  plink2Rd(pfile, pfileRd); gc()
} # pfile.traw.gz & .raw.gz should exist

# load genotype data and scale
load(pfileRd)
loginfo("genotype scaling ...")
dat$G <- scaleRcpp(dat$G)
gc()

# simulate phenotype
outname <- args[4]
loginfo("generating imputed expression file: %s",  paste0(outname, "-cis-expr.Rd"))

weight <- args[2]
exprres <- cis_expr(dat, weight, method = "best", checksnps = F)
save(exprres, file = paste0(outname, "-cis-expr.Rd"))

# load(paste0(outname, "-cis-expr.Rd"))
loginfo("expr scaling ...")
dat$expr <-  scaleRcpp(exprres$expr)

source(args[3])
loginfo("generating phenotype file: %s", paste0(outname, "-pheno.Rd"))
phenores <- simulate_phenotype(dat, mode = "snp-expr",
                                pve.expr = pve.expr,
                                pve.snp = pve.snp,
                                pi_beta =  pi_beta,
                                pi_theta = pi_theta,
                                tau = tau)

# write output
save(phenores, file = paste0(outname, "-pheno.Rd"))

logwarn(str(summary(warnings())))
loginfo("simulation done.")
