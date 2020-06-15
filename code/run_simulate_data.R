
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("4 arguments must be supplied:
       * genotype file name
       * weight dir
       * param R file
       * out file name", call.=FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir,"input_reformat.R"))
source(paste0(codedir,"simulate_phenotype.R"))

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

# simulate phenotype
weight <- args[2]
exprres <- cis_expr(dat, weight, method = "best", checksnps = F)

outname <- args[4]
save(exprres, file = paste0(outname, "-cis-expr.Rd"))
# load(paste0(outname, "-cis-expr.Rd"))
dat$expr <-  scaleRcpp(exprres$expr)

source(args[3])
phenores <- simulate_phenotype(dat, mode = "snp-expr",
                                pve.expr = pve.expr,
                                pve.snp = pve.snp,
                                pi_beta =  pi_beta,
                                pi_theta = pi_theta,
                                tau = tau)

# write output
save(phenores, file = paste0(outname, "-pheno.Rd"))
print(warnings())
