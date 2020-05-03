library(mr.ash.alpha)
library(VEB.Boost)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("3 arguments must be supplied:
       * genotype file name
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
source(args[2])
phenolist <- simulate_phenotype(dat, mode = "snp-only",
                                pve.expr = pve.expr,
                                pve.snp = pve.snp,
                                pi_beta =  pi_beta,
                                pi_theta = pi_theta,
                                tau = tau)

# write output
outname <- args[3]
save(phenolist, file = paste0(outname, "-pheno.Rd"))
print(warnings())

# run mr.ash
mr.ash.fit <- mr.ash(dat$G, phenolist$Y)
mr.ash.fit$data <- NULL
save(mr.ash.fit, file = paste0(outname,"-mr.ash-res.Rd"))



