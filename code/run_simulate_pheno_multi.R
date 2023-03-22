library(ctwas)
library(logging)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  stop("5 arguments must be supplied:
       * genotype file name (multiple file names in a .txt file)
       * imputed expr file (multiple file names in a .txt file)
       * param R file
       * out file name
       * outputdir", call.=FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
#source(paste0(codedir, "stats_func.R"))
source(paste0(codedir, "simulate_phenotype.R"))

pgenfs <- read.table(args[1], header = F, stringsAsFactors = F)[,1]
exprfs <- read.table(args[2], header = F, stringsAsFactors = F)[,1]

source(args[3])
outname <- args[4]
outputdir <- args[5]

loginfo("simulating phenotype for %s", outname)
set.seed(SED)

phenores <- simulate_phenotype_multi(pgenfs = pgenfs,
                                     exprfs = exprfs,
                                     weight_1 = weight_1,
                                     pve_expr_1 = pve_expr_1,
                                     pi_beta_1 = pi_beta_1,
                                     weight_2 = weight_2,
                                     pve_expr_2 = pve_expr_2,
                                     pi_beta_2 = pi_beta_2,
                                     pve_snp = pve_snp,
                                     pi_theta = pi_theta)

save(phenores, file = paste0(outputdir,"/", outname, "-pheno.Rd"))

