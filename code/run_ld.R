library(logging)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("3 arguments must be supplied:
       * genotype file name
       * expr Rd
       * out file name", call.=FALSE)
}

addHandler(writeToFile, file="run_ld.R.log", level='DEBUG')
loginfo('script started ... ')

loginfo("input arg 1 (genotype): %s ", args[1])
loginfo("input arg 2 (expr): %s ", args[2])
loginfo("input arg 3 (outname): %s ", args[3])


codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir,"gen_mr.ash2_output.R"))

# load genotype data
pfile <- args[1]
pfileRd <- paste0(pfile, ".Rd")
load(pfileRd)

# load expression Rd file, variable: exprres
load(args[2])
gc()

outname <- args[3]

# get ld information before genotype scaling
get_ld(outname)

loginfo("ld for genes and causal genes/SNPs done.")



