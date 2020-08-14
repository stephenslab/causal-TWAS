library(logging)

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

# load genotype data
pfile <- args[1]
pfileRd <- paste0(pfile, ".unscaled.FBM.Rd")
if (file.exists(pfileRd)){
  load(pfileRd)
} else {
  load(paste0(pfile, ".Rd"))
}


# load expression Rd file, variable: exprres
load(args[2])
gc()

# load pheno Rd file, variable: pheno
load(args[3])

outname <- args[4]

# get ld information before genotype scaling
get_ld(outname)

loginfo("ld for genes and causal genes/SNPs done.")



