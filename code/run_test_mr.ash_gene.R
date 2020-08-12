library(logging)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("3 arguments must be supplied:
       * expr Rd
       * pheno Rd
       * out file name", call.=FALSE)
}

addHandler(writeToFile, file="run_mr.ash2s_gene.R.log", level='DEBUG')
loginfo('script started ... ')

loginfo("input arg 1 (expr): %s ", args[1])
loginfo("input arg 2 (pheno): %s ", args[2])
loginfo("input arg 3 (outname): %s ", args[3])


codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir, "fit_mr.ash.R"))

# load expression Rd file, variable: exprres
load(args[1])
loginfo("expr scaling ...")
dat$expr <-  scaleRcpp(exprres$expr)
gc()

# load phenotype Rd file, variable: phenores
load(args[2])
phenores$Y <-  phenores$Y - mean(phenores$Y)

outname <- args[3]

# run mr.ash2s
loginfo("run mr.ash for gene expression ...")
mr.ash.fit <- mr.ash( X = dat$expr,
                      y = phenores$Y)
save(mr.ash.fit, file = paste0(outname, ".mr.ashGene.Rd"))

