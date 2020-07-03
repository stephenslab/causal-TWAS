library(logging)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("4 arguments must be supplied:
       * genotype file name
       * expr Rd
       * pheno Rd
       * out file name", call.=FALSE)
}

addHandler(writeToFile, file="run_mr.ash2s.R.log", level='DEBUG')
loginfo('script started ... ')

loginfo("input arg 1 (genotype): %s ", args[1])
loginfo("input arg 2 (expr): %s ", args[2])
loginfo("input arg 3 (pheno): %s ", args[3])
loginfo("input arg 4 (outname): %s ", args[4])


codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir,"input_reformat.R"))
source(paste0(codedir,"simulate_phenotype.R"))
source(paste0(codedir,"mr.ash2.R"))
source(paste0(codedir,"fit_mr.ash.R"))
source(paste0(codedir,"gen_mr.ash2_output.R"))

# genotype plink to R data file
pfile <- args[1]
pfileRd <- paste0(pfile, ".Rd")

if (!file.exists(pfileRd)) {
  print("format converting ...")
  plink2Rd(pfile, pfileRd); gc()
} # pfile.traw.gz & .raw.gz should exist

# load genotype data
load(pfileRd)

# load expression Rd file, variable: exprres
load(args[2])
loginfo("expr scaling ...")
dat$expr <-  scaleRcpp(exprres$expr)
gc()

# load phenotype Rd file, variable: phenores
load(args[3])
phenores$Y <-  phenores$Y - mean(phenores$Y)

outname <- args[4]

# get ld information before genotype scaling
get_ld(outname)

# scale genotype
loginfo("genotype scaling ...")
dat$G <- scaleRcpp(dat$G)
gc()

# run mr.ash2s (a simplified version of veb_boost)

## mr.ash.init = lasso for gene and snp
loginfo("mr.ash2s init from lasso ...")
mr.ash2s.fit <- mr.ash2s(dat$expr, dat$G, phenores$Y, iter = 30, mr.ash.init = "lasso")

## save output
g.fit <- mr.ash2s.fit$fit1
s.fit <-  mr.ash2s.fit$fit2

gen_mr.ash2_output(g.fit, s.fit, paste0(outname,"-mr.ash2s.lasso"))

mr.ash2s.fit$fit2$data$X <- NULL
save(mr.ash2s.fit, file = paste0(outname,"-mr.ash2s.lasso.Rd"))

logwarn(str(summary(warnings())))
gc()
loginfo("mr.ash2s init from lasso done.")


## mr.ash.init = lasso for snp
loginfo("mr.ash2s init from lasso SNP ...")
mr.ash2s.fit <- mr.ash2s(dat$expr, dat$G, phenores$Y, iter = 30, mr.ash.init = "lassoX2")

## save output
g.fit <- mr.ash2s.fit$fit1
s.fit <-  mr.ash2s.fit$fit2

gen_mr.ash2_output(g.fit, s.fit, paste0(outname,"-mr.ash2s.lassoSNP"))

mr.ash2s.fit$fit2$data$X <- NULL
save(mr.ash2s.fit, file = paste0(outname,"-mr.ash2s.lassoSNP.Rd"))

logwarn(str(summary(warnings())))
gc()
loginfo("mr.ash2s init from lasso SNP done.")



