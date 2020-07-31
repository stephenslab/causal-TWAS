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
source(paste0(codedir,"mr.ash2_FBM.R"))
source(paste0(codedir,"fit_mr.ash.R"))
source(paste0(codedir,"gen_mr.ash2_output.R"))

# load genotype data
pfile <- args[1]
pfileRd <- paste0(drop_ext(pfile), ".FBM.Rd")
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

# run mr.ash2s (a simplified version of veb_boost)
# run mr.ash2s
loginfo("mr.ash2s initiation (beta: NULL, update order: expr-snp)")
mr.ash2s.fit <- mr.ash2s(snp = dat$G,
                         expr = dat$expr,
                         y= phenores$Y,
                         mr.ash.init = NULL,
                         init.order = "es", iter = 30, outname = outname, ncores = 5)
mr.ash2s_save(mr.ash2s.fit, outname)
rm(mr.ash2s.fit); gc()

loginfo("mr.ash2s initiation (beta: NULL, update order: snp-expr)")
mr.ash2s.fit <- mr.ash2s(snp = dat$G,
                         expr = dat$expr,
                         y= phenores$Y,
                         mr.ash.init = NULL,
                         init.order = "se", iter = 30, outname = outname, ncores = 5)
mr.ash2s_save(mr.ash2s.fit, outname)
rm(mr.ash2s.fit); gc()


loginfo("mr.ash2s initiation (beta: lasso, update order: expr-snp)")
mr.ash2s.fit <- mr.ash2s(snp = dat$G,
                         expr = dat$expr,
                         y= phenores$Y,
                         mr.ash.init = "lasso",
                         init.order = "es", iter = 30, outname = outname, ncores = 5)
mr.ash2s_save(mr.ash2s.fit, outname)
rm(mr.ash2s.fit); gc()

loginfo("mr.ash2s initiation (beta: lasso, iter order: snp-expr)")
mr.ash2s.fit <- mr.ash2s(snp = dat$G,
                         expr = dat$expr,
                         y= phenores$Y,
                         mr.ash.init = "lasso",
                         init.order = "es", iter.order = "se", iter = 30, outname = outname, ncores = 5)
mr.ash2s_save(mr.ash2s.fit, outname)
rm(mr.ash2s.fit); gc()


loginfo("mr.ash2s initiation (beta: lassoSNP, update order: expr-snp)")
mr.ash2s.fit <- mr.ash2s(snp = dat$G,
                         expr = dat$expr,
                         y= phenores$Y,
                         mr.ash.init = "lassoSNP",
                         init.order = "es", iter = 30, outname = outname, ncores = 5)
mr.ash2s_save(mr.ash2s.fit, outname)
rm(mr.ash2s.fit); gc()

# clean up
unlink(paste0(outname, c(".bk", ".rds")))
logwarn(str(summary(warnings())))



