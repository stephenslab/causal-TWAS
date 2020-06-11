library(mr.ash.alpha)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("4 arguments must be supplied:
       * genotype file name
       * expr Rd
       * pheno Rd
       * out file name", call.=FALSE)
}

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

# load genotype data and scale
load(pfileRd)
print("scaling ...")
dat$G <- scaleRcpp(dat$G)

# load expression Rd file, variable: exprres
load(args[2])
dat$expr <-  scaleRcpp(exprres$expr)

# load phenotype Rd file, variable: phenores
load(args[3])

# run mr.ash2s (a simplified version of veb_boost)

outname <- args[4]

## start with expression
mr.ash2s.fit <- mr.ash2s(dat$expr, dat$G, phenores$Y, iter = 30)

## save output
g.fit <- mr.ash2s.fit$fit1
s.fit <-  mr.ash2s.fit$fit2

gen_mr.ash2_output(g.fit, s.fit, paste0(outname,"-mr.ash2s.expr-res"))

mr.ash2s.fit$fit2$data$X <- NULL
save(mr.ash2s.fit, file = paste0(outname,"-mr.ash2s.expr-res.Rd"))

## start with snp
mr.ash2s.fit <- mr.ash2s(dat$G, dat$expr, phenores$Y, iter = 30)

## save output
g.fit <- mr.ash2s.fit$fit2
s.fit <-  mr.ash2s.fit$fit1

gen_mr.ash2_output(g.fit, s.fit, paste0(outname,"-mr.ash2s.snp-res"))

mr.ash2s.fit$fit1$data$X <- NULL
save(mr.ash2s.fit, file = paste0(outname,"-mr.ash2s.snp-res.Rd"))





