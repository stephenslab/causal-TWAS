args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("5 argument must be supplied:  M.c, J.c, PVE.snp, PVE.expr, outname", call.=FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir,"input_reformat.R"))

simulate_phenotype<- function(G.scaled, weightslist, M.c, J.c, PVE.snp, PVE.expr, outname) {
  # G.scaled: scaled genotype matrix, row: sample, column:SNP. column mean= 0, sd=1
  # weightslist: each gene is a list item, providing a dataframe with columns: chr, pos, major minor, labels, alpha and rows eQTLs.
  # M.c: number of causal SNPs
  # J.c: number of causal genes

  N <- dim(G.scaled)[1]
  M <- dim(G.scaled)[2]

  # generate expr: cis expression of gene (generated using true weights), matrix, row: sample, column: gene
  exprlist <- list()
  for (gname in names(weightslist)){
     weights <- weightslist[[gname]]
     exprlist[[gname]] <- as.matrix(G.scaled[ , match(weights$labels, labels)]) %*% weights$alpha
     if (is.na(exprlist[[gname]])) exprlist[[gname]] <- 0 #TODO: eQTL not genotyped should be omited, instead of simply put gene expression to
  }

  expr <- do.call(cbind, exprlist)
  J <- dim(expr)[2]
  expr.meanvar <- mean(apply(expr[,1:1000],2,var))
  save(expr, file = paste0(paste(outname, M.c, J.c, PVE.snp, PVE.expr,sep="-"),".trueexpr.Rd"))

  set.seed(999)

  idx.cgene <- sample(1:J, J.c)
  idx.cSNP <- sample(1:M, M.c)

  sigma_theta =sqrt(PVE.snp/(M.c * (1-PVE.snp-PVE.expr)))
  sigma_gamma =sqrt(PVE.expr/(J.c * expr.meanvar * (1-PVE.snp-PVE.expr)))
  s.theta <- rnorm(M.c, mean = 0, sd = sigma_theta)
  e.gamma <- rnorm(J.c, mean = 0, sd = sigma_gamma)

  Y <- expr[,idx.cgene] %*% e.gamma + G.scaled[,idx.cSNP] %*% s.theta + rnorm(N)
  write.gemma.pheno(paste0(paste(outname, M.c, J.c, PVE.snp, PVE.expr,sep="-"), ".simulated.pheno.txt"), Y)

  PVE.snp.truth <- var(s.theta)/var(Y)
  PVE.expr.truth <- var(expr[,idx.cgene] %*% e.gamma)


  save(Y, sigma_theta, sigma_gamma, s.theta, e.gamma, idx.cSNP, idx.cgene, M.c, J.c, expr.meanvar, file=paste0(paste(outname,M.c, J.c, PVE.snp, PVE.expr,sep="-"), ".simulated_phenotype.Rd"))
  return(list("Y"=Y,"theta"=s.theta, "gamma"=e.gamma))
}

# load genotype data and scale
load("/home/simingz/causalTWAS/WTCCC/bd.RData")
X.scaled <- scaleRcpp(X)
rm(X);gc()

# load expression data
exprfile <-"/project2/mstephens/causalTWAS/simulations/simulation_WTCCC_20191127/simulated_gene_cis_expr.Rd"
load(exprfile)

# simulate phenotype
M.c <- as.numeric(args[1])
J.c <- as.numeric(args[2])
PVE.snp <- as.numeric(args[3])
PVE.expr <- as.numeric(args[4])
outname <- args[5]

phenolist <- simulate_phenotype(X.scaled, weightslist, M.c, J.c, PVE.snp, PVE.expr, outname)

print(warnings())
