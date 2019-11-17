args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  M.c, J.c, PVE.snp, PVE.expr
  stop("5 argument must be supplied:  M.c, J.c, PVE.snp, PVE.expr, outname", call.=FALSE)
}

source("stats_func.R")

simulate_phenotype<- function(G.scaled, chr, pos, expr, M.c, J.c, PVE.snp, PVE.expr) {
  # G.scaled: scaled genotype matrix, row: sample, column:SNP. column mean= 0, sd=1
  # expr: cis expression of gene, matrix, row: sample, column: imputed gene
  # M.c: number of causal SNPs
  # J.c: number of causal genes

  N <- dim(G.scaled)[1]
  M <- dim(G.scaled)[2]
  J <- dim(expr)[2]
  expr.meanvar <- mean(apply(expr[,1:100],2,var))

  idx.gene <- sample(1:J, J.c)
  idx.SNP <- sample(1:M, M.c)


  Y <- expr * gamma + G.scaled[, idx.SNP] %*% theta + rnorm(N)
  n.eQTL.mean= mean(n.eQTL.range)
  sigma_alpha <- sqrt(h2.eQTL/(n.eQTL.mean - n.eQTL.mean*h2.eQTL))

  N <- dim(G)[1]
  M <- dim(G)[2]
  ga <- read.table(gafile,header = F,stringsAsFactors = F)
  simulate_eQTL_geno <- function(i){
    print(i);gc()
    idx <- c(1:M)[paste0("chr",chr) == ga[i,2] & pos > ga[i,4] & pos < ga[i,5]]
    if (length(idx) < max(n.eQTL.range)) {
      cat(paste0(ga[i,6], " not included due to too few SNPs\n"))
      return(NULL) # will not simulate if too few SNPs in region
    } else {
      idx.eQTL <- sample(idx, sample(n.eQTL.range,1))
      return(G[,idx.eQTL])
    }
  }
  ga.eQTL.geno <-lapply(1:dim(ga)[1], simulate_eQTL_geno)
  names(ga.eQTL.geno) <- ga$V6
  save(ga.eQTL.geno, file="eQTL_genotype.Rd")

  # simulate gene expression data
  exprlist <- list()
  for (i in 1:length(ga.eQTL.geno)){
    if (!is.null(ga.eQTL.geno[[i]])){
      gname <- names(ga.eQTL.geno[i])
      geno <- scale(as.matrix(ga.eQTL.geno[[i]]))
      alpha <- rnorm(dim(geno)[2], mean=0, sd=sigma_alpha)
      expr <- geno %*% alpha + rnorm(N, sd=1)
      exprlist[[gname]] <-expr
    }
  }
  gnames <- names(exprlist)
  expr <- do.call(cbind, exprlist) # rows are samples, columns are genes.
  J <- length(exprlist)
  save(gnames, expr, sigma_alpha, J ,file="simulated_gene_cis_expr.Rd")
}

# load genotype data and scale
load("/home/simingz/causalTWAS/WTCCC/bd.RData")
X.scaled <- scaleRcpp(X); gc()

# load expression data
exprfile <-"/project2/mstephens/causalTWAS/simulations/simulation_WTCCC_20191111/simulated_gene_cis_expr.Rd"
load(exprfile)

# simulate phenotype
M.c <-
J.c <-
PVE.snp <-
PVE.expr <-
simulate_phenotype(X.scaled, chr, pos, expr, M.c, J.c, PVE.snp, PVE.expr)



