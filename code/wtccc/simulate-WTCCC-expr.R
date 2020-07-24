# Simulation function
codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir,"input_reformat.R"))
library(data.table)

simulate_expr<- function(G, gafile, n.eQTL.range=1:5, h2.eQTL=0.2) {
  # G: genotype matrix, row: sample, column:SNP
  # gafile: gene annotation file, define the start and end of gene
  # n.eQTL.range: ranges of number of eQTL SNPs, will be sampled randomly for each gene
  # h2.eQTL: h2 of eQTL.
    n.eQTL.mean= mean(n.eQTL.range)
    sigma_alpha <- sqrt(h2.eQTL/(n.eQTL.mean - n.eQTL.mean*h2.eQTL))

    N <- dim(G)[1]
    M <- dim(G)[2]

    ga <- read.table(gafile, header = F, stringsAsFactors = F)

    set.seed(345)

    simulate_eQTL <- function(i){
      print(i);gc()
      idx <- c(1:M)[paste0("chr", chr) == ga[i,2] & pos > ga[i,4] & pos < ga[i,5]]
      if (length(idx) < max(n.eQTL.range)) {
        cat(paste0(ga[i,6], " not included due to too few SNPs\n"))
        return(NULL) # will not simulate if too few SNPs in region
      } else {
        idx.eQTL <- sample(idx, sample(n.eQTL.range,1))
        return(idx.eQTL)
      }
    }

    #ga.eQTL.idx <- lapply(1:dim(ga)[1], simulate_eQTL) # this step takes around 2 hours
    # names(ga.eQTL.idx) <- ga$V6
    # save(ga.eQTL.idx, file="eQTL_idx.Rd")
    load("eQTL_idx.Rd")

    # simulate gene expression data
    exprlist <- list()
    weightslist <- list()
    for (i in 1:length(ga.eQTL.idx)){
      print(i)
      idx <- ga.eQTL.idx[[i]]
      if (!is.null(idx)){
        geno <- G[,idx]
        alpha <- rnorm(length(idx), mean=0, sd=sigma_alpha)
        expr <- scale(as.matrix(geno)) %*% alpha + rnorm(N, sd=1)
        if (!is.na(expr[1])){
          gname <- names(ga.eQTL.idx[i])
          exprlist[[gname]] <- expr
          weightslist[[gname]] <- data.frame("chr" = chr[idx], "pos" = pos[idx], "major" = major[idx], "minor" = minor[idx], "labels" = labels[idx], "alpha" = alpha)
        }
      }
    }
    gnames <- names(exprlist)
    expr <- do.call(cbind, exprlist) # rows are samples, columns are genes.
    J <- length(exprlist)
    save(G, labels, chr, pos, gnames, expr, sigma_alpha, J, weightslist, file="simulated_gene_cis_expr2.Rd")
}

load("/home/simingz/causalTWAS/WTCCC/cad.RData")
X <- X[y==0,][1:500,]
gc()
# simulate expression
gafile <- "/home/simingz/causalTWAS/annotation/refGene_processed_5e+05.txt"
simulate_expr(X, gafile, n.eQTL.range=1:5, h2.eQTL=0.2)










