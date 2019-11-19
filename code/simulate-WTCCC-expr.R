# Simulation function
codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir,"input_reformat.R"))

simulate_expr<- function(G, chr, pos, gafile, n.eQTL.range=1:5, h2.eQTL=0.2) {
  # G: genotype matrix, row: sample, column:SNP
  # gafile: gene annotation file, define the start and end of gene
  # n.eQTL.range: ranges of number of eQTL SNPs, will be sampled randomly for each gene
  # h2.eQTL: h2 of eQTL.
    n.eQTL.mean= mean(n.eQTL.range)
    sigma_alpha <- sqrt(h2.eQTL/(n.eQTL.mean - n.eQTL.mean*h2.eQTL))

    N <- dim(G)[1]
    M <- dim(G)[2]

    ga <- read.table(gafile,header = F,stringsAsFactors = F)

    set.seed(345)

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

    ga.eQTL.geno <-lapply(1:dim(ga)[1], simulate_eQTL_geno) # this step takes around 2 hours
    names(ga.eQTL.geno) <- ga$V6
    save(ga.eQTL.geno, file="eQTL_genotype.Rd")

    # simulate gene expression data
    exprlist <- list()
    for (i in 1:length(ga.eQTL.geno)){
      print(i)
      if (!is.null(ga.eQTL.geno[[i]])){
        gname <- names(ga.eQTL.geno[i])
        alpha <- rnorm(dim(ga.eQTL.geno[[i]])[2], mean=0, sd=sigma_alpha)
        expr <- scale(ga.eQTL.geno[[i]]) %*% alpha + rnorm(N, sd=1)
        exprlist[[gname]] <-expr
      }
    }
    gnames <- names(exprlist)
    expr <- do.call(cbind, exprlist) # rows are samples, columns are genes.
    J <- length(exprlist)
    save(gnames, expr, sigma_alpha, J ,file="simulated_gene_cis_expr.Rd")

}

# get expression heritability, this can be used to verify simulation.
expr_h2 <- function(){

}

load("/home/simingz/causalTWAS/WTCCC/bd.RData")
# simulate expression
gafile <- "/home/simingz/causalTWAS/annotation/refGene_processed_5e+05.txt"
simulate_expr(X, chr, pos, gafile, n.eQTL.range=1:5, h2.eQTL=0.2)









