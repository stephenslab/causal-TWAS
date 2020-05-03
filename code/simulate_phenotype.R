library(data.table)
# requires tibble

# generate cis expression of gene given weights
# weightslist: each gene is a list item, providing a dataframe with columns: chr, pos, major minor, labels, alpha and rows eQTLs.
# TODO
cis_expr <- function(dat, weightslist){
  exprlist <- list()
  for (gname in names(weightslist)){
    weights <- weightslist[[gname]]
    exprlist[[gname]] <- as.matrix(dat$G[ , match(weights$labels, labels)]) %*% weights$alpha
    if (is.na(exprlist[[gname]])) exprlist[[gname]] <- 0 #TODO: eQTL not genotyped should be omited, instead of simply put gene expression to 0
  }

  expr <- do.call(cbind, exprlist)
  return(expr)
}

# dat: a list with the following items:
#   - G: matrix, N x M, genotype, need to be scaled
#   - snp: matrix, N x 1, SNP name
#   - ref
#   - alt
#   - chr
#   - pos
#   - expr: matrix, N x J
# mode: string, "snp-only", "snp-expr"
# pve.expr: scalar [0,1]
# pve.snp: scalar [0,1]
# pi_beta: scalar [0,1]
# pi_theta: scalar [0,1]
# tau: 1 (to add)
simulate_phenotype<- function(dat,
                              mode = "snp-expr",
                              pve.expr = 0.1,
                              pve.snp = 0.1,
                              pi_beta = 0.01,
                              pi_theta = 0.0001,
                              tau = 1) {

  N <- dim(dat$G)[1]
  M <- dim(dat$G)[2]
  J <- dim(dat$expr)[2]

  J.c <- round(J * pi_beta)
  M.c <- round(M * pi_theta)

  set.seed(999)

  if (mode == "snp-only"){
    expr.meanvar <- NULL
    idx.cgene <- NULL
    idx.cSNP <- sample(1:M, M.c)
    sigma_beta <- NULL
    sigma_theta <- sqrt(pve.snp / (M.c * (1 - pve.snp - pve.expr)))
    e.beta <- NULL
    s.theta <- rnorm(M.c, mean = 0, sd = sigma_theta)

    Y <- dat$G[, idx.cSNP, drop = F] %*% s.theta +
          rnorm(N)

    pve.expr.truth <- 0
    pve.snp.truth <- var(dat$G[ , idx.cSNP, drop = F] %*% s.theta)/var(Y)
  }
  if (mode == "snp-expr"){
    expr.meanvar <- mean(apply(dat$expr, 2, var))
    idx.cgene <- sample(1:J, J.c)
    idx.cSNP <- sample(1:M, M.c)
    sigma_beta <- sqrt(pve.expr / (J.c * expr.meanvar * (1 - pve.snp - pve.expr)))
    sigma_theta <- sqrt(pve.snp / (M.c * (1 - pve.snp - pve.expr)))
    e.beta <- rnorm(J.c, mean = 0, sd = sigma_beta)
    s.theta <- rnorm(M.c, mean = 0, sd = sigma_theta)

    Y <- dat$expr[, idx.cgene, drop = F] %*% e.beta +
        dat$G[, idx.cSNP, drop = F] %*% s.theta +
        rnorm(N)

    pve.expr.truth <- var(dat$expr[ , idx.cgene, drop = F] %*% e.beta)/var(Y)
    pve.snp.truth <- var(dat$G[ , idx.cSNP, drop = F] %*% s.theta)/var(Y)
  }

  param <- tibble::lst(s.theta, e.beta,
                       sigma_theta, sigma_beta,
                       idx.cSNP, idx.cgene, M.c, J.c,
                       expr.meanvar,
                       pve.snp.truth, pve.expr.truth)

  return(list("Y"=Y,"param"= param))
}
