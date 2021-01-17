#' Simulate phenotype
#' @description simulate Y under the sparse model (spike and slab prior).
#' will standardize gene expression and SNP genotype when simulate.
#' @param pve.expr scalar [0,1], not used when if `sigma_beta` is given
#' @param pve.snp scalar [0,1], not used if `sigma_theta` is given
#' @param pi_beta scalar [0,1]
#' @param pi_theta scalar [0,1]
#' @param sigma_beta default NULL, if not given, will calculate from pve.
#' @param sigma_theta default NULL, if not given, will calculate from pve.
#' @return a list, Y.g is phenotype determined by genetic, will need to add
#'  error term ( pve calculation is based on N(0,1)).
#' requires tibble, pgenlibr, scaleRcpp
simulate_phenotype<- function(pgenfs,
                              exprfs,
                              pve.expr = 0.1,
                              pve.snp = 0.1,
                              pi_beta = 0.01,
                              pi_theta = 0.0001,
                              sigma_beta = NULL,
                              sigma_theta = NULL
) {
  pvarfs <- sapply(pgenfs, prep_pvar, outputdir = outputdir)
  exprvarfs <- sapply(exprfs, prep_exprvar)

  pgen1 <- prep_pgen(pgenf = pgenfs[1], pvarfs[1])
  N <- pgenlibr::GetRawSampleCt(pgen1)

  M.b <- sapply(pvarfs, function(x) nrow(read_pvar(x)))
  J.b <- sapply(exprvarfs, function(x) nrow(read_exprvar(x)))

  M <- sum(M.b)
  J <- sum(J.b)

  if (is.null(sigma_beta)){
    expr.meanvar <- 1 # always scale before simulate
    J.c <- round(J * pi_beta)
    sigma_beta <- sqrt(pve.expr / (J.c * expr.meanvar * (1 - pve.snp - pve.expr)))
    if (is.infinite(sigma_beta)) sigma_beta <- 0
  }

  if (is.null(sigma_theta)){
    M.c <- round(M * pi_theta)
    sigma_theta <- sqrt(pve.snp / (M.c * (1 - pve.snp - pve.expr)))
  }

  phenores <- list("batch" = list())
  for ( b in 1:length(pgenfs)){

    idx.cgene <- which(rbinom(J.b[b], 1, pi_beta) == 1)
    idx.cSNP <- which(rbinom(M.b[b], 1, pi_theta) == 1)

    X.g <- read_expr(exprfs[b], variantidx = idx.cgene)
    if (is.null(X.g)){
      X.g <- matrix(, nrow= N, ncol =0)
    }
    X.g <- scaleRcpp(X.g) # always scale expr

    pgen <- prep_pgen(pgenf = pgenfs[b], pvarfs[b])

    X.s <- read_pgen(pgen, variantidx = idx.cSNP)
    X.s <- scaleRcpp(X.s) # always scale genotype

    e.beta <- rnorm(length(idx.cgene), mean = 0, sd = sigma_beta)
    s.theta <- rnorm(length(idx.cSNP), mean = 0, sd = sigma_theta)

    Y.g <- X.g %*% e.beta + X.s %*% s.theta

    var.gene <- var(X.g %*% e.beta)
    var.snp <- var(X.s %*% s.theta)

    phenores[["batch"]][[b]] <- tibble::lst(Y.g, s.theta, e.beta,
                                             sigma_theta, sigma_beta,
                                             idx.cSNP, idx.cgene,
                                             J = J.b[b], M = M.b[b],
                                             N, var.gene, var.snp)
  }

  Y <- matrix(rowSums(do.call(cbind, lapply(phenores[["batch"]], '[[', "Y.g"))) + rnorm(N), ncol = 1)
  phenores$Y <- Y

  var.y <- var(Y)

  phenores$param$pve.snp.truth <-
    sum(unlist(lapply(phenores[["batch"]], '[[', "var.snp")))/var.y
  phenores$param$pve.gene.truth <-
    sum(unlist(lapply(phenores[["batch"]], '[[', "var.gene")))/var.y

  return(phenores)

}


