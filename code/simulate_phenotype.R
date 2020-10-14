library(data.table)
# requires tibble
#' @param dat a list with the following items:
#'   - G: matrix, N x M, genotype, need to be scaled
#'   - snp: matrix, N x 1, SNP name
#'   - ref
#'   - alt
#'   - chr
#'   - pos
#'   - expr: matrix, N x J
#' @param mode string, "snp-only", "snp-expr"
#' @param pve.expr scalar [0,1], not used when if `sigma_beta` is given
#' @param pve.snp scalar [0,1], not used if `sigma_theta` is given
#' @param pi_beta scalar [0,1]
#' @param pi_theta scalar [0,1]
#' @param sigma_beta default NULL, if not given, will calculate from pve.
#' @param sigma_theta default NULL, if not given, will calculate from pve.
#' @return a list, Y.g is phenotype determined by genetic, will need to add
#'  error term ( pve calculation is based on N(0,1)).
simulate_phenotype<- function(dat,
                              mode = "snp-expr",
                              pve.expr = 0.1,
                              pve.snp = 0.1,
                              pi_beta = 0.01,
                              pi_theta = 0.0001,
                              sigma_beta = NULL,
                              sigma_theta = NULL
                              ) {

  N <- dim(dat$G)[1]
  M <- dim(dat$G)[2]
  J <- dim(dat$expr)[2]

  J.c <- round(J * pi_beta)
  M.c <- round(M * pi_theta)

  if (mode == "snp-only"){
    expr.meanvar <- NULL
    idx.cgene <- NULL
    idx.cSNP <- sample(1:M, M.c)
    sigma_beta <- NULL
    if (is.null(sigma_theta)){
      sigma_theta <- sqrt(pve.snp / (M.c * (1 - pve.snp - pve.expr)))
    }

    e.beta <- NULL
    s.theta <- rnorm(M.c, mean = 0, sd = sigma_theta)

    Y.g <- dat$G[, idx.cSNP, drop = F] %*% s.theta

    var.gene <- NULL
    var.snp <- var(dat$G[ , idx.cSNP, drop = F] %*% s.theta)
  }
  if (mode == "snp-expr"){
    idx.cgene <- sample(1:J, J.c)
    idx.cSNP <- sample(1:M, M.c)

    if (is.null(sigma_beta)){
      expr.meanvar <- mean(apply(dat$expr, 2, var))
      sigma_beta <- sqrt(pve.expr / (J.c * expr.meanvar * (1 - pve.snp - pve.expr)))
    }

    if (is.null(sigma_theta)){
      sigma_theta <- sqrt(pve.snp / (M.c * (1 - pve.snp - pve.expr)))
    }
    e.beta <- rnorm(J.c, mean = 0, sd = sigma_beta)
    s.theta <- rnorm(M.c, mean = 0, sd = sigma_theta)

    Y.g <- dat$expr[, idx.cgene, drop = F] %*% e.beta +
        dat$G[, idx.cSNP, drop = F] %*% s.theta

    var.gene <- var(dat$expr[ , idx.cgene, drop = F] %*% e.beta)
    var.snp <- var(dat$G[ , idx.cSNP, drop = F] %*% s.theta)
  }

  param <- tibble::lst(s.theta, e.beta,
                       sigma_theta, sigma_beta,
                       idx.cSNP, idx.cgene, M.c, J.c,
                       J, M, N, var.gene, var.snp)

  return(list("Y.g" = Y.g, "param" = param))
}
