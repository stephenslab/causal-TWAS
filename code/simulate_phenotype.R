library(data.table)
# requires tibble

#' @description  generate cis expression of gene given weights
#' @param dat: a list with the following items:
#'   - G: matrix, N x M, genotype, need to be scaled
#'   - snp: matrix, N x 1, SNP name
#' @param weight a string, pointing to the fusion/twas format of weights
#' @param method a string,  blup/bslmm/lasso/top1/enet/best
#'   "best" means the method giving the best cross validation R^2
#' @param checksnps T/F, if need to check SNP consistency between weights
#'    file and genotype. fusion format of weights gives snp information from plink .bim file. checking include chr, pos, A1, A2.
#' @return a matrix, N x J. Missing genotypes will not contribute to cis-expr and give warnings.
cis_expr <- function(dat, weight, method = "bslmm", checksnps = F){
  exprlist <- list()
  wgtfiles <- Sys.glob(file.path(weight, "*.wgt.RDat"))

  print(paste0("Number of genes with weights provided: ", length(wgtfiles)))

  for (wf in wgtfiles){
    load(wf)
    gname <- strsplit(basename(wf), "." , fixed = T)[[1]][2]
    print(gname)

    g.method = method
    if (g.method == "best"){
      g.method = names(which.max(cv.performance["rsq", ]))
    }

    wgt.idx <- match(rownames(wgt.matrix), dat$snp)

    if (checksnps){
      datsnps <- data.frame("V1" = dat$chr,
                            "V4" = dat$pos,
                            "V5" = dat$counted,
                            "V6" = dat$alt)

      if (!identical(datsnps, snps)){
        stop("didn't pass snp consistency test, stopped ...")
      }
    }

    g <- as.matrix(dat$G[, wgt.idx])
    g[, is.na(wgt.idx)] <- 0
    exprlist[[gname]][["expr"]] <- g %*% wgt.matrix[, g.method]

    nmiss <- length(wgt.idx[is.na(wgt.idx)])
    n <- length(wgt.idx)
    exprlist[[gname]][["n"]] <- n
    exprlist[[gname]][["nmiss"]] <- nmiss
    exprlist[[gname]][["missrate"]] <- nmiss/n
  }

  expr <- do.call(cbind, lapply(exprlist, '[[', "expr"))
  colnames(expr) <- names(exprlist)

  # filter gene based on genotype missing rate
  gfilter <- unlist(lapply(exprlist, function(x) x[["missrate"]] < 1))
  expr <- expr[ , gfilter]
  print(paste0("Number of genes with imputed expression: ", dim(expr)[2]))

  return(list("expr"= expr, "exprlist" = exprlist))
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

  return(list("Y" = Y,"param" = param))
}
