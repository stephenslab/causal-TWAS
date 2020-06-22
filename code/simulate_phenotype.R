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
  wgtdir <- dirname(weight)
  wgtposfile <- file.path(wgtdir, paste0(basename(weight), ".pos"))

  wgtpos <- read.table(wgtposfile, header = T, stringsAsFactors = F)
  
  wgtpos <- transform(wgtpos, 
            "ID" = ifelse(duplicated(ID) | duplicated(ID, fromLast=TRUE), 
                                 paste(ID, ave(ID, ID, FUN=seq_along), sep='_ID'), 
                                 ID))
  
  write.table(wgtpos , file= paste0(wgtposfile, ".IDfixed") , row.names=F, col.names=T, sep="\t", quote = F)
  
  print(paste0("Number of genes with weights provided: ", dim(wgtpos)[1]))

  for (i in 1: dim(wgtpos)[1]){
    wf <- file.path(wgtdir, wgtpos[i, "WGT"])
    load(wf)
    gname <- wgtpos[i, "ID"]
    print(gname)

    g.method = method
    if (g.method == "best"){
      g.method = names(which.max(cv.performance["rsq", ]))
    }

    wgt.idx <- match(rownames(wgt.matrix), dat$snp)

    if (checksnps){
      datsnps <- data.frame("V1" = dat$chr[wgt.idx, 1],
                            "V4" = dat$pos[wgt.idx, 1],
                            "V5" = dat$counted[wgt.idx, 1],
                            "V6" = dat$alt[wgt.idx, 1])

      if (!identical(datsnps, snps)){
        stop("didn't pass snp consistency test, stopped ...")
      }
    }

    g <- as.matrix(dat$G[, wgt.idx])
    g[, is.na(wgt.idx)] <- 0
    
    exprlist[[gname]][["expr"]] <- g %*% wgt.matrix[, g.method]
    exprlist[[gname]][["chr"]] <- wgtpos[wgtpos$ID == gname, "CHR"]
    exprlist[[gname]][["p0"]] <- wgtpos[wgtpos$ID == gname, "P0"]# position boundary 1
    exprlist[[gname]][["p1"]] <- wgtpos[wgtpos$ID == gname, "P1"]# position boundary 2

    nmiss <- length(wgt.idx[is.na(wgt.idx)])
    n <- length(wgt.idx)
    exprlist[[gname]][["n"]] <- n
    exprlist[[gname]][["nmiss"]] <- nmiss
    exprlist[[gname]][["missrate"]] <- nmiss/n
    
    gc()
  }

  expr <- do.call(cbind, lapply(exprlist, '[[', "expr"))
  colnames(expr) <- names(exprlist)
  chrom <- unlist(lapply(exprlist,'[[', "chr"))
  p0 <- unlist(lapply(exprlist,'[[', "p0"))
  p1 <- unlist(lapply(exprlist,'[[', "p1"))

  # filter gene based on genotype missing rate and also expression among samples are not the same.
  gfilter <- unlist(lapply(exprlist, function(x) x[["missrate"]] < 1 & abs(max(x[["expr"]]) - min(x[["expr"]])) > 1e-8))
  expr <- expr[ , gfilter]
  chrom <- chrom[gfilter]
  p0 <- p0[gfilter]
  p1 <- p1[gfilter]

  print(paste0("Number of genes with imputed expression: ", dim(expr)[2]))
  gc()
  
  return(list("expr"= expr, "exprlist" = exprlist, "chrom" = chrom, "p0" = p0, "p1" = p1))
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

  set.seed(SED)

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
                       pve.snp.truth, pve.expr.truth, J, M)

  return(list("Y" = Y,"param" = param))
}
