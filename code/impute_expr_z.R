library(data.table)
# requires tibble

#' @description  generate cis expression of gene given weights
#' @param dat: a list with the following items:
#'   - G: matrix, N x M, genotype, need to be scaled
#'   - snp: matrix, N x 1, SNP name
#' @param weight a string, pointing to the fusion/twas format of weights
#' @param method a string,  blup/bslmm/lasso/top1/enet/best
#'   "best" means the method giving the best cross validation R^2
#' @param zdf a data frame, with columns "snp", "chr", "pos", "ref", "alt", "z"
#' @param checksnps T/F, if need to check SNP consistency between weights
#'    file and genotype. fusion format of weights gives snp information from plink .bim file. checking include chr, pos, A1, A2.
#' @return a matrix, N x J. Missing genotypes will not contribute to cis-expr and give warnings.
impute_expr_z <- function(dat, weight, zdf, method = "bslmm", checksnps = F){
  exprlist <- list()
  qclist <- list()
  wgtdir <- dirname(weight)
  wgtposfile <- file.path(wgtdir, paste0(basename(weight), ".pos"))

  if (!file.exists(paste0(wgtposfile, ".IDfixed"))){
    wgtpos <- read.table(wgtposfile, header = T, stringsAsFactors = F)
    wgtpos <- transform(wgtpos,
                        "ID" = ifelse(duplicated(ID) | duplicated(ID, fromLast=TRUE),
                                      paste(ID, ave(ID, ID, FUN=seq_along), sep='_ID'),
                                      ID))

    write.table(wgtpos, file = paste0(wgtposfile, ".IDfixed") , row.names=F, col.names=T, sep="\t", quote = F)
  }
  wgtpos <- read.table(paste0(wgtposfile, ".IDfixed"), header = T, stringsAsFactors = F)

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
    if (!(g.method %in% names(cv.performance[1,]))) next

    wgt.matrix <- wgt.matrix[abs(wgt.matrix[, g.method]) > 0, , drop = F]
    wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), ,drop = F]

    if (dim(wgt.matrix)[1] == 0) next

    snpnames <- Reduce(intersect, list(rownames(wgt.matrix), dat$snp, zdf$snp))

    if (length(snpnames) == 0) next

    wgt.idx <- match(snpnames, rownames(wgt.matrix))
    dat.idx <- match(snpnames, dat$snp)
    zdf.idx <- match(snpnames, zdf$snp)

    if (checksnps){
      datsnps <- data.frame("V1" = dat$chr[dat.idx, 1],
                            "V4" = dat$pos[dat.idx, 1],
                            "V5" = dat$counted[dat.idx, 1],
                            "V6" = dat$alt[dat.idx, 1])
      zdfsnps <- zdf[zdf.idx, c("chr", "pos", "ref", "alt"), drop = F]
      colnames(zdfsnps) <- c("V1", "V4", "V5", "V6")
      wgtsnps <- snps[wgt.idx, c("V1", "V4", "V5", "V6"), drop = F ]

      if (all(sapply(list(zdfsnps, datsnps), FUN = identical, wgtsnps))){
        stop("didn't pass snp consistency test, stopped ...")
      }
    }

    wgt <-  wgt.matrix[wgt.idx, g.method, drop = F]
    X.g <- dat$G[, dat.idx, drop = F]

    gexpr <- X.g %*% wgt
    if (abs(max(gexpr) - min(gexpr)) < 1e-8) next

    z.s <- as.matrix(zdf[zdf.idx, "z"])
    var.s <- apply(X.g, 2, var)
    Gamma.g <- cov(X.g)

    z.g <-  (t(wgt) * var.s) %*% z.s/ sqrt(t(wgt) %*% Gamma.g %*% wgt)

    exprlist[[gname]][["expr"]] <- gexpr
    exprlist[[gname]][["z.g"]] <- z.g
    exprlist[[gname]][["chr"]] <- wgtpos[wgtpos$ID == gname, "CHR"]
    exprlist[[gname]][["p0"]] <- wgtpos[wgtpos$ID == gname, "P0"]# position boundary 1
    exprlist[[gname]][["p1"]] <- wgtpos[wgtpos$ID == gname, "P1"]# position boundary 2
    exprlist[[gname]][["wgt"]] <- wgt

    qclist[[gname]][["n"]] <- nrow(wgt.matrix)
    qclist[[gname]][["nmiss"]] <- nrow(wgt.matrix) - length(snpnames)
    qclist[[gname]][["missrate"]] <-  qclist[[gname]][["nmiss"]]/ qclist[[gname]][["n"]]

    gc()
  }

  z.g <- unlist(lapply(exprlist,'[[', "z.g"))
  expr <- do.call(cbind, lapply(exprlist, '[[', "expr"))
  colnames(expr) <- names(exprlist)
  chrom <- unlist(lapply(exprlist,'[[', "chr"))
  p0 <- unlist(lapply(exprlist,'[[', "p0"))
  p1 <- unlist(lapply(exprlist,'[[', "p1"))
  wgtlist <- lapply(exprlist, '[[', "wgt")

  print(paste0("Number of genes with imputed expression: ", dim(expr)[2]))
  gc()

  return(list("z.g" = z.g,
              "expr"= expr,
              "wgtlist" = wgtlist,
              "qclist" = qclist,
              "chrom" = chrom,
              "p0" = p0,
              "p1" = p1,
              "gnames" = names(exprlist)))
}
