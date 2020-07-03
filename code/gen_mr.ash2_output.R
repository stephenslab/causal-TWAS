# global variables: phenores, exprres, dat
library(snpStats)

para_mr.ash <- function(fit){
  pi1 <-  1-fit$pi[[1]]
  pve <- get_pve(fit)
  p <- c(pi1, pve)
  names(p) <- c("pi1", "pve")
  return(p)
}

para_mr.ash2 <- function(g.fit, s.fit){

  gpara <- para_mr.ash(g.fit)
  names(gpara) <- paste0("gene.", names(gpara))

  spara <- para_mr.ash(s.fit)
  names(spara) <- paste0("snp.", names(spara))

  outdf <- data.frame("estimated" = c(gpara,spara))
  outdf$truth <- c(phenores$param$J.c/phenores$param$J,
                   phenores$param$pve.expr.truth,
                   phenores$param$M.c/phenores$param$M,
                   phenores$param$pve.snp.truth)

  return(outdf)
}

#' @description calculates measures of linkage disequilibrium
#' between each column of x and each column of y
#' @return the maximum of ld measures of these pairs
#'
ld_max <- function(x,  y, stats = "R.squared"){
  x <- as(round(x),"SnpMatrix")
  y <- as(round(y),"SnpMatrix")
  s <- ld(x, y =y, stats = stats)
  inds <- arrayInd(which.max(s), dim(s))
  s.max <- max(s)
  names(s.max) <- paste(rownames(s)[inds[,1]], colnames(s)[inds[,2]], sep = ':')
  return(s.max)
}

.eqtl_geno <- function(gname){
  x <- exprres$exprlist[[gname]]
  eqtlgeno <- dat$G[ , x$wgt[abs(x$wgt[, "wgt"]) > 0, "idx"], drop = F]
  colnames(eqtlgeno) <- paste(gname, dat$snp[x$wgt[abs(x$wgt[, "wgt"]) > 0, "idx"], 1], sep = ".")
  return(eqtlgeno)
}

get_ld <- function(outname){
  g.eqtl <- lapply(colnames(exprres$expr), .eqtl_geno)

  causalg.eqtl <-lapply(colnames(exprres$expr)[phenores$param$idx.cgene], .eqtl_geno)
  cg <- do.call(cbind, causalg.eqtl)

  csnp <- dat$G[ ,phenores$param$idx.cSNP]
  colnames(csnp) <- dat$snp[phenores$param$idx.cSNP ,1]

  rsq <- unlist(lapply(g.eqtl, ld_max, y = cbind(cg,csnp), stats = "R.squared"))

  outdf <- cbind(colnames(exprres$expr), rsq, names(rsq))
  colnames(outdf) <- c("name", "max_r2", "max_pair")

  write.table(outdf, file= paste0(outname, ".LD.txt") , row.names=F, col.names=T, sep="\t", quote = F)
}


gen_mr.ash2_output <- function(g.fit, s.fit, outname){
  e.b <- rep(0, length(g.fit$beta))
  e.b[phenores$param$idx.cgene] <- phenores$param$e.beta

  indi <- rep(0, length(g.fit$beta))
  indi[phenores$param$idx.cgene] <- 1
  pip <- get_pip(g.fit)

  g.outdf <- data.frame("chr" = exprres$chrom,
                        "p0" = exprres$p0,
                        "p1" = exprres$p1,
                        "ifcausal"= indi,
                        "name" = colnames(exprres$expr),
                        "TrueBeta" = e.b,
                        "EstBeta" = g.fit$beta,
                        "PIP" = pip)

  g.outdf.sorted <- g.outdf[with(g.outdf, order(chr, p0)), ]
  write.table(g.outdf.sorted , file= paste0(outname, ".expr.txt") , row.names=F, col.names=T, sep="\t", quote = F)

  s.b <- rep(0, length(s.fit$beta))
  s.b[phenores$param$idx.cSNP] <- phenores$param$s.theta

  indi <- rep(0, length(s.fit$beta))
  indi[phenores$param$idx.cSNP] <- 1
  pip <- get_pip(s.fit)

  s.outdf <- data.frame("chr" = dat$chr[,1],
                        "p0" = dat$pos[,1],
                        "p1" = dat$pos[,1],
                        "ifcausal"= indi,
                        "name" = dat$snp[,1],
                        "TrueBeta" = s.b,
                        "EstBeta" = s.fit$beta,
                        "PIP" = pip)

  s.outdf.sorted <- s.outdf[with(s.outdf, order(chr, p0)), ]
  write.table(s.outdf.sorted , file= paste0(outname, ".snp.txt") , row.names=F, col.names=T, sep="\t", quote = F)

  cg.outdf <- g.outdf[g.outdf[, "ifcausal"] == 1, ]
  cg.outdf[, "ifcausal"] <- "causal gene"
  cs.outdf <- s.outdf[s.outdf[, "ifcausal"] == 1, ]
  cs.outdf[, "ifcausal"] <- "causal snp"

  c.outdf <- rbind(cg.outdf, cs.outdf)
  c.outdf[, "p0"] <- c.outdf[, "p0"] - 1 # bed format, for locuszoom

  write.table(c.outdf, file= paste0(outname, ".causal-snpexpr.bed") , row.names=F, col.names=F, sep="\t", quote = F)

  para <- para_mr.ash2(g.fit, s.fit)
  write.table(para , file= paste0(outname, ".param.txt") , row.names=T, col.names=T, sep="\t", quote = F)
}
