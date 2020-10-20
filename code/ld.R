# global variables: phenores, exprres, dat
library(snpStats)

#' @description calculates measures of linkage disequilibrium
#' between each column of x and each column of y
#' need global variables: phenores, exprres, dat
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
  x <- exprres$wgtlist[[gname]]
  if (is.null(dim(x))){
    dim(x) <- c(1,2)
    colnames(x) <- c("idx", "wgt")
  }
  eqtlgeno <- dat$G[ , x[abs(x[, "wgt"]) > 0, "idx"], drop = F]
  colnames(eqtlgeno) <- paste(gname, dat$snp[x[abs(x[, "wgt"]) > 0, "idx"], 1], sep = ".")
  return(eqtlgeno)
}


#' @description get ld for genes and all other causal genes or snps
#' and report max r^2 pairs
#'
#'
get_ld <- function(outname){
  g.eqtl <- lapply(colnames(exprres$expr), .eqtl_geno)

  causalg.eqtl <-lapply(colnames(exprres$expr)[phenores$param$idx.cgene], .eqtl_geno)
  cg <- do.call(cbind, causalg.eqtl)

  csnp <- dat$G[ ,phenores$param$idx.cSNP]
  colnames(csnp) <- dat$snp[phenores$param$idx.cSNP, 1]

  rsq <- lapply(g.eqtl, ld_max, y = cbind(cg,csnp), stats = "R.squared")

  outdf <- cbind(colnames(exprres$expr), rsq, names(rsq))
  colnames(outdf) <- c("name", "max_r2", "max_pair")

  write.table(outdf, file= paste0(outname, ".LD.txt") , row.names=F, col.names=T, sep="\t", quote = F)
}

#' @description get ld for gene pairs, select gene pairs that have r^2 more than
#' some value `r2cut`
#'
get_ld_gene <- function(r2cut, outname){
  g.eqtl <- lapply(colnames(exprres$expr), .eqtl_geno)
  names(g.eqtl) <- colnames(exprres$expr)

  ggr2 <- list()
  for (i in 1:length(g.eqtl)){
    gname <- names(g.eqtl)[i]
    rsq <- lapply(g.eqtl, ld_max, y = g.eqtl[[i]], stats = "R.squared")
    rsq[[i]] <- NULL
    ggr2[[gname]] <-  cbind(unlist(rsq), colnames(exprres$expr), gname)
  }

  outdf <- do.call(rbind, ggr2)
  colnames(outdf) <- c("g_max_r2", "gene1", "gene2")
  outdf[, "max_pair"] <- rownames(outdf)
  rownames(outdf) <- NULL
  save(outdf, file = paste0(outname, ".geneLD.txt") , row.names=F, col.names=T, sep="\t", quote = F)

  outdf2 <- outdf[outdf[,"g_max_r2"] > r2cut, ]
  return(outdf2)
}
