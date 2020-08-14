
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

gen_mr.ash_gene_output <- function(g.fit, outname){
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
  return(g.outdf.sorted)
}

gen_mr.ash_snp_output <- function(s.fit, outname){
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
  return(s.outdf.sorted)
}

gen_mr.ash2_output <- function(g.fit, s.fit, outname){
  g.outdf <- gen_mr.ash_gene_output(g.fit, outname)
  s.outdf <- gen_mr.ash_snp_output(s.fit, outname)

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


mr.ash2s_save <- function(mr.ash2s.fit, outname){

  g.fit <- mr.ash2s.fit$g.fit
  s.fit <-  mr.ash2s.fit$s.fit

  mr.ash.init <- mr.ash2s.fit$mr.ash.init
  init.order <- mr.ash2s.fit$init.order
  iter.order <- mr.ash2s.fit$iter.order

  outprefix <- paste0(outname,"-mr.ash2s.", mr.ash.init, init.order,"-", iter.order)

  save(mr.ash2s.fit, file = paste0(outprefix, ".Rd"))
  gen_mr.ash2_output(g.fit, s.fit, outprefix)

}
