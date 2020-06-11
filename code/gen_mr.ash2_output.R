
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

}
