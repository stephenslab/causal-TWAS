weight <- "/home/simingz/causalTWAS/fusion_weights/YFS.BLOOD.RNAARR"
wgtfiles <- Sys.glob(file.path(weight, "*.wgt.RDat"))
wgtposfile <- file.path(dirname(weight), paste0(basename(weight), ".pos"))
wgtpos <- read.table(wgtposfile, header = T, stringsAsFactors = F)

f <- "~/causalTWAS/simulations/simulation_ashtest_20200503/mr.ash2_20200503-4-cis-expr.Rd"
load(f)

for (wf in wgtfiles){
  gname <- strsplit(basename(wf), "." , fixed = T)[[1]][2]
  #print(gname)
  exprres$exprlist[[gname]][["chr"]] <- wgtpos[wgtpos$ID == gname, "CHR"]
  exprres$exprlist[[gname]][["p0"]] <- wgtpos[wgtpos$ID == gname, "P0"] # position boundary 1
  exprres$exprlist[[gname]][["p1"]] <- wgtpos[wgtpos$ID == gname, "P1"] # position boundary 2
}

gfilter <- unlist(lapply(exprres$exprlist, function(x) x[["missrate"]] < 1))

chrom <- unlist(lapply(exprres$exprlist,'[[', "chr"))
p0 <- unlist(lapply(exprres$exprlist,'[[', "p0"))
p1 <- unlist(lapply(exprres$exprlist,'[[', "p1"))

chrom <- chrom[gfilter]
p0 <- p0[gfilter]
p1 <- p1[gfilter]

exprres$chrom <- chrom
exprres$p0 <- p0
exprres$p1 <- p1

save(exprres, file = f)

# ====================================================
for (i in 1:4){
  outname0 <- paste0("mr.ash2_20200503-",i)

  load(paste0(outname0, "-cis-expr.Rd"))
  load(paste0(outname0, "-pheno.Rd"))

  load(paste0(outname0,"-mr.ash2.expr-res.Rd"))

  g.fit <- mr.ash2.fit$fit1
  s.fit <-  mr.ash2.fit$fit2
  s.fit$data$X <- dat$G
  outname <- paste0(outname0, "-mr.ash2.expr-res")
  gen_mr.ash2_output(g.fit, s.fit, outname)

  load(paste0(outname0,"-mr.ash2.snp-res.Rd"))

  g.fit <- mr.ash2.fit$fit2
  s.fit <-  mr.ash2.fit$fit1
  s.fit$data$X <- dat$G
  outname <- paste0(outname0, "-mr.ash2.snp-res")
  gen_mr.ash2_output(g.fit, s.fit, outname)
}


fi0 <- "20200616-3-mr.ash2s.expr-res.expr.txt"
fi <- read.table(fi0, header = T)

fo <- merge(fi, wgtpos, by = "name")
outc <- colnames(fi)
outc[1:3] <- c("CHR", "P0", "P1")

out <- fo[,outc]
colnames(out)[1:3] <- colnames(fi)[1:3]
write.table( out , file= fi0 , row.names=F, col.names=T, sep="\t", quote = F)
#====================================================================================
load("/project2/mstephens/causalTWAS/ukbiobank/ukb_chr17to22_s20000.FBM.Rd")
for (i in 1:50) {
  load(paste0(simdatadir, "simu_20200721-1-",i, "-pheno.Rd"))
  load(paste0(simdatadir, "simu_20200721-1-",i, "-cis-expr.Rd"))
  cau <- c(dat$snp[phenores$param$idx.cSNP,], colnames(exprres$expr)[phenores$param$idx.cgene])
  tf <- paste0(susiedir,"20200721-1-", i, ".fixedprior2.L1.susieres.expr.txt")
  outgdf <- read.table(tf, header = T)
  outgdf$ifcausal <- ifelse(outgdf$name %in% cau, 1, 0)
  write.table(outgdf, file = tf , row.names=F, col.names= T, sep="\t", quote = F)
}
