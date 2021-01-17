#-------------------------------------------------------------------------------------------
cl <- makeCluster(Ncore, outfile="")
registerDoParallel(cl)


iter <- 1
b <- 1
# test 1 --> not working, low CPU
a <- dat$G$bm()
bm.desc <- describe(a)


susiefunc <- function(rn, bm.desc){
  print(c(b, rn))
  
  sidx <- regionlist[[b]][[rn]][["sidx"]]
  
  a <- bigmemory::attach.big.matrix(bm.desc)
  sidx.n <-c(1:dim(a)[2])[sidx]
  x <- a[, sidx.n]
  1:5
  
}

parallel::clusterCall(cl, function() {
  
  library(bigmemory)
  library(susieR)
})

parallel::clusterExport(cl, c("regionlist", "bm.desc", "b"),
                        envir=environment())

fold.results <- parallel::parLapply(cl = cl, X = seq_len(length(regionlist[[b]][1:100])), fun = susiefunc, bm.desc = bm.desc)


# test 2 --> not working, low CPU
outdf.b <- foreach(r = 1: length(regionlist[[b]]), .combine = "rbind",
                   .packages = c("susieR", "bigstatsr")) %dopar% {
                     
                     r.s <- rs[r]
                     print(r.s)
                     x <- dat$G[, (1+r.s):(50+r.s)]
                     5
                   }

# test 3 --> not working, low CPU
sidx <- regionlist[[b]][["1"]][["sidx"]]
outdf.b <- foreach(rn = names(regionlist[[b]]), .combine = "rbind",
                   .packages = c("susieR", "bigstatsr")) %dopar% {
                     print(c(b, rn))
                     
                     X <- dat$G[, sidx]
                     1: 6
                   }

#test 4 --> slow, when assigning jobs
outdf.b.list <- list()
for (b in 1:22){
  
  # load genotype data
  load(pfileRds[b])
  
  # load expression Rd file, variable: exprres
  load(efiles[b])
  
  cau <- c(dat$snp[phenores$batch[[b]]$param$idx.cSNP,],
           exprres$gnames[phenores$batch[[b]]$param$idx.cgene])
  
  rns <- names(regionlist[[b]])
  
  outdf.b.b.list <- list()
  chunk <- Ncore * 3
  nb.b <- floor(length(regionlist[[b]])/chunk) + 1
  print("here")
  for (b.b in 1:nb.b) {
    regionlist.temp <- list()
    for (r in 1:chunk){
      rn.b <- rns[(b.b - 1) * chunk + r]
      gidx <- regionlist[[b]][[rn.b]][["gidx"]]
      sidx <- regionlist[[b]][[rn.b]][["sidx"]]
      
      regionlist.temp[[rn.b]][["X"]] <- cbind(exprres$expr[, gidx], dat$G[, sidx])
      
      anno.gene <- cbind(exprres$gnames[gidx], exprres$chrom[gidx],  exprres$p0[gidx],
                         rep("gene", length(gidx[gidx])))
      anno.SNP <- cbind(dat$snp[sidx,], dat$chr[sidx,], dat$pos[sidx,], "SNP")
      anno <- rbind(anno.gene, anno.SNP)
      colnames(anno) <-  c("name", "chr", "pos", "type")
      anno <- as.data.frame(anno)
      anno$ifcausal <- ifelse(anno$name %in% cau, 1, 0)
      
      regionlist.temp[[rn.b]][["anno"]] <- anno
    }
    print("here2")
    y <- phenores$Y
    # run susie for each region (in parallel)
    outdf.b.b.list[[b.b]] <- foreach(rn = names(regionlist.temp), .combine = "rbind",
                                     .packages = c("susieR")) %dopar% {
                                       print(c(b, b.b, rn))
                                       
                                       susieres <- susie(regionlist.temp[[rn]][["X"]], y , L = L)
                                       anno <- cbind(regionlist.temp[[rn]][["anno"]], susieres$pip)
                                       anno
                                       print("done")
                                     }
    print(b.b)
  }
  
  outdf.b.list[[b]] <- do.call(rbind, outdf.b.b.list)
  
}
outdf <- do.call(rbind, outdf.b.list)

stopCluster(cl)

loginfo("susie done for %s", outname)

#test 5

outdf <- foreach (b =  1:22, .combine = "rbind",.packages = c("susieR")) %dopar% {
  
  # load genotype data
  load(pfileRds[b])
  
  # load expression Rd file, variable: exprres
  load(efiles[b])
  
  cau <- c(dat$snp[phenores$batch[[b]]$param$idx.cSNP,],
           exprres$gnames[phenores$batch[[b]]$param$idx.cgene])
  
  # run susie for each region (in parallel)
  outdf.b.list <- list()
  for(rn in names(regionlist[[b]])){
    print(c(b, rn))
    
    gidx <- regionlist[[b]][[rn]][["gidx"]]
    sidx <- regionlist[[b]][[rn]][["sidx"]]
    
    print("here1")
    X <- cbind(exprres$expr[, gidx], dat$G[, sidx])
    
    anno.gene <- cbind(exprres$gnames[gidx], exprres$chrom[gidx],  exprres$p0[gidx],
                       rep("gene", length(gidx[gidx])))
    anno.SNP <- cbind(dat$snp[sidx,], dat$chr[sidx,], dat$pos[sidx,], "SNP")
    anno <- rbind(anno.gene, anno.SNP)
    colnames(anno) <-  c("name", "chr", "pos", "type")
    anno <- as.data.frame(anno)
    anno$ifcausal <- ifelse(anno$name %in% cau, 1, 0)
    print("here2")
    # susieres <- susie(X, phenores$Y , L = L)
    print("here3")
    # anno <- cbind(anno, susieres$pip)
    outdf.b.list[[rn]] <- anno
  }
  
  outdf.b <- do.call(rbind, outdf.b.list)
  
  outdf.b
  
}

stopCluster(cl)

loginfo("susie done for %s", outname)

