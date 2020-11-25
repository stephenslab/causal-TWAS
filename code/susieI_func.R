#' require `pfileRds`,`efiles`, `phenores`, `codedir`
#' `plugin_prior_variance`, `V.g`, `V.s`
#' `L`
#' `ifnullweight`
#' `coverage`

susieI <- function(prior.gene_init =NULL,
                   prior.SNP_init =NULL,
                   regionlist = NULL,
                   outname = "susieI",
                   niter = 1,
                   Ncore = 15) {

  varY <- var(phenores$Y)

  cl <- makeCluster(Ncore,outfile="")
  registerDoParallel(cl)

  prior.SNP_rec <- rep(0, Niter)
  prior.gene_rec <- rep(0, Niter)

  for (iter in 1:niter){

    loginfo("run iteration %s", iter)

    snp.rpiplist <- list()
    gene.rpiplist <- list()
    outdf <- foreach (b = 1:22, .combine = "rbind",.packages = c("susieR", "bigstatsr"),.export = ls(globalenv())) %dopar% {
      source(paste0(codedir, "susie_func.R"))

      # load genotype data
      load(pfileRds[b])

      # load expression Rd file, variable: exprres
      load(efiles[b])

      cau <- c(dat$snp[phenores$batch[[b]]$param$idx.cSNP,],
               exprres$gnames[phenores$batch[[b]]$param$idx.cgene])

      # run susie for each region (in parallel by batch)
      outdf.b.list <- list()
      for (rn in names(regionlist[[b]])) {
        print(c(iter, b, rn))

        gidx <- regionlist[[b]][[rn]][["gidx"]]
        sidx <- regionlist[[b]][[rn]][["sidx"]]
        p <- length(gidx[gidx]) + length(sidx[sidx])
        prior <- rep(0, p)

        if (is.null(prior.gene_init) | is.null(prior.SNP_init)){
          prior.gene_init <- 1/p
          prior.SNP_init <- 1/p
        }

        if (iter == 1) {
          prior <- c(rep(prior.gene_init, length(gidx[gidx])), rep(prior.SNP_init, length(sidx[sidx])))
        } else {
          prior <- c(rep(prior.gene, length(gidx[gidx])), rep(prior.SNP, length(sidx[sidx])))
        }

        if (isTRUE(plugin_prior_variance)){
          V.scaled <- c(rep(V.g/varY, length(gidx[gidx])), rep(V.s/varY, length(sidx[sidx])))
          V.scaled <- matrix(rep(V.scaled, each = L), nrow=L)
        } else{
          V.scaled <- matrix(rep(0.2, L * p), nrow = L) # following the default in susieR
        }

        if (isTRUE(ifnullweight)){
          nw <- max(0, 1 - sum(prior))
          prior <- prior/(1-nw)
        } else {
          nw <- NULL
        }

        X <- cbind(exprres$expr[, gidx], dat$G[, sidx])

        susieres <- susie(X, phenores$Y, L = L, prior_weights = prior,
                          null_weight = nw, scaled_prior_variance = V.scaled,
                          estimate_prior_variance = !plugin_prior_variance, coverage = coverage)

        anno.gene <- cbind(exprres$gnames[gidx], exprres$chrom[gidx],  exprres$p0[gidx],
                           rep("gene", length(gidx[gidx])))
        anno.SNP <- cbind(dat$snp[sidx,], dat$chr[sidx,], dat$pos[sidx,],
                          rep("SNP", length(sidx[sidx])))

        anno <- rbind(anno.gene, anno.SNP)
        colnames(anno) <-  c("name", "chr", "pos", "type")
        anno <- as.data.frame(anno)
        anno$ifcausal <- ifelse(anno$name %in% cau, 1, 0)
        anno$b <- b
        anno$rn <- rn
        anno$cs_index <- 0
        if (!is.null(susieres$sets$cs)){
          for (cs_i in susieres$sets$cs_index){
            X.idx <- susieres$sets$cs[[paste0("L", cs_i)]]
            anno$cs_index[X.idx] <- cs_i #TODO: note this ignore the fact that some variant can belong to multiple CS
          }
        }

        outdf.rn <- cbind(anno, susieres$pip)
        colnames(outdf.rn)[9] <- "susie_pip"

        outdf.b.list[[rn]] <- outdf.rn
      }

      outdf.b <- do.call(rbind, outdf.b.list)
      outdf.b
    }

    prior.SNP <- mean(outdf[outdf[ , "type"] == "SNP", "susie_pip"])
    prior.gene <- mean(outdf[outdf[ , "type"] == "gene", "susie_pip"])
    print(c(prior.SNP, prior.gene))

    prior.SNP_rec[iter] <- prior.SNP
    prior.gene_rec[iter] <- prior.gene
    save(prior.gene_rec, prior.SNP_rec, file = paste(outname, "susieIres.Rd", sep = "."))
    gc()
    write.table(outdf, file= paste0(outname, ".susieI.txt" ) , row.names=F, col.names=T, sep="\t", quote = F)
  }

  stopCluster(cl)

  return(list("prior.gene"= prior.gene,
              "prior.SNP" = prior.SNP))

}

