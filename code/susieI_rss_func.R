#' require `pfileRds`, `pfileRds2`, `efiles`, `phenores`, `codedir`, `zdf`
#' `plugin_prior_variance`, `V.gene`, `V.SNP`, `est_prior_variance`
#' `L`,  `z_ld_weight`
#' `ifnullweight`
#' `coverage`

susieI_rss <- function(prior.gene_init =NULL,
                   prior.SNP_init =NULL,
                   regionlist = NULL,
                   outname = "susieI",
                   niter = 1,
                   Ncore = 15) {

  cl <- makeCluster(Ncore,outfile="")
  registerDoParallel(cl)

  prior.SNP_rec <- rep(0, Niter)
  prior.gene_rec <- rep(0, Niter)
  V.gene_rec <- rep(0, Niter)
  V.SNP_rec <- rep(0, Niter)

  for (iter in 1:niter){

    loginfo("run iteration %s", iter)

    snp.rpiplist <- list()
    gene.rpiplist <- list()
    outdf <- foreach (b = 1:22, .combine = "rbind",.packages = c("susieR", "bigstatsr"),.export = ls(globalenv())) %dopar% {
      source(paste0(codedir, "susie_func.R"))

      # load GWAS genotype data (only used for getting true parameter)
      load(pfileRds2[b])

      # load expression Rd file, variable: exprres
      load(efiles[b])

      cau <- c(dat$snp[phenores$batch[[b]]$param$idx.cSNP,],
               exprres$gnames[phenores$batch[[b]]$param$idx.cgene])

      # load LD genotype data
      load(pfileRds[b])

      # run susie for each region (in parallel by batch)
      outdf.b.list <- list()
      for (rn in names(regionlist[[b]])) {
        print(c(iter, b, rn))

        gidx <- regionlist[[b]][[rn]][["gidx"]]
        sidx <- regionlist[[b]][[rn]][["sidx"]]
        ldsidx <- regionlist[[b]][[rn]][["ldsidx"]]
        p <- length(gidx) + length(sidx)
        prior <- rep(0, p)

        if (is.null(prior.gene_init) | is.null(prior.SNP_init)){
          prior.gene_init <- 1/p
          prior.SNP_init <- 1/p
        }

        if (iter == 1) {
          prior <- c(rep(prior.gene_init, length(gidx)), rep(prior.SNP_init, length(sidx)))
        } else {
          prior <- c(rep(prior.gene, length(gidx)), rep(prior.SNP, length(sidx)))
        }

        if (isTRUE(plugin_prior_variance)){
          V <- c(rep(V.gene, length(gidx)), rep(V.SNP, length(sidx)))
          V <- matrix(rep(V, each = L), nrow=L)
        } else{
          V <- matrix(rep(50, L * p), nrow = L) # following the default in susieR
        }

        if (isTRUE(ifnullweight)){
          nw <- max(0, 1 - sum(prior))
          prior <- prior/(1-nw)
        } else {
          nw <- NULL
        }

        z <- c(exprres$z.g[gidx], zdf[sidx, "z"])
        X <- cbind(exprres$expr[, gidx], dat$G[, ldsidx])
        R <- cor(X)
        susieres <- susie_rss(z, R,
                              z_ld_weight = z_ld_weight,
                              L = L, prior_weights = prior,
                              null_weight = nw,
                              prior_variance = V,
                              estimate_prior_variance = !plugin_prior_variance,
                              coverage = coverage)

        # note: susie_rss use prior variance instead of scaled prior variance

        anno.gene <- cbind(exprres$gnames[gidx], exprres$chrom[gidx],  exprres$p0[gidx],
                           rep("gene", length(gidx)))
        anno.SNP <- cbind(dat$snp[ldsidx,], dat$chr[ldsidx,], dat$pos[ldsidx,],
                          rep("SNP", length(sidx)))

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
            X.idx <- X.idx[X.idx != susieres$null_index] # susie_rss' bug
            anno$cs_index[X.idx] <- cs_i #TODO: note this ignore the fact that some variant can belong to multiple CS
          }
        }

        outdf.rn <- cbind(anno, susieres$pip)
        colnames(outdf.rn)[9] <- "susie_pip"
        outdf.rn$mu2 <- colSums(susieres$mu2[,  1:ncol(X), drop = F]) #WARN: not sure for L>1, mu2 has an extra column (last) given null_weight.
        outdf.b.list[[rn]] <- outdf.rn
      }

      outdf.b <- do.call(rbind, outdf.b.list)
      outdf.b
    }

    prior.SNP <- mean(outdf[outdf[ , "type"] == "SNP", "susie_pip"])
    prior.gene <- mean(outdf[outdf[ , "type"] == "gene", "susie_pip"])
    print(c(prior.SNP, prior.gene))

    if (isTRUE(est_prior_variance)){
      outdf.g <- outdf[outdf[ , "type"] == "gene", ]
      outdf.s <- outdf[outdf[ , "type"] == "SNP", ]

      V.gene <- sum(outdf.g$susie_pip * outdf.g$mu2)/sum(outdf.g$susie_pip)
      V.SNP <- sum(outdf.s$susie_pip * outdf.s$mu2)/sum(outdf.s$susie_pip)
    }

    prior.SNP_rec[iter] <- prior.SNP
    prior.gene_rec[iter] <- prior.gene
    V.gene_rec[iter] <- V.gene
    V.SNP_rec[iter] <- V.SNP
    save(prior.gene_rec, prior.SNP_rec, V.gene_rec, V.SNP_rec, file = paste(outname, "susieIrssres.Rd", sep = "."))
    gc()
    write.table(outdf, file= paste0(outname, ".susieIrss.txt" ) , row.names=F, col.names=T, sep="\t", quote = F)
  }

  stopCluster(cl)

  return(list("prior.gene"= prior.gene,
              "prior.SNP" = prior.SNP,
              "V.gene" = V.gene,
              "V.SNP" = V.SNP))

}

