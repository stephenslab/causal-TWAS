library(susieR)
library(logging)
library(tools)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("at least 5 arguments must be supplied:
       * genotype file name (or multiple file names in a .txt file)
       * expr Rd (or multiple file names in a .txt file, matching geno)
       * pheno Rd
       * filter file names in txt (put any string here if no filter, or match geno )
       * out file name
       * config file (optional)
       ", call.=FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir,"input_reformat.R"))
source(paste0(codedir,"susie_filter.R"))
source(paste0(codedir,"SER.R"))

addHandler(writeToFile, file="run_test_SER.R.log", level='DEBUG')
loginfo('script started ... ')

# default parameters
filtertype <- "mrash2s" # "nofilter": keep all regions
                        # or mrash: mr.ash for genes only
                        # or pgene: twas for genes only
                        # or p2: twas for genes and snps
                        # mrash2s: mr.ash2s for genes and snps
                        # rBF: regional BF for genes and snps
                        # rfdr: regional min fdr for genes and snps

PIPfilter <- 0.3 # for "mrash2s" or "mrash", regions with sum of PIP < PIPfilter will be filtered
prankfilter.gene <- 1 # for pgene or p2, lowest rank of gene/SNP pvalues > prankfilter that will be filtered.
prankfilter.snp <- 1 # for pgene or p2, lowest rank of gene/SNP pvalues > prankfilter that will be filtered.
rBFfilter <- 1 # for rBF, regional BF of genes & SNPs < rBFfilter will be filtered
rfdrfilter <- 1 # for rfdr, if region contains genes or snp with fdr < fdrfilter, it will be kept.


ifgene <- F # if limiting to regions with at least one gene.
ifca <- T # if limiting to regions with at least one causal signal
ifsingleca <- F # if limiting to regions with at most 1 causal signal
regionsfile <- "/home/simingz/causalTWAS/simulations/shared_files/chr1to22-500kb_fn.txt" # need to match input

Niter <- 30
ifnullweight <- F # if include/update null weight in iterations.
prior.gene_init <- NULL
prior.SNP_init <- NULL
L <- 1

Ncore <- 1

# user provided parameters
if (length(args) == 6){
  source(args[6])
}

loginfo("filter type: %s", filtertype)
loginfo("region file: %s", regionsfile)

outname <- args[5]

# prepare input files
pfile <- args[1]
if (file_ext(pfile) == "txt"){
  pfiles <- read.table(pfile, header =F, stringsAsFactors = F)[,1]
  efiles <- read.table(args[2], header =F, stringsAsFactors = F)[,1]
} else {
  pfiles <- pfile
  efiles <- args[2]
}

pfileRds <- paste0(pfiles, ".FBM.Rd")
phenofile <- args[3]

# read regions files
reg <- read.table(regionsfile, header = F, stringsAsFactors = F)
reglist <- list()
for (b in 1: nrow(reg)){
  reglist[[b]] <- read.table(reg[b,], header = T, stringsAsFactors = F)
}

loginfo("No. intervals: %s", sum(unlist(lapply(reglist, nrow))))

# load phenotype Rd file, variable: phenores
load(phenofile)
phenores$Y <-  phenores$Y - mean(phenores$Y)

regionlist <- list()
regionsall <- NULL
for (b in 1:length(pfiles)){
  # load genotype data
  load(pfileRds[b])

  # load expression Rd file, variable: exprres
  load(efiles[b])

  loginfo("No. intervals before filtering %s for batch %s: %s", filtertype, b, nrow(reglist[[b]]))

  # filter regions based on filtertype
  cau <- c(dat$snp[phenores$batch[[b]]$param$idx.cSNP,],
           exprres$gnames[phenores$batch[[b]]$param$idx.cgene])
  regout <- filter_region(reg = reglist[[b]], filtertype = filtertype, filterfile = args[4])
  regions <- regout[["filtered"]]
  regionsall <- rbind(regionsall, regout[["all"]])

  # filter regions based on ifca
  if (isTRUE(ifca)){
    loginfo("keep only causal regions")
    regions <- regions[regions[, "nCausal"] > 0,  ] # select regions
  }

  if (isTRUE(ifsingleca)){
    loginfo("keep only causal regions")
    regions <- regions[regions[, "nCausal"] <=1,  ] # select regions
  }

  loginfo("No. intervals after filtering %s for batch %s: %s", filtertype, b, nrow(regions))

  #TODO: join regions based on LD genes, add a column indicating which region should be linked
  regions$rname <- as.character(1:nrow(regions))

  regionlist[[b]] <- list()
  # prep geno for each region
  loginfo("prepare geno for each region ...")
  for (rn in unique(regions[, "rname"])){
    regions.n <- regions[regions[, "rname"] == rn, ,drop = F]
    for (i in 1:nrow(regions.n)){
      chr <- regions.n[i, "chrom"]
      p0 <- regions.n[i, "p0"]
      p1 <- regions.n[i, "p1"]

      idx.gene <- exprres$chrom == chr & exprres$p0 > p0 & exprres$p0 < p1
      idx.SNP <- dat$chr == chr & dat$pos > p0 & dat$pos < p1

      if (length(idx.SNP[idx.SNP]) == 0 & length(idx.gene[idx.gene]) == 0) {next}

      regionlist[[b]][[rn]][["gidx"]] <- idx.gene
      regionlist[[b]][[rn]][["sidx"]] <- idx.SNP

    }
  }

  loginfo("No. intervals with genotype for batch %s: %s",  b, length(regionlist[[b]]))
  gc()

  if (isTRUE(ifgene)){
    loginfo("limiting to regions with at least one gene for batch %s", b)

    for (rn in names(regionlist[[b]])){
      idx.gene <- regionlist[[b]][[rn]][["gidx"]]
      if (length(idx.gene[idx.gene]) == 0) {
        regionlist[[b]][[rn]] <- NULL
      }
    }
  }
}

write.table(regionsall, file= paste0(outname,".", filtertype, ".r.txt" ) , row.names=F, col.names=T, sep="\t", quote = F)

# start SERI
loginfo("get log BF for %s", outname)

cl <- makeCluster(Ncore,outfile="")
registerDoParallel(cl)

V.g <- phenores$batch[[1]]$param$sigma_beta ** 2
V.s <- phenores$batch[[1]]$param$sigma_theta ** 2

lbfdf <- foreach (b = 1:22, .combine = "rbind",.packages = c("susieR", "bigstatsr")) %dopar% {

    # load genotype data
    load(pfileRds[b])

    # load expression Rd file, variable: exprres
    load(efiles[b])

    cau <- c(dat$snp[phenores$batch[[b]]$param$idx.cSNP,],
             exprres$gnames[phenores$batch[[b]]$param$idx.cgene])

    # get lbf for each region (in parallel by batch)
    lbfdflist.b <- list()
    for (rn in names(regionlist[[b]])) {

      gidx <- regionlist[[b]][[rn]][["gidx"]]
      sidx <- regionlist[[b]][[rn]][["sidx"]]

      X.g <- exprres$expr[, gidx, drop =F]
      X.s <- dat$G[, sidx, drop=F]

      X <- cbind(X.g, X.s)

      V <- c(rep(V.g, dim(X.g)[2]), rep(V.s, dim(X.s)[2]))

      lbfout <- get_lbf(Y =  phenores$Y, X = X, V = V, residual_variance = var(phenores$Y)[1])

      anno.gene <- cbind(exprres$gnames[gidx], exprres$chrom[gidx],  exprres$p0[gidx],
                         rep("gene", length(gidx[gidx])))
      anno.SNP <- cbind(dat$snp[sidx,], dat$chr[sidx,], dat$pos[sidx,],
                        rep("SNP", length(sidx[sidx])))

      anno <- rbind(anno.gene, anno.SNP)
      colnames(anno) <-  c("name", "chr", "pos", "type")
      anno <- as.data.frame(anno)
      anno <- cbind(anno, lbfout)
      anno$ifcausal <- ifelse(anno$name %in% cau, 1, 0)
      anno$b <- b
      anno$rn <- rn

      lbfdflist.b[[rn]] <- anno
    }
    lbfdf.b <- do.call(rbind, lbfdflist.b)
    lbfdf.b
}

save(lbfdf, file = paste0(outname, ".lbfdf.Rd"))
loginfo("log BF done for %s", outname)


loginfo("SER started for %s", outname)

prior.SNP_rec <- rep(0, Niter)
prior.gene_rec <- rep(0, Niter)
loglik_rec <- rep(0, Niter)

for (iter in 1:Niter){

  loginfo("run iteration %s", iter)

  SERpip <- foreach (b = 1:22, .combine = rbind) %dopar% {

    pip.b.list <- list()

    lbfdf.b <- lbfdf[lbfdf$b ==b, ]

    # run SER for each region
    for (rn in unique(lbfdf.b[, "rn"])) {

      df <- lbfdf.b[lbfdf.b$rn == rn, ]

      n.g <- nrow(df[df$type == "gene", , drop =F])
      n.s <- nrow(df[df$type == "SNP", , drop =F])
      p <- n.g + n.s
      prior <- rep(0, p)

      if (is.null(prior.gene_init) | is.null(prior.SNP_init)){
        prior.gene_init <- 1/p
        prior.SNP_init <- 1/p
      }

      if (iter == 1) {
        prior <- c(rep(prior.gene_init, n.g), rep(prior.SNP_init, n.s))
      } else {
        prior <- c(rep(prior.gene, n.g), rep(prior.SNP, n.s))
      }

      nw <- max(0, 1 - sum(prior))

      lbf <- as.numeric(df$lbf)

      if (isTRUE(ifnullweight)){
        prior = prior/(1-nw)
      } else {
        nw = 0
      }

      SERres <- SER(lbf, prior_weights = prior, null_weight = nw)

      pip.b.list[[rn]] <- cbind(SERres$alpha, prior, nw,  SERres$loglik)
      colnames(pip.b.list[[rn]]) <- c("SERpip", "prior", "region_nw", "region_loglik")
    }

    pip.b <- do.call(rbind, pip.b.list)
    pip.b
  }

  outdf <- cbind(lbfdf, SERpip)
  prior.gene <- mean(outdf[outdf[ , "type"] == "gene", "SERpip"])
  prior.SNP <- mean(outdf[outdf[ , "type"] == "SNP", "SERpip"])
  loglikall <- sum(unique(outdf[, c("b", "rn", "region_loglik")])[, "region_loglik"])
  print(c(prior.gene, prior.SNP, loglikall))

  prior.SNP_rec[iter] <- prior.SNP
  prior.gene_rec[iter] <- prior.gene
  loglik_rec[iter] <- loglikall
  save(prior.gene_rec, prior.SNP_rec, loglik_rec, file = paste(outname, "SERIres.Rd", sep = "."))
  gc()
  write.table(outdf, file= paste0(outname, ".SERI.txt" ) , row.names=F, col.names=T, sep="\t", quote = F)
}

stopCluster(cl)

loginfo("susie done for %s", outname)





