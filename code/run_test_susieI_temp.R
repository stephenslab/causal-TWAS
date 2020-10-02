library(susieR)
library(dplyr)
library(logging)
library(tools)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("at least 5 arguments must be supplied:
       * genotype file name (or multiple file names in a .txt file)
       * expr Rd (or multiple file names in a .txt file, matching genotype file)
       * pheno Rd (or multiple file names in a .txt file, matching genotype file)
       * filter file (put any string here if no filter)
       * out file name
       * config file (optional)
       ", call.=FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir,"input_reformat.R"))
source(paste0(codedir,"susie_filter.R"))

addHandler(writeToFile, file="run_test_susie.R.log", level='DEBUG')
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
regionsfile <- "/home/simingz/causalTWAS/simulations/shared_files/chr17-22-500kb_fn.txt" # need to match input

Niter <- 30
ifnullweight <- F # if include/update null weight in iterations.
nullweight_init <- NULL # initializing value
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

pfile <- args[1]
if (file_ext(pfile) == "txt"){
  pfiles <- read.table(pfile, header =F, stringsAsFactors = F)[,1]
  efiles <- read.table(args[2], header =F, stringsAsFactors = F)[,1]
  phenofiles <- read.table(args[3], header =F, stringsAsFactors = F)[,1]
} else {
  pfiles <- pfile
  efiles <- args[2]
  phenofiles <- args[3]
}
pfileRds <- paste0(drop_ext(pfiles), ".FBM.Rd")
efileRds <- paste0(drop_ext(efiles), ".FBM.scaled.rds")

outname <- args[5]

# read regions files
reg <- read.table(regionsfile, header = F, stringsAsFactors = F)
reglist <- list()
for (b in 1: nrow(reg)){
  reglist[[b]] <- read.table(reg[b,], header = T, stringsAsFactors = F)
}

loginfo("No. intervals: %s", sum(unlist(lapply(reglist, nrow))))

regionlist <- list()


for (b in 1:length(pfiles)){
  # load genotype data
  load(pfileRds[b])

  # load expression Rd file, variable: exprres
  load(efiles[b])
  if (!file.exists(efileRds[b])) {
    expr <- scaleRcpp(exprres$expr)
    as_FBM(expr, backingfile = paste0(drop_ext(efiles[b]), ".FBM.scaled"), is_read_only = T)$save()
  } else {
    expr <- readRDS(efileRds[b])
  }
  exprres$gnames <- colnames(exprres$expr)

  # load phenotype Rd file, variable: phenores
  load(phenofiles[b])
  phenores$Y <-  phenores$Y - mean(phenores$Y) # Y is the same, but param is different for each b.

  loginfo("No. intervals before filtering %s for batch %s: %s", filtertype, b, nrow(reglist[[b]]))

  # filter regions based on filtertype
  regions <- filter_region(reg = reglist[[b]], filtertype = filtertype, filterfile = paste0(args[4], "-B", b), outname = paste0(outname, ".B", b))

  # filter regions based on ifca
  if (isTRUE(ifca)){
    loginfo("keep only causal regions")
    regions <- regions[regions[, "nCausal"] > 0,  ] # select regions
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

# start susieI
loginfo("susie started for %s", outname)

prior.SNP_rec <- rep(0, Niter)
prior.gene_rec <- rep(0, Niter)

for (iter in 1:Niter){
  loginfo("run iteration ", iter)

  snp.rpiplist <- list()
  gene.rpiplist <- list()
  outdf <- NULL
  for (b in 1:length(regionlist)){

    # load genotype data
    load(pfileRds[b])

    # load expression Rd file, variable: exprres
    load(efiles[b])
    exprres$gnames <- colnames(exprres$expr)
    exprres$expr <- NULL; gc()
    exprres$exprlist <- NULL
    exprres$qclist <- NULL
    expr <- readRDS(efileRds[b])

    # load phenotype Rd file, variable: phenores
    load(phenofiles[b])
    phenores$Y <-  phenores$Y - mean(phenores$Y)

    cau <- c(dat$snp[phenores$param$idx.cSNP,], exprres$gnames[phenores$param$idx.cgene])

    snp.rpiplist[[b]] <- list()
    gene.rpiplist[[b]] <- list()

    # run susie for each region (in parallel)
    cl <- makeCluster(Ncore)
    registerDoParallel(cl)

    outdf.b <- foreach(rn = names(regionlist[[b]]), .combine = "rbind") %dopar% {
      print(c(b, rn))
      gidx <- regionlist[[b]][[rn]][["gidx"]]
      sidx <- regionlist[[b]][[rn]][["sidx"]]
      p <- length(gidx[gidx]) + length(sidx[sidx])
      prior <- rep(0, p)

      if (is.null(prior.gene_init) | is.null(prior.SNP_init)){
        prior.gene_init <- 1/p
        prior.SNP_init <- 1/p
      }

      if (iter == 1) {
        prior[1: length(gidx[gidx])] <- prior.gene_init
        prior[(length(gidx[gidx])+ 1): p ] <- prior.SNP_init
        nw <- nullweight_init
      } else {
        prior[1: length(gidx[gidx])] <- prior.gene
        prior[ (length(gidx[gidx])+ 1): p] <- prior.SNP
        nw <- max(0, 1 - sum(prior))
      }

      X <- cbind(expr[, gidx], dat$G[, sidx])

      if (isTRUE(ifnullweight)){
        susieres <- susie( X , phenores$Y, L = L, prior_weights = prior, null_weight = nw)
      } else {
        susieres <- susie(X, phenores$Y, L = L, prior_weights = prior)
      }

      gene.rpiplist[[b]][[rn]] <- susieres$pip[1: length(gidx[gidx])]
      snp.rpiplist[[b]][[rn]] <- susieres$pip[(length(gidx[gidx])+ 1): p]


      anno.gene <- cbind(exprres$gnames[gidx], exprres$chrom[gidx],  exprres$p0[gidx],
                         rep("gene", length(gidx[gidx])))
      anno.SNP <- cbind(dat$snp[sidx,], dat$chr[sidx,], dat$pos[sidx,], "SNP")

      anno <- rbind(anno.gene, anno.SNP)
      colnames(anno) <-  c("name", "chr", "pos", "type")
      anno <- as.data.frame(anno)
      anno$ifcausal <- ifelse(anno$name %in% cau, 1, 0)

      outdf.rn <- cbind(anno, susieres$pip)
      colnames(outdf.rn)[6] <- "susie_pip"

      outdf.rn
    }

    outdf <- rbind(outdf, outdf.b)
  }

  prior.SNP <- mean(unlist(snp.rpiplist))
  prior.gene <- mean(unlist(gene.rpiplist))
  print(c(prior.SNP, prior.gene))

  prior.SNP_rec[iter] <- prior.SNP
  prior.gene_rec[iter] <- prior.gene
  save(prior.gene_rec, prior.SNP_rec, file = paste(outname, "susieIres.Rd", sep = "."))
  gc()
  write.table(outdf, file= paste0(outname, ".susieI.txt" ) , row.names=F, col.names=T, sep="\t", quote = F)
}

loginfo("susie done for %s", outname)




