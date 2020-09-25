library(susieR)
library(dplyr)
library(logging)
library(tools)

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

PIPfilter <- 0.3 # for "mrash2s" or "mrash", regions with sum of PIP < PIPfilter will be filtered
prankfilter.gene <- 1 # for pgene or p2, lowest rank of gene/SNP pvalues > prankfilter that will be filtered.
prankfilter.snp <- 1 # for pgene or p2, lowest rank of gene/SNP pvalues > prankfilter that will be filtered.
rBFfilter <- 1 # for rBF, regional BF of genes & SNPs < rBFfilter will be filtered

ifgene <- F # if limiting to regions with at least one gene.
ifca <- T # if limiting to regions with at least one causal signal
regionsfile <- "/home/simingz/causalTWAS/simulations/shared_files/chr17-22-500kb_B_bins.txt"

Niter <- 30
ifnullweight <- F # if include/update null weight in iterations.
nullweight_init <- NULL # initializing value
prior.gene_init <- NULL
prior.SNP_init <- NULL
L <- 1

# user provided parameters
if (length(args) == 6){
  source(args[6])
}

loginfo("regional PIP cut: %s ", PIPfilter)
loginfo("region size: %s", chunksize)


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

m <- as_FBM(scaleRcpp(exprres$expr), backingfile = , is_read_only = T)$save()


# read regions files
reg <- read.table(regionsfile, header = T)
loginfo("No. intervals: %s", nrow(reg))
regionlist <- list()

for (b in 1:length(pfiles)){
  # load genotype data
  load(pfileRds[b])

  # load expression Rd file, variable: exprres
  if (!file.exists(efileRds[b])) {
    load(efiles[b])
    expr <- scaleRcpp(exprres$expr)
    as_FBM(expr, backingfile = paste0(drop_ext(efiles[b]), ".FBM.scaled"), is_read_only = T)$save()
  }
  dat$expr <- readRDS(efilesRds[b])

  # load phenotype Rd file, variable: phenores
  load(phenofiles[b])
  phenores$Y <-  phenores$Y - mean(phenores$Y)

  outname <- args[5]

  # filter regions based on filtertype
  regions <- filter_region(reg = reg, filtertype = filtertype, filterfile = args[4])

  # filter regions based on ifca
  if (isTRUE(ifca)){
    loginfo("keep only causal regions")
    regions <- regions[regions[, "nCausal"] > 0,  ] # select regions
  }

  loginfo("No. intervals after filtering %s for batch %s: %s", filtertype, b, nrow(regions))

  #TODO: join regions based on LD genes, add a column indicating which region should be linked
  regions <- as.data.frame(regions)
  regions$rname <- as.character(1:nrow(regions))

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

      regionlist[[b]][[rn]][["gidx"]] <- idx.gene
      regionlist[[b]][[rn]][["sidx"]] <- idx.SNP
    }
  }

  rm(dat); gc()

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
outname <- args[5]
loginfo("susie started for %s", outname)

prior.SNP_rec <- rep(0, Niter)
prior.gene_rec <- rep(0, Niter)
elbo_rec <- list()
for (iter in 1:Niter){
  print(iter)
  susieres <- list()
  snp.rpiplist <- list()
  gene.rpiplist <- list()

  outdf <- NULL
  for (b in 1:length(regionlist)){

    for (rn in names(regionlist[[b]])){

      p <- ncol(regionlist[[b]][[rn]][["X"]])
      prior <- rep(0, p)

      if (is.null(prior.gene_init) | is.null(prior.SNP_init)){
        prior.gene_init <- 1/p
        prior.SNP_init <- 1/p
      }

      if (iter == 1) {
        prior[regionlist[[b]][[rn]][["anno"]][, "type"] == "SNP"] <- prior.SNP_init
        prior[regionlist[[b]][[rn]][["anno"]][, "type"] == "gene"] <- prior.gene_init
        nw <- nullweight_init
      } else {
        prior[regionlist[[rn]][["anno"]][, "type"] == "SNP"] <- prior.SNP
        prior[regionlist[[rn]][["anno"]][, "type"] == "gene"] <- prior.gene
        nw <- max(0, 1 - sum(prior))
      }

      if (isTRUE(ifnullweight)){
        susieres[[rn]] <- susie(regionlist[[rn]][["X"]], phenores$Y, L = L, prior_weights = prior, null_weight = nw)
        if (susieres[[rn]]$null_index !=0 ){
          susieres[[rn]]$alpha <- susieres[[rn]]$alpha[, - susieres[[rn]]$null_index, drop = F]
        }
      } else {
        susieres[[rn]] <- susie(regionlist[[rn]][["X"]], phenores$Y, L = L, prior_weights = prior)
      }

      snp.rpiplist[[rn]] <- susieres[[rn]]$pip[regionlist[[rn]][["anno"]][,"type"] == "SNP"]
      gene.rpiplist[[rn]] <- susieres[[rn]]$pip[regionlist[[rn]][["anno"]][,"type"] == "gene"]

      outdf.rn <- cbind(regionlist[[rn]]$anno, t(susieres[[rn]]$alpha), susieres[[rn]]$pip)
      colnames(outdf.rn)[6:ncol(outdf.rn)] <- c(paste0("susie_alpha", 1:nrow(susieres[[rn]]$alpha)), "susie_pip")

      outdf <- dplyr::bind_rows(outdf, outdf.rn)
    }

    prior.SNP <- mean(unlist(snp.rpiplist))
    prior.gene <- mean(unlist(gene.rpiplist))
    print(c(prior.SNP, prior.gene))

    prior.SNP_rec[iter] <- prior.SNP
    prior.gene_rec[iter] <- prior.gene
    elbo_rec[[iter]] <- lapply(susieres, '[[', "elbo")
    save(prior.gene_rec, prior.SNP_rec, elbo_rec, file = paste(outname, "susieIres.Rd", sep = "."))
    gc()
  }
}

loginfo("susie done for %s", outname)

write.table(outdf, file= paste0(outname, ".susieI.txt" ) , row.names=F, col.names=T, sep="\t", quote = F)


