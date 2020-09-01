library(susieR)
library(dplyr)
library(logging)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("at least 5 arguments must be supplied:
       * genotype file name
       * expr Rd
       * pheno Rd
       * mr.ash.res file
       * out file name
       * config file (optional)
       ", call.=FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir,"input_reformat.R"))

addHandler(writeToFile, file="run_test_susie.R.log", level='DEBUG')
loginfo('script started ... ')

# default parameters
PIPfilter <- 0.3
chunksize = 500000
ifshuffle <- F # if shuffle genes to different regions to break LD
ifgene <- F # if limiting to regions with at least one gene.
ifca <- T # if limiting to regions with at least one causal signal

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

# load genotype data
pfile <- args[1]
pfileRd <- paste0(drop_ext(pfile), ".FBM.Rd")
load(pfileRd)

dat$G <- dat$G[]

# load expression Rd file, variable: exprres
load(args[2])
dat$expr <-  scaleRcpp(exprres$expr)

# load phenotype Rd file, variable: phenores
load(args[3])
phenores$Y <-  phenores$Y - mean(phenores$Y)

# prep regions
regf <- paste0(args[4], ".rPIP.txt")
if (!file.exists(regf)) {
  gres <- paste0(args[4], ".expr.txt")
  sres <- paste0(args[4], ".snp.txt") # may not need sres
  pipres <- rbind(read.table(gres, header =T), read.table(sres, header =T))

  regions <- NULL
  chroms <- unique(pipres$chr)
  for (chrom in chroms){
    pipres.chr <- pipres[pipres$chr == chrom, ]
    stpos <- min(pipres.chr$p0) - chunksize/2
    edpos <- max(pipres.chr$p1) + chunksize/2
    p0 <- stpos + 0: floor((edpos - stpos)/chunksize) * chunksize
    p1 <- stpos + 1: ceiling((edpos - stpos)/chunksize) * chunksize
    itv <- cbind(p0, p1)
    rPIP <- apply(itv, 1, function(x) sum(pipres.chr[x[1] < pipres.chr$p0 & pipres.chr$p0 < x[2], "PIP"]))
    nCausal <- apply(itv, 1, function(x) sum(pipres.chr[x[1] < pipres.chr$p0 & pipres.chr$p0 < x[2], "ifcausal"]))
    itv <- cbind(chrom, itv, rPIP, nCausal)
    regions <- rbind(regions, itv)
  }
  write.table(regions , file = paste0(args[4], ".rPIP.txt")  , row.names=F, col.names=T, sep="\t", quote = F)
} else{
  regions <- read.table(regf, header = T)
}
loginfo("No. intervals: %s", nrow(regions))

# filter regions based on rPIP and others
if (isTRUE(ifca)){
  regions <- regions[regions[, "rPIP"] > PIPfilter & regions[, "nCausal"] > 0,  ] # select regions
} else {
  regions <- regions[regions[, "rPIP"] > PIPfilter,  ] # select regions
}
loginfo("No. intervals after PIP filter: %s", nrow(regions))


#TODO: join regions based on LD genes, add a column indicating which region should be linked
regions <- as.data.frame(regions)
regions$rname <- as.character(1:nrow(regions))

# prep geno for each region
loginfo("prepare geno for each region ...")
cau <- c(dat$snp[phenores$param$idx.cSNP,], colnames(exprres$expr)[phenores$param$idx.cgene])
regionlist <- list()
for (rn in unique(regions[, "rname"])){
  regions.n <- regions[regions[, "rname"] == rn, ,drop = F]
  regionlist[[rn]] <- list(X = NULL, anno = NULL)
  for (i in 1:nrow(regions.n)){
    chr <- regions.n[i, "chrom"]
    p0 <- regions.n[i, "p0"]
    p1 <- regions.n[i, "p1"]

    name <- paste(chr, p0, p1, sep = "-")

    idx.gene <- exprres$chrom == chr & exprres$p0 > p0 & exprres$p0 < p1
    idx.SNP <- dat$chr == chr & dat$pos > p0 & dat$pos < p1

    X.gene <- exprres$expr[ , idx.gene, drop = F]
    X.SNP <-  dat$G[ ,idx.SNP, drop = F]
    X <- cbind(X.gene, X.SNP)

    anno.gene <- cbind(colnames(X.gene), exprres$chrom[idx.gene],  exprres$p0[idx.gene],
                       rep("gene", length(idx.gene[idx.gene])))
    anno.SNP <- cbind(dat$snp[idx.SNP,], dat$chr[idx.SNP,], dat$pos[idx.SNP,], "SNP")
    anno <- rbind(anno.gene, anno.SNP)
    colnames(anno) <-  c("name", "chr", "pos", "type")
    anno <- as.data.frame(anno)
    anno$ifcausal <- ifelse(anno$name %in% cau, 1, 0)

    regionlist[[rn]][["X"]] <- cbind(regionlist[[rn]][["X"]], X)
    regionlist[[rn]][["anno"]] <- rbind(regionlist[[rn]][["anno"]], anno)
    gc()
  }
}

rm(dat); gc()

if (isTRUE(ifgene)){
  loginfo("limiting to regions with at least one gene: %s", length(regionlist))

  for (rn in names(regionlist)){
    anno <- regionlist[[rn]][["anno"]]
    if (nrow(anno[anno$type == "gene",]) == 0) {
      regionlist[[rn]] <- NULL
    }
  }

}

if (isTRUE(ifshuffle)){
  loginfo("shuffling genes to different regions")

  regionlist.sh <- list()
  rn.sh <- sample(names(regionlist))

  for (i in 1:length(regionlist)){
    rn.snp <- names(regionlist)[i]
    rn.gene <- rn.sh[i]
    anno.snp <- regionlist[[rn.snp]][["anno"]]
    anno.gene <- regionlist[[rn.gene]][["anno"]]
    rn <- rn.snp
    regionlist.sh[[rn]][['X']] <- cbind(regionlist[[rn.gene]][['X']][, anno.gene$type == "gene"],
                                        regionlist[[rn.snp]][['X']][, anno.snp$type == "SNP"])

    regionlist.sh[[rn]][['anno']] <- rbind(anno.gene[anno.gene$type == "gene", ],
                                           anno.snp[anno.snp$type == "SNP", ])
  }
  regionlist <- regionlist.sh
}


# start susieI
outname <- args[5]
loginfo("susie started for %s", outname)

prior.SNP_rec <- rep(0, Niter)
prior.gene_rec <- rep(0, Niter)
for (iter in 1:Niter){
  susieres <- list()
  snp.rpiplist <- list()
  gene.rpiplist <- list()

  outdf <- NULL
  for (rn in names(regionlist)){
    p <- ncol(regionlist[[rn]][["X"]])
    prior <- rep(0, p)

    if (is.null(prior.gene_init) | is.null(prior.SNP_init)){
      prior.gene_init <- 1/p
      prior.SNP_init <- 1/p
    }

    if (iter == 1) {
      prior[regionlist[[rn]][["anno"]][, "type"] == "SNP"] <- prior.SNP_init
      prior[regionlist[[rn]][["anno"]][, "type"] == "gene"] <- prior.gene_init
      nw <- nullweight_init
    } else {
      prior[regionlist[[rn]][["anno"]][, "type"] == "SNP"] <- prior.SNP
      prior[regionlist[[rn]][["anno"]][, "type"] == "gene"] <- prior.gene
      nw <- max(0, prod(1 - prior))
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
  save(prior.gene_rec, prior.SNP_rec, file = paste(outname, "susieIres.Rd", sep = "."))
  gc()
}

loginfo("susie done for %s", outname)

write.table(outdf , file= paste0(outname, ".susieI.txt" ) , row.names=F, col.names=T, sep="\t", quote = F)

