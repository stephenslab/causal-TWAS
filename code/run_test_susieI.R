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
       ", call.=FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir,"input_reformat.R"))

addHandler(writeToFile, file="run_test_susie.R.log", level='DEBUG')
loginfo('script started ... ')

PIPfilter <- 0.3
chunksize = 500000
Niter <- 10
ifshuffle <- T

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


# gres <- paste0(args[4], ".expr.txt")
# pipres <- read.table(gres, header =T)
#
# regions <- NULL
# chroms <- unique(pipres$chr)
# for (chrom in chroms){
#   pipres.chr <- pipres[pipres$chr == chrom, ]
#   stpos <- min(pipres.chr$p0) - chunksize/2
#   edpos <- max(pipres.chr$p1) + chunksize/2
#   p0 <- stpos + 0: floor((edpos - stpos)/chunksize) * chunksize
#   p1 <- stpos + 1: ceiling((edpos - stpos)/chunksize) * chunksize
#   itv <- cbind(p0, p1)
#   rPIP <- apply(itv, 1, function(x) sum(pipres.chr[x[1] < pipres.chr$p0 & pipres.chr$p0 < x[2], "PIP"]))
#   nCausal <- apply(itv, 1, function(x) sum(pipres.chr[x[1] < pipres.chr$p0 & pipres.chr$p0 < x[2], "ifcausal"]))
#   itv <- cbind(chrom, itv, rPIP, nCausal)
#   regions <- rbind(regions, itv)
# }
#
# loginfo("No. intervals: %s", nrow(regions))
# write.table(regions , file = paste0(args[4], ".rPIP.txt")  , row.names=F, col.names=T, sep="\t", quote = F)

regions <- read.table(paste0(args[4], ".rPIP.txt"), header = T)
print(nrow(regions))

gres <- paste0(args[4], ".expr.txt")
sres <- paste0(args[4], ".snp.txt")
pipres <- rbind(read.table(gres, header =T), read.table(sres, header =T))
cau <- pipres[pipres$ifcausal == 1, "name"]


regions <- regions[regions[, "rPIP"] > PIPfilter & regions[, "nCausal"] > 0,  ] # select regions
loginfo("No. intervals after PIP filter: %s", nrow(regions))

#TODO: join regions based on LD genes, add a column indicating which region should be linked
regions <- as.data.frame(regions)
regions$rname <- as.character(1:nrow(regions))

loginfo("prepare geno for each region ...")

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
    X.SNP <-  dat$G[ ,idx.SNP]
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
  }
}


for (rn in names(regionlist)){
  anno <- regionlist[[rn]][["anno"]]
  if (nrow(anno[anno$type == "gene",]) == 0) {
    regionlist[[rn]] <- NULL
  }
}
loginfo("limiting to regions with at least one gene: %s", length(regionlist))


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


outname <- args[5]
loginfo("susie started for %s", outname)

rm(dat); gc()
save.image(file = "temp.Rd")
saveRDS(regionlist, file = paste0(outname, "_regionlist.rds"))


prior.SNP <- NULL
prior.gene <- NULL
prior.SNP_rec <- rep(0, Niter)
prior.gene_rec <- rep(0, Niter)
for (iter in 1:Niter){
  print(iter)
  susieres <- list()
  snp.rpiplist <- list()
  gene.rpiplist <- list()
  outdf <- NULL
  for (rn in names(regionlist)){
    if (iter == 1) {
      prior <- NULL
    } else {
      prior <- rep(0, dim(regionlist[[rn]][["X"]])[2])
      prior[regionlist[[rn]][["anno"]][, "type"] == "SNP"] <- prior.SNP
      prior[regionlist[[rn]][["anno"]][, "type"] == "gene"] <- prior.gene
    }
    susieres[[rn]] <- susie(regionlist[[rn]][["X"]], phenores$Y, L=1, prior_weights = prior)
    snp.rpiplist[[rn]] <- sum(susieres[[rn]]$alpha[regionlist[[rn]][["anno"]][,"type"] == "SNP"])
    gene.rpiplist[[rn]] <- sum(susieres[[rn]]$alpha[regionlist[[rn]][["anno"]][,"type"] == "gene"])
    outdf <- rbind(outdf, cbind(regionlist[[rn]]$anno, t(susieres[[rn]]$alpha), susieres[[rn]]$pip))
  }
  colnames(outdf)[6:7] <- c("susie_alpha", "susie_pip")
  write.table(outdf , file= paste0(outname, ".", iter, ".txt" ) , row.names=F, col.names=T, sep="\t", quote = F)
  save( snp.rpiplist, gene.rpiplist, file = paste(outname, iter, "susieIres.Rd", sep = "."))
  prior.SNP <- mean(unlist(snp.rpiplist))
  prior.gene <- mean(unlist(gene.rpiplist))
  print(c(prior.SNP, prior.gene))
  prior.SNP_rec[iter] <- prior.SNP
  prior.gene_rec[iter] <- prior.gene
}

loginfo("susie done for %s", outname)

save(prior.gene_rec, prior.SNP_rec, file = paste(outname, "susieIres.Rd", sep = "."))




