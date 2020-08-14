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

PIPfilter <- 0.02
chunksize = 500000
Niter <- 10

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


gres <- paste0(args[4], ".expr.txt")
pipres <- read.table(gres, header =T)

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

loginfo("No. intervals: %s", nrow(regions))
write.table(regions , file = paste0(args[4], ".rPIP.txt")  , row.names=F, col.names=T, sep="\t", quote = F)

regions <- regions[regions[, "rPIP"] > PIPfilter,  ]
loginfo("No. intervals for after PIP filter: %s", nrow(regions))

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

    idx.gene <- exprres$chrom == chr & exprres$p0 > p0 & exprres$p1 < p1
    idx.SNP <- dat$chr == chr & dat$pos > p0 & dat$pos < p1

    X.gene <- exprres$expr[ , idx.gene, drop = F]
    X.SNP <-  dat$G[ ,idx.SNP]
    X <- cbind(X.gene, X.SNP)

    anno.gene <- cbind(colnames(X.gene), exprres$chrom[idx.gene],  exprres$p0[idx.gene], "gene")
    anno.SNP <- cbind(dat$snp[idx.SNP,], dat$chr[idx.SNP,], dat$pos[idx.SNP,], "SNP")
    anno <- rbind(anno.gene, anno.SNP)
    colnames(anno) <-  c("name", "chr", "pos", "type")
    regionlist[[rn]][["X"]] <- cbind(regionlist[[rn]][["X"]], X)
    regionlist[[rn]][["anno"]] <- rbind(regionlist[[rn]][["anno"]], anno)
  }
}

rm(dat);gc()
save.image(file = "temp.Rd")

outname <- args[5]
loginfo("susie started for %s", outname)

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

  write.table( outdf , file= paste0(outname, ".", iter, ".txt" ) , row.names=F, col.names=T, sep="\t", quote = F)
  prior.SNP <- mean(unlist(snp.rpiplist))
  prior.gene <- mean(unlist(gene.rpiplist))
  print(c(prior.SNP, prior.gene))
  prior.SNP_rec[iter] <- prior.SNP
  prior.gene_rec[iter] <- prior.gene
}

loginfo("susie done for %s", outname)

save(prior.gene_rec, prior.SNP_rec, file = paste(outname, "susieIres.Rd", sep = "."))




