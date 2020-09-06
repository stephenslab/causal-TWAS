library(susieR)
library(dplyr)
library(logging)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
  stop("at least 6 arguments must be supplied:
       * genotype file name
       * expr Rd
       * pheno Rd
       * param txt
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
ifcausal <- T
L = 1

# user provided parameters
if (length(args) == 7){
  source(args[7])
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

# causal signal names
cau <- c(dat$snp[phenores$param$idx.cSNP,], colnames(exprres$expr)[phenores$param$idx.cgene])

# run susie using given priors
param <- read.table(args[4], header = T, row.names = 1)

prior.gene <- param["gene.pi1", "estimated"]
prior.SNP <- param["snp.pi1", "estimated"]

gres <- paste0(args[5], ".expr.txt")
sres <- paste0(args[5], ".snp.txt")
pipres <- rbind(read.table(gres, header =T),
                 read.table(sres, header =T))

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
  nCausal <- apply(itv, 1, function(x) {gnames <- pipres.chr[x[1] < pipres.chr$p0 & pipres.chr$p0 < x[2], "name"]; sum(ifelse(gnames %in% cau, 1, 0))})
  itv <- cbind(chrom, itv, rPIP, nCausal)
  regions <- rbind(regions, itv)
}

loginfo("No. intervals for chr %s : %s", chrom, nrow(regions))
write.table(regions, file = paste0(args[5], ".rPIP.txt"), row.names=F, col.names=T, sep="\t", quote = F)

regions <- regions[regions[, "rPIP"] > PIPfilter,  ]

if (isTRUE(ifcausal)){
  regions <- regions[regions[, "nCausal"] > 0,  ]
}

loginfo("No. intervals for chr %s after PIP filter: %s", chrom, nrow(regions))

outname <- args[6]


loginfo("L = %s for susie run", L)

loginfo("susie started for %s", outname)

outlist <- list()

for (i in 1:nrow(regions)){

  chr <- regions[i, "chrom"]
  p0 <- regions[i, "p0"]
  p1 <- regions[i, "p1"]

  name <- paste(chr,p0,p1, sep="-")

  idx.gene <- exprres$chrom == chr & exprres$p0 > p0 & exprres$p1 < p1
  idx.SNP <- dat$chr == chr & dat$pos > p0 & dat$pos < p1

  X.gene <- exprres$expr[ , idx.gene, drop = F]
  X.SNP <-  dat$G[ , idx.SNP, drop = F]

  if (ncol(X.gene) + ncol(X.SNP) == 0) next

  prior <- c(rep(prior.gene, dim(X.gene)[2]), rep(prior.SNP, dim(X.SNP)[2]))
  wgt_null <- max(0, prod(1 - prior))

  #-----------------run susie------------
  susieres <- susie(cbind(X.gene, X.SNP), phenores$Y, L=L, prior_weights = prior)
  susieres.null <- susie(cbind(X.gene, X.SNP), phenores$Y, L=L)
  susieres.w0 <- susie(cbind(X.gene, X.SNP), phenores$Y, L=L, null_weight = wgt_null, prior_weights = prior)
  #--------------------------------------

  anno.gene <- cbind(colnames(X.gene), exprres$chrom[idx.gene],  exprres$p0[idx.gene])
  anno.SNP <- cbind(dat$snp[idx.SNP,], dat$chr[idx.SNP,], dat$pos[idx.SNP,])
  anno <- rbind(anno.gene, anno.SNP)

  outdf <- cbind(anno, susieres$pip, susieres.null$pip, susieres.w0$pip, "SNP")
  colnames(outdf) <- c("name", "chr", "pos", "pip", "pip.null", "pip.w0", "type")
  outdf[1:nrow(anno.gene), "type"] <- "gene"

  outlist[[name]] <- outdf
  rm(susieres,  susieres.null,  susieres.w0);gc()
}

loginfo("susie done for %s", outname)

saveRDS(outlist, file = paste(outname, "susieres.rds", sep = "."))

outgdf <- do.call(rbind, outlist)
outgdf <- as.data.frame(outgdf)

outgdf$ifcausal <- ifelse(outgdf$name %in% cau, 1, 0)

write.table(outgdf[outgdf[, "type"] == "gene", ], file= paste(outname, "susieres.expr.txt", sep = ".")  , row.names=F, col.names= T, sep="\t", quote = F)
write.table(outgdf[outgdf[, "type"] == "SNP", ], file= paste(outname, "susieres.snp.txt", sep = ".")  , row.names=F, col.names= T, sep="\t", quote = F)




