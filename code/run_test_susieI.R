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
source(paste0(codedir,"susieI_func.R"))

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
regionsfile <- "/home/simingz/causalTWAS/simulations/shared_files/chr1to22-500kb_fn.txt" # need to match input

Niter <- 30
ifnullweight <- F # if include/update null weight in iterations.
prior.gene_init <- NULL
prior.SNP_init <- NULL
L <- 1
plugin_prior_variance <- F # using truth when T for now
coverage <- 0.95 # credible set coverage

Niter2 <- 20
filterandrerun <- F
P2cut <- 0.2

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
  cau <- #c(dat$snp[phenores$batch[[b]]$param$idx.cSNP,],
           exprres$gnames[phenores$batch[[b]]$param$idx.cgene] #)
  regout <- filter_region(reg = reglist[[b]], filtertype = filtertype, filterfile = args[4])
  regions <- regout[["filtered"]]
  regionsall <- rbind(regionsall, regout[["all"]])

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

write.table(regionsall, file= paste0(outname,".", filtertype, ".r.txt" ) , row.names=F, col.names=T, sep="\t", quote = F)

# start susieI
loginfo("susie started for %s", outname)

V.g <- phenores$batch[[1]]$param$sigma_beta ** 2 # prior variance for gene
V.s <- phenores$batch[[1]]$param$sigma_theta ** 2 # prior variance for SNP

pars <- susieI(prior.gene_init = prior.gene_init,
                   prior.SNP_init = prior.SNP_init,
                   regionlist = regionlist,
                   outname = outname,
                   niter = Niter,
                   Ncore = Ncore)



if (isTRUE(filterandrerun)){

  prior.gene_init2 <- pars[["prior.gene"]]
  prior.SNP_init2 <- pars[["prior.SNP"]]

  regionlist2 <- regionlist

  for (b in 1:length(pfiles)){

    for (rn in names(regionlist[[b]])) {

      gidx <- regionlist[[b]][[rn]][["gidx"]]
      sidx <- regionlist[[b]][[rn]][["sidx"]]
      p.g <- length(gidx[gidx])
      p.s <- length(sidx[sidx])
      P2 <- 1 - (1- prior.gene_init2)**p.g * (1- prior.SNP_init2)**p.s -
         p.g * prior.gene_init2 * (1 - prior.gene_init2) ** (p.g - 1) * (1- prior.SNP_init2) ** p.s -
         p.s * (1 - prior.gene_init2) ** p.g * prior.SNP_init2 * (1- prior.SNP_init2) ** (p.s - 1)

      if (P2 >= P2cut){
        regionlist2[[b]][[rn]] <- NULL
      }
    }
  }

  loginfo("Blocks are filtered for parameter estimation: %s blocks left",  sum(unlist(lapply(lapply(regionlist2, names), length))))

  pars2 <- susieI(prior.gene_init = prior.gene_init2,
                 prior.SNP_init = prior.SNP_init2,
                 regionlist = regionlist2,
                 outname = paste0(outname, ".fl"),
                 niter = Niter2,
                 Ncore = Ncore)

  loginfo("Parameter estimation done for %s", outname)

  susieI(prior.gene_init = pars2[["prior.gene"]],
         prior.SNP_init = pars2[["prior.SNP"]],
         regionlist = regionlist,
         outname = paste0(outname, ".flrerun"),
         niter = 1,
         Ncore = Ncore)
}


loginfo("susie done for %s", outname)

# temp <- lapply(regionlist2, names)
# alist <- list()
# for (i in 1:22){
#   dt.b <- dt[b == i,]
#   alist[[i]] <- dt.b[rn %in% as.numeric(temp[[i]]),]
# }
# a <- do.call(rbind, alist)
# nrow(a[type == "SNP" & ifcausal ==1])/nrow(a[type == "SNP"])
# nrow(a[type == "gene" & ifcausal ==1])/nrow(a[type == "gene"])
# a2 <- a[ ,sum(ifcausal), by = list(b, rn)]
# dim(a2)
# a0 <- dt[,sum(ifcausal), by = list(b, rn)]
# nrow(a2[a2$V1<=1,])/nrow(a2)

