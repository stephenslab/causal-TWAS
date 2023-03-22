
library(tools)
library(ctwas)

args = commandArgs(trailingOnly=TRUE)

args <- c("/project2/mstephens/causalTWAS/ukbiobank/ukb_pgen_s80.45/ukb-s80.45_pgenfs.txt",
          "/project2/mstephens/causalTWAS/simulations/simulation_ctwas_rss_20220127_wc/ukb-s80.45-liv_wb_simu1-1-pheno.Rd",
          "ukb-s80.45-liv_wb_simu1-1",
          "/project2/mstephens/causalTWAS/simulations/simulation_ctwas_rss_20220127_wc/",
          "27",
          "1000")


if (length(args) < 4 ) {
  stop(" 6 arguments (last two optional):
       * geno file (txt file)
       * phenotype file (Rd file)
       * out file name
       * outputdir
       * ncore
       * nsplits", call.= FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "gwas_wc.R"))

outputdir <- args[4]

pgenfs <- read.table(args[1], header = F, stringsAsFactors = F)[,1]
pvarfs <- sapply(pgenfs, prep_pvar, outputdir = outputdir)

phenofile <- args[2]
load(phenofile)
pheno <- phenores$Y

outnames <- paste0(args[3], "-B", 1:length(pgenfs), ".snpgwas.txt")

ncore <- as.numeric(args[5])
if (is.na(args[5])) ncore <- 3

nsplits <- as.numeric(args[6])
if (is.na(args[6])) nsplits <- 100

for (b in 1:length(pgenfs)){

  outname <- outnames[b]

  # load genotype data

  snpinfo <- read_pvar(pvarfs[b])
  anno <- snpinfo[, - "id"] # chrom pos alt ref

  GWAA(pgenfs[b], mode = "snp", pheno, snpname = snpinfo$id, anno = anno,
       outname, outputdir, family = gaussian,
       ncore = ncore, nSplits = nsplits, compress = F)
}

# combine all results
outdflist <- list()

# for (b in 1:length(pgenfs)){
#   outdflist[[b]] <- read.table(paste0(outnames[b], ".gz"),
#                                header = T, stringsAsFactors = F,
#                                comment.char = "")
# }
for (b in 1:length(pgenfs)){
  outdflist[[b]] <- read.table(paste0(outputdir, outnames[b]),
                               header = F, stringsAsFactors = F,
                               comment.char = "")
}
outdf <- do.call(rbind, outdflist)

outname <-  paste0(outputdir, "/", args[3], ".snpgwas.txt")
write.table(outdf, file = outname , row.names=F, col.names=F, sep="\t", quote = F)
system(paste0("gzip ", outname))

# clean
# for (b in 1:length(pgenfs)){
#   system(paste0("rm ", outputdir, outnames[b]))
# }




