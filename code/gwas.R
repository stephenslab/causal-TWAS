library(foreach)
library(doParallel)


#' @description this is a function to run GWAS analysis for small genotype matrix (using glm),for faster results use plink (start from plink files). Credit to http://www.stat-gen.org/tut/tut_analysis.html, requires tabix and bgzip if compress is true.
#'
#' @param geno genotype matrix, columns: snps, rows: individuals.
#' @param pheno phenotype, matrix, ncol = 1
#' @param snpname a vector of N, correpond to geno type matrix, will be written to output.
#' @param anno a data frame/matrix, anno column with N rows.
#' @param compress T/F, if T will be compressed and indexed by tabix. (TODO)
#'
GWAA <- function(geno, pheno, snpname = NULL, anno = NULL, outname, family = gaussian, ncore = 1, nSplits = 10, compress = F){

  cl <- makeCluster(ncore)
  show(cl)
  registerDoParallel(cl)

  colnames(pheno) <- "pheno"

  if (is.null(snpname)){
    snpname <- 1: ncol(geno)
  }

  nSNPs <- ncol(geno)
  genosplit <- ceiling(nSNPs/nSplits) # number of SNPs in each subset

  snp.start <- seq(1, nSNPs, genosplit) # index of first SNP in group
  snp.stop <- pmin(snp.start+genosplit-1, nSNPs) # index of last SNP in group

  columns<-c("#CHROM",	"BEGIN",	"END", "MARKER_ID", "Estimate", "Std.Error", "t-value", "PVALUE") # header for LOCUSZOOM.
  write.table(t(columns), outname, row.names=FALSE, col.names=FALSE, quote=FALSE, sep = "\t")

  foreach (part=1:nSplits) %do% {

    geno.i <- geno[ ,snp.start[part]:snp.stop[part]]
    snpname.i <- snpname[snp.start[part]:snp.stop[part]]
    if (!is.null(anno)) {
      anno.i <- anno[snp.start[part]:snp.stop[part], ]
    } else {
      anno.i <- NULL
    }

    res <- foreach(snp.idx=1: length(snpname.i), .combine='rbind') %dopar% {
      snp <- geno.i[, snp.idx, drop = F]
      a <- summary(glm(pheno ~ snp,
                    family=family,
                    data=data.frame("pheno" = pheno, "snp" = snp)))
      a$coefficients['snp',]
    }

    # write results so far to a file
    write.table(cbind(anno.i, snpname.i, res), outname, append=TRUE,
                quote=FALSE, col.names=FALSE, row.names=FALSE,
                sep = "\t")

    cat(sprintf("GWAS SNPs %s%% finished\n", 100 * part/nSplits))
  }

  stopCluster(cl)

  # compress
  if (isTRUE(compress)){
    system(paste0("bgzip ", outname))
    system(paste0("tabix -p bed ", outname, ".gz"))
  }

  return(print("Done."))
}

