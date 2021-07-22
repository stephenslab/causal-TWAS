library(data.table)
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(" 5 arguments:
       * prune ld score file
       * gwas file
       * eqtl summary statistics
       * TWAS results files
       * out file name", call.=FALSE)
}

mrjti <- "/home/simingz/causalTWAS/MR-JTI/mr/MR-JTI.r"


ldfile <- args[1]
gwas <- args[2]
eqtl <- args[3]
twas <- args[4]
outname <- args[5]

harmonize_gwas_eqtl <- function(gwas, eqtl){
  snpnames <- intersect(gwas$id, eqtl$id)
  if (length(snpnames) != 0){
    g.idx <- match(snpnames, gwas$id)
    e.idx <- match(snpnames, eqtl$id)
    qc <- ctwas:::allele.qc(gwas[g.idx, ]$alt, gwas[g.idx, ]$ref,
                    eqtl[e.idx, ]$alt, eqtl[e.idx, ]$ref)

    ifflip <- qc[["flip"]]
    ifremove <- !qc[["keep"]]

    flip.idx <- g.idx[ifflip]
    gwas[flip.idx, c("alt", "ref")] <- gwas[flip.idx, c("ref", "alt")]
    gwas[flip.idx, "Estimate"] <- - gwas[flip.idx, "Estimate"]

    remove.idx <- g.idx[ifremove]
    if (length(remove.idx) != 0) {
      gwas <- gwas[-remove.idx, ,drop = F]
    }
  }
  return(gwas)
}

# save.image("temp.Rd")
# GWAS
gwas <- fread(gwas, header = T)

#  GTEx: the eQTL effect allele is the ALT allele
eqtl <- fread(eqtl, header = T)
eqtl <- eqtl[, c("rs_id_dbSNP147_GRCh37p13", "gene_name", "ref", "alt", "slope", "slope_se", "pval_nominal")]

# harmonize gwas and eqtl
gwas <- harmonize_gwas_eqtl(gwas, eqtl)

# merge gwas, eqtl and ldinfo

ldfiles <- fread(ldfile, header = F)

ldinfolist <- list()
for (i in 1:22){
  prune <- fread(ldfiles[[i,1]], header = F)
  colnames(prune) <- "SNP"
  ldsc <- fread(ldfiles[[i,2]], header = T)
  ldinfolist[[i]] <- merge(prune, ldsc, by = "SNP")
}
ldinfo <- do.call(rbind, ldinfolist)

snpall <- merge(merge(gwas, eqtl, by.x = "id", by.y = "rs_id_dbSNP147_GRCh37p13"),
      ldinfo, by.x = "id", by.y = "SNP")

# get genes that pass TWAS FDR cut off
a <- fread(twas, header = T)
a$TWAS.P.FDR <- p.adjust(a$TWAS.P, method = "fdr")
genes <- a[a$TWAS.P.FDR <0.05, ID]
ngene <- length(genes)

# run MR-JTI
gnlist <- list()
gene.test <- NULL
for (gene in genes){
  gout <- snpall[snpall$gene_name == gene, c("id", "ldscore", "slope", "slope_se", "pval_nominal", "Estimate", "PVALUE")]
  colnames(gout) <- c("rsid", "ldscore", "eqtl_beta", "eqtl_se", "eqtl_p", "gwas_beta", "gwas_p")
  gnlist[[gene]] <- nrow(gout)
  if (nrow(gout) < 20) next
  gfile <- paste0(outname, "_", gene,"_input.txt")
  resfile <- paste0(outname, "_", gene, "_res.txt")
  write.table(gout, file= gfile, row.names=F, col.names=T, sep="\t", quote = F)
  gene.test<- c(gene.test, gene)
  # system(paste("Rscript", mrjti,  "--df_path", gfile, "--n_genes", ngene, "--out_path", resfile))
  # res <- read.table(resfile, header = T, stringsAsFactors = F, sep = ",")
  # res[1, 1] <- gene
  # reslist[[gene]] <- res[1, ]
  # file.remove(resfile)
  # file.remove(gfile)
}

outdf <- do.call(rbind, reslist)
write.table(outdf, file=outname, row.names=F, col.names=T, sep="\t", quote = F)

