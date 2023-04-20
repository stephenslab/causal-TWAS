library(data.table)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 6) {
  stop(" 4 arguments:
       * gwas file
       * eqtl summary statistics
       * TWAS results files
       * out file name", call.=FALSE)
}

# gwas <- as.character(args[1])
# eqtl <- as.character(args[2])
# twas <- as.character(args[3])
# ld_bedfs <- as.character(args[4])
# outname <- as.character(args[5])

gwas <- "/home/wcrouse/causalTWAS/simulations/simulation_ctwas_rss_20210416/ukb-s80.45-adi_simu1-1.snpgwas.txt.gz"
eqtl <- "/project2/mstephens/causalTWAS/GTEx_v7_all/Adipose_Subcutaneous.allpairs_processed.txt.gz"
weightf <- "/project2/mstephens/causalTWAS/fusion_weights/Adipose_Subcutaneous.pos"
ld_bedfs <- "/home/wcrouse/causalTWAS/simulations/simulation_ctwas_rss_20210416_compare/ukb_chr22_s80.45.2.FUSION22.bed"
outname <- "test_mrlocus"

####################

harmonize_eqtl_gwas <- function(eqtl_sumstats, gwas_sumstats){
  snpnames <- intersect(eqtl_sumstats$id, gwas_sumstats$id)

  if (length(snpnames) != 0){
    ss.idx <- which(eqtl_sumstats$id %in% snpnames)
    ld.idx <- match(eqtl_sumstats$id[ss.idx], gwas_sumstats$id)
    qc <- ctwas:::allele.qc(eqtl_sumstats$alt[ss.idx], eqtl_sumstats$ref[ss.idx],
                            gwas_sumstats$alt[ld.idx], gwas_sumstats$ref[ld.idx])
    ifflip <- qc[["flip"]]
    ifremove <- !qc[["keep"]]

    flip.idx <- ss.idx[ifflip]
    eqtl_sumstats[flip.idx, c("alt", "ref")] <- eqtl_sumstats[flip.idx, c("ref", "alt")]
    eqtl_sumstats[flip.idx, "slope"] <- -eqtl_sumstats[flip.idx, "slope"]

    remove.idx <- ss.idx[ifremove]
    if (length(remove.idx) != 0) {
      eqtl_sumstats <- eqtl_sumstats[-remove.idx,,drop = F]
    }
  }
  return(eqtl_sumstats)
}

####################
#load weight names
weights <- as.data.frame(fread(weightf, header = T))
genes <- sapply(weights$WGT, function(x){unlist(strsplit(unlist(strsplit(x,"/"))[2], "[.]"))[2]})
names(genes) <- NULL
rm(weights)

####################
# GWAS
gwas <- as.data.frame(fread(gwas, header = T))

#drop variants that are non uniquely identified by ID
gwas <- gwas[!(gwas$id %in% gwas$id[duplicated(gwas$id)]),]

####################
#load eQTL; in GTEx, the eQTL effect allele is the ALT allele
eqtl <- as.data.frame(fread(eqtl, header = T))

#trim version number from ensembl IDs
eqtl_gene_id_crosswalk <- unique(eqtl$gene_id)
eqtl_gene_id_crosswalk <- data.frame(original=eqtl_gene_id_crosswalk,
                                     trimmed=sapply(eqtl_gene_id_crosswalk, function(x){unlist(strsplit(x, "[.]"))[1]}))
eqtl$gene_id <- eqtl_gene_id_crosswalk[eqtl$gene_id, "trimmed"]
eqtl <- eqtl[, c("rs_id_dbSNP147_GRCh37p13", "gene_id", "gene_name", "ref", "alt", "slope", "slope_se", "pval_nominal")]
eqtl <- dplyr::rename(eqtl, id="rs_id_dbSNP147_GRCh37p13")

#keep genes analyzed by TWAS
eqtl <- eqtl[eqtl$gene_id %in% genes,,drop=F]

#subset to SNPs in GWAS
eqtl <- eqtl[eqtl$id %in% gwas$id,,drop=F]

#drop entries not uniquely identified by gene_id and variant id (variant not biallelic)
eqtl_unique_id_gene <- paste0(eqtl$id, eqtl$gene_id)
eqtl <- eqtl[!(eqtl_unique_id_gene %in% eqtl_unique_id_gene[duplicated(eqtl_unique_id_gene)]),]

#harmonize gwas and eqtl
eqtl <- harmonize_eqtl_gwas(eqtl, gwas)













