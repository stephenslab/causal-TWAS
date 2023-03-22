library(PMR)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 8) {
  stop(" 7 arguments:
       * LD_R directory
       * gwas file
       * eqtl summary statistics
       * weight dir
       * gwas sample size
       * eqtl sample size
       * out file name", call.=FALSE)
}

ld_R_dir = as.character(args[1])
gwas <- as.character(args[2])
eqtl <- as.character(args[3])
weightf <- as.character(args[4])
gwasn = as.numeric(args[5])
eqtln = as.numeric(args[6])
outname <- as.character(args[7])
batch <- as.numeric(args[8])

batch_size <- 100

#window_size_kb <- NULL
window_size_kb <- 100

#lambda <- 0
lambda <- 0.25




# simutag <- "4-1"
#
# ld_R_dir = "/project2/mstephens/causalTWAS/ukbiobank/ukb_LDR_s80.45/s80.45.2_LDR"
# gwas <- paste0("/project2/mstephens/causalTWAS/simulations/simulation_ctwas_rss_20210416/ukb-s80.45-adi_simu", simutag, ".snpgwas.txt.gz")
# eqtl <- "/project2/mstephens/causalTWAS/GTEx_v7_all/Adipose_Subcutaneous.allpairs_processed.txt.gz"
# weightf <- "/project2/mstephens/causalTWAS/fusion_weights/Adipose_Subcutaneous.pos"
# gwasn = 45000
# eqtln = 380
# outname <- paste0("/project2/mstephens/causalTWAS/simulations/simulation_ctwas_rss_20210416_compare/ukb-s80.45-adi_simu", simutag, ".Adipose_Subcutaneous.pmr.result.batch_2")
# batch <- 2

####################
#define function for harmonization
harmonize_sumstat_ld <- function(sumstats, ldref){
  snpnames <- intersect(sumstats$id, ldref$id)

  if (length(snpnames) != 0){
    ss.idx <- which(sumstats$id %in% snpnames)
    ld.idx <- match(sumstats$id[ss.idx], ldref$id)
    qc <- ctwas:::allele.qc(sumstats$alt[ss.idx], sumstats$ref[ss.idx],
                            ldref$alt[ld.idx], ldref$ref[ld.idx])
    ifflip <- qc[["flip"]]
    ifremove <- !qc[["keep"]]

    flip.idx <- ss.idx[ifflip]
    sumstats[flip.idx, c("alt", "ref")] <- sumstats[flip.idx, c("ref", "alt")]
    sumstats[flip.idx, "z"] <- -sumstats[flip.idx, "z"]

    remove.idx <- ss.idx[ifremove]
    if (length(remove.idx) != 0) {
      sumstats <- sumstats[-remove.idx,,drop = F]
    }
  }
  return(sumstats)
}

####################
#load GWAS
gwas <- as.data.frame(fread(gwas, header = T))

#rename z score
gwas <- dplyr::rename(gwas, z="t.value")
gwas <- gwas[,!(colnames(gwas) %in% c("Estimate", "Std.Error", "PVALUE"))]

#drop variants that are non uniquely identified by ID
gwas <- gwas[!(gwas$id %in% gwas$id[duplicated(gwas$id)]),]

####################
#load weight names
weights <- as.data.frame(fread(weightf, header = T))
weights$ENSEMBL_ID <- sapply(weights$WGT, function(x){unlist(strsplit(unlist(strsplit(x,"/"))[2], "[.]"))[2]})

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

#drop genes without FUSION weights
eqtl <- eqtl[eqtl$gene_id %in% weights$ENSEMBL_ID,]

#drop entries not uniquely identified by gene_id and variant id (variant not biallelic)
eqtl_unique_id_gene <- paste0(eqtl$id, eqtl$gene_id)
eqtl <- eqtl[!(eqtl_unique_id_gene %in% eqtl_unique_id_gene[duplicated(eqtl_unique_id_gene)]),]

#compute z score
eqtl$z <- eqtl$slope/eqtl$slope_se
eqtl <- eqtl[,!(colnames(eqtl) %in% c("slope", "slope_se", "pval_nominal"))]

####################
#LD files
ld_R_files <- paste0(ld_R_dir, "/", list.files(ld_R_dir))
ld_R_files <- ld_R_files[grep(".RDS", ld_R_files)]

ld_R_info_files <- paste0(ld_R_dir, "/", list.files(ld_R_dir))
ld_R_info_files <- ld_R_info_files[grep(".Rvar", ld_R_info_files)]

ld_R_info <- lapply(1:length(ld_R_info_files), function(x){y <- fread(ld_R_info_files[x]); y$index <- x; as.data.frame(y)})
ld_R_info <- do.call(rbind, ld_R_info)

####################
#subset to SNPs in all three datasets
snplist <- intersect(intersect(gwas$id, ld_R_info$id), eqtl$id)

gwas <- gwas[gwas$id %in% snplist,]
eqtl <- eqtl[eqtl$id %in% snplist,]
ld_R_info <- ld_R_info[ld_R_info$id %in% snplist,]

rm(snplist)

####################
#harmonize GWAS to LD reference
gwas <- harmonize_sumstat_ld(gwas, ld_R_info)

#harmonize eQTL to LD reference
eqtl <- harmonize_sumstat_ld(eqtl, ld_R_info)

####################
#load gene location information

# gtf <- rtracklayer::import("/project2/mstephens/causalTWAS/GTEx_v7_all/gencode.v19.genes.v7.patched_contigs.gtf")
# gtf_df <- as.data.frame(gtf)
# save(gtf_df, file="/project2/mstephens/causalTWAS/GTEx_v7_all/gencode.v19.genes.v7.patched_contigs.RData")
load("/project2/mstephens/causalTWAS/GTEx_v7_all/gencode.v19.genes.v7.patched_contigs.RData")

#keep only gene-level information
gtf_df <- gtf_df[gtf_df$type=="gene",]

#trim ensembl ID version number
gtf_df$gene_id <- sapply(gtf_df$gene_id, function(x){unlist(strsplit(x, "[.]"))[1]})

#drop genes without FUSION weights
gtf_df <- gtf_df[gtf_df$gene_id %in% weights$ENSEMBL_ID,]

####################
#gene list from eqtl
genes <- unique(eqtl$gene_id)

if (file.exists(paste0(outname, "_temp"))){
  outdf <- as.data.frame(fread(paste0(outname, "_temp")))
} else {
  outdf <- data.frame(gene_id=as.character(),
                      causal_effect=as.numeric(),
                      causal_pvalue=as.numeric(),
                      pleiotropy_effect=as.numeric(),
                      pleiotropy_pvalue=as.numeric(),
                      sigma_cisSNP=as.numeric(),
                      sigma_error_1=as.numeric(),
                      sigma_error_2=as.numeric(),
                      nsnps=as.numeric(),
                      pmr_message=as.character())
}

# outdf <- data.frame(gene_id=as.character(),
#                     causal_effect=as.numeric(),
#                     causal_pvalue=as.numeric(),
#                     pleiotropy_effect=as.numeric(),
#                     pleiotropy_pvalue=as.numeric(),
#                     sigma_cisSNP=as.numeric(),
#                     sigma_error_1=as.numeric(),
#                     sigma_error_2=as.numeric(),
#                     nsnps=as.numeric(),
#                     pmr_message=as.character())

j <- nrow(outdf) + 1

batch_start <- (batch-1)*batch_size+j
batch_end <- min((batch)*batch_size, length(genes))

for (i in batch_start:batch_end){

  print(i)

  gene <- genes[i]

  #subset eQTL data to current gene
  eqtl_current <- eqtl[eqtl$gene_id==gene,]

  #all SNPs for current gene
  snplist <- eqtl_current$id

  if (!is.null(window_size_kb)){
    #positions +/- window_size of gene start and end
    gtf_df_idx <- which(gtf_df$gene_id==gene)
    gene_start <- gtf_df$start[gtf_df_idx] - window_size_kb*1000
    gene_end <- gtf_df$end[gtf_df_idx] + window_size_kb*1000

    #subset snplist to SNPs within window
    ld_R_info_current <- ld_R_info[match(snplist, ld_R_info$id),]
    ld_R_info_current <- ld_R_info_current[ld_R_info_current$pos >= gene_start  & ld_R_info_current$pos <= gene_end,]
    snplist <- ld_R_info_current$id
  }

  if (length(snplist)>0){
    #construct PMR inputs
    Zscore_1 <- eqtl_current$z[match(snplist, eqtl_current$id)]
    Zscore_2 <- gwas$z[match(snplist, gwas$id)]

    ld_R_idx <- unique(ld_R_info$index[match(snplist, ld_R_info$id)])

    R_snp <- lapply(ld_R_files[ld_R_idx], readRDS)
    R_snp <- as.matrix(Matrix::bdiag(R_snp))

    R_snp_info <- lapply(ld_R_info_files[ld_R_idx], fread)
    R_snp_info <- as.data.frame(do.call(rbind, R_snp_info))

    colnames(R_snp) <- R_snp_info$id
    rownames(R_snp) <- R_snp_info$id

    R_snp <- R_snp[snplist,snplist,drop=F]

    error_flag <- FALSE
    pmr_message <- ""

    tryCatch(PMR_results <- PMR_summary_Egger(Zscore_1=Zscore_1,
                                              Zscore_2=Zscore_2,
                                              Sigma1sin=R_snp,
                                              Sigma2sin=R_snp,
                                              samplen1=eqtln,
                                              samplen2=gwasn,
                                              lambda=lambda),
             error = function(cond){error_flag <<- TRUE
                                    pmr_message <<- cond$message})
  } else {
    error_flag <- T
    pmr_message <- ""
  }

  if(error_flag){
    PMR_results <- list(causal_effect=NA,
                        causal_pvalue=NA,
                        pleiotropy_effect=NA,
                        pleiotropy_pvalue=NA,
                        sigma_cisSNP=NA,
                        sigma_error_1=NA,
                        sigma_error_2=NA)
  }

  outdf <- rbind(outdf,
                 data.frame(gene_id=gene,
                            causal_effect=PMR_results$causal_effect,
                            causal_pvalue=PMR_results$causal_pvalue,
                            pleiotropy_effect=PMR_results$pleiotropy_effect,
                            pleiotropy_pvalue=PMR_results$pleiotropy_pvalue,
                            sigma_cisSNP=PMR_results$sigma_cisSNP,
                            sigma_error_1=PMR_results$sigma_error_1,
                            sigma_error_2=PMR_results$sigma_error_2,
                            nsnps=length(snplist),
                            pmr_message=pmr_message))

  write.table(outdf, file=paste0(outname, "_temp"), row.names=F, col.names=T, sep="\t", quote = F)
}

write.table(outdf, file=outname, row.names=F, col.names=T, sep="\t", quote = F)

#file.remove(paste0(outname, "_temp"))
