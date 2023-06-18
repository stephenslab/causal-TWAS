library(data.table)
library(mrlocus)
library(ctwas)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 8) {
  stop(" 8 arguments:
       * LD_R directory
       * gwas file
       * eqtl summary statistics
       * weight .pos file
       * output directory
       * out file name
       * batch
       * prune clumps", call.=FALSE)
}

print(args)

ld_R_dir = as.character(args[1])
gwas <- as.character(args[2])
eqtl <- as.character(args[3])
weightf <- as.character(args[4])
outputdir <- as.character(args[5])
outname <- as.character(args[6])
batch <- as.numeric(args[7])
prune_clumps <- as.logical(args[8])

batch_size <- 250
ncores <- 2

####################

# ld_R_dir = "/project2/mstephens/causalTWAS/ukbiobank/ukb_LDR_s80.45/s80.45.2_LDR"
# gwas <- "/project2/mstephens/causalTWAS/simulations/simulation_ctwas_rss_20210416/ukb-s80.45-adi_simu4-1.snpgwas.txt.gz"
# eqtl <- "/project2/mstephens/causalTWAS/GTEx_v7_all/Adipose_Subcutaneous.allpairs_processed.txt.gz"
# weightf <- "/project2/mstephens/causalTWAS/fusion_weights/Adipose_Subcutaneous.pos"
# outputdir <- "/home/wcrouse/causalTWAS/simulations/simulation_ctwas_rss_20210416_compare/"
# outname <- "/project2/mstephens/causalTWAS/simulations/simulation_ctwas_rss_20210416_compare/ukb-s80.45-adi_simu4-1.Adipose_Subcutaneous.mrlocus.result.batch_1"
# batch <- 1
# prune_clumps <- TRUE

####################

ld_bedfs <- list.files(outputdir, pattern="*.bed")
ld_bedfs <- data.frame(chr=sapply(ld_bedfs, function(x){as.numeric(unlist(strsplit(unlist(strsplit(x, "_"))[2], "chr"))[2])}),
                       file=ld_bedfs)
rownames(ld_bedfs) <- NULL
ld_bedfs <- ld_bedfs[order(ld_bedfs$chr),]
ld_bedfs$file <- paste0(outputdir, tools::file_path_sans_ext(ld_bedfs$file))

temp_dir <- paste0(outputdir, "temp_mrlocus/")
outname_base <- rev(unlist(strsplit(outname, "/")))[1]

####################
#define function for harmonization
harmonize_sumstat_ld <- function(sumstats, ldref, effect_var="z"){
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
    sumstats[flip.idx, effect_var] <- -sumstats[flip.idx, effect_var]

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

#format GWAS data
gwas <- gwas[,!(colnames(gwas) %in% c("t.value", "PVALUE"))]

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
rm(eqtl_unique_id_gene)

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
gwas <- harmonize_sumstat_ld(gwas, ld_R_info, effect_var="Estimate")

#harmonize eQTL to LD reference
eqtl <- harmonize_sumstat_ld(eqtl, ld_R_info, effect_var="slope")

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
#load variant information for genotype files
ldref_files <- readLines("/project2/mstephens/causalTWAS/ukbiobank/ukb_pgen_s80.45/ukb-s80.45.2_pgenfs.txt")
ld_pvarfs <- sapply(ldref_files, prep_pvar, outputdir = temp_dir)

####################
#run MR-Locus for genes in this batch

genes <- weights$ENSEMBL_ID

#save.image("/home/wcrouse/scratch-midway2/temp_image.RData")
#save.image("/home/wcrouse/temp_image.RData")

if (file.exists(paste0(outname, "_temp"))){
  outdf <- as.data.frame(fread(paste0(outname, "_temp")))
  #outlist <- load(paste0(outname, "_samples.RData"))
} else {
  outdf <- data.frame(chromosome=as.numeric(),
                      gene_id=as.character(),
                      n_clumps=as.numeric(),
                      alpha=as.numeric(),
                      CI_10=as.numeric(),
                      CI_90=as.numeric())
  #outlist <- list()
}

j <- nrow(outdf) + 1

batch_start <- (batch-1)*batch_size+j
batch_end <- min((batch)*batch_size, length(genes))

#drop eqtl not in this batch to free memory
eqtl <- eqtl[eqtl$gene_id %in% genes[batch_start:batch_end],]

#subset to SNPs in all three datasets again - catch discordant variants between GWAS and eQTL after harmonization
snplist <- intersect(intersect(gwas$id, ld_R_info$id), eqtl$id)
gwas <- gwas[gwas$id %in% snplist,]
eqtl <- eqtl[eqtl$id %in% snplist,]
ld_R_info <- ld_R_info[ld_R_info$id %in% snplist,]
rm(snplist)

for (i in batch_start:batch_end){

  gene <- genes[i]

  chr <- as.numeric(gtf_df$seqnames[gtf_df$gene_id==gene])

  #prepare files for plink
  plink_exclude_file <- paste0(temp_dir, "ukb_chr", chr, "_s80.45.2.exclude")

  if (!file.exists(plink_exclude_file)){
    #store list of SNPs in .bed file using plink
    system_cmd <- paste0("plink --bfile ",
                         ld_bedfs$file[chr],
                         " --write-snplist --out ",
                         temp_dir, "ukb_chr", chr, "_s80.45.2")
    system(system_cmd)

    #store list of duplicate SNPs in .bed file to exclude
    system_cmd <- paste0("sort ",
                         temp_dir,
                         "ukb_chr", chr, "_s80.45.2.snplist|uniq -d > ",
                         temp_dir,
                         "ukb_chr", chr, "_s80.45.2.exclude")

    system(system_cmd)
  }

  #clumping using plink
  eqtl_current <- eqtl[eqtl$gene_id==gene,]
  eqtl_current <- dplyr::rename(eqtl_current, SNP="id", P="pval_nominal")

  eqtl_current_file <- paste0(temp_dir, outname_base, ".eqtl_sumstats.", gene, ".temp")
  write.table(eqtl_current, file=eqtl_current_file, sep="\t", col.names=T, row.names=F, quote=F)

  eqtl_clump_file <- paste0(temp_dir, outname_base, ".eqtl_clumps.", gene, ".temp")
  system_cmd <- paste0("plink --bfile ",
                       ld_bedfs$file[chr],
                       " --clump ",
                       eqtl_current_file,
                       " --clump-p1 0.001 --clump-p2 1 --clump-r2 0.1 --clump-kb 500 --out ",
                       eqtl_clump_file,
                       " --exclude ",
                       plink_exclude_file)

  system(system_cmd)

  eqtl_clump_file <- paste0(eqtl_clump_file, ".clumped")

  #load and parse clumps
  if (!file.exists(eqtl_clump_file)){
    outdf_current <- data.frame(chromosome=chr,
                                gene_id=gene,
                                n_clumps=0,
                                alpha=NA,
                                CI_10=NA,
                                CI_90=NA,
                                row.names = NULL)
  } else {
    clumps <- read.table(eqtl_clump_file, header=T)
    clumps$SP2 <- sapply(1:nrow(clumps), function(x){paste0(clumps$SNP[x], "(1),", clumps$SP2[x])})
    clumps <- lapply(clumps$SP2, function(y){unname(sapply(unlist(strsplit(y, ",")), function(x){unlist(strsplit(x, "[(]"))[1]}))})
    clumps <- lapply(clumps, function(x){x[x!="NONE"]})

    # #load LD matrix for SNPs in all clumps
    # snplist_all <- unlist(clumps)
    #
    # ld_R_idx <- unique(ld_R_info$index[match(snplist_all, ld_R_info$id)])
    #
    # R_snp <- lapply(ld_R_files[ld_R_idx], readRDS)
    #
    # R_snp_info <- lapply(ld_R_info_files[ld_R_idx], fread)
    # R_snp_info <- lapply(R_snp_info, as.data.frame)
    #
    # R_snp_index <- lapply(R_snp_info, function(x){x$id %in% snplist_all})
    # rm(snplist_all)
    #
    # for (j in 1:length(R_snp)){
    #   R_snp[[j]] <- R_snp[[j]][R_snp_index[[j]],R_snp_index[[j]]]
    #   R_snp_info[[j]] <- R_snp_info[[j]][R_snp_index[[j]],,drop=F]
    # }
    #
    # R_snp_info <- as.data.frame(do.call(rbind, R_snp_info))
    # R_snp <- as.matrix(Matrix::bdiag(R_snp))
    #
    # colnames(R_snp) <- R_snp_info$id
    # rownames(R_snp) <- R_snp_info$id

    #load genotypes and compute LD matrix for SNPs in all clumps
    snplist_all <- unlist(clumps)

    ld_pvarf <- ld_pvarfs[chr]
    ldref_file <- ldref_files[chr]

    R_snp_info <- read_pvar(ld_pvarf)
    ld_pgen <- prep_pgen(pgenf = ldref_file, ld_pvarf)

    sidx <- match(snplist_all, R_snp_info$id)
    X.g <- read_pgen(ld_pgen, variantidx = sidx)

    R_snp <- Rfast::cora(X.g)
    colnames(R_snp) <- snplist_all
    rownames(R_snp) <- snplist_all

    rm(snplist_all, X.g)

    #prepare data object
    data <- list(sum_stats=list(),
                 ld_mat=list())

    for (j in 1:length(clumps)){
      snplist <- clumps[[j]]

      R_snp_current <- R_snp[snplist,snplist,drop=F]
      dimnames(R_snp_current) <- NULL

      R_snp_info_current <- R_snp_info[match(snplist, R_snp_info$id),c("id", "ref", "alt")]

      sumstats_current <- data.frame(id=snplist,
                                     ref=R_snp_info_current$ref,
                                     eff=R_snp_info_current$alt,
                                     beta_hat_eqtl=eqtl_current$slope[match(snplist, eqtl_current$SNP)],
                                     beta_hat_gwas=gwas$Estimate[match(snplist, gwas$id)],
                                     se_eqtl=eqtl_current$slope_se[match(snplist, eqtl_current$SNP)],
                                     se_gwas=gwas$Std.Error[match(snplist, gwas$id)])
      sumstats_current$abs_z <- abs(sumstats_current$beta_hat_eqtl/sumstats_current$se_eqtl)

      data$sum_stat[[j]] <- sumstats_current
      data$ld_mat[[j]] <- R_snp_current
    }

    #collapse high correlation SNP keeping most extreme z score (most significant)
    data <- collapseHighCorSNPs(data$sum_stat,
                                data$ld_mat,
                                score="abs_z",
                                plot=F)

    data <- flipAllelesAndGather(data$sum_stat,
                                 data$ld_mat,
                                 snp_id="id",
                                 ref="ref",
                                 a="eqtl",
                                 sep="_",
                                 b="gwas",
                                 eff="eff",
                                 beta="beta_hat",
                                 se="se",
                                 a2_plink="ref_eqtl",
                                 alleles_same=T,
                                 plot=F)

    #colocalization
    coloc_fit <- list()
    nclust <- length(data$beta_hat_a)

    options(mc.cores=ncores)
    for (j in 1:nclust) {
      if (length(data$beta_hat_a[[j]])>1){
        coloc_fit[[j]] <- with(data,
                               fitBetaColoc(
                                 beta_hat_a = beta_hat_a[[j]],
                                 beta_hat_b = beta_hat_b[[j]],
                                 se_a = se_a[[j]],
                                 se_b = se_b[[j]],
                                 Sigma_a = Sigma[[j]],
                                 Sigma_b = Sigma[[j]]
                               ))
      } else {
        coloc_fit[[j]] <- list(beta_hat_a=data$beta_hat_a[[j]],
                               beta_hat_b=data$beta_hat_b[[j]])
      }
    }

    res <- list(beta_hat_a = lapply(coloc_fit, `[[`, "beta_hat_a"),
                beta_hat_b = lapply(coloc_fit, `[[`, "beta_hat_b"),
                sd_a = data$se_a,
                sd_b = data$se_b,
                alleles = data$alleles)

    res <- extractForSlope(res, plot=F)

    #sort clusters by significance
    res$abs_z <- abs(res$beta_hat_a/res$sd_a)
    res_sig_order <- order(-res$abs_z)

    res$beta_hat_a <- res$beta_hat_a[res_sig_order]
    res$beta_hat_b <- res$beta_hat_b[res_sig_order]
    res$sd_a <- res$sd_a[res_sig_order]
    res$sd_b <- res$sd_b[res_sig_order]
    res$alleles <- res$alleles[res_sig_order,,drop=F]
    res$abs_z <- NULL

    #prune clumps if (absolute) correlation greater than 0.05
    if (length(res$alleles$id)>1 & prune_clumps){
      trim_index <- trimClusters(r2=abs(R_snp[res$alleles$id,res$alleles$id,drop=F]),
                                 r2_threshold=0.05)

      if (length(trim_index)>0){
        res$beta_hat_a <- res$beta_hat_a[-trim_index]
        res$beta_hat_b <- res$beta_hat_b[-trim_index]
        res$sd_a <- res$sd_a[-trim_index]
        res$sd_b <- res$sd_b[-trim_index]
        res$alleles <- res$alleles[-trim_index,,drop=F]
      }
    }

    if (length(res$beta_hat_a)==1){
      #parametric sample from MRLocus code
      m <- 1e+05
      post_samples <- list(alpha=rnorm(m, res$beta_hat_b, res$sd_b)/rnorm(m, res$beta_hat_a, res$sd_a))
    } else {
      #fit MRLocus slope model
      res <- fitSlope(res, iter=10000)
      post_samples <- rstan::extract(res$stanfit)

      #display results
      #print(res$stanfit, pars=c("alpha","sigma"), probs=c(.1,.9), digits=3)
    }

    outdf_current <- data.frame(chromosome=chr,
                                gene_id=gene,
                                n_clumps=length(res$beta_hat_a),
                                alpha=mean(post_samples$alpha),
                                CI_10=quantile(post_samples$alpha, 0.1),
                                CI_90=quantile(post_samples$alpha, 0.9),
                                row.names = NULL)

    attributes(post_samples$alpha) <- NULL
    #outlist[[gene]] <- post_samples$alpha
  }

  outdf <- rbind(outdf,
                 outdf_current)

  cleanup_files <- list.files(temp_dir, outname_base)

  if (length(cleanup_files)>0){
    file.remove(paste0(temp_dir, list.files(temp_dir, outname_base)))
  }

  write.table(outdf, file=paste0(outname, "_temp"), row.names=F, col.names=T, sep="\t", quote = F)

  #save(outlist, file=paste0(outname, "_samples.RData"))
  #save(outlist, file=paste0(outname, "_samples_", i, ".RData"))
}

write.table(outdf, file=paste0(outname), row.names=F, col.names=T, sep="\t", quote = F)

#file.remove(paste0(outname, "_temp"))
