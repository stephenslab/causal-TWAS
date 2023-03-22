library(ctwas)

args = commandArgs(trailingOnly=TRUE)

print(args)

if (length(args) <9) {
  stop(" 9 arguments:
       * zscore file name
       * ld R directory
       * weight_1
       * weight_2
       * weight_null
       * config file name (.R)
       * out expr z file name
       * out file name
       * outputdir", call.=FALSE)
}

sessionInfo()

####################

outputdir <- args[9]

dir.create(outputdir, showWarnings=F, recursive=T)

z_snp <- data.table::fread(args[1], header = T)
data.table::setnames(z_snp, old = c("alt", "ref", "t.value"), new = c("A1", "A2", "z"))
z_snp <- z_snp[, c("id", "A1", "A2", "z")]

ld_R_dir <- args[2]

weight_1 <- args[3]
weight_2 <- args[4]
weight_null <- args[5]
weight <- c(weight_1, weight_2, weight_null)

outname.e <- args[7]
outname <- args[8]

source(args[6]) # config

#set other parameters, previously specified in config
thin <- 0.1
ncore <- 10
ncore.rerun <- 1
prob_single <- 0.8
max_snp_region <- 20000
ld_regions <- "EUR"
ld_regions_version <- "b38"

# get gene z score
if (file.exists(paste0(outputdir, "/", outname.e, "_z_gene.Rd"))){
  ld_exprfs <- paste0(outputdir, "/", outname.e, "_chr", 1:22, ".expr.gz")
  load(file = paste0(outputdir, "/", outname.e, "_z_gene.Rd"))
  load(file = paste0(outputdir, "/", outname.e, "_z_snp.Rd"))
} else {
  res <- impute_expr_z(z_snp, weight = weight, ld_R_dir = ld_R_dir,
                       method = NULL, outputdir = outputdir, outname = outname.e,
                       harmonize_z = T, harmonize_wgt = T,
                       strand_ambig_action_z = "drop", recover_strand_ambig_wgt = F)
  z_gene <- res$z_gene
  ld_exprfs <- res$ld_exprfs
  z_snp <- res$z_snp

  save(z_gene, file = paste0(outputdir, "/", outname.e, "_z_gene.Rd"))
  save(z_snp, file = paste0(outputdir, "/", outname.e, "_z_snp.Rd"))
}

warnings()

#set the type for each gene (i.e. tissue)
z_gene$type <- sapply(z_gene$id, function(x){paste(unlist(strsplit(unlist(strsplit(x, "[|]"))[2],"_"))[-1], collapse="_") })

#run ctwas_rss
ctwas_rss(z_gene, z_snp, ld_exprfs, ld_pgenfs = NULL, ld_R_dir = ld_R_dir, ld_regions = ld_regions, ld_regions_version = ld_regions_version, thin = thin, max_snp_region = max_snp_region, outputdir = outputdir, outname = outname, ncore = ncore, ncore.rerun = ncore.rerun, prob_single = prob_single,
          group_prior_var_structure = group_prior_var_structure)

warnings()

sessionInfo()
