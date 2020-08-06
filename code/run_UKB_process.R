codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir,"input_reformat.R"))

#----------------------------chr17-22 combined-------------------
weight <- "/home/simingz/causalTWAS/fusion_weights/Adipose_Subcutaneous"
method <- "lasso"
wgtcut <- 1e-8 # abs(wight) > wgtcut
prune <- c(T, rep(F,9)) # T will be kept

wgtfs <- read.table(paste0(weight,".pos"), header = T)
wgtdir <-dirname(weight)
wgtfs <- paste0(wgtdir, "/", wgtfs$WGT)
snpall <- NULL
snpannoall <- NULL
for (wgtf in wgtfs){
  load(wgtf)
  snp <- rownames(wgt.matrix[abs(wgt.matrix[, method]) > wgtcut,])
  snpanno <- snps[snps[,2] %in% snp, ]
  snpanno <- cbind(snpanno, wgt.matrix[snpanno[,2],])
  snpall <- c(snpall, snp)
  snpannoall <- rbind(snpannoall, snpanno)
}

write.table(snpannoall, file= paste(weight, method, "eqtl.txt", sep= ".") , row.names=F, col.names=T, sep="\t", quote = F)

pfiles <- paste0("/home/simingz/causalTWAS/ukbiobank/ukb_chr", 17:22, "_s20000.pgen")

mtxall <- NULL
mtxall.unscaled <- NULL
varall <- NULL
for (pfile in pfiles){
  pfile <- drop_ext(pfile)
  pvar <- paste0(pfile, ".pvar")
  var <- read.table(pvar, header = F)
  eqtl <- var[,3] %in% snpall
  prunetag <- rep_len(prune, nrow(var))
  select <- (1:nrow(var))["|"(eqtl, prunetag)]
  varall <- rbind(varall, var[select,])
  mtx <- pgen2mtx(pfile, nb_parts = 1, select=select, progress = T)
  mtxall.unscaled <- cbind(mtxall.unscaled, mtx)
  # mtx.scaled <- scaleRcpp(mtx)
  # mtxall <- cbind(mtxall, mtx.scaled)
}

# dat <- list("G"       = mtxall,
#             "chr"     = as.matrix(varall[,1]),
#             "pos"     = as.matrix(varall[,2]),
#             "snp"     = as.matrix(varall[,3]),
#             "counted" = as.matrix(varall[,4]),
#             "alt"     = as.matrix(varall[,5]))
#
#
# m <- as_FBM(mtxall, backingfile = "/home/simingz/causalTWAS/ukbiobank/ukb_chr17to22_s20000")$save()
# dat$G <- m
# save(dat, file = "/home/simingz/causalTWAS/ukbiobank/ukb_chr17to22_s20000.FBM.Rd")

dat <- list("G"       = mtxall.unscaled,
            "chr"     = as.matrix(varall[,1]),
            "pos"     = as.matrix(varall[,2]),
            "snp"     = as.matrix(varall[,3]),
            "counted" = as.matrix(varall[,4]),
            "alt"     = as.matrix(varall[,5]))

m <- as_FBM(mtxall.unscaled, backingfile = "/home/simingz/causalTWAS/ukbiobank/ukb_chr17to22_s20000.unscaled")$save()
dat$G <- m
save(dat, file = "/home/simingz/causalTWAS/ukbiobank/ukb_chr17to22_s20000.unscaled.FBM.Rd")

#----------------------------chr22 only-------------------
pfile <- "/home/simingz/causalTWAS/ukbiobank/ukb_chr22_s20000"
pgen2fbm(pfile, select = NULL, scale = T, type = "double")
