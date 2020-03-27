library(foreach)
library(doParallel)
library(glmnet)

Ncore <- 2

run_ctwas <- function(phenof, genof, exprtrf, gafile, geno, chr, pos, labels, outname) {
  # phenof: phenotype file, in BIMBAM phenotype format
  # genof:  genotype file in BIMBAM mean genotype format
  # exprtrf: .Rdata file for expr trained model has variable trlist. trlist[[genename]]: $fit, $labels: SNP labels for SNPs used for this gene.
  # geno: genotype matrix, column SNPs row samples
  # pos: vector, positions for each SNP
  # labels: vector, SNP labels for each SNP

  cat(" * Computing kinship matrix for snps.\n")
  genodir <-  normalizePath(dirname(genof))
  genokdir <- paste0(genodir, "/", basename(genof), ".kdir/")
  if (!dir.exists(genokdir)){
    dir.create(genokdir)
    system(paste("cd", genokdir, "; gemma -g", genof, "-p", phenof, "-gk 1 -o", paste0(genof,"-geno"))) # currently this step requires 80G memory for 3 hours
  }

  exprf <-  normalizePath(paste0(outname, "_imputed_expr.txt"))
  if (!file.exists(exprf)){
    cat(" * loading trained expression model.\n")
    load(exprtrf)

    cat(" * imputing expression.\n")
    gnames <- names(trlist)
    N <- dim(geno)[1] # number of samples
    M <- dim(geno)[2] # number of SNPs
    ga <- read.table(gafile, header = F, stringsAsFactors = F)

    impute_expr_g <- function(i){
      # i is the gene index
      gname <- gnames[i]
      if (is.null(trlist[[gname]])){
        return(NULL)
      }

      gai <- ga[ga$V6 == gname,]
      idx <- c(1:M)[paste0("chr", chr) == gai[[2]] & pos > gai[[4]] & pos < gai[[5]]]
      ggeno <- as.matrix(geno[,idx])
      glabels <- labels[idx]
      colnames(ggeno) <- glabels
      sharedlabels <- as.character(intersect(glabels, trlist[[gname]]$SNPlabel))
      if (length(sharedlabels) == 0) {
        return(NULL)
      }

      ggeno.test <- matrix(NA, nrow = N , ncol = length(trlist[[gname]]$SNPlabel))
      colnames(ggeno.test) <- as.character(trlist[[gname]]$SNPlabel)
      ggeno.test[ ,sharedlabels] <- ggeno[ ,sharedlabels]
      if (length(trlist[[gname]]$SNPlabel) == 1){
        X.tilde <- predict(trlist[[gname]]$fit, ggeno.test)
        return(X.tilde)
      } else {
        X.tilde <- predict(trlist[[gname]]$fit, ggeno.test, s = "lambda.1se")
        return(X.tilde)
      }
    }

    # ncore <- Ncore # note this is not working, too much memory usage
    # cl <- makeCluster(ncore,outfile="")
    # registerDoParallel(cl)
    # print(paste0("start parallel computing using ", ncore, " cores ..."))

    implist <- list()
    for (i in 1:length(gnames)){ #
      print(i)
      implist[[i]] <- impute_expr_g(i)
    }

    # implist <- foreach(i= 1:length(gnames),.packages = c("glmnet")) %dopar% {
    #   print(i)
    #   impute_expr_g(i)
    # }
    # print("end parallel computing...")
    # stopCluster(cl)
    names(implist) <- gnames
    save(implist, file = paste0(outname, ".expr_imputed.Rdata"))

    outdf <- round(t(do.call(cbind, implist)),3)
    outdf2 <- cbind(gnames, "A", "B", outdf)
    write.table(outdf2 , file= exprf , row.names=F, col.names=F, sep=" ", quote = F)
  }

  cat(" * Computing correlation matrix for expression.\n")
  exprdir <-  normalizePath(dirname(exprf))
  exprkdir <- paste0(exprdir, "/", basename(exprf), ".kdir/")
  exprfn <- basename(exprf)
  if (!dir.exists(exprkdir)){
    dir.create(exprkdir)
    print(paste("cd", exprkdir, "; gemma -g", exprf, "-p", phenof, "-gk 1 -notsnp -o", paste0(exprfn,"-expr")))
    system(paste("cd", exprkdir, "; gemma -g", exprf, "-p", phenof, "-gk 1 -notsnp -o", paste0(exprfn,"-expr")))
  }

  # Estimate variance component by gemma
  cat(" * Estimating variance components.\n")
  vcdir <- paste0(outname,"_VC/")
  if (!dir.exists(vcdir)){
    dir.create(vcdir)
    vcdir <- normalizePath(vcdir)
    mfile <- paste0(vcdir,"/", outname,"-Kfiles.txt")
    file1 <- paste0(exprkdir,"/output/", exprfn, "-expr.cXX.txt")
    file2 <- paste0(genokdir,"/output/", basename(genof), "-geno.cXX.txt")
    fileConn <- file(mfile)
    writeLines(c(file1, file2), fileConn)
    close(fileConn)
    system(paste("cd", vcdir, "; gemma -p", phenof, "-mk", mfile,  "-n 1 -vc 1 -o", outname))
  }
}


args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("4 arguments must be supplied:  phenofile, genofilebase, exprRdfile, outname\n
       phenofile: phenotype file in BIMBAM format\n
       genofilebase: genofilebase.geno.txt.gz is genotype file in BIMBAM format,        genofilebase.RData is genotype file in R data format\n
       exprRdfile: expression file in Rdata format, data used to train expression model, have the following the variables: G, labels, chr, pos, gnames, expr.
      ", call.=FALSE)
}

phenof <-  normalizePath(args[1])
genofbase <- args[2]
exprfR <- args[3]
outname <- args[4]

gafile <- "/home/simingz/causalTWAS/annotation/refGene_processed_5e+05.txt"

exprtrf <- paste0(exprfR, ".expr_train_model.Rdata")
if (!file.exists(exprtrf)){
  print(" * training expression model")
  load(exprfR)
  geno.tr <- G; rm(G)
  labels.tr <- labels; rm(labels)
  chr.tr <- chr; rm(chr)
  pos.tr <- pos; rm(pos)
  gnames.tr <- gnames; rm(gnames)
  expr.tr <- expr; rm(expr)
  train_expr(expr.tr, gnames.tr, geno.tr, chr.tr, pos.tr, labels.tr, gafile, exprfR)
  gc()
}

genof <- normalizePath(paste0(genofbase, ".geno.txt.gz"))
genofR <- paste0(genofbase, ".RData")
load(genofR)
geno <- X; rm(X); gc()
run_ctwas(phenof, genof, exprtrf, gafile, geno, chr, pos, labels, outname)
