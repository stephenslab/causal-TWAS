train_expr <- function(expr.tr, gnames.tr, geno.tr, chr.tr, pos.tr, labels.tr, gafile, outname){
  # expr.tr:  expression matrix for training data, column genes row samples
  # gnames.tr: vector of gene names, they should match the columns of expr.tr
  # geno.tr: genotype matrix for training data, column SNPs row samples
  # chr.tr: vector, chroms for each SNP for training data
  # pos.tr: vector, positions for each SNP for training data
  # labels.tr: vector, SNP labels for each SNP for training data
  # gafile: Gene annotation file: NM_130786 chr19  - 63048360 64056670    A1BG

  N <- dim(geno.tr)[1]
  M <- dim(geno.tr)[2]

  ga <- read.table(gafile, header = F, stringsAsFactors = F)

  train_expr_g <- function(i){
    # i is the gene index.
    gname <- gnames.tr[i]
    gai <- ga[ga$V6 == gname,]
    idx <- c(1:M)[paste0("chr", chr.tr) == gai[[2]] & pos.tr > gai[[4]] & pos.tr < gai[[5]]]
    ggeno <- as.matrix(geno.tr[,idx])
    if (dim(ggeno)[2] == 0) {
      return(NULL)
    } else if (dim(ggeno)[2] == 1){
      lmfit <- lm(expr.tr[,i] ~ as.matrix(ggeno))
      res <- list("fit" = lmfit, "SNPlabel" = labels.tr[idx])
      return(res)
    } else {
      cvfit <- cv.glmnet(ggeno, expr.tr[,i])
      res <- list("fit" = cvfit, "SNPlabel" = labels.tr[idx])
      return(res)
    }
  }

  # trlist <- lapply(1:length(gnames.tr), train_expr_g)

  ncore <- Ncore
  cl <- makeCluster(ncore,outfile="")
  registerDoParallel(cl)
  print(paste0("start parallel computing using ", ncore, " cores ..."))
  trlist <- foreach(i= 1:length(gnames.tr),.packages = c("glmnet", "data.table")) %dopar% {
    print(i)
    train_expr_g(i)
  }
  print("end parallel computing...")
  stopCluster(cl)
  names(trlist) <- gnames.tr
  save(trlist, file = paste0(outname, ".expr_train_model.Rdata"))
}
