get_files <- function(tag, tag2){
  par <- paste0(outputdir, tag, "-mr.ash2s.", tag2, ".param.txt")
  rpip <- paste0(outputdir, tag, "-mr.ash2s.", tag2, ".rPIP.txt")

  gmrash <- paste0(outputdir, tag, "-mr.ash2s.", tag2, ".expr.txt")
  smrash <- paste0(outputdir, tag, "-mr.ash2s.", tag2, ".snp.txt")

  ggwas <- paste0(outputdir, tag, ".exprgwas.txt.gz")
  sgwas <- paste0(outputdir, tag, ".snpgwas.txt.gz")

  gsusie <- paste0(susiedir, tag, ".", tag2, ".L3.susieres.expr.txt")
  ssusie <- paste0(susiedir, tag, ".", tag2, ".L3.susieres.snp.txt")

  return(tibble::lst(par, rpip, gmrash, ggwas, smrash, sgwas, gsusie, ssusie))
}

get_tags <- function(globpattern, extrpattern, tag2){
  lapply(lapply(get_files(globpattern, tag2), Sys.glob), str_extract, pattern = extrpattern)
}

show_param <- function(tags, tag2){
  f <- lapply(tags, get_files, tag2 = tag2)
  parf <- lapply(f, '[[', "par")
  param <- do.call(rbind, lapply(parf, function(x) t(read.table(x))[2:1,]))
  truth <- param[1:(nrow(param)/2)*2-1,]
  est <- param[1:(nrow(param)/2)*2,]
  outdt <- matrix(0, ncol = 2*ncol(param), nrow = nrow(param)/2)
  outdt[,c(1,3,5,7)] <- truth
  outdt[,c(2,4,6,8)] <- est
  outdt <-cbind(1:nrow(outdt),outdt)
  colnames(outdt) <- c("Simulation#", paste0(rep(c("Truth","Est."),4)))
  knitr::kable(outdt) %>%
    kable_styling("striped") %>%
    add_header_above(c(" " = 1, "Gene.pi1" = 2, "Gene.PVE" = 2, "SNP.pi1" = 2, "SNP.PVE" =2))
}

cp_plot <- function(s, ifcausal, mode = c("PIP", "FDR"), main = mode[1], ...){
  # ifcausal:0,1
  a_bin <- cut(s, breaks= seq(0,1, by=0.1))
  if (mode == "PIP") {
    Expected = c(by(s, a_bin, FUN = mean))
    Observed = c(by(ifcausal, a_bin, FUN = mean))
  } else if (mode == "FDR"){
    Expected = c(by(s, a_bin, FUN = mean))
    Observed = 1 - c(by(ifcausal, a_bin, FUN = mean))
  }
  plot(Expected, Observed, xlim= c(0,1), ylim=c(0,1), pch =19, main = main, ...)
  lines(x = c(0,1), y = c(0,1), col ="grey", lty = 2)
}

caliPIP_plot <- function(tags, tag2){
  f <- lapply(tags, get_files, tag2 = tag2)
  mrashf <- lapply(f, '[[', "gmrash")
  names(mrashf) <- tags

  susief <- lapply(f, '[[', "gsusie")
  names(susief) <- tags

  .tagname <- function(x, flist){
    a <- read.table(flist[[x]], header =T)
    a[, "name"] <- paste0(x, ":", a[, "name"])
    a
  }
  mrashres <- do.call(rbind, lapply(tags, .tagname, flist = mrashf))
  susieres <- do.call(rbind, lapply(tags, .tagname, flist = susief))

  res <- merge(mrashres, susieres, by = "name", all = T)

  res <- res[complete.cases(res),]
  res <- rename(res, c("PIP" = "mr.ash_PIP", "pip" = "SUSIE.w_PIP", "pip.null" = "SUSIE.u_PIP", "pip.w0" = "SUSIE.w0_PIP"))
  par(mfrow=c(1,4), mar=c(5, 6, 4, 1))
  cp_plot(res$mr.ash_PIP, res$ifcausal, main = "Mr.ash PIP")
  cp_plot(res$SUSIE.w_PIP, res$ifcausal, main = "SUSIE.weighted PIP")
  cp_plot(res$SUSIE.u_PIP, res$ifcausal, main = "SUSIE.uniform PIP")
  cp_plot(res$SUSIE.w0_PIP, res$ifcausal, main = "SUSIE.weighted-w-null PIP")
  return(res)
}


#' s is pip or fdr.
caliFDR_plot <- function(tags, tag2){

  f <- lapply(tags, get_files, tag2 = tag2)
  gwasf <- lapply(f, '[[', "ggwas")
  names(gwasf) <- tags

  .tagname <- function(x, flist, colnames = NULL){
    a <- read.table(flist[[x]], header =T)
    if (!is.null(colnames)){
      colnames(a) <- colnames
    }
    a[, "name"] <- paste0(x, ":", a[, "name"])
    a
  }

  .addFDR <- function(res){
    res$FDR <- p.adjust(res$PVALUE, method = "fdr")
    res
  }


  gwasres <- do.call(rbind, lapply(lapply(tags, .tagname, flist = gwasf,
                                          colnames =  c("chr",	"p0",	"p1", "name",
                                                        "Estimate", "Std.Error", "t-value", "PVALUE")), .addFDR))

  f <- lapply(tags, get_files, tag2 = tag2)
  mrashf <- lapply(f, '[[', "gmrash")
  names(mrashf) <- tags

  susief <- lapply(f, '[[', "gsusie")
  names(susief) <- tags

  .tagname <- function(x, flist){
    a <- read.table(flist[[x]], header =T)
    a[, "name"] <- paste0(x, ":", a[, "name"])
    a
  }
  mrashres <- do.call(rbind, lapply(tags, .tagname, flist = mrashf))
  susieres <- do.call(rbind, lapply(tags, .tagname, flist = susief))

  res <- merge(mrashres, susieres, by = "name", all = T)

  res <- merge(res, gwasres, by = "name", all = T)

  res <- res[complete.cases(res),]

  cp_plot(res$FDR, res$ifcausal, mode ="FDR", main = "TWAS FDP")
  cat("FDP at bonferroni corrected p = 0.05: ", 1 - mean(res[res$PVALUE < 0.05 /dim(res)[1], "ifcausal"]))
  return(res)
}

scatter_plot_PIP<- function(tags, tag2){
  f <- lapply(tags, get_files, tag2 = tag2)
  mrashf <- lapply(f, '[[', "gmrash")
  names(mrashf) <- tags

  susief <- lapply(f, '[[', "gsusie")
  names(susief) <- tags

  .tagname <- function(x, flist){
    a <- read.table(flist[[x]], header =T)
    a[, "name"] <- paste0(x, ":", a[, "name"])
    a
  }
  mrashres <- do.call(rbind, lapply(tags, .tagname, flist = mrashf))
  susieres <- do.call(rbind, lapply(tags, .tagname, flist = susief))

  res <- merge(mrashres, susieres, by = "name", all = T)

  res <- res[complete.cases(res),]
  res <- rename(res, c("PIP" = "mr.ash_PIP", "pip" = "SUSIE.w_PIP", "pip.null" = "SUSIE.u_PIP", "pip.w0" = "SUSIE.w0_PIP") )
  res$ifcausal <- mapvalues(res$ifcausal,
                            from=c(0,1),
                            to=c("Non causal", "Causal"))

  fig1 <- plot_ly(data = res, x = ~ mr.ash_PIP, y = ~ SUSIE.w_PIP, color = ~ ifcausal,
                  colors = c( "salmon", "darkgreen"), text = ~ paste("Name: ", name,
                                                                     "\nChr: ", chr.x,
                                                                     "\nPos:", p0))

  fig2 <- plot_ly(data = res, x = ~ mr.ash_PIP, y = ~ SUSIE.u_PIP, color = ~ ifcausal,
                  colors = c( "salmon", "darkgreen"), text = ~ paste("Name: ", name,
                                                                     "\nChr: ", chr.x,
                                                                     "\nPos:", p0))

  fig3 <- plot_ly(data = res, x = ~ mr.ash_PIP, y = ~ SUSIE.w0_PIP, color = ~ ifcausal,
                  colors = c( "salmon", "darkgreen"), text = ~ paste("Name: ", name,
                                                                     "\nChr: ", chr.x,
                                                                     "\nPos:", p0))

  fig <- subplot(fig1, fig2, fig3, titleX = TRUE, titleY = T, margin = 0.05)
  fig
}


scatter_plot_PIP_p <- function(tags, tag2){
  f <- lapply(tags, get_files, tag2 = tag2)
  mrashf <- lapply(f, '[[', "gmrash")
  names(mrashf) <- tags

  susief <- lapply(f, '[[', "gsusie")
  names(susief) <- tags

  gwasf <- lapply(f, '[[', "ggwas")
  names(gwasf) <- tags

  .tagname <- function(x, flist, colnames = NULL){
    a <- read.table(flist[[x]], header =T)
    if (!is.null(colnames)){
      colnames(a) <- colnames
    }
    a[, "name"] <- paste0(x, ":", a[, "name"])
    a
  }
  mrashres <- do.call(rbind, lapply(tags, .tagname, flist = mrashf))
  susieres <- do.call(rbind, lapply(tags, .tagname, flist = susief))
  gwasres <- do.call(rbind, lapply(tags, .tagname, flist = gwasf,
                                   colnames =  c("chr",	"p0",	"p1", "name", "Estimate", "Std.Error", "t-value", "PVALUE")))

  res <- merge(mrashres, susieres, by = "name", all = T)
  res <- merge(res, gwasres, by = "name", all = T)

  res <- res[complete.cases(res),]
  res <- rename(res, c("PIP" = "mr.ash", "pip" = "SUSIE.w", "pip.null"= "SUSIE.u", "pip.w0" = "SUSIE.w0", "PVALUE" = "TWAS") )
  res[,"TWAS"] <- -log10(res[, "TWAS"])
  res$ifcausal <- mapvalues(res$ifcausal,
                            from=c(0,1),
                            to=c("Non causal", "Causal"))

  fig1 <- plot_ly(data = res, x = ~ TWAS, y = ~ mr.ash, color = ~ ifcausal,
                  colors = c( "salmon", "darkgreen"), text = ~ paste("Name: ", name,
                                                                       "\nChr: ", chr.x,
                                                                       "\nPos:", p0.x))

  fig2 <- plot_ly(data = res, x = ~ TWAS, y = ~ SUSIE.w, color = ~ ifcausal,
                  colors = c( "salmon", "darkgreen"), text = ~ paste("Name: ", name,
                                                                     "\nChr: ", chr.x,
                                                                     "\nPos:", p0.x))

  fig3 <- plot_ly(data = res, x = ~ TWAS, y = ~ SUSIE.u, color = ~ ifcausal,
                  colors = c( "salmon", "darkgreen"), text = ~ paste("Name: ", name,
                                                                     "\nChr: ", chr.x,
                                                                     "\nPos:", p0.x))

  fig4 <- plot_ly(data = res, x = ~ TWAS, y = ~ SUSIE.w0, color = ~ ifcausal,
                  colors = c( "salmon", "darkgreen"), text = ~ paste("Name: ", name,
                                                                     "\nChr: ", chr.x,
                                                                     "\nPos:", p0.x))

  fig <- subplot(fig1, fig2, fig3, fig4, nrows=2, titleY = T, margin = 0.05)
  print(fig)
  return(res)
}

ROC_plot<- function(tags, tag2){
  f <- lapply(tags, get_files, tag2 = tag2)
  mrashf <- lapply(f, '[[', "gmrash")
  names(mrashf) <- tags

  susief <- lapply(f, '[[', "gsusie")
  names(susief) <- tags

  gwasf <- lapply(f, '[[', "ggwas")
  names(gwasf) <- tags

  .tagname <- function(x, flist, colnames = NULL){
    a <- read.table(flist[[x]], header =T)
    if (!is.null(colnames)){
      colnames(a) <- colnames
    }
    a[, "name"] <- paste0(x, ":", a[, "name"])
    a
  }
  mrashres <- do.call(rbind, lapply(tags, .tagname, flist = mrashf))
  susieres <- do.call(rbind, lapply(tags, .tagname, flist = susief))
  gwasres <- do.call(rbind, lapply(tags, .tagname, flist = gwasf,
                                   colnames =  c("chr",	"p0",	"p1", "name", "Estimate", "Std.Error", "t-value", "PVALUE")))

  res <- merge(mrashres, susieres, by = "name", all = T)
  res <- merge(res, gwasres, by = "name", all = T)

  res <- res[complete.cases(res),]
  res <- rename(res, c("PIP" = "mr.ash", "pip" = "SUSIE.w", "pip.null"= "SUSIE.u", "pip.w0" = "SUSIE.w0", "PVALUE" = "TWAS") )
  res[,"TWAS"] <- -log10(res[, "TWAS"])

  roccolors <-  c("red", "green", "orange", "pink", "blue")
  methods <- c("mr.ash", "SUSIE.w", "SUSIE.u", "SUSIE.w0", "TWAS")
  plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab = "FPR", ylab = "TPR")
  for (i in 1:length(methods)){
    method <- methods[i]
    bordered <- res[order(res[,method]),]
    actuals <- bordered$ifcausal == 1
    sens <- (sum(actuals) - cumsum(actuals))/sum(actuals)
    spec <- cumsum(!actuals)/sum(!actuals)
    lines(1 - spec, sens, type = "l", col = roccolors[i])
    abline(c(0,0),c(1,1))
    auc <- sum(spec*diff(c(0, 1 - sens)))
    cat("AUC for ", method, ": ", auc)
  }
  legend(0.6,0.3, legend= methods, col=roccolors, lty=1, cex=0.5 )
  grid()
}
