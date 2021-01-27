library(ggplot2)
library(cowplot)
library(plotrix)


# Calculate power in a GWAS study
pow <- function(total, n, beta, cutp){
  rec <- rep(0, total)
  for (i in 1:total){
    x <- rnorm(n)
    y <- x * rnorm(1, sd = beta) + rnorm(n, sd = sqrt(2.5))
    lm.s <- lm(y~x)
    pv <- summary(lm.s)$coefficients[2,4]
    rec[i] <- pv
  }
  length(rec[rec < cutp])/length(rec)
}

# p value distribution

pdist_plot <- function(gwasf, chr, cau, thin = 1){

  a <- data.table::fread(gwasf, header = T)
  a <- a[a$chrom == chr]

  a$ifcausal <- ifelse(a$id %in% cau, 1, 0)
  a.causal <- a[a$ifcausal == 1, ]
  a <- a[sample(nrow(a), nrow(a)* thin), ]
  a.nonc <- a[a$ifcausal == 0, ]
  a <- rbind(a.causal, a.nonc)

  ax <- pretty(0:max(-log10(a$PVALUE)), n = 30)

  if (!("p0" %in% colnames(a))) a$p0 <- a$pos

  par(mfrow=c(3,1))

  h1 <- hist(-log10(a$PVALUE), breaks = 100, xlab = "-log10(p)", main = "P value distribution-all", col = "grey", xlim= c(3,20), ylim =c(0,50)); grid()

  h2 <- hist(-log10(a[a$ifcausal == 1, ]$PVALUE), breaks = h1$breaks, xlab = "-log10(p)", main = "P value distribution-causal", col = "salmon", xlim= c(3,20), ylim =c(0,50));grid()

  plot(a$p0, -log10(a$PVALUE), col = a$ifcausal + 1, xlab = paste0("chr", chrom), ylab = "-log10(pvalue)")
  points(a[a$ifcausal ==1, ]$p0, -log10(a[a$ifcausal ==1, ]$PVALUE), col = "red", pch =19)
  grid()

  # cat("number of p < ", p.cut, ": ",nrow(a[a$PVALUE < p.cut, ]))
  # cat("number of causal p < ", p.cut, ": ", nrow(a[a$PVALUE < p.cut & a$ifcausal ==1, ]))

}


# used by show_param, caliPIP_plot
get_causal_id <- function(phenores){
  gene.causal <- unlist(sapply(1:22, function(x) phenores$batch[[x]]$id.cgene))
  snp.causal <- unlist(sapply(1:22, function(x) phenores$batch[[x]]$id.cSNP))
  return(c(gene.causal, snp.causal))
}

# used by nca_plot
.obn <- function(pips, ifcausal, mode = c("PIP", "FDR")){
  a_bin <- cut(pips, breaks= seq(0, 1, by=0.1))
  if (mode == "PIP") {
    ob = c(by(ifcausal, a_bin, FUN = sum))
  } else if (mode == "FDR"){
    ob = c(by((1-ifcausal), a_bin, FUN = sum))
  }
  return(ob)
}

# used by ncausal plot
nca_plot <- function(pips, ifcausal, runtag = NULL, mode = c("PIP", "FDR"), main = mode[1], ...){

  # ifcausal:0,1, runtag: for adding std.
  if (is.null(runtag)){
    se = 0
  } else{
    dflist <- list()
    for (rt in unique(runtag)){
      pips.rt <- pips[runtag == rt]
      ifcausal.rt <- ifcausal[runtag == rt]
      dflist[[rt]] <- .obn(pips.rt, ifcausal.rt, mode = mode)
    }
    nca_mean <- colMeans(do.call(rbind, dflist))
    se <- apply(do.call(rbind, dflist), 2, plotrix::std.error)
  }

  df <- data.frame("ncausal" = nca_mean, "se" = se)

  fig <- ggplot(df) +
    geom_bar( aes(x=seq(0, 1, by=0.1)[1:10]+0.05, y=ncausal), color ="black", stat="identity", fill="salmon", alpha=0.7, width = 0.1) +
    geom_errorbar( aes(x= seq(0, 1, by=0.1)[1:10] + 0.05, ymin= ncausal-se, ymax=ncausal+se), width=0.05, colour="black", alpha=0.9, size=0.5) +
    ggtitle(main) +
    xlab("PIP") + ylab("No. causal genes")
  theme_cowplot()
  return(fig)
}

# Power plot, plot the number of causal genes in PIP bins.
ncausal_plot <- function(phenofs, pipfs, main = "PIP"){

  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})

  df <- NULL
  for (i in 1:length(pipfs)) {
    res <- fread(pipfs[i], header = T)
    res <- data.frame(res[res$type  =="gene", ])
    res$ifcausal <- ifelse(res$id %in% cau[[i]], 1, 0)
    res$runtag <- i
    res <- res[complete.cases(res),]
    df <- rbind(df, res)
  }

  fig <- nca_plot(df$susie_pip, df$ifcausal, df$runtag, mode ="PIP", main = main)
  return(fig)
}

# used by cp_plot
.exob <- function(pips, ifcausal, mode = c("PIP", "FDR")){
  a_bin <- cut(pips, breaks= seq(0, 1, by=0.1))
  if (mode == "PIP") {
    ex = c(by(pips, a_bin, FUN = mean))
    ob = c(by(ifcausal, a_bin, FUN = mean))
  } else if (mode == "FDR"){
    ex = c(by(pips, a_bin, FUN = mean))
    ob = 1 - c(by(ifcausal, a_bin, FUN = mean))
  }
  return(list("expected" = ex, "observed" = ob))
}

# used by caliPIP_plot, caliFDP_plot
cp_plot <- function(pips, ifcausal, runtag = NULL, mode = c("PIP", "FDR"), main = mode[1], ...){

  # ifcausal:0,1, runtag: for adding std.
  if (is.null(runtag)){
    se = 0
  } else{
    dflist <- list()
    for (rt in unique(runtag)){
      pips.rt <- pips[runtag == rt]
      ifcausal.rt <- ifcausal[runtag == rt]
      dflist[[rt]] <- .exob(pips.rt, ifcausal.rt, mode = mode)
    }
    mean_pip <- colMeans(do.call(rbind, lapply(dflist, '[[', "expected")))
    observed_freq <- colMeans(do.call(rbind, lapply(dflist, '[[', "observed")))
    se <- apply(do.call(rbind, lapply(dflist, '[[', "observed")), 2, plotrix::std.error)
  }

  df <- data.frame("mean_pip" = mean_pip, "observed_freq"= observed_freq, "se" = se)

  ggplot(df, aes(x=mean_pip, y=observed_freq)) +
    geom_errorbar(aes(ymin=observed_freq-se, ymax=observed_freq+se), colour="black", size = 0.5, width=.01) +
    geom_point(size=1.5, shape=21, fill="#002b36") + # 21 is filled circle
    xlab("Expected") +
    ylab("Observed") +
    coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
    geom_abline(slope=1,intercept=0,colour='red', size=0.2) +
    ggtitle(main) +
    expand_limits(y=0) +                        # Expand y range
    theme_cowplot()

  #plot(Expected, Observed, xlim= c(0,1), ylim=c(0,1), pch =19, main = main, ...)
  #lines(x = c(0,1), y = c(0,1), col ="grey", lty = 2)
}

# PIP calibration plot
caliPIP_plot <- function(phenofs, pipfs, main = "PIP Calibration"){
  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})
  df <- NULL
  for (i in 1:length(pipfs)) {
    res <- fread(pipfs[i], header = T)
    res <- data.frame(res[res$type  =="gene", ])
    res$runtag <- i
    res <- res[complete.cases(res),]
    res$ifcausal <- ifelse(res$id %in% cau[[i]], 1, 0)
    df <- rbind(df, res)
  }

  fig <- cp_plot(df$susie_pip, df$ifcausal, df$runtag, mode ="PIP", main = main)
  return(fig)
}

# TWAS FDP calibration
caliFDP_plot <- function(phenofs, gwasfs, main = "FDP"){
  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})
  df <- NULL
  for (i in 1:length(gwasfs)) {
    res <- read.table(gwasfs[i], header = T)
    res$FDR <- p.adjust(gwasres$PVALUE, method = "fdr")
    res$ifcausal <- ifelse(res$id %in% cau[[i]], 1, 0)
    res <- res[complete.cases(res),]
    df <- rbind(df, res)
  }
  fig <- cp_plot(df$FDR, df$ifcausal, df$runtag, mode ="FDR", main = main)
  cat("FDP at bonferroni corrected p = 0.05: ", 1 - mean(df[df$PVALUE < 0.05 /dim(df)[1], "ifcausal"]))
  return(fig)
}

scatter_plot_PIP_p <- function(phenofs, pipfs, gwasfs, main ="PIP-p"){
  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})
  df <- NULL
  for (i in 1:length(pipfs)) {
    pipres <- fread(pipfs[i], header = T)
    pipres <- data.frame(pipres[pipres$type == "gene", ])
    pipres$runtag <- i
    res$ifcausal <- ifelse(res$id %in% cau[[i]], 1, 0)
    gwasres <- read.table(gwasfs[i], header = T)
    res <- merge(gwasres, pipres, by = "id", all = T)
    res <- res[complete.cases(res),]
    df <- rbind(df, res)
  }

  df <- rename(df, c( "PVALUE" = "TWAS.p"))
  df[,"TWAS.p"] <- -log10(df[, "TWAS.p"])

  df$ifcausal <- mapvalues(df$ifcausal, from=c(0,1),to=c("darkgreen", "salmon"))
  plot(df$TWAS.p, df$susie_pip, col = df$ifcausal, main = main, xlab = "-log10(TWAS p value)", ylab = "PIP")

  # df$ifcausal <- mapvalues(df$ifcausal, from=c(0,1), to=c("Non causal", "Causal"))
  # fig <- plot_ly(data = df, x = ~ TWAS.p, y = ~ pip, color = ~ ifcausal,
  #                colors = c( "salmon", "darkgreen"), type ="scatter", text = ~ paste("Name: ", paste0(runtag,":",name),
  #                                                                  "\nChr: ", chr,  "\nPos:", pos))

}

# used by show_param
get_pi1 <- function(phenores){
  gene.causal <- sum(sapply(1:22, function(x) length(phenores$batch[[x]]$idx.cgene)))
  gene.total <- sum(sapply(1:22, function(x) phenores$batch[[x]]$J))

  snp.causal <- sum(sapply(1:22, function(x) length(phenores$batch[[x]]$idx.cSNP)))
  snp.total <- sum(sapply(1:22, function(x) phenores$batch[[x]]$M))

  return(c(gene.causal/gene.total, snp.causal/snp.total))
}

# return data.frame, rows are each simu run, columns are different parameters.
show_param <- function(phenofs, susieIfs, susieIfs2, thin = 1){

  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})

  # overall truth
  truth <- do.call(rbind, lapply(phenofs, function(x) {load(x);
    c(phenores$param$pve.gene.truth, phenores$param$pve.snp.truth,
      get_pi1(phenores))}))
  colnames(truth) <- c("PVE.gene_truth", "PVE.SNP_truth", "pi1.gene_truth", "pi1.SNP_truth")

  # truth in selected regions for param estimation
  truth.se <- do.call(rbind, lapply(1: length(susieIfs2), function(x) {
    a <- fread(susieIfs2[x], header = T)
    a$ifcausal <- ifelse(a$id %in% cau[[x]], 1, 0)
    c(nrow(a[a$ifcausal == 1 & a$type == "gene" ])/ nrow(a[a$type == "gene"]),
      nrow(a[a$ifcausal == 1 & a$type == "SNP"])/ nrow(a[a$type == "SNP"]))
  }))
  colnames(truth.se) <- c("pi1.gene_se_truth", "pi1.SNP_se_truth")

  # estimated parameters
  est <- do.call(rbind, lapply(susieIfs, function(x) {load(x); group_prior_rec[, ncol(group_prior_rec)]}))
  colnames(est) <- c("pi1.gene_est", "pi1.SNP_est")

  est[, "pi1.SNP_est"] <- est[, "pi1.SNP_est"] * thin

  mtx <- cbind(truth, truth.se, est)
  return(mtx)

  # df %>%
  # kable("html", escape = F) %>%
  # kable_styling("striped", full_width = F) %>%
  #  row_spec(c(1:5, 11:15), background = "#FEF3B9") %>%
  # scroll_box(width = "100%", height = "600px", fixed_thead = T)
}

plot_param <- function(mtx){
  n <- nrow(mtx)
  # plot gene pi1
  y1 <- mtx[,"pi1.gene_truth"]
  y2 <- mtx[,"pi1.gene_se_truth"]
  y3 <- mtx[,"pi1.gene_est"]
  plot(jitter(rep(1, n)), y1, col = 1:n, pch = 19, ylab = "gene pi1", xaxt = "n", xlab="", xlim = c(0.8, 3.5), ylim = c(0, max(c(y1, y2, y3)) * 1.05) ,frame.plot=FALSE)
  points(jitter(rep(2, n)), y2, col = 1:n, pch =19)
  points(jitter(rep(3, n)), y3, col = 1:length(y1), pch =19)
  axis(side=1, at=1:2, labels = FALSE, tick = F)
  text(x=1:3, 0, labels = c( "truth", "truth(selected)", "ctwas"), xpd = T, pos =1)
  grid()

  # plot SNP pi1
  y1 <- mtx[,"pi1.SNP_truth"]
  y2 <- mtx[,"pi1.SNP_se_truth"]
  y3 <- mtx[,"pi1.SNP_est"]
  plot(jitter(rep(1, n)), y1, col = 1:n, pch = 19, ylab = "SNP pi1", xaxt = "n", xlab="", xlim = c(0.8, 3.5), ylim = c(0, max(c(y1, y2, y3)) * 1.05) ,frame.plot=FALSE)
  points(jitter(rep(2, n)), y2, col = 1:n, pch =19)
  points(jitter(rep(3, n)), y3, col = 1:length(y1), pch =19)
  axis(side=1, at=1:2, labels = FALSE, tick = F)
  text(x=1:3, 0, labels = c( "truth", "truth(selected)", "ctwas"), xpd = T, pos =1)
  grid()

  #plot enrichment

  y1 <- mtx[,"pi1.gene_truth"]/mtx[,"pi1.SNP_truth"]
  y2 <- mtx[,"pi1.gene_se_truth"]/mtx[,"pi1.SNP_se_truth"]
  y3 <- mtx[,"pi1.gene_est"]/mtx[,"pi1.SNP_est"]
  plot(jitter(rep(1, n)), y1, col = 1:n, pch = 19, ylab = "Enrich", xaxt = "n", xlab="", xlim = c(0.8, 3.5), ylim = c(0, max(c(y1, y2, y3)) * 1.05) ,frame.plot=FALSE)
  points(jitter(rep(2, n)), y2, col = 1:n, pch =19)
  points(jitter(rep(3, n)), y3, col = 1:length(y1), pch =19)
  axis(side=1, at=1:2, labels = FALSE, tick = F)
  text(x=1:3, 0, labels = c( "truth", "truth(selected)", "ctwas"), xpd = T, pos =1)
  grid()
}


