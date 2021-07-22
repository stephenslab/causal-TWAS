source("~/causalTWAS/causal-TWAS/analysis/summarize_basic_plots.R")

#
caliMR_plot <- function(phenofs, mrfs, main = "MR-JTI Bonferroni 0.05"){
  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})
  df <- NULL
  for (i in 1:length(mrfs)) {
    res <- read.table(mrfs[i], header = T)
    res <- res[complete.cases(res),]
    res$FDR <- 1
    res[res$CI_significance== "sig", "FDR"] <- 0.05
    res$ifcausal <- ifelse(res$variable %in% cau[[i]], 1, 0)
    res$runtag <- i
    df <- rbind(df, res)
  }
  fig <- cp_plot(df$FDR, df$ifcausal, df$runtag, mode ="FDR", main = main)
  return(fig)
}


ncausalMR_plot <- function(phenofs, mrfs, main = "MR-JTI Bonferroni 0.05"){

  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})

  df <- NULL
  for (i in 1:length(mrfs)) {
    res <- fread(mrfs[i], header = T)
    res <- res[complete.cases(res),]
    res$FDR <- 1
    res[res$CI_significance== "sig", "FDR"] <- 0.05
    res$ifcausal <- ifelse(res$variable %in% cau[[i]], 1, 0)
    res$runtag <- i
    df <- rbind(df, res)
  }

  fig <- nca_plot(df$FDR, df$ifcausal, df$runtag, mode ="FDR", xmin = 0.2, main = main)
  return(fig)
}
