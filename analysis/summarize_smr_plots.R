source("~/causalTWAS/causal-TWAS/analysis/summarize_basic_plots.R")

#
caliSMRp_plot <- function(phenofs, smrfs, main = "SMR-HEIDI FDP"){
  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})
  df <- NULL
  for (i in 1:length(smrfs)) {
    res <- read.table(smrfs[i], header = T)
    res <- res[complete.cases(res),]
    res <- res[res$p_HEIDI > 0.05,]
    res$FDR <- p.adjust(res$p_SMR, method = "fdr")
    res$ifcausal <- ifelse(res$probeID %in% cau[[i]], 1, 0)
    res$runtag <- i
    df <- rbind(df, res)
  }
  fig <- cp_plot(df$FDR, df$ifcausal, df$runtag, mode ="FDR", main = main)
  return(fig)
}


ncausalSMRp_plot <- function(phenofs, smrfs, main = "SMR-HEIDI FDP"){

  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})

  df <- NULL
  for (i in 1:length(smrfs)) {
    res <- fread(smrfs[i], header = T)
    res <- res[complete.cases(res),]
    res <- res[res$p_HEIDI > 0.05,]
    res$FDR <- p.adjust(res$p_SMR, method = "fdr")
    res$ifcausal <- ifelse(res$probeID %in% cau[[i]], 1, 0)
    res$runtag <- i
    df <- rbind(df, res)
  }

  fig <- nca_plot(df$FDR, df$ifcausal, df$runtag, mode ="FDR", xmin = 0.2, main = main)
  return(fig)
}
