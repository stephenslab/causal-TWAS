source("~/causalTWAS/causal-TWAS/analysis/summarize_basic_plots.R")

# TWAS bonferroni calibration
caliFUSIONbon_plot <- function(phenofs, fusionfs, main = "FUSION bonferroni"){
  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})
  df <- NULL
  for (i in 1:length(fusionfs)) {
    res <- read.table(fusionfs[i], header = T)
    res$FDR <- p.adjust(res$TWAS.P, method = "bonferroni")
    res$ifcausal <- ifelse(res$ID %in% cau[[i]], 1, 0)
    res$runtag <- i
    df <- rbind(df, res)
  }
  fig <- cp_plot(df$FDR, df$ifcausal, df$runtag, mode ="FDR", main = main)
  return(fig)
}

# TWAS FDP calibration
caliFUSIONp_plot <- function(phenofs, fusionfs, main = "FUSION FDR"){
  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})
  df <- NULL
  for (i in 1:length(fusionfs)) {
    res <- read.table(fusionfs[i], header = T)
    res$FDR <- p.adjust(res$TWAS.P, method = "fdr")
    res$ifcausal <- ifelse(res$ID %in% cau[[i]], 1, 0)
    res$runtag <- i
    df <- rbind(df, res)
  }
  fig <- cp_plot(df$FDR, df$ifcausal, df$runtag, mode ="FDR", main = main)
  return(fig)
}


# coloc PP4 calibration plot
caliPP4_plot <- function(phenofs, colocfs, twas.p = 1e-4, main = "coloc PP4 Calibration"){
  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})
  df <- NULL
  for (i in 1:length(colocfs)) {
    res <- fread(colocfs[i], header = T)
    res$runtag <- i
    res <- res[complete.cases(res),]
    res$ifcausal <- ifelse(res$ID %in% cau[[i]], 1, 0)
    res <- res[res$TWAS.P < twas.p,]
    df <- rbind(df, res)
  }

  fig <- cp_plot(df$COLOC.PP4, df$ifcausal, df$runtag, mode ="PIP", main = main)
  return(fig)
}

# power plot FUSION bonferroni
ncausalFUSIONbon_plot <- function(phenofs, fusionfs, main = "FUSION Bonferroni"){

  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})

  df <- NULL
  for (i in 1:length(fusionfs)) {
    res <- fread(fusionfs[i], header = T)
    res$FDR <- p.adjust(res$TWAS.P, method = "bonferroni")
    res$ifcausal <- ifelse(res$ID %in% cau[[i]], 1, 0)
    res$runtag <- i
    df <- rbind(df, res)
  }

  fig <- nca_plot(df$FDR, df$ifcausal, df$runtag, mode ="FDR", xmin = 0.2, main = main)
  return(fig)
}

# power plot FUSION FDR
ncausalFUSIONp_plot <- function(phenofs, fusionfs, main = "FUSION FDR"){

  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})

  df <- NULL
  for (i in 1:length(fusionfs)) {
    res <- fread(fusionfs[i], header = T)
    res$FDR <- p.adjust(res$TWAS.P, method = "fdr")
    res$ifcausal <- ifelse(res$ID %in% cau[[i]], 1, 0)
    res$runtag <- i
    df <- rbind(df, res)
  }

  fig <- nca_plot(df$FDR, df$ifcausal, df$runtag, mode ="FDR", xmin = 0.2, main = main)
  return(fig)
}


# power plot coloc PP4
ncausalPP4_plot <- function(phenofs, colocfs, twas.p = 1e-4, main = "coloc PP4"){

  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})

  df <- NULL
  for (i in 1:length(colocfs)) {
    res <- fread(colocfs[i], header = T)
    res$ifcausal <- ifelse(res$ID %in% cau[[i]], 1, 0)
    res$runtag <- i
    res <- res[complete.cases(res),]
    res <- res[res$TWAS.P < twas.p,]
    df <- rbind(df, res)
  }

  fig <- nca_plot(df$COLOC.PP4, df$ifcausal, df$runtag, mode ="PIP", xmin = 0.5, main = main)
  return(fig)
}
