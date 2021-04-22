source("~/causalTWAS/causal-TWAS/analysis/summarize_basic_plots.R")

# coloc PP4 calibration plot
califocusPIP_plot <- function(phenofs, focusfs, main = "FOCUS PIP Calibration"){
  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})
  df <- NULL
  for (i in 1:length(focusfs)) {
    res <- fread(focusfs[i], header = T)
    res$runtag <- i
    res <- res[res$mol_name != "NULL",]
    res$ifcausal <- ifelse(res$mol_name %in% cau[[i]], 1, 0)
    df <- rbind(df, res)
  }

  fig <- cp_plot(df$pip, df$ifcausal, df$runtag, mode ="PIP", main = main)
  return(fig)
}


ncausalfocusPIP_plot <- function(phenofs, focusfs, main = "FOCUS PIP"){

  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})

  df <- NULL
  for (i in 1:length(focusfs)) {
    res <- fread(focusfs[i], header = T)
    res$ifcausal <- ifelse(res$mol_name %in% cau[[i]], 1, 0)
    res$runtag <- i
    res <- res[res$mol_name != "NULL",]
    df <- rbind(df, res)
  }

  fig <- nca_plot(df$pip, df$ifcausal, df$runtag, mode ="PIP", xmin = 0.5, main = main)
  return(fig)
}
