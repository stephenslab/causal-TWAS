#' Filter regions
#' require global variables: `phenores`, `exprres`, `dat`, `cau`
#' and cutoff for different filter types.
#' @param reg a matrix/data.frame, with columns chr, p0, p1
#' @param filtertype "nofilter", "mr.ash", "mr.ash2s", "p2", "pgene", "rBF"
#' @param filterfile, for "nofilter", filterfile is not used
filter_region <- function(reg, filtertype, filterfile){

  # prep/filter regions
  anno <- rbind(cbind(dat$chr, dat$pos, dat$snp, "SNP"),
                cbind(exprres$chrom, exprres$p0, exprres$gnames, "gene"))
  colnames(anno) <- c("chrom", "p0", "name", "type")

  regions <- NULL
  chroms <- unique(reg$chrom)

  if (filtertype == 'nofilter'){
    for (chrom in chroms){
      itv <- cbind(reg[reg$chrom == chrom, ]$p0, reg[reg$chrom == chrom, ]$p1)
      nCausal <- apply(itv, 1, function(x) {gnames <- anno[x[1] < anno[, "p0"] & anno[, "p0"] < x[2], "name"]; sum(ifelse(gnames %in% cau, 1, 0))})
      itv <- data.frame(chrom, itv, nCausal, stringsAsFactors = F)
      regions <- rbind(regions, itv)
    }
    colnames(regions) <- c("chrom", "p0", "p1", "nCausal")
    regions.f <- regions

  } else if (filtertype %in% c("mrash2s", "mrash")){

    gres <- paste0(filterfile, ".expr.txt")
    pipres <- read.table(gres, header =T)
    if (filtertype == "mrash2s"){
      sres <- paste0(filterfile, ".snp.txt")
      pipres <- rbind(read.table(gres, header =T), read.table(sres, header =T))
    }
    for (chrom in chroms){
      pipres.chr <- pipres[pipres$chr == chrom, ]
      itv <- cbind(reg[reg$chrom == chrom, ]$p0, reg[reg$chrom == chrom, ]$p1)
      rPIP <- apply(itv, 1, function(x) sum(pipres.chr[x[1] < pipres.chr$p0 & pipres.chr$p0 < x[2], "PIP"]))
      nCausal <- apply(itv, 1, function(x) {gnames <- pipres.chr[x[1] < pipres.chr$p0 & pipres.chr$p0 < x[2], "name"]; sum(ifelse(gnames %in% cau, 1, 0))})
      itv <- data.frame(chrom, itv, rPIP, nCausal, stringsAsFactors = F)
      regions <- rbind(regions, itv)
    }
    colnames(regions) <- c("chrom", "p0", "p1", "rPIP", "nCausal")
    regions.f <- regions[regions[, "rPIP"] > PIPfilter,  ]

  } else if (filtertype %in% c("pgene", "p2")){

    gres <- paste0(filterfile,".exprgwas.txt.gz")
    gres <- read.table(gres, header = T, comment.char = "")
    gres[, "prank"] <- rank(gres$PVALUE)/nrow(gres)
    pres <- gres
    if (filtertype =="p2"){
      sres <- paste0(filterfile,".snpgwas.txt.gz")
      sres <- read.table(sres, header =T, comment.char = "")
      sres[, "prank"] <- apply(cbind(rank(sres$PVALUE)/nrow(sres) * prankfilter.gene/prankfilter.snp, 1), 1, min)
      pres <- rbind(gres, sres)
    }
    for (chrom in chroms){
      pres.chr <- pres[pres$X.CHROM == chrom, ]
      itv <- cbind(reg[reg$chrom == chrom, ]$p0, reg[reg$chrom == chrom, ]$p1)
      rprank <- apply(itv, 1, function(x) min(c(pres.chr[x[1] < pres.chr$BEGIN & pres.chr$BEGIN < x[2], "prank"],1)))
      nCausal <- apply(itv, 1, function(x) {gnames <- pres.chr[x[1] < pres.chr$BEGIN & pres.chr$BEGIN < x[2], "MARKER_ID"]; sum(ifelse(gnames %in% cau, 1, 0))})
      rpmin <- apply(itv, 1, function(x) min(pres.chr[x[1] < pres.chr$BEGIN & pres.chr$BEGIN < x[2], "PVALUE"]))
      itv <- data.frame(chrom, itv, rpmin, rprank, nCausal, stringsAsFactors = F)
      regions <- rbind(regions, itv)
    }
    colnames(regions) <- c("chrom", "p0", "p1", "rpmin", "rprank", "nCausal")
    regions.f <- regions[regions[, "rprank"] <= prankfilter.gene,  ]

  } else if (filtertype == "rBF") {

    gres <- paste0(filterfile,".exprgwas.txt.gz")
    gres <- read.table(gres, header = T, comment.char = "")
    gres[, "BF"] <- apply(gres[, c("Estimate", "Std.Error")], 1, function(x) gwasbf(x[1], x[2], w = 0.02**2*0.1))
    sres <- paste0(filterfile,".snpgwas.txt.gz")
    sres <- read.table(sres, header =T, comment.char = "")
    sres[, "BF"] <- apply(sres[, c("Estimate", "Std.Error")], 1, function(x) gwasbf(x[1], x[2], w = 0.02**2*0.1))
    pres <- rbind(gres, sres)
    for (chrom in chroms){
      pres.chr <- pres[pres$X.CHROM == chrom, ]
      itv <- cbind(reg[reg$chrom == chrom, ]$p0, reg[reg$chrom == chrom, ]$p1)
      rBF <- apply(itv, 1, function(x) mean(pres.chr[x[1] < pres.chr$BEGIN & pres.chr$BEGIN < x[2], "BF"]))
      nCausal <- apply(itv, 1, function(x) {gnames <- pres.chr[x[1] < pres.chr$BEGIN & pres.chr$BEGIN < x[2], "MARKER_ID"]; sum(ifelse(gnames %in% cau, 1, 0))})
      itv <- data.frame(chrom, itv, rBF, nCausal, stringsAsFactors = F)
      regions <- rbind(regions, itv)
    }
    colnames(regions) <- c("chrom", "p0", "p1", "rBF", "nCausal")
    regions.f <- regions[regions[, "rBF"] > rBFfilter,  ]
    regions.f <- regions.f[complete.cases(regions.f),]

  } else if (filtertype == "rfdr") {

    gres <- paste0(filterfile,".exprgwas.txt.gz")
    gres <- read.table(gres, header = T, comment.char = "")
    sres <- paste0(filterfile,".snpgwas.txt.gz")
    sres <- read.table(sres, header =T, comment.char = "")
    pres <- rbind(gres, sres)
    pres[, "fdr"] <- p.adjust(pres$PVALUE, "BH")

    for (chrom in chroms){
      pres.chr <- pres[pres$X.CHROM == chrom, ]
      itv <- cbind(reg[reg$chrom == chrom, ]$p0, reg[reg$chrom == chrom, ]$p1)
      rfdr <- apply(itv, 1, function(x) min(c(pres.chr[x[1] < pres.chr$BEGIN & pres.chr$BEGIN < x[2], "fdr"],1)))
      nCausal <- apply(itv, 1, function(x) {gnames <- pres.chr[x[1] < pres.chr$BEGIN & pres.chr$BEGIN < x[2], "MARKER_ID"]; sum(ifelse(gnames %in% cau, 1, 0))})
      itv <- data.frame(chrom, itv, rfdr, nCausal, stringsAsFactors = F)
      regions <- rbind(regions, itv)
    }
    colnames(regions) <- c("chrom", "p0",  "p1",  "rfdr", "nCausal")
    regions.f <- regions[regions[, "rfdr"] <= rfdrfilter,  ]

  } else {

    stop("unknown filter type")
  }
  out <- list("all" = regions, "filtered"= regions.f)
  return(out)
}
