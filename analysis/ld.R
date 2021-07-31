
get_geno <- function(id, pgen, pvar, exprf, exprvar){
  if (id %in% exprvar$id){
    idx <- which(exprvar$id == id)
    return(read_expr(exprf, variantidx = idx))
  } else if (id %in% pvar$id){
    idx <- which(pvar$id == id)
    return(read_pgen(pgen, variantidx = idx))
  } else{
    stop("no such id")
  }
}

get_ld <- function(ids, phenores, pgenfs, exprfs, chrom){
  pvarf <- prep_pvar(pgenfs[chrom])
  pgen <- prep_pgen(pgenf = pgenfs[chrom],pvarf = pvarf)
  pvar <- read_pvar(pvarf)
  exprf <- exprfs[chrom]
  exprvarf <- prep_exprvar(exprfs[chrom])
  exprvar <- read_exprvar(exprvarf)

  geno.t <- lapply(ids, get_geno, pgen = pgen, pvar = pvar, exprf = exprf, exprvar = exprvar)
  ids.cau <- c(phenores$batch[[chrom]]$id.cSNP, phenores$batch[[chrom]]$id.cgene)
  geno.cau <-do.call(cbind, lapply(ids.cau, get_geno, pgen = pgen, pvar = pvar, exprf = exprf, exprvar = exprvar))
  r2max <- unlist(lapply(geno.t, function(x) max(cor(x=x, y=geno.cau))))
  return(r2max)
}

