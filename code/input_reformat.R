library(data.table) # requires v.1.12.7 or higher
library(bigstatsr)
library(bigreadr)
library(pgenlibr) #install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr"), Notice that the installation of pgenlibr requires zstd(>=1.4.4). It can be built from source or simply available from conda, pip or brew.

source(paste0(codedir, "stats_func.R"))

write.gemma.pheno <- function(outfile, y) {
  if (is.numeric(y))
    y <- round(y,digits = 6)
  write.table(y, outfile,quote = FALSE,row.names = FALSE,col.names = FALSE)
}

write.gemma.geno <- function(outfile, geno, labels, append=F) {
  geno <- t(geno)
  geno <- data.table(geno)
  geno <- round(geno,digits = 3)
  geno <- cbind(paste0("rs",labels),"A","B",geno) # fake allele for now
  fwrite(geno, outfile, col.names = FALSE, sep = " ", quote = FALSE, row.names = FALSE, append=append)
  rm(geno);gc()
}

write.gemma.geno.chunks <- function(outfile, geno, labels, n_per_chunk){

	l = split(1:ncol(geno), ceiling(seq_along(1:ncol(geno))/n_per_chunk))

	# Write large data frame to csv in chunks
	count = 1
	for(i in l){
	  cat(paste0("writing chunks: ",count ,"/",length(l),"...\n"))
	  if (count==1){
		  write.gemma.geno(outfile, geno[,i], labels[i], append=F)
	  }
	  else{
		  write.gemma.geno(outfile, geno[,i], labels[i], append=T)
	  }
	  count = count + 1
	}
}

gemma.prepare <- function(Rdatafile,
                       outnamebase, chromosomes = NULL) {
  # read in Rdata file
  cat("* loading data.\n")
  load(Rdatafile)
  nsample <- dim(X)[1]
  nSNP <- dim(X)[2]
  save(nsample,nSNP,file= paste0(outnamebase,".QC.RData"))

  # Write the phenotype data
  cat("* Writing phenotype to file.\n")
  phenofile <- paste0(outnamebase,".pheno.txt")
  write.gemma.pheno(phenofile, y)

  # Write genotype data
  cat("* Writing genotype to file.\n")
  genofile <- paste0(outnamebase,".geno.txt.gz")
  write.gemma.geno.chunks(genofile, X, labels, N_per_chunk)
}

plink2Rd <- function(plinkfile, Rdfile){
  # plinkfile: in .traw.gz format and .raw.gz format
  # Note: this way is more memory efficient and faster than plink2R package.
  dat <- list()

  print("reading .traw.gz file")
  dat0t <- fread(paste0(plinkfile, ".traw.gz"), select = c(1:6))
  dat$chr <- as.matrix(dat0t[ , "CHR", with=F])
  dat$pos <- as.matrix(dat0t[ , "POS", with=F])
  dat$snp <- as.matrix(dat0t[ , "SNP", with=F])
  dat$counted <- as.matrix(dat0t[ , "COUNTED" , with=F])
  dat$alt <- as.matrix(dat0t[ , "ALT", with=F])

  print("reading .raw.gz file")
  dat0 <- fread(paste0(plinkfile, ".raw.gz"))
  dat0[ , 1:6] <- NULL
  colnames(dat0) <- NULL
  dat$G <- as.matrix(dat0)

  save(dat, file = Rdfile)
}

################################################################################

cut_in_nb <- function(x, nb) {
  split(x, sort(rep_len(seq_len(nb), length(x))))
}

drop_ext <- function(file) tools::file_path_sans_ext(file)

cbind_df <- function(list_df) {
  do.call(cbind, list_df)
}

################################################################################

#' Read from pgen to FBM/BM
#' Adapted from https://github.com/privefl/bigstatsr/blob/master/R/read-write.R
#' Read a file as a Filebacked Big Matrix by using package {pgenlibr} {bigreadr}.
#'
#' @param pfile plink files to read. Assume presence of pfile.pgen/pvar/psam
#' @param select Indices of columns to read (sorted).
#'   The length of `select` will be the number of columns of the resulting FBM.
#' @param scale T/F whether to scale genotype matrix
#' @param filter Vector used to subset the rows of each data frame.
#' @param backingfile Path to the file storing the FBM data on disk.
#'   An extension ".bk" will be automatically added.
#'   Default uses `file` without its extension.
#' @inheritDotParams bigreadr::big_fread2 -file -select -.transform -.combine
#' @inheritParams FBM
#'
#' @return A Filebacked Big Matrix of type `type` with `length(select)` columns.

pgen2fbm <- function(pfile, select=NULL, scale = F,
                     filter = NULL,
                     type = c("double", "float", "integer",
                              "unsigned short", "unsigned char", "raw"),
                     backingfile = drop_ext(pfile),
                     ...) {

  # Prepare reading
  X <- NULL
  offset <- 0

  pfile <- drop_ext(pfile)
  vars <- data.table::fread(paste0(pfile, ".pvar"))

  # Selected columns

  if (is.null(select)) {
    nb_cols <- dim(vars)[1]
    select <- seq_len(nb_cols)
  } else {
    bigassertr::assert_int(select); bigassertr::assert_pos(select)
    if (is.unsorted(select, strictly = TRUE))
      stop2("Argument 'select' should be sorted.")
  }

  fill_FBM <- function(mtx) {

    # Filter rows
    if (!is.null(filter)) mtx <- mtx[filter, , drop = FALSE]

    # Scale option
    if (isTRUE(scale)) {
      mtx <- scaleRcpp(mtx)
    }

    # Initialize FBM on first round
    if (offset == 0) {
      # Resulting FBM
      X <<- FBM(nrow(mtx), length(select), type = type, init = NULL,
                backingfile = backingfile, create_bk = TRUE)$save()
    }

    # Fill part of the FBM
    ind <- cols_along(mtx)
    X[, offset + ind] <- mtx
    offset <<- offset + length(ind)

    return(NULL)
  }

  # Read and fill by parts
  qc <- pgen2mtx(
    pfile, select = select, .transform = fill_FBM, progress = T, part_size = 50 * 1024^2, ... = ...)
  rm(X);gc()

  X <- readRDS(paste0(backingfile, ".rds"))


  # Verify scaling results
  if (isTRUE(scale)){
    cat("After scaling, mean range:", range(big_scale()(X)[,"center"]), "\n")
    cat("After scaling, SD range:", range(big_scale()(X)[,"scale"]), "\n")
  }

  # Save
  if (!is.null(select)){
    vars <- vars[select]
  }

  dat <- list("X"       = X,
              "chr"     = as.matrix(vars[,"#CHROM"]),
              "pos"     = as.matrix(vars[,"POS"]),
              "snp"     = as.matrix(vars[,"ID"]),
              "counted" = as.matrix(vars[,"REF"]),
              "alt"     = as.matrix(vars[,"ALT"]))

  save(dat, file = paste0(backingfile, ".FBM.Rd"))
}


#' Read from pgen to R matrix format. read by columns. can specify chunks.
pgen2mtx <- function(pfile, nb_parts = NULL,
                       .transform = identity,
                       .combine = cbind_df,
                       select = NULL,
                       progress = FALSE,
                       part_size = 500 * 1024^2,  ## 500 MB
                       meanimpute = F,
                       ...) {

  pfile <- drop_ext(pfile)
  pvar <- pgenlibr::NewPvar(paste0(pfile, ".pvar"))
  pgen <- pgenlibr::NewPgen(paste0(pfile, ".pgen"), pvar = pvar)

  ## Split selected columns in nb_parts
  if (is.null(select)) {
    nb_cols <- pgenlibr::GetVariantCt(pgen)
    select <- seq_len(nb_cols)
  } else {
    bigassertr::assert_int(select); bigassertr::assert_pos(select)
    if (is.unsorted(select, strictly = TRUE))
      stop2("Argument 'select' should be sorted.")
  }
  # Number of parts
  if (is.null(nb_parts)) {
    nb_parts <- ceiling(file.size(paste0(pfile, ".pgen")) / part_size)
    if (progress) bigassertr::message2("Will read the file in %d parts.", nb_parts)
  }
  split_cols <- cut_in_nb(select, nb_parts)

  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = length(select), style = 3)
    on.exit(close(pb), add = TRUE)
  }

  ## Read + transform other parts
  already_read <- 0
  all_parts <- lapply(split_cols, function(cols) {
    part <- .transform(
      pgenlibr::ReadList(pgen, cols, meanimpute = meanimpute)
    )
    already_read <<- already_read + length(cols)
    if (progress) utils::setTxtProgressBar(pb, already_read)
    part
  })
  all_parts <- unname(all_parts)

  ## Combine
  tryCatch(.combine(all_parts), error = function(e) {
    warning2("Combining failed. Returning list of parts instead..")
    all_parts
  })
}
