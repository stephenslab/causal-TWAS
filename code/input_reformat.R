library(data.table) # requires v.1.12.7 or higher

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


