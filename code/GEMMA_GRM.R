
write.gemma.pheno <- function (outfile, y) {
  if (is.numeric(y))
    y <- round(y,digits = 6)
  write.table(y,file,quote = FALSE,row.names = FALSE,col.names = FALSE)
}

write.gemma.geno <- function (outfile, geno, labels) {
  geno <- t(geno)
  geno <- as.data.frame(geno,check.names = FALSE)
  geno <- round(geno,digits = 3)
  geno <- cbind(paste0("rs",labels),"A","B",geno) # fake allele for now
  write.table(geno,file,sep = " ",quote = FALSE,row.names = FALSE,
              col.names = FALSE)
}

run.gemma <- function (Rdatafile,
                       outnamebase, gemma.exe, chromosomes = NULL) {
  # read in Rdata file
  load(Rdatafile)
  nsample <- dim(X)[1]
  nSNP <- dim(X)[2]
  save(nsample,nSNP,file= paste0(outnamebase,".QC.txt"))
	
  # Write the phenotype data 
  cat("* Writing phenotype to file.\n")
  phenofile <- paste0(outnamebase,".pheno.txt")
  write.gemma.pheno(phenofile, y) 

  # Write genotype data
  cat("* Writing genotype to file.\n")
  genofile <- paste0(outnamebase,".geno.txt")
  write.gemma.geno(genofile), X)
	
    # Now we are finally ready to run GEMMA for all markers on the
    # chromosome using the kinship matrix computed using all the
    # markers *not* on the chromosome.
    cat(" * GRM\n")
    system(paste(gemma.exe,"-g geno.txt -a map.txt -p pheno.txt",
                 "-c covariates.txt -k kinship.txt -notsnp -lmm 2",
                 "-lmin 0.01 -lmax 100"),
           ignore.stdout = TRUE)
      
    # Load the results of the GEMMA association analysis.
    scans[[chr]] <-
      read.gemma.assoc(paste0(gemmadir,"/output/result.assoc.txt"))
  }

  # Restore the working directory.
  setwd(srcdir)
  
  # Merge the mapping results from all chromosomes into a single table.
  gwscan           <- do.call(rbind,scans)
  rownames(gwscan) <- do.call(c,lapply(scans,rownames))
  class(gwscan)    <- c("scanone","data.frame")

  # Return the genome-wide scan.
  return()
}
