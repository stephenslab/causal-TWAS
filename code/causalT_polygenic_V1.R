run.gemma <- function (phenotype, covariates, pheno, geno, map,
                       gemmadir, gemma.exe, chromosomes = NULL) {

  # Set the local directory to the location of the GEMMA files.
  srcdir <- getwd()
  setwd(gemmadir)


    # Compute the kinship matrix using all markers that are *not* on
    # the chromosome.
    cat("Mapping QTLs on chromosome ",chr,".\n",sep="")
    cat(" * Computing kinship matrix.\n")
    markers <- which(map$chr != chr)
    K <- tcrossprod(center.columns(geno[,markers])) / length(markers)

    # Save the kinship matrix to a text file.
    cat(" * Writing kinship matrix to file.\n")
    write.table(round(K,digits = 6),paste0(gemmadir,"/kinship.txt"),sep = " ",
                quote = FALSE,row.names = FALSE,col.names = FALSE)

    # Write out the mean genotypes and map information for all markers
    # on the chromosome.
    markers <- which(map$chr == chr)
    cat(" * Writing to file genotypes for ",length(markers),
        " markers on chromosome ",chr,".\n",sep="")
    write.gemma.geno(paste0(gemmadir,"/geno.txt"),geno[,markers],map[markers,])
    cat(" * Writing genetic map for",length(markers),
        "markers on chromosome",chr,"to file.\n")
    write.gemma.map(paste0(gemmadir,"/map.txt"),map[markers,])

    # Now we are finally ready to run GEMMA for all markers on the
    # chromosome using the kinship matrix computed using all the
    # markers *not* on the chromosome.
    cat(" * Computing p-values for ",length(markers),
        " markers on chromosome ",chr,".\n",sep="")
    system(paste(gemma.exe,"-g geno.txt -a map.txt -p pheno.txt",
                 "-c covariates.txt -k kinship.txt -notsnp -lmm 2",
                 "-lmin 0.01 -lmax 100"),
           ignore.stdout = TRUE)

    # Load the results of the GEMMA association analysis.
    scans[[chr]] <-
      read.gemma.assoc(paste0(gemmadir,"/output/result.assoc.txt"))

  # Return the genome-wide scan.
  return(gwscan)
}
