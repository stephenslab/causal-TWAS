GEMMA = "gemma"

run.ctwas <- function (pheno, geno, expr) {

  cat(" * Computing kinship matrix for snps.\n")
  GEMMA

  cat(" * Computing correlation matrix for expression.\n")


  # Estimate variance component by gemma
  system(paste(GEMMA,"-g geno.txt -a map.txt -p pheno.txt",
               "-c covariates.txt -k kinship.txt -notsnp -lmm 2",
               "-lmin 0.01 -lmax 100"),
         ignore.stdout = TRUE)

}
