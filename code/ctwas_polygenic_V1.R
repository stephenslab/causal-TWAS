
run.ctwas <- function (pheno, geno, expr, G.expr, outname) {
  # pheno: phenotype file, in BIMBAM phenotype format
  # geno:  genotype file, in BIMBAM mean genotype format
  # expr:  expression matrix, column genes row samples
  # G.expr: genotype matrix for expresion data, will be used for training

  cat(" * Computing kinship matrix for snps.\n")
  genodir <-  normalizePath(dirname(geno))
  genokdir <- paste0(genodir, "/", basename(geno), ".kdir/")
  genof <- basename(geno)
  if (!dir.exists(genokdir)){
    dir.create(genokdir)
    system(paste("cd", genokdir, "; gemma -g", geno, "-p", pheno, "-gk 1 -o", paste0(genof,"-geno"))) # currently this step requires 80G memory for 3 hours.
  }

  cat(" * Training expression model.\n")

  cat(" * Imputing cis-expression (X tilde).\n")


  cat(" * Computing correlation matrix for expression.\n")
  exprdir <-  normalizePath(dirname(expr))
  exprkdir <- paste0(exprdir, "/", basename(expr), ".kdir/")
  exprf <- basename(expr)
  if (!dir.exists(exprkdir)){
    dir.create(exprkdir)
    system(paste("cd", exprkdir, "; gemma -g", expr, "-p", pheno, "-gk 1 -notsnp -o", paste0(exprf,"-expr")))
  }

  # Estimate variance component by gemma
  cat(" * Estimating variance components.\n")
  vcdir <- paste0(outname,"_VC/")
  if (!dir.exists(vcdir)){
    dir.create(vcdir)
    vcdir <- normalizePath(vcdir)
    mfile <- paste0(vcdir,"/", outname,"-Kfiles.txt")
    file1 <- paste0(exprkdir,"/output/", exprf, "-expr.cXX.txt")
    file2 <- paste0(genokdir,"/output/", genof, "-geno.cXX.txt")
    fileConn <- file(mfile)
    writeLines(c(file1, file2), fileConn)
    close(fileConn)
    system(paste("cd", vcdir, "; gemma -p", pheno, "-mk", mfile,  "-n 1 -vc 1 -o", outname))
  }
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("4 argument must be supplied:  pheno, geno, expr, outname", call.=FALSE)
}

phenofile <-  normalizePath(args[1])
genofile <- normalizePath(args[2])
exprfile <- normalizePath(args[3])

run.ctwas(phenofile, genofile, exprfile, args[4])
