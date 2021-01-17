library(data.table)
library(logging)
library(tools)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("3 arguments must be supplied:
       * genotype file name (or multiple file names in a .txt file)
       * weight dir
       * out file name", call.=FALSE)
}

addHandler(writeToFile, file="run_impute_expr.R.log", level='DEBUG')
loginfo('script started ... ')

loginfo("input arg 1 (genotype): %s ", args[1])
loginfo("input arg 2 (weight): %s ", args[2])
loginfo("input arg 3 (outname): %s ", args[3])


codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir,"input_reformat.R"))
source(paste0(codedir,"cis_expr.R"))


outname <- args[3]
exprtxt <- file(paste0(outname, "-cis-expr.txt"), 'w')

# load genotype data
pfile <- args[1]

if (file_ext(pfile) == "txt"){
  pfiles <- read.table(pfile, header = F, stringsAsFactors = F)[,1]
  outnames <- paste0(outname, "-B", 1:length(pfiles))
} else {
  outnames <- outname
  pfiles <- pfile
}
pfileRds <- paste0(pfiles, ".FBM.Rd")

loginfo("start to impute expression ...")
for ( b in 1:length(pfileRds)){
  # load genotype
  pfileRd <- pfileRds[b]
  outname <- outnames[b]
  load(pfileRd)

  # impute expression
  loginfo("generating imputed expression file: %s",  paste0(outname, "-cis-expr.Rd"))
  weight <- args[2]
  exprres <- cis_expr(dat, weight, method = "lasso", checksnps = F)

  # normalize and covert to FBM
  expr <- scaleRcpp(exprres$expr)
  exprres$expr <- as_FBM(expr, backingfile = paste0(outname, "-cis-expr.FBM.scaled"), is_read_only = T)$save()

  #save
  save(exprres, file = paste0(outname, "-cis-expr.Rd"))
  write(file_path_as_absolute(paste0(outname, "-cis-expr.Rd")), file = exprtxt, append = T)
 }

