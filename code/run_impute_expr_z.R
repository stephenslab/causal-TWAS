library(data.table)
library(logging)
library(tools)

# impute z score using summary statistics

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("4 arguments must be supplied:
       * LD reference genotype file name (or multiple file names in a .txt file)
       * weight dir
       * sumstats (snp z scores)
       * out file name", call.=FALSE)
}

addHandler(writeToFile, file="run_impute_expr_z.R.log", level='DEBUG')
loginfo('script started ... ')

loginfo("input arg 1 (genotype): %s ", args[1])
loginfo("input arg 2 (weight): %s ", args[2])
loginfo("input arg 3 (sumstats): %s ", args[3])
loginfo("input arg 4 (outname): %s ", args[4])


codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))
source(paste0(codedir,"input_reformat.R"))
source(paste0(codedir,"impute_expr_z.R"))

weight <- args[2]

zdf <- fread(args[3], header =T)
setnames(zdf, old = c('t.value','MARKER_ID'), new = c('z','snp'))
zdf <- data.frame(zdf)

outname <- args[4]
exprtxt <- file(paste0(outname, "-expr-z.txt"), 'w')

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

  # impute expression z scores
  loginfo("generating imputed expression file: %s",  paste0(outname, "-expr-z.Rd"))
  exprres <- impute_expr_z(dat, weight, zdf, method = "lasso", checksnps = F)

  # covert to FBM
  exprres$expr <- as_FBM(exprres$expr, backingfile = paste0(outname, "-expr-z.FBM"), is_read_only = T)$save()

  # save
  save(exprres, file = paste0(outname, "-expr-z.Rd"))
  write(file_path_as_absolute(paste0(outname, "-expr-z.Rd")), file = exprtxt, append = T)
 }

